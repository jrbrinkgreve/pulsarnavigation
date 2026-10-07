function testChannelTOA()
%TESTCHANNELTOA  Unit tests for estimateTOA (A2a) and detectPulsar (A2b) on channelized folds.
%{
Run from the PulsarSimMatlab folder: run('tests/testChannelTOA.m'). ~5 s (measured).
Needs data/mc/toaRef_preA2.mat (tests/makeToaReference.m), the frozen folds
of data/mc/foldRef_pre3c.mat and the 128-channel power file
data/chan/test_dedisp_chan_power.dat (tests/testDetectChannels.m).

  1. Regression: on single-channel (full-band) folds estimateTOA reproduces
     the pre-A2 outputs bit for bit (linear, 'offpulse', MinCoverage 0.3 with
     gap filling, nearest, 3 turns per sub-int).
  2. Algebra: phaseErr and ampErr equal sqrt(d' * Sigma * d) / |C''| and
     sqrt(c' * Sigma * c) / sum|S|^2 with Sigma the full covariance matrix of
     the combined profile, built here channel by channel from the fold
     weights (Sigma = sum_c D_c Q_c D_c, Q_c banded from weight2 / weightX,
     D_c = diag(m_c / sqrt(Bnoise*binDt) / W_c)) and the channels chosen by
     the exclusion rule. Cases: (i) 128 channels with NoiseCoeffs (shared
     weights, 3 lags); (ii) per-channel data weights with fake blanking:
     5 % random per channel, and channels 1-10 blanked at phase 0.40-0.42 in
     every other turn (so excluded from those sub-ints). Plus: Bnoise as one
     value per channel gives the same result as the scalar.
  3. Physics, seed-43 data:
     (a) 128-channel fold vs the fold of the summed channel power (the old
         way, nChan = 1, Bnoise of the whole band): the same profile, so the
         same TOAs (pass < 1e-4 sigma; the float32 summed file alone moves
         them by ~1e-6 sigma); error bars and red. chi^2 reported (the old
         noise model is ~1.3 % too high per phase bin, ~0 for TOAs).
     (b) blanked vs unblanked: TOA differences consistent with zero (from
         their own scatter); error bars and channels used reported.
     Honest error bars over many noise realizations: A3.
  4. detectPulsar regression: single-channel folds bit-identical to the pre-A2
     outputs (linear, known phase 0.01, MinCoverage 0.3, nearest, 3 turns).
  5. detectPulsar algebra: T0 (known phase), normProfile and Tmax (at the best
     phase) equal the explicit formulas with the full H0 covariance matrix of
     the combined profile (each channel at its baseline a_c), for the folds of
     test 2.
  6. Physics: the H0 noise model per bin is honest on the channel path:
     noiseRatio (off-pulse variance / model) = 1 within 4 sigma (sigma from
     the scatter between sub-ints), unblanked and blanked; per-channel vs
     summed power: same best phase and detections (noise ratio of the old
     path ~1.3 % lower, as the fold level showed).
Errors at the end if any check fails. Writes ~0.5 GB of temporary files to
tempdir (deleted at the end).
%}

root = fullfile(fileparts(mfilename('fullpath')), '..');
addpath(fullfile(root, 'functions'));
oldDir = cd(root); restoreDir = onCleanup(@() cd(oldDir));
pipelineParams;
tmp = fullfile(tempdir, 'testChannelTOA');
if ~isfolder(tmp), mkdir(tmp); end
cleanTmp = onCleanup(@() delete(fullfile(tmp, '*.dat')));
fails = {};
template = gaussianTemplate(nBin, ephem.profileFWHM);

% ---------------------------------------------------------------------------------
% 1. Regression on single-channel folds
% ---------------------------------------------------------------------------------
R = load(fullfile(dataDir, 'mc', 'toaRef_preA2.mat')); tr = R.ref;
F = load(fullfile(dataDir, 'mc', 'foldRef_pre3c.mat')); fr = F.ref;
cs = {'A', fr.A, {}; 'Aoff', fr.A, {'NoiseModel', 'offpulse'}; ...
      'Acov', fr.A, {'MinCoverage', 0.3}; 'B', fr.B, {}; 'C', fr.C, {}};
for i = 1:size(cs, 1)
    [t1, i1] = estimateTOA(cs{i, 2}.fold, cs{i, 2}.info, template, 'Bnoise', tr.Bnoise, ...
        'Verbose', false, cs{i, 3}{:});
    same = sameFields(t1, tr.toa.(cs{i, 1})) && sameFields(i1, tr.toaInfo.(cs{i, 1}));
    used = all(t1.nChanUsed(t1.valid) == 1) && t1.total.nChanUsed == 1;
    fprintf('1. regression %-4s (%d TOAs): bit-identical %d, nChanUsed 1 %d: %s\n', cs{i, 1}, ...
        nnz(t1.valid), same, used, passStr(same && used));
    if ~(same && used), fails{end+1} = ['regression ' cs{i, 1}]; end %#ok<AGROW>
end

% ---------------------------------------------------------------------------------
% 2. Error bars vs the explicit covariance matrix of the combined profile
% ---------------------------------------------------------------------------------
detD = loadInfo(fullfile(dataDir, 'chan', 'test_dedisp_chan_power.dat'));
nc = detD.noise; nCh = detD.nChan; Nd = detD.N;
Bc = nc.Bnoise;                                      % noise bandwidth of one channel
base = {'F1', ephem.F1, 'TRef', ephem.TRef, 'NBin', nBin, 'SubintPeriods', subintPeriods, ...
        'SaveFile', false, 'Verbose', false};
q = {'Bnoise', Bc, 'Verbose', false};
[ifC, foldC] = foldProfile(detD, 'x', ephem.f0, base{:}, 'NoiseCoeffs', [nc.V, nc.X]);
toaC = estimateTOA(foldC, ifC, template, q{:});

% fake blanking (whole detected bins): 5 % random, channels 1-10 at phase 0.40-0.42 in odd turns
Pall = reshape(readF32(detD.file), nCh, Nd);
kAll = 1:Nd;
phAll = (detD.binTime0 + (kAll - 1)*detD.binDt - ephem.TRef) * ephem.f0;
turnK = round(phAll); phk = phAll - turnK;
rng(5);
keep = rand(nCh, Nd) >= 0.05;
keep(1:10, phk >= 0.40 & phk < 0.42 & mod(turnK, 2) == 1) = false;
detB = detD; detB.file = fullfile(tmp, 'p_blank.dat');
writeF32(detB.file, Pall .* single(keep));
coef = [nc.V, nc.X];
infoW = struct('file', fullfile(tmp, 'w_blank.dat'), 'nChan', nCh, 'N', Nd, 'Lmax', nc.Lmax, ...
    'byteOrder', 'ieee-le');
writeWeightBlocks(infoW.file, Nd, @(kk) blankStream(kk, keep, coef));
[ifB, foldB] = foldProfile(detB, 'x', ephem.f0, base{:}, 'DataWeights', infoW);
toaB = estimateTOA(foldB, ifB, template, q{:});

T = templateStuff(template);
cases = {'128 chan, NoiseCoeffs', foldC, ifC, toaC; '128 chan, blanked', foldB, ifB, toaB};
for i = 1:2
    [fo, ifo, to] = cases{i, 2:4};
    errT = 0; errA = 0; usedOk = true; nUsed = [];
    for s = find(to.valid & to.coverage == 1).'
        [eT, eA, use] = explicitErrors(fo, s, to.phase(s), Bc * ones(1, nCh), ifo.binDt, T);
        errT = max(errT, abs(eT - to.phaseErr(s)) / to.phaseErr(s));
        errA = max(errA, abs(eA - to.ampErr(s)) / to.ampErr(s));
        usedOk = usedOk && to.nChanUsed(s) == nnz(use);
        nUsed(end+1) = nnz(use); %#ok<AGROW>
    end
    pass = errT < 1e-9 && errA < 1e-9 && usedOk;
    fprintf(['2. %-22s %d sub-ints, %d lags: phaseErr vs explicit max rel %.1e, ampErr %.1e; ' ...
             'channels used %s (= rule %d): %s\n'], cases{i, 1}, numel(nUsed), size(fo.weightX, 4), ...
        errT, errA, mat2str(nUsed), usedOk, passStr(pass));
    if ~pass, fails{end+1} = ['explicit errors: ' cases{i, 1}]; end %#ok<AGROW>
end
oddTurnSubs = mod(foldB.subint.turnRef, 2) == 1;
exclOk = all(toaB.nChanUsed(toaB.valid & oddTurnSubs) <= nCh - 10) && ...
         all(toaB.nChanUsed(toaB.valid & ~oddTurnSubs) == nCh);
toaBv = estimateTOA(foldB, ifB, template, 'Bnoise', Bc * ones(nCh, 1), 'Verbose', false);
vecSame = isequaln(toaBv, toaB);
fprintf(['   exclusion: channels 1-10 left out in odd turns, all in even turns %d; ' ...
         'Bnoise per channel = scalar %d: %s\n'], exclOk, vecSame, passStr(exclOk && vecSame));
if ~(exclOk && vecSame), fails{end+1} = 'exclusion / Bnoise vector'; end

% ---------------------------------------------------------------------------------
% 3a. Per-channel fold vs the fold of the summed channel power (old way)
% ---------------------------------------------------------------------------------
detS = detD; detS.nChan = 1; detS.file = fullfile(tmp, 'p_sum.dat');
writeF32(detS.file, sum(double(Pall), 1));
[ifS, foldS] = foldProfile(detS, 'x', ephem.f0, base{:});
toaS = estimateTOA(foldS, ifS, template, 'Bnoise', detD.BnoiseTotal, 'Verbose', false);
v = toaC.valid & toaS.valid;
dTau = max(abs(toaC.phase(v) - toaS.phase(v)) ./ toaC.phaseErr(v));
fprintf(['3a. per-channel (%d chan) vs summed power: %d common TOAs, max |dTOA| %.1e sigma; ' ...
         'error ratio median %.4f, red. chi2 %.3f vs %.3f (Bnoise chan %.6g vs %.6g MHz): %s\n'], ...
    nCh, nnz(v), dTau, median(toaC.phaseErr(v) ./ toaS.phaseErr(v)), median(toaC.redChi2(v)), ...
    median(toaS.redChi2(v)), Bc/1e6, detD.BnoiseChan/1e6, passStr(dTau < 1e-4));
if ~(dTau < 1e-4), fails{end+1} = 'per-channel vs summed'; end

% ---------------------------------------------------------------------------------
% 3b. Blanked vs unblanked TOAs
% ---------------------------------------------------------------------------------
v = toaC.valid & toaB.valid;
dT = (toaB.toa(v) - toaC.toa(v)) * 1e6;
nV = nnz(v);
okb = abs(mean(dT)) < 4 * std(dT) / sqrt(nV);
fprintf(['3b. blanked (valid fraction %.3f) vs unblanked: %d TOAs, difference mean %+.3f +- %.3f us ' ...
         '(rms %.3f us), error ratio median %.3f (1/sqrt(valid fraction) %.3f), ' ...
         'channels used %s: %s\n'], ifB.validFraction, nV, mean(dT), std(dT)/sqrt(nV), rms(dT), ...
    median(toaB.toaErr(v) ./ toaC.toaErr(v)), 1/sqrt(ifB.validFraction), ...
    mat2str(toaB.nChanUsed(v).'), passStr(okb));
if ~okb, fails{end+1} = 'blanked TOAs'; end

% ---------------------------------------------------------------------------------
% 4. detectPulsar: regression on single-channel folds
% ---------------------------------------------------------------------------------
dcs = {'A', fr.A, {}; 'Aph', fr.A, {'Phase', 0.01}; 'Acov', fr.A, {'MinCoverage', 0.3}; ...
       'B', fr.B, {}; 'C', fr.C, {}};
for i = 1:size(dcs, 1)
    [d1, di1] = detectPulsar(dcs{i, 2}.fold, dcs{i, 2}.info, template, 'Bnoise', tr.Bnoise, ...
        'Verbose', false, dcs{i, 3}{:});
    same = sameFields(d1, tr.det.(dcs{i, 1})) && sameFields(di1, tr.detInfo.(dcs{i, 1}));
    used = all(d1.nChanUsed(d1.tested) == 1) && d1.total.nChanUsed == 1;
    fprintf('4. detectPulsar regression %-4s (%d tested): bit-identical %d, nChanUsed 1 %d: %s\n', ...
        dcs{i, 1}, nnz(d1.tested), same, used, passStr(same && used));
    if ~(same && used), fails{end+1} = ['detect regression ' dcs{i, 1}]; end %#ok<AGROW>
end

% ---------------------------------------------------------------------------------
% 5. detectPulsar vs the explicit H0 covariance matrix of the combined profile
% ---------------------------------------------------------------------------------
detC = detectPulsar(foldC, ifC, template, q{:});
detB = detectPulsar(foldB, ifB, template, q{:});
dcases = {'128 chan, NoiseCoeffs', foldC, ifC, detC; '128 chan, blanked', foldB, ifB, detB};
for i = 1:2
    [fo, ifo, de] = dcases{i, 2:4};
    e0 = 0; eZ = 0; eM = 0; usedOk = true; nT = 0;
    for s = find(de.tested & de.coverage == 1).'
        im = mod(round(de.phaseMax(s) * nBin), nBin) + 1;
        [T0e, ze, Tme, use] = explicitDetect(fo, s, Bc * ones(1, nCh), ifo.binDt, template, im);
        e0 = max(e0, abs(T0e - de.T0(s)) / abs(de.T0(s)));
        eZ = max(eZ, max(abs(ze - de.normProfile(:, s))) / max(abs(ze)));
        eM = max(eM, abs(Tme - de.Tmax(s)) / abs(de.Tmax(s)));
        usedOk = usedOk && de.nChanUsed(s) == nnz(use);
        nT = nT + 1;
    end
    pass = e0 < 1e-9 && eZ < 1e-9 && eM < 1e-9 && usedOk;
    fprintf(['5. %-22s %d sub-ints: T0 vs explicit max rel %.1e, normProfile %.1e, ' ...
             'Tmax %.1e; channels used = rule %d: %s\n'], dcases{i, 1}, nT, e0, eZ, eM, usedOk, ...
        passStr(pass));
    if ~pass, fails{end+1} = ['explicit detection: ' dcases{i, 1}]; end %#ok<AGROW>
end

% ---------------------------------------------------------------------------------
% 6. H0 noise model per bin honest (noiseRatio = 1); channel path vs summed power
% ---------------------------------------------------------------------------------
detS = detectPulsar(foldS, ifS, template, 'Bnoise', detD.BnoiseTotal, 'Verbose', false);
for i = 1:2
    de = dcases{i, 4};
    nr = de.noiseRatio(de.tested);
    m = mean(nr); sg = std(nr) / sqrt(numel(nr));
    ok6 = abs(m - 1) < 4*sg;
    fprintf('6. noise ratio, %-22s mean over %d sub-ints %.4f +- %.4f (total fold %.4f): %s\n', ...
        dcases{i, 1}, numel(nr), m, sg, de.total.noiseRatio, passStr(ok6));
    if ~ok6, fails{end+1} = ['noise ratio: ' dcases{i, 1}]; end %#ok<AGROW>
end
v = detC.tested & detS.tested;
samePh = all(abs(detC.phaseMax(v) - detS.phaseMax(v)) < 0.5/nBin) && ...
    isequal(detC.detectedKnown(v), detS.detectedKnown(v)) && ...
    isequal(detC.detectedUnknown(v), detS.detectedUnknown(v));
fprintf(['   per-channel vs summed power (%d sub-ints): same phaseMax and detections %d; T0 ratio ' ...
         'median %.4f; noise ratio %.4f vs %.4f (ratio %.4f; fold level A1a: 1.0134): %s\n'], ...
    nnz(v), samePh, median(detC.T0(v) ./ detS.T0(v)), mean(detC.noiseRatio(v)), ...
    mean(detS.noiseRatio(v)), mean(detC.noiseRatio(v) ./ detS.noiseRatio(v)), passStr(samePh));
if ~samePh, fails{end+1} = 'detection per-channel vs summed'; end

% ---------------------------------------------------------------------------------
if isempty(fails)
    fprintf('testChannelTOA: ALL PASSED\n');
else
    error('testChannelTOA:failed', 'FAILED: %s', strjoin(fails, ', '));
end
end


% =================================================================================
function T = templateStuff(template)
% Template quantities as estimateTOA defines them.
t = template(:) / max(template);
N = numel(t);
T.N = N; T.kk = (1:floor(N/2) - 1).'; T.Sx = fft(t); T.Sk = T.Sx(T.kk + 1);
T.sumS2 = sum(abs(T.Sk).^2); T.ksg = [0:ceil(N/2)-1, -floor(N/2):-1].'; T.w1 = 2*pi*T.kk;
end

function [eTau, eAmp, use] = explicitErrors(fold, s, tau, B, binDt, T)
% Error bars of phase and amplitude from the full covariance matrix of the
% combined profile of sub-int s (channels chosen by the exclusion rule).
N = T.N; nChan = size(fold.prof, 3);
W  = reshape(fold.weight(:, s, :), N, []);
W2 = reshape(fold.weight2(:, s, :), N, []);
WX = reshape(fold.weightX(:, s, :, :), N, size(W, 2), []);
if size(W, 2) == 1                                   % shared weights
    W = repmat(W, 1, nChan); W2 = repmat(W2, 1, nChan); WX = repmat(WX, 1, nChan, 1);
end
use = all(W > 0, 1);                                 % full coverage: data in every bin
Pc = reshape(fold.prof(:, s, :), N, nChan);
Pc = Pc(:, use);
p  = sum(Pc, 2);
e  = exp(1i * T.w1 * tau);
PC = fft(Pc);
bC = real(sum(PC(T.kk + 1, :) .* conj(T.Sk) .* e, 1)) / T.sumS2;
aC = (real(PC(1, :)) - bC * real(T.Sx(1))) / N;
shape = real(ifft(T.Sx .* exp(-2i*pi*T.ksg*tau)));
j = (1:N).';
Sig = sparse(N, N);
cs = find(use);
for ic = 1:numel(cs)
    c = cs(ic);
    Q = sparse(j, j, W2(:, c), N, N);
    for d = 1:size(WX, 3)
        jd = mod(j - 1 + d, N) + 1;
        Q = Q + sparse([j; jd], [jd; j], [WX(:, c, d); WX(:, c, d)], N, N);
    end
    sc = abs(aC(ic) + bC(ic) * shape) / sqrt(B(c) * binDt) ./ W(:, c);
    Dg = sparse(j, j, sc, N, N);
    Sig = Sig + Dg * Q * Dg;
end
P  = fft(p);
X  = P(T.kk + 1) .* conj(T.Sk);
c2 = real(sum(-(T.w1.^2) .* X .* e));
Dv = zeros(N, 1); Dv(T.kk + 1) = 1i * T.w1 .* conj(T.Sk) .* e; dvec = real(fft(Dv));
Cv = zeros(N, 1); Cv(T.kk + 1) = conj(T.Sk) .* e;              cvec = real(fft(Cv));
eTau = sqrt(dvec.' * Sig * dvec) / abs(c2);
eAmp = sqrt(cvec.' * Sig * cvec) / T.sumS2;
end

function [T0e, ze, Tme, use] = explicitDetect(fold, s, B, binDt, template, im)
% detectPulsar's statistics of sub-int s from the full H0 covariance matrix of
% the combined profile (each channel at its own baseline a_c). im: index of the
% best phase shift (template shifted by im-1 bins).
t = template(:) / max(template); N = numel(t); c = t - mean(t);
nChan = size(fold.prof, 3);
W  = reshape(fold.weight(:, s, :), N, []);
W2 = reshape(fold.weight2(:, s, :), N, []);
WX = reshape(fold.weightX(:, s, :, :), N, size(W, 2), []);
if size(W, 2) == 1
    W = repmat(W, 1, nChan); W2 = repmat(W2, 1, nChan); WX = repmat(WX, 1, nChan, 1);
end
use = all(W > 0, 1);
Pc = reshape(fold.prof(:, s, :), N, nChan);
Pc = Pc(:, use);
p  = sum(Pc, 2);
a  = mean(p);
aC = mean(Pc, 1);
j = (1:N).';
Sig = sparse(N, N);
cs = find(use);
for ic = 1:numel(cs)
    k = cs(ic);
    Q = sparse(j, j, W2(:, k), N, N);
    for d = 1:size(WX, 3)
        jd = mod(j - 1 + d, N) + 1;
        Q = Q + sparse([j; jd], [jd; j], [WX(:, k, d); WX(:, k, d)], N, N);
    end
    Dg = sparse(j, j, aC(ic) / sqrt(B(k) * binDt) ./ W(:, k), N, N);
    Sig = Sig + Dg * Q * Dg;
end
dd  = p - a;
T0e = (c.' * dd) / sqrt(c.' * Sig * c);              % known phase 0
ze  = dd ./ sqrt(full(diag(Sig)));
cm  = circshift(c, im - 1);
Tme = (cm.' * dd) / sqrt(cm.' * Sig * cm);
end

function ok = sameFields(new, old)
% Every field of old is bit-identical in new (new may have extra fields).
ok = true;
fn = fieldnames(old);
for i = 1:numel(fn)
    if ~isfield(new, fn{i}), ok = false; return; end
    if isstruct(old.(fn{i}))
        ok = ok && sameFields(new.(fn{i}), old.(fn{i}));
    else
        ok = ok && isequaln(new.(fn{i}), old.(fn{i}));
    end
end
end

function s = passStr(p)
if p, s = 'PASS'; else, s = 'FAIL'; end
end

function writeF32(file, p)
fid = fopen(file, 'w', 'ieee-le');
fwrite(fid, single(p), 'single');
fclose(fid);
end

function p = readF32(file)
fid = fopen(file, 'r', 'ieee-le');
c = onCleanup(@() fclose(fid));
p = fread(fid, Inf, 'single=>single');
end

function writeWeightBlocks(file, N, blockFun)
% Data-weight file [nChan x (2+Lmax) x N] (float32), written in blocks of bins.
fid = fopen(file, 'w', 'ieee-le');
c = onCleanup(@() fclose(fid));
blk = 8192;
for k0 = 1:blk:N
    fwrite(fid, blockFun(k0:min(N, k0 + blk - 1)), 'single');
end
end

function Y = blankStream(kk, keep, coef)
% Exact W, V, X(L) of time bins kk when whole detected bins are zeroed
% (as in tests/testFoldWeights.m).
[nCh, N] = size(keep);
Lmax = numel(coef) - 1;
kp = keep(:, kk);
Y = zeros(nCh, 2 + Lmax, numel(kk), 'single');
Y(:, 1, :) = reshape(kp, nCh, 1, []);
Y(:, 2, :) = reshape(coef(1) * kp, nCh, 1, []);
for L = 1:Lmax
    kn = false(nCh, numel(kk));
    inside = kk + L <= N;
    kn(:, inside) = keep(:, kk(inside) + L);
    Y(:, 2 + L, :) = reshape(coef(1 + L) * (kp & kn), nCh, 1, []);
end
end
