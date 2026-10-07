function testChannelTOA()
%TESTCHANNELTOA  Unit tests for estimateTOA on channelized folds (A2a).
%{
Run from the PulsarSimMatlab folder: run('tests/testChannelTOA.m'). ~4 s (measured).
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
