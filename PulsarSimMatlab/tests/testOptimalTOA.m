function testOptimalTOA()
%TESTOPTIMALTOA  Unit tests for estimateTOA 'Weighting', 'optimal' (A2c-1).
%{
Run from the PulsarSimMatlab folder: run('tests/testOptimalTOA.m').
Needs data/mc/foldRef_pre3c.mat, data/mc/toaRef_preA2.mat and the
128-channel power file data/chan/test_dedisp_chan_power.dat.

  1. Default unchanged: 'equal' (default) bit-identical to the pre-A2 outputs.
  2. Equivalence: one channel with NoiseModel 'offpulse' has uniform weights,
     so 'optimal' must give FFTFIT's tau, amplitude and baseline (to
     rounding); error bars agree on average (model vs data curvature).
  3. Algebra, 128 channels (unblanked with NoiseCoeffs; fake-blanked with
     per-channel weights, including channels with empty phase bins): from the
     outputs (tau, b, a_c) rebuild the model, the weights 1/var and the full
     covariance of every channel as sparse matrices; tau is a maximum of
     N/sqrt(Dn); b and a_c are the weighted least-squares values; phaseErr and
     ampErr = the sandwich A^-1 B A^-1 with these matrices.
  4. Monte Carlo: 800 profiles per SNR level drawn with the exact noise
     covariance (radiometer model, all lags) of 64 channels of a fake-blanked
     fold (5 of them with empty phase bins), true tau 0.0123 turns, optimal
     SNR ~10 and ~40. For 'equal' and 'optimal': bias within 4 sigma; pulls
     (tau - tau0)/sigma with std 1 within 4/sqrt(2R); 'optimal' scatter vs
     'equal' scatter = the predicted ratio of their error bars.
  5. The fake-blanked data: 'optimal' uses all channels (also those with empty
     phase bins); TOAs consistent with 'equal'.
Errors at the end if any check fails.
%}

root = fullfile(fileparts(mfilename('fullpath')), '..');
addpath(fullfile(root, 'functions'));
oldDir = cd(root); restoreDir = onCleanup(@() cd(oldDir));
pipelineParams;
tmp = fullfile(tempdir, 'testOptimalTOA');
if ~isfolder(tmp), mkdir(tmp); end
cleanTmp = onCleanup(@() delete(fullfile(tmp, '*.dat')));
fails = {};
template = gaussianTemplate(nBin, ephem.profileFWHM);
R = load(fullfile(dataDir, 'mc', 'toaRef_preA2.mat')); tr = R.ref;
F = load(fullfile(dataDir, 'mc', 'foldRef_pre3c.mat')); fr = F.ref;

% ---------------------------------------------------------------------------------
% 1. Default unchanged
% ---------------------------------------------------------------------------------
t1 = estimateTOA(fr.A.fold, fr.A.info, template, 'Bnoise', tr.Bnoise, 'Verbose', false);
same = sameFields(t1, tr.toa.A);
fprintf('1. default (equal) bit-identical to pre-A2: %d: %s\n', same, passStr(same));
if ~same, fails{end+1} = 'default'; end

% ---------------------------------------------------------------------------------
% 2. One channel, uniform weights: 'optimal' = FFTFIT
% ---------------------------------------------------------------------------------
q = {'Bnoise', tr.Bnoise, 'Verbose', false, 'NoiseModel', 'offpulse'};
tE = estimateTOA(fr.A.fold, fr.A.info, template, q{:});
tO = estimateTOA(fr.A.fold, fr.A.info, template, q{:}, 'Weighting', 'optimal');
v = tE.valid & tO.valid;
dTau = max(abs(tO.phase(v) - tE.phase(v)));
dB = max(abs(tO.amp(v) ./ tE.amp(v) - 1)); dA = max(abs(tO.baseline(v) ./ tE.baseline(v) - 1));
er = tO.phaseErr(v) ./ tE.phaseErr(v);
okE = abs(mean(er) - 1) < 4 * std(er) / sqrt(nnz(v));
pass = isequal(v, tE.valid) && dTau < 1e-12 && dB < 1e-10 && dA < 1e-10 && okE;
fprintf(['2. one channel, offpulse (uniform weights): %d TOAs; tau diff %.1e turns, amp %.1e, ' ...
         'baseline %.1e; error ratio optimal/FFTFIT %.4f +- %.4f: %s\n'], nnz(v), dTau, dB, dA, ...
    mean(er), std(er)/sqrt(nnz(v)), passStr(pass));
if ~pass, fails{end+1} = 'equivalence'; end

% ---------------------------------------------------------------------------------
% 3. Algebra on 128 channels (unblanked, and fake-blanked with empty phase bins)
% ---------------------------------------------------------------------------------
detD = loadInfo(fullfile(dataDir, 'chan', 'test_dedisp_chan_power.dat'));
nc = detD.noise; nCh = detD.nChan; Nd = detD.N;
base = {'F1', ephem.F1, 'TRef', ephem.TRef, 'NBin', nBin, 'SubintPeriods', subintPeriods, ...
        'SaveFile', false, 'Verbose', false};
qc = {'Bnoise', nc.Bnoise, 'Verbose', false};
[ifC, fC] = foldProfile(detD, 'x', ephem.f0, base{:}, 'NoiseCoeffs', [nc.V, nc.X]);
% fake blanking of whole detected bins: 5 % random, channels 1-10 at phase 0.40-0.42 in odd turns
Pall = reshape(readF32(detD.file), nCh, Nd);
kAll = 1:Nd;
phAll = (detD.binTime0 + (kAll - 1)*detD.binDt - ephem.TRef) * ephem.f0;
turnK = round(phAll); phk = phAll - turnK;
rng(5);
keep = rand(nCh, Nd) >= 0.05;
keep(1:10, phk >= 0.40 & phk < 0.42 & mod(turnK, 2) == 1) = false;
detB = detD; detB.file = fullfile(tmp, 'p_blank.dat');
writeF32(detB.file, Pall .* single(keep));
infoW = struct('file', fullfile(tmp, 'w_blank.dat'), 'nChan', nCh, 'N', Nd, 'Lmax', nc.Lmax, ...
    'byteOrder', 'ieee-le');
writeWeightBlocks(infoW.file, Nd, @(kk) blankStream(kk, keep, [nc.V, nc.X]));
[ifB, fB] = foldProfile(detB, 'x', ephem.f0, base{:}, 'DataWeights', infoW);
clear Pall keep

T = templateStuff(template);
cases = {'unblanked', fC, ifC; 'blanked', fB, ifB};
for i = 1:2
    [fo, ifo] = cases{i, 2:3};
    to = estimateTOA(fo, ifo, template, qc{:}, 'Weighting', 'optimal');
    eT = 0; eB = 0; eFit = 0; isMax = true; nS = 0;
    for s = find(to.valid & to.coverage == 1).'
        ex = explicitOptimal(fo, s, to, nc.Bnoise, ifo.binDt, T);
        eT = max(eT, abs(ex.tauErr - to.phaseErr(s)) / to.phaseErr(s));
        eB = max(eB, abs(ex.bErr - to.ampErr(s)) / to.ampErr(s));
        eFit = max([eFit, abs(ex.b - to.amp(s)) / to.amp(s), ...
            max(abs(ex.aC - to.chanBaseline(s, ex.use))) / mean(ex.aC)]);
        isMax = isMax && ex.isMax;
        nS = nS + 1;
    end
    pass = eT < 1e-8 && eB < 1e-8 && eFit < 1e-9 && isMax;
    fprintf(['3. algebra, %-9s %d sub-ints, channels used %s: tau a maximum of N/sqrt(Dn) %d; b, a_c ' ...
             'vs WLS %.1e; phaseErr vs sandwich %.1e, ampErr %.1e: %s\n'], cases{i, 1}, nS, ...
        mat2str(unique(to.nChanUsed(to.valid)).'), isMax, eFit, eT, eB, passStr(pass));
    if ~pass, fails{end+1} = ['algebra ' cases{i, 1}]; end %#ok<AGROW>
end

% ---------------------------------------------------------------------------------
% 4. Monte Carlo with the exact noise covariance of a blanked fold
% ---------------------------------------------------------------------------------
tE_B = estimateTOA(fB, ifB, template, qc{:});
sMC = find(mod(fB.subint.turnRef, 2) == 1 & tE_B.valid, 1);   % channels 1-10 have empty bins
chMC = 1:2:nCh; nM = numel(chMC); nLag = size(fB.weightX, 4); N = nBin;
Wm  = reshape(fB.weight(:, sMC, chMC), N, nM);
W2m = reshape(fB.weight2(:, sMC, chMC), N, nM);
WXm = reshape(fB.weightX(:, sMC, chMC, :), N, nM, nLag);
tau0 = 0.0123;
T0 = T.shape(tau0);                                  % band-limited template at tau0
sTB2 = 1 / (nc.Bnoise * ifB.binDt);
nEmpty = nnz(any(Wm == 0, 1));
Rtot = 800; Rb = 100;
fprintf('4. Monte Carlo: %d channels (%d with empty phase bins), %d profiles per SNR level\n', ...
    nM, nEmpty, Rtot);
for snrT = [10 40]
    % amplitude for the target optimal SNR: error of b from a noise-free fold
    b = 0.05;
    for it = 1:2
        [foldMC, infoMC] = mcFold(ones(N, nM) + b * T0, zeros(N, nM, 1), Wm, W2m, WXm, ifB, 1);
        t0 = estimateTOA(foldMC, infoMC, template, qc{:}, 'Weighting', 'optimal');
        b = b * snrT / t0.snr(1);
    end
    m = ones(N, nM) + b * T0;                        % true power model per channel
    L = noiseFactors(m, Wm, W2m, WXm, sTB2);
    rng(100 + snrT);
    est = struct('equal', [], 'optimal', []);
    for bt = 1:Rtot/Rb
        noise = zeros(N, nM, Rb);
        for c = 1:nM
            h = Wm(:, c) > 0;
            noise(h, c, :) = reshape(L{c}.' * randn(nnz(h), Rb), nnz(h), 1, Rb);
        end
        [foldMC, infoMC] = mcFold(m, noise, Wm, W2m, WXm, ifB, Rb);
        for wn = {'equal', 'optimal'}
            t = estimateTOA(foldMC, infoMC, template, qc{:}, 'Weighting', wn{1});
            est.(wn{1}) = [est.(wn{1}); t.phase, t.phaseErr, t.snr, t.redChi2, t.valid];
        end
    end
    sd = struct();
    for wn = {'equal', 'optimal'}
        e = est.(wn{1}); e = e(e(:, 5) == 1, :); nR = size(e, 1);
        d = mod(e(:, 1) - tau0 + 0.5, 1) - 0.5;
        z = d ./ e(:, 2);
        okBias = abs(mean(d)) < 4 * std(d) / sqrt(nR);
        okPull = abs(std(z) - 1) < 4 / sqrt(2*nR);
        sd.(wn{1}) = [std(d), median(e(:, 2))];
        fprintf(['   SNR %2d %-7s %4d fits: bias %+.2e +- %.2e turns (%+.1f sigma); pull std %.3f ' ...
                 '(+- %.3f); scatter %.3e, median error bar %.3e; median SNR %.1f, red. chi2 %.3f: %s\n'], ...
            snrT, wn{1}, nR, mean(d), std(d)/sqrt(nR), mean(d)/(std(d)/sqrt(nR)), std(z), ...
            1/sqrt(2*nR), std(d), median(e(:, 2)), median(e(:, 3)), median(e(:, 4)), ...
            passStr(okBias && okPull));
        if ~(okBias && okPull), fails{end+1} = sprintf('MC SNR %d %s', snrT, wn{1}); end %#ok<AGROW>
    end
    rMeas = sd.optimal(1) / sd.equal(1); rPred = sd.optimal(2) / sd.equal(2);
    sR = rMeas * sqrt(2 / Rtot);                    % generous: ignores the pairing
    okR = abs(rMeas - rPred) < 4 * sR;
    fprintf('   SNR %2d scatter optimal/equal %.3f (predicted from the error bars %.3f, +- %.3f): %s\n', ...
        snrT, rMeas, rPred, sR, passStr(okR));
    if ~okR, fails{end+1} = sprintf('MC ratio SNR %d', snrT); end %#ok<AGROW>
end

% ---------------------------------------------------------------------------------
% 5. Fake-blanked data: all channels used, TOAs consistent with 'equal'
% ---------------------------------------------------------------------------------
tO_B = estimateTOA(fB, ifB, template, qc{:}, 'Weighting', 'optimal');
v = tE_B.valid & tO_B.valid;
dT = (tO_B.toa(v) - tE_B.toa(v)) * 1e6;
okT = abs(mean(dT)) < 4 * std(dT) / sqrt(nnz(v));
okU = all(tO_B.nChanUsed(v) == nCh);
fprintf(['5. fake-blanked data: channels used optimal %s vs equal %s; TOA optimal - equal mean ' ...
         '%+.3f +- %.3f us; errors %.3f vs %.3f us: %s\n'], mat2str(tO_B.nChanUsed(v).'), ...
    mat2str(tE_B.nChanUsed(v).'), mean(dT), std(dT)/sqrt(nnz(v)), median(tO_B.toaErr(v))*1e6, ...
    median(tE_B.toaErr(v))*1e6, passStr(okT && okU));
if ~(okT && okU), fails{end+1} = 'blanked data'; end

% ---------------------------------------------------------------------------------
if isempty(fails)
    fprintf('testOptimalTOA: ALL PASSED\n');
else
    error('testOptimalTOA:failed', 'FAILED: %s', strjoin(fails, ', '));
end
end


% =================================================================================
function T = templateStuff(template)
% Band-limited template (harmonics |k| <= K and DC) as estimateTOA 'optimal' uses it.
t = template(:) / max(template); N = numel(t);
K = floor(N/2) - 1; k = (1:K).';
S = fft(t); Sbl = zeros(N, 1); Sbl([1; k + 1; N - k + 1]) = S([1; k + 1; N - k + 1]);
ksg = [0:ceil(N/2)-1, -floor(N/2):-1].';
T.N = N; T.Sbl = Sbl; T.ksg = ksg;
T.shape = @(tau) real(ifft(Sbl .* exp(-2i*pi*ksg*tau)));
T.dshape = @(tau) real(ifft(Sbl .* (-2i*pi*ksg) .* exp(-2i*pi*ksg*tau)));
end

function ex = explicitOptimal(fold, s, to, B, binDt, T)
% Rebuild the 'optimal' fit of sub-int s from its outputs with explicit matrices.
N = T.N; nChan = size(fold.prof, 3);
W  = reshape(fold.weight(:, s, :), N, []);
W2 = reshape(fold.weight2(:, s, :), N, []);
WX = reshape(fold.weightX(:, s, :, :), N, size(W, 2), []);
if size(W, 2) == 1
    W = repmat(W, 1, nChan); W2 = repmat(W2, 1, nChan); WX = repmat(WX, 1, nChan, 1);
end
use = any(W > 0, 1); nU = nnz(use);
P = reshape(fold.prof(:, s, :), N, nChan); P = P(:, use);
W = W(:, use); W2 = W2(:, use); WX = WX(:, use, :); has = W > 0; P(~has) = 0;
tau = to.phase(s); b = to.amp(s) / nU; aC = to.chanBaseline(s, use);
Tt = T.shape(tau); dT = T.dshape(tau);
m = aC + b * Tt;
sTB = abs(m) / sqrt(B * binDt);
j = (1:N).';
Sig = cell(1, nU); w = zeros(N, nU);
for c = 1:nU
    h = has(:, c);
    vd = zeros(N, 1); vd(h) = sTB(h, c).^2 .* W2(h, c) ./ W(h, c).^2;
    Sg = sparse(j, j, vd, N, N);
    for d = 1:size(WX, 3)
        jd = mod(j - 1 + d, N) + 1;
        ok = h & h(jd);
        cvd = zeros(N, 1);
        cvd(ok) = sTB(ok, c) .* sTB(jd(ok), c) .* WX(ok, c, d) ./ (W(ok, c) .* W(jd(ok), c));
        Sg = Sg + sparse([j; jd], [jd; j], [cvd; cvd], N, N);
    end
    Sig{c} = Sg;
    w(h, c) = 1 ./ vd(h);
end
Wc = sum(w, 1);
% WLS for this tau: b and a_c, and the objective around tau
F = @(t) wlsObjective(T.shape(t), P, w, Wc);
[~, bW, aW] = F(tau);
ex.b = bW * nU; ex.aC = aW; ex.use = use;
ex.isMax = F(tau) >= F(tau + 1e-7) && F(tau) >= F(tau - 1e-7);
% sandwich
G1 = b * (dT - sum(w .* dT, 1) ./ Wc);
G2 = Tt - sum(w .* Tt, 1) ./ Wc;
A = zeros(2); Bm = zeros(2);
for c = 1:nU
    J = [G1(:, c), G2(:, c)];
    A = A + J.' * (w(:, c) .* J);
    X = w(:, c) .* J;
    Bm = Bm + X.' * Sig{c} * X;
end
Cov = (A \ Bm) / A;
ex.tauErr = sqrt(Cov(1, 1));
ex.bErr = sqrt(Cov(2, 2)) * nU;
end

function [F, b, aC] = wlsObjective(Tt, P, w, Wc)
% N/sqrt(Dn) and the WLS b, a_c for template values Tt (s_c = 1).
pbar = sum(w .* P, 1) ./ Wc;
u = sum(w .* (P - pbar), 2);
WT = sum(w .* Tt, 1);
Dn = sum(w, 2).' * Tt.^2 - sum(WT.^2 ./ Wc);
F = (u.' * Tt) / sqrt(Dn);
b = (u.' * Tt) / Dn;
aC = pbar - b * WT ./ Wc;
end

function L = noiseFactors(m, W, W2, WX, sTB2)
% Cholesky factors of every channel's profile covariance (bins with data only).
[N, nM] = size(m); j = (1:N).'; L = cell(1, nM);
for c = 1:nM
    h = W(:, c) > 0;
    s = sqrt(sTB2) * m(:, c);
    vd = zeros(N, 1); vd(h) = s(h).^2 .* W2(h, c) ./ W(h, c).^2;
    Sg = sparse(j, j, vd, N, N);
    for d = 1:size(WX, 3)
        jd = mod(j - 1 + d, N) + 1;
        ok = h & h(jd);
        cvd = zeros(N, 1);
        cvd(ok) = s(ok) .* s(jd(ok)) .* WX(ok, c, d) ./ (W(ok, c) .* W(jd(ok), c));
        Sg = Sg + sparse([j; jd], [jd; j], [cvd; cvd], N, N);
    end
    L{c} = chol(Sg(h, h));                           % Sg(h,h) = L' * L
end
end

function [fold, info] = mcFold(m, noise, W, W2, WX, ifB, R)
% A fold struct with R sub-ints: profile = m + noise, the same weights in each.
[N, nM] = size(m); nLag = size(WX, 3);
prof = repmat(reshape(m, N, 1, nM), 1, R, 1);
if size(noise, 3) == R
    prof = prof + permute(noise, [1 3 2]);
end
Wf = repmat(reshape(W, N, 1, nM), 1, R, 1);
prof(Wf == 0) = NaN;
fold.prof    = prof;
fold.weight  = Wf;
fold.weight2 = repmat(reshape(W2, N, 1, nM), 1, R, 1);
fold.weightX = repmat(reshape(WX, N, 1, nM, nLag), 1, R, 1, 1);
fold.sum     = prof .* Wf;
pt = mean(prof, 2, 'omitnan');
fold.profTotal = reshape(pt, N, nM);
fold.subint = struct('tRef', zeros(R, 1), 'fRef', ones(R, 1), 'turnRef', (1:R).');
info = struct('NBin', N, 'nSub', R, 'binDt', ifB.binDt, 'f0', ifB.f0, 'nChan', nM);
end

function ok = sameFields(new, old)
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
fid = fopen(file, 'w', 'ieee-le');
c = onCleanup(@() fclose(fid));
for k0 = 1:8192:N
    fwrite(fid, blockFun(k0:min(N, k0 + 8191)), 'single');
end
end

function Y = blankStream(kk, keep, coef)
% Exact W, V, X(L) when whole detected bins are zeroed (as tests/testFoldWeights.m).
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
