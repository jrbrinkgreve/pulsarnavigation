function testRFIModel()
%TESTRFIMODEL  Unit tests for the RFI model of the fast simulator (block C4b): rfiProfile and
%simulateFold('RFI') against foldProfile on synthetic RFI data with known statistics.
%{
Run from the PulsarSimMatlab folder: run('tests/testRFIModel.m'), or runAllTests("testRFIModel").
Needs no data files (synthetic power and weights in tempdir, ~0.9 GB, deleted at the end); does
not use pipelineParams.

The RFI model (simulator.html, C4b): rfiProfile measures per channel what RFI leaves in the
detected power after excision - level, variance V and lag covariances X(L) of the time bins,
valid fraction - and simulateFold('RFI', prof) draws folds with those statistics. Here the
detected power is synthetic with exactly known statistics: one complex sample a_k per time bin
(binDt 1 us, Bnoise 1 MHz, so Bnoise*binDt = 1), power p_k = |a_k|^2. Receiver noise n_k:
complex Gaussian filtered by h = [0.8 0.5 0.33]/|h| (autocorrelation R(L), R(0) = 1), so p has
mean 1, V = 1 and X(L) = R(L)^2, L = 1, 2 (the clean model, Lmax 2). Seven channels:
  1, 6, 7  noise only
  2        constant-envelope carrier a = n + A e^{iwk}, I = |A|^2 = 4: level 1 + I,
           V = (1 + 2I)/(1 + I)^2, X(L) = (R(L)^2 + 2I R(L) cos(wL))/(1 + I)^2 - quieter
           than the model (B5e)
  3        noise-like narrowband RFI a = n + r, r complex Gaussian filtered by [1 1 1]/sqrt(3)
           (R_r(L) = (3 - L)/3), power J = 3: level 1 + J, V = 1, X(L) = R_a(L)^2/R_a(0)^2
           with R_a = R + J*R_r - noisier than the model (narrower spectrum)
  4        impulses (bursts of 5-50 bins, ~10 % of the time, 30x the noise amplitude),
           blanked by a mask: W = 0 in the bursts, V and X of the clean model x keep (as
           testFoldWeights 6); valid fraction = the mask's mean
  5        carrier (as 2) plus impulses (as 4)
The pulsar (-6 dB, unit-peak periodic Gaussian, FWHM 0.05 turns, 100 Hz, peak at phase 0.123)
scales the noise amplitude, a = sqrt(1 + rho*q(t))*n + RFI (the pulsar is noise, as in the
generator).

  1. Estimator: rfiProfile on a calibration stretch without pulsar (2.1e6 bins = 210 turns):
     level, V, X(1), X(2) per channel against the closed forms above within 4 sigma (sigma
     from the scatter of 32 blocks), valid fraction = the mask's mean (1e-6).
  2. Folded profiles: a second stretch with the pulsar (4.2e6 bins, another seed) folded by
     foldProfile with its weights file (NBin 256, 2 turns per sub-int: 210 sub-ints), against
     simulateFold with the profile of test 1 (1000 sub-ints). Per channel: off-pulse mean
     level; mean weight; noise relative to the observer's model (mean u^2, u = (prof -
     level)/s_j with the fold's own weight2 - the noise ratio the error bars miss); neighbour-
     bin correlation. Real and simulated within 4 sigma (scatter between sub-ints, plus the
     calibration error of the profile from test 1's blocks).
  3. TOAs: foldToTOA ('optimal') on both folds (SNR ~10 per TOA): at least 95 % valid; median
     SNR, median sigma, mean reduced chi^2, pull mean and std, within 4 sigma.
  4. Negative control: one channel with noise-like RFI correlated far beyond the measured lags
     (AR(1), r_k = 0.98 r_k-1 + ..., power 1). The profile cannot hold it, so the real fold's
     noise must exceed the simulated one by more than 4 sigma: the fold comparison of test 2
     (and of testRFIChain on the voltage chain) does see correlations the profile misses.
  5. Per-channel noise calibration on the observation itself (20 Oct; processIQ proc.noiseCal):
     rfiProfile on the pulsar stretch of test 2 with 'Ephem' and 'OnPulse' (1.5 FWHM left out)
     recovers the closed-form coefficients (4 sigma); foldProfile with 'ChannelNoise' (those
     coefficients) and the weights file: noise against the fold's own model 1 in every channel
     (4 sigma; test 2 with the clean model: carrier 0.49, narrowband 1.14), TOAs red. chi^2 1 and
     pull std 1 (4 sigma); the model coefficients as ChannelNoise give the uncalibrated fold bit
     for bit, and without DataWeights 'ChannelNoise' equals 'NoiseCoeffs' (1e-12).
  6. More lags (rfiProfile 'Lmax', 21 Oct): the negative control's channel measured with 300 lags
     (0.98^300 ~ 0.2 %; the closed form X(L) = (R(L) + 0.98^L)^2 / 4, sum rule V + 2 sum X =
     14.13, against 2.95 measured with lags 1-2 only): measured sum within 10 % of the closed form, and
     the fold comparison of test 4 now agrees (real and simulated noise/model within 4 sigma).
Errors at the end if any check fails. (20 Oct 2026, block C4b.)
%}

root = fullfile(fileparts(mfilename('fullpath')), '..');
addpath(fullfile(root, 'functions'));
oldDir = cd(root); restoreDir = onCleanup(@() cd(oldDir)); %#ok<NASGU>
tmp = tempname(tempdir);                         % own folder per run (21 Oct): two runs at once
mkdir(tmp);                                      %   (e.g. Jasper's and Claude's) cannot collide
cleanTmp = onCleanup(@() rmdir(tmp, 's')); %#ok<NASGU>
fails = {};

f0 = 100; T = 1/f0; dt = 1e-6; Bn = 1/dt; NB = 256; fwhm = 0.05; dphi = 0.123;
TRef = 0.5/f0; rhoDB = -6; per = 2; nPerTurn = round(T/dt);
nPhase = per / (NB * f0 * dt);                       % time bins per phase bin per sub-int
h  = [0.8 0.5 0.33]; h = h / norm(h);
R  = [sum(h.^2), sum(h(1:2) .* h(2:3)), h(1)*h(3)];  % R(0..2) of the receiver noise
sp = struct('h', h, 'I', 4, 'w', 2*pi*0.05, 'J', 3, 'g', ones(1, 3)/sqrt(3), 'beta', 0.98, ...
    'rho', 0, 'f0', f0, 'dt', dt, 'TRef', TRef, 'dphi', dphi, 'fwhm', fwhm);
model = R.^2 / R(1)^2;                               % clean model [V X(1) X(2)]
kind  = ["noise", "carrier", "narrowband", "impulses", "carrier+impulses", "noise", "noise"];
nCh   = numel(kind);

% closed-form truth per channel
levT = ones(1, nCh); coefT = repmat(model, nCh, 1);
for c = find(contains(kind, "carrier"))
    levT(c) = R(1) + sp.I;
    coefT(c, :) = (R.^2 + 2*sp.I*R .* cos(sp.w*(0:2))) / levT(c)^2;
end
c = find(kind == "narrowband");
Ra = R + sp.J * (3 - (0:2))/3;
levT(c) = Ra(1); coefT(c, :) = Ra.^2 / Ra(1)^2;

% ---------------------------------------------------------------------------------
% 1. Estimator: rfiProfile on a calibration stretch without pulsar
% ---------------------------------------------------------------------------------
N1 = 210 * nPerTurn;
[det1, w1, keep1] = makeRun(fullfile(tmp, 'cal'), kind, N1, sp, 1, model, Bn);
prof = rfiProfile(det1, w1, 'RefLevel', 1, 'Verbose', false);
K = 32; e = round(linspace(1, N1 + 1, K + 1));
bLev = zeros(nCh, K); bCo = zeros(nCh, 3, K); bVf = zeros(nCh, K);
for b = 1:K
    d = det1; d.fullySupportedBins = [e(b), e(b+1) - 1];
    pb = rfiProfile(d, w1, 'RefLevel', 1, 'Verbose', false);
    bLev(:, b) = pb.level.'; bCo(:, :, b) = pb.noiseCoeffs; bVf(:, b) = pb.validFraction.';
end
sLev = std(bLev, 0, 2) / sqrt(K);                    % calibration errors (also used in 2)
sCo  = std(bCo, 0, 3) / sqrt(K);
sVf  = std(bVf, 0, 2) / sqrt(K);
sumB = squeeze(bCo(:, 1, :) + 2*sum(bCo(:, 2:3, :), 2));
sumP = prof.noiseCoeffs(:, 1) + 2*sum(prof.noiseCoeffs(:, 2:3), 2);
sSumRel = std(sumB, 0, 2) / sqrt(K) ./ sumP;          % relative error of the sum rule
zLev = (prof.level.' - levT.') ./ sLev;
zCo  = (prof.noiseCoeffs - coefT) ./ sCo;
vfT  = mean(keep1, 2).';
dVf  = max(abs(prof.validFraction - vfT));
fprintf('1. rfiProfile vs closed forms (%d bins, %d blocks for sigma):\n', N1, K);
fprintf('   %-17s | level est / true  | V est / true    | X(1) est / true | X(2) est / true | valid    | max |z|\n', 'channel');
for c = 1:nCh
    fprintf('   %d %-15s | %7.4f / %7.4f | %.4f / %.4f | %.4f / %.4f | %.4f / %.4f | %.5f | %.2f\n', ...
        c, kind(c), prof.level(c), levT(c), prof.noiseCoeffs(c, 1), coefT(c, 1), ...
        prof.noiseCoeffs(c, 2), coefT(c, 2), prof.noiseCoeffs(c, 3), coefT(c, 3), ...
        prof.validFraction(c), max(abs([zLev(c), zCo(c, :)])));
end
pass = all(abs(zLev) < 4) && all(abs(zCo) < 4, 'all') && dVf < 1e-6;
fprintf('   valid fraction vs the mask: max |diff| %.1e; noise / model (sum rule) %s: %s\n', dVf, ...
    mat2str(round(sumP.' / sum([1 2 2] .* model), 3)), passStr(pass));
if ~pass, fails{end+1} = 'estimator'; end

% ---------------------------------------------------------------------------------
% 2. Folded profiles: foldProfile on RFI data vs simulateFold with the profile
% ---------------------------------------------------------------------------------
N2 = 420 * nPerTurn;
sp2 = sp; sp2.rho = 10^(rhoDB/10);
[det2, w2] = makeRun(fullfile(tmp, 'psr'), kind, N2, sp2, 2, model, Bn);
[ifr, fr] = foldProfile(det2, 'x', f0, 'TRef', TRef, 'NBin', NB, 'SubintPeriods', per, ...
    'DataWeights', w2, 'SaveFile', false, 'Verbose', false);
complete = all(fr.subint.nTimeBins == per * nPerTurn);
[ifs, fs, ts] = simulateFold(f0, rhoDB, 2000, 'NChan', nCh, 'RFI', prof, 'FWHM', fwhm, ...
    'NBin', NB, 'SubintPeriods', per, 'PhaseOffset', dphi, 'Seed', 3, 'Verbose', false);
sr = foldStats(fr, dphi, Bn, dt);
ss = foldStats(fs, dphi, Bn, dt);
zL = (sr.lev - ss.lev) ./ sqrt(sr.sLev.^2 + ss.sLev.^2 + sLev.^2);
zR = (sr.ratio - ss.ratio) ./ sqrt(sr.sRatio.^2 + ss.sRatio.^2 + (ss.ratio .* sSumRel).^2);
zC = (sr.rho1 - ss.rho1) ./ sqrt(sr.sRho1.^2 + ss.sRho1.^2);
dW = sr.wt - ss.wt;
sW = sqrt(sr.sWt.^2 + ss.sWt.^2 + (nPhase * sVf).^2);
zW = dW ./ max(sW, 1e-9 * ss.wt);                     % noise-only channels: exact, 0/0
fprintf(['2. folded profiles: foldProfile on the RFI data (%d sub-ints of %d turns, complete %d) ' ...
         'vs simulateFold with the profile (%d sub-ints), off-pulse:\n'], ifr.nSub, per, complete, ifs.nSub);
fprintf('   %-17s | level real / sim  | weight real / sim | noise/model real / sim | nb. corr real / sim | max |z|\n', 'channel');
for c = 1:nCh
    fprintf('   %d %-15s | %7.4f / %7.4f | %6.3f / %6.3f   | %.4f / %.4f          | %.4f / %.4f     | %.2f\n', ...
        c, kind(c), sr.lev(c), ss.lev(c), sr.wt(c), ss.wt(c), sr.ratio(c), ss.ratio(c), ...
        sr.rho1(c), ss.rho1(c), max(abs([zL(c), zR(c), zC(c), zW(c)])));
end
pass = complete && all(abs([zL; zR; zC; zW]) < 4, 'all');
fprintf('   %s\n', passStr(pass));
if ~pass, fails{end+1} = 'folded profiles'; end

% ---------------------------------------------------------------------------------
% 3. TOAs on both folds
% ---------------------------------------------------------------------------------
tr = foldToTOA(fr, ifr, fwhm, Bn, 'Weighting', 'optimal', 'Verbose', false);
tm = foldToTOA(fs, ifs, fwhm, Bn, 'Weighting', 'optimal', 'Verbose', false);
qr = toaStats(tr, fr.subint.tRef(:) + dphi/f0);
qm = toaStats(tm, ts.pulseCenters(:));
fn = {'snr', 'err', 'chi2', 'pullMean', 'pullStd'};
lb = {'median SNR', 'median sigma [us]', 'mean red. chi2', 'pull mean', 'pull std'};
sc = [1 1e6 1 1 1];
fprintf('3. TOAs, foldToTOA optimal: real %d of %d valid, simulated %d of %d\n', qr.n, numel(tr.valid), ...
    qm.n, numel(tm.valid));
pass = qr.n >= 0.95 * numel(tr.valid) && qm.n >= 0.95 * numel(tm.valid);
for i = 1:numel(fn)
    z = (qr.(fn{i}) - qm.(fn{i})) / sqrt(qr.(['s' fn{i}])^2 + qm.(['s' fn{i}])^2);
    fprintf('   %-18s real %8.4f +- %.4f, simulated %8.4f +- %.4f, z %+.2f\n', lb{i}, ...
        sc(i)*qr.(fn{i}), sc(i)*qr.(['s' fn{i}]), sc(i)*qm.(fn{i}), sc(i)*qm.(['s' fn{i}]), z);
    pass = pass && abs(z) < 4;
end
fprintf('   %s\n', passStr(pass));
if ~pass, fails{end+1} = 'TOAs'; end

% ---------------------------------------------------------------------------------
% 4. Negative control: correlations beyond the measured lags
% ---------------------------------------------------------------------------------
[detC, ~] = makeRun(fullfile(tmp, 'ctrlCal'), "longcorr", 100 * nPerTurn, sp, 4, model, Bn);
profC = rfiProfile(detC, [], 'RefLevel', 1, 'Verbose', false);
[detF, ~] = makeRun(fullfile(tmp, 'ctrlFold'), "longcorr", 200 * nPerTurn, sp, 5, model, Bn);
[~, fc] = foldProfile(detF, 'x', f0, 'TRef', TRef, 'NBin', NB, 'SubintPeriods', per, ...
    'NoiseCoeffs', model, 'SaveFile', false, 'Verbose', false);
[~, fcs] = simulateFold(f0, -Inf, 800, 'NChan', 1, 'RFI', profC, 'NBin', NB, ...
    'SubintPeriods', per, 'Seed', 6, 'Verbose', false);
cr = foldStats(fc, NaN, Bn, dt);
cs = foldStats(fcs, NaN, Bn, dt);
z = (cr.ratio - cs.ratio) / sqrt(cr.sRatio^2 + cs.sRatio^2);
pass = z > 4;
fprintf(['4. negative control (AR(1) 0.98 RFI, profile lags 1-2 only): noise/model real %.3f, ' ...
         'simulated %.3f (x%.2f, z %+.1f): %s\n'], cr.ratio, cs.ratio, cr.ratio / cs.ratio, z, ...
    passStr(pass));
if ~pass, fails{end+1} = 'negative control'; end

% ---------------------------------------------------------------------------------
% 5. Per-channel noise calibration on the observation itself (rfiProfile 'OnPulse',
%    foldProfile 'ChannelNoise'; processIQ proc.noiseCal)
% ---------------------------------------------------------------------------------
ephC = struct('f0', f0, 'F1', 0, 'TRef', TRef, 'Phi0', -dphi);   % the observer's phase 0 at the pulse
pOff = rfiProfile(det2, w2, 'RefLevel', 1, 'Ephem', ephC, 'OnPulse', 1.5*fwhm, 'Verbose', false);
sOff = sCo * sqrt(N1 / (N2 * (1 - pOff.onPulseFraction)));   % errors scaled from test 1's blocks
zOff = (pOff.noiseCoeffs - coefT) ./ sOff;
[i5, f5] = foldProfile(det2, 'x', f0, 'TRef', TRef, 'NBin', NB, 'SubintPeriods', per, ...
    'DataWeights', w2, 'ChannelNoise', pOff.noiseCoeffs, 'SaveFile', false, 'Verbose', false);
s5 = foldStats(f5, dphi, Bn, dt);
z5 = (s5.ratio - 1) ./ sqrt(s5.sRatio.^2 + (sSumRel * sqrt(N1 / N2)).^2);
q5 = toaStats(foldToTOA(f5, i5, fwhm, Bn, 'Weighting', 'optimal', 'Verbose', false), f5.subint.tRef(:) + dphi/f0);
[~, fM] = foldProfile(det2, 'x', f0, 'TRef', TRef, 'NBin', NB, 'SubintPeriods', per, ...
    'DataWeights', w2, 'ChannelNoise', repmat(model, nCh, 1), 'SaveFile', false, 'Verbose', false);
okId = isequal(fM.weight2, fr.weight2) && isequal(fM.weightX, fr.weightX) && isequal(fM.prof, fr.prof);
[~, fA] = foldProfile(detF, 'x', f0, 'TRef', TRef, 'NBin', NB, 'SubintPeriods', per, ...
    'NoiseCoeffs', model, 'SaveFile', false, 'Verbose', false);
[~, fB] = foldProfile(detF, 'x', f0, 'TRef', TRef, 'NBin', NB, 'SubintPeriods', per, ...
    'ChannelNoise', model, 'SaveFile', false, 'Verbose', false);
okNoW = max(abs(fA.weight2(:) - fB.weight2(:))) <= 1e-12 * max(fA.weight2(:)) && ...
    max(abs(fA.weightX(:) - fB.weightX(:))) <= 1e-12 * max(fA.weight2(:)) && isequal(fA.prof, fB.prof);
okChi = abs(q5.chi2 - 1) < 4 * q5.schi2 + 0.005 && abs(q5.pullStd - 1) < 4 * q5.spullStd;
pass = all(abs(zOff) < 4, 'all') && all(abs(z5) < 4) && okChi && okId && okNoW;
fprintf(['5. per-channel noise calibration: rfiProfile on the pulsar stretch with the pulse left out ' ...
         '(%.1f %% of the bins): coefficients vs the closed forms max |z| %.2f; fold with ''ChannelNoise'': ' ...
         'noise / own model %s (clean model in test 2: %s); TOAs red. chi2 %.4f +- %.4f (test 3: %.4f), ' ...
         'pull std %.3f; model coefficients as ChannelNoise = no calibration, bit for bit %d; without ' ...
         'DataWeights = NoiseCoeffs %d: %s\n'], 100 * pOff.onPulseFraction, max(abs(zOff), [], 'all'), ...
    mat2str(round(s5.ratio.', 3)), mat2str(round(sr.ratio.', 3)), q5.chi2, q5.schi2, qr.chi2, q5.pullStd, ...
    okId, okNoW, passStr(pass));
if ~pass, fails{end+1} = 'noise calibration'; end

% ---------------------------------------------------------------------------------
% 6. More lags (rfiProfile 'Lmax'): the negative control's channel with 300 lags
% ---------------------------------------------------------------------------------
L6 = 300;
prof6 = rfiProfile(detC, [], 'RefLevel', 1, 'Lmax', L6, 'Verbose', false);
Ra6 = [R, zeros(1, L6 - 2)] + sp.beta.^(0:L6);       % R_a(L) = R_n(L) + beta^L, power 2
sumT6 = 1 + 2*sum(Ra6(2:end).^2) / Ra6(1)^2;         % closed form of V + 2 sum X
sum6  = prof6.noiseCoeffs(1) + 2*sum(prof6.noiseCoeffs(2:end));
[~, f6s] = simulateFold(f0, -Inf, 800, 'NChan', 1, 'RFI', prof6, 'NBin', NB, ...
    'SubintPeriods', per, 'Seed', 6, 'Verbose', false);
c6 = foldStats(f6s, NaN, Bn, dt);
z6 = (cr.ratio - c6.ratio) / sqrt(cr.sRatio^2 + c6.sRatio^2);
pass = abs(sum6 / sumT6 - 1) < 0.10 && abs(z6) < 4;
fprintf(['6. more lags (''Lmax'' %d) on the negative control: sum rule %.2f (closed form %.2f, 2 lags %.2f); ' ...
         'noise/model real %.3f, simulated %.3f (z %+.1f; test 4 with 2 lags: %.3f): %s\n'], L6, sum6, sumT6, ...
    profC.noiseCoeffs(1) + 2*sum(profC.noiseCoeffs(2:end)), cr.ratio, c6.ratio, z6, cs.ratio, passStr(pass));
if ~pass, fails{end+1} = 'more lags'; end

if isempty(fails)
    fprintf('testRFIModel: all passed\n');
else
    error('testRFIModel:fail', 'testRFIModel: failed: %s', strjoin(fails, ', '));
end
end


% =================================================================================
function [det, infoW, keep] = makeRun(base, kind, N, sp, seed, model, Bn)
%MAKERUN  Write a synthetic power file (and a weights file if any channel is blanked).
rs = RandStream('mt19937ar', 'Seed', seed);
nCh = numel(kind);
P = zeros(nCh, N, 'single');
keep = true(nCh, N);
q = zeros(N, 1);
if sp.rho > 0                                       % unit-peak periodic Gaussian at phase dphi
    t  = ((1:N).' - 0.5) * sp.dt;
    ph = mod((t - sp.TRef) * sp.f0 - sp.dphi + 0.5, 1) - 0.5;
    sg = sp.fwhm / (2*sqrt(2*log(2)));
    for m = -1:1
        q = q + exp(-0.5 * ((ph + m) / sg).^2);
    end
    q = q / sum(exp(-0.5 * ((-1:1) / sg).^2));
end
cn = @(n) (randn(rs, n, 1) + 1i*randn(rs, n, 1)) / sqrt(2);
k = (1:N).';
for c = 1:nCh
    x = filter(sp.h, 1, cn(N + 2));
    a = sqrt(1 + sp.rho * q) .* x(3:end);
    if contains(kind(c), "carrier")
        a = a + sqrt(sp.I) * exp(1i * (sp.w * k + 2*pi*rand(rs)));
    end
    if kind(c) == "narrowband"
        r = filter(sp.g, 1, cn(N + 2));
        a = a + sqrt(sp.J) * r(3:end);
    end
    if kind(c) == "longcorr"
        r = filter(sqrt(1 - sp.beta^2), [1, -sp.beta], cn(N + 500));
        a = a + r(501:end);
    end
    if contains(kind(c), "impulses")
        kp = true(N, 1);
        nB = round(0.1 * N / 27.5);
        st = randi(rs, N, nB, 1); ln = randi(rs, [5 50], nB, 1);
        for b = 1:nB
            kp(st(b) : min(N, st(b) + ln(b) - 1)) = false;
        end
        a(~kp) = a(~kp) + 30 * cn(nnz(~kp));
        keep(c, :) = kp.';
    end
    P(c, :) = single(abs(a).^2).';
end
det = struct('file', [base '_power.dat'], 'nChan', nCh, 'N', N, 'binTime0', 0.5*sp.dt, ...
    'binDt', sp.dt, 'fullySupportedBins', [1 N], 'byteOrder', 'ieee-le', ...
    'chanFreqs', (1:nCh) * 1e6, 'chanWidth', 1e6, ...
    'noise', struct('V', model(1), 'X', model(2:end), 'Bnoise', Bn));
fid = fopen(det.file, 'w', 'ieee-le');
fwrite(fid, P, 'single');
fclose(fid);
infoW = [];
if ~all(keep, 'all')
    infoW = struct('file', [base '_weights.dat'], 'nChan', nCh, 'N', N, 'Lmax', numel(model) - 1, ...
        'byteOrder', 'ieee-le');
    writeWeights(infoW.file, keep, model);
end
end

function writeWeights(file, keep, coef)
%WRITEWEIGHTS  Data-weight file [nChan x (2+Lmax) x N] of whole blanked time bins:
% W = keep, V = V0*keep_k, X(L) = X0(L)*keep_k*keep_k+L (as blankingWeights' format).
[nCh, N] = size(keep);
Lmax = numel(coef) - 1;
fid = fopen(file, 'w', 'ieee-le');
cl = onCleanup(@() fclose(fid));
blk = 2^18;
for k0 = 1:blk:N
    kk = k0 : min(N, k0 + blk - 1);
    kp = keep(:, kk);
    Y = zeros(nCh, 2 + Lmax, numel(kk), 'single');
    Y(:, 1, :) = reshape(kp, nCh, 1, []);
    Y(:, 2, :) = reshape(coef(1) * kp, nCh, 1, []);
    for L = 1:Lmax
        nb = keep(:, min(kk + L, N)) & (kk + L <= N);
        Y(:, 2 + L, :) = reshape(coef(1 + L) * (kp & nb), nCh, 1, []);
    end
    fwrite(fid, Y, 'single');
end
end

function st = foldStats(fold, dphi, Bn, dt)
%FOLDSTATS  Per channel, off-pulse (|phase - dphi| > 0.15; dphi NaN: all bins): level, mean
% weight, noise relative to the observer's model (mean u^2, u = (prof - level)/s_j, s_j =
% level*sqrt(weight2/(Bn*dt))/weight) and the neighbour-bin correlation of u; each with its
% error from the scatter between sub-ints.
[NBin, nSub, nCh] = size(fold.prof);
ph = fold.phase(:);
off = true(NBin, 1);
if ~isnan(dphi), off = abs(mod(ph - dphi + 0.5, 1) - 0.5) > 0.15; end
nW = size(fold.weight, 3);
z0 = zeros(nCh, 1);
st = struct('lev', z0, 'sLev', z0, 'wt', z0, 'sWt', z0, 'ratio', z0, 'sRatio', z0, ...
    'rho1', z0, 'sRho1', z0);
for c = 1:nCh
    W  = fold.weight(:, :, min(c, nW));
    W2 = fold.weight2(:, :, min(c, nW));
    P  = fold.prof(:, :, c);
    ok = off & W > 0;
    Pz = P; Pz(~ok) = 0;
    levS = sum(Pz, 1) ./ sum(ok, 1);
    L  = mean(levS);
    u  = (P - L) ./ (L * sqrt(W2 / (Bn*dt)) ./ W);
    u(~ok) = NaN;
    rS = mean(u.^2, 1, 'omitnan');
    cS = mean(u .* circshift(u, -1, 1), 1, 'omitnan') ./ rS;
    wS = mean(W, 1);                                % all bins: sum = valid time bins / NBin
    st.lev(c)  = L;          st.sLev(c)   = std(levS) / sqrt(nSub);
    st.wt(c)   = mean(wS);   st.sWt(c)    = std(wS) / sqrt(nSub);
    st.ratio(c) = mean(rS);  st.sRatio(c) = std(rS) / sqrt(nSub);
    st.rho1(c) = mean(cS);   st.sRho1(c)  = std(cS) / sqrt(nSub);
end
end

function q = toaStats(toa, truth)
%TOASTATS  Medians / means of the valid TOAs with their errors.
v = toa.valid(:);
z = (toa.toa(v) - truth(v)) ./ toa.toaErr(v);
n = nnz(v);
q = struct('n', n, ...
    'snr', median(toa.snr(v)), 'ssnr', 1.2533 * std(toa.snr(v)) / sqrt(n), ...
    'err', median(toa.toaErr(v)), 'serr', 1.2533 * std(toa.toaErr(v)) / sqrt(n), ...
    'chi2', mean(toa.redChi2(v)), 'schi2', std(toa.redChi2(v)) / sqrt(n), ...
    'pullMean', mean(z), 'spullMean', std(z) / sqrt(n), ...
    'pullStd', std(z), 'spullStd', std(z) / sqrt(2*n));
end

function s = passStr(p)
if p, s = 'PASS'; else, s = 'FAIL'; end
end
