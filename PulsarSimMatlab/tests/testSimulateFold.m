function testSimulateFold()
%TESTSIMULATEFOLD  Unit tests for simulateFold (block C1: fast simulator, full-band equivalent).
%{
Run from the PulsarSimMatlab folder: run('tests/testSimulateFold.m'). Needs no data files
(synthetic power in tempdir, deleted at the end); does not use pipelineParams.

  1. Conventions against foldProfile: a noise-free synthetic power series (the generator's
     unit-peak periodic Gaussian, FWHM 0.05 turns, averaged over 1 us time bins with erf,
     pulses at phase 0.0123, P = 10 ms) folded by foldProfile (linear, NBin 2048, 1 turn per
     sub-int) versus simulateFold with 'Bnoise' Inf (the mean profile alone): every complete
     sub-int within 5e-4 of the pulse peak (a one-bin shift would give ~1e-2). Also the
     field names of fold, fold.subint and info (simulateFold has all of foldProfile's).
  2. Noise statistics, asymptotic weights: 4000 sub-ints without pulsar; u = (prof - 1)/s_j
     must have variance 1 and neighbour correlation 0.25 (lag 2: 0); with a pulse (0 dB) the
     residuals (prof - m)/s_j in the brightest bins also variance 1 (s_j ~ m_j: self-noise).
     Pass within 4 sigma.
  3. Copied weights: the foldProfile weights of test 1 (non-uniform), 400 draws: pooled
     variance 1 and per-bin correlation equal to WX/sqrt(W2*W2') (mean difference within
     4 sigma).
  4. Estimator: foldToTOA on simulated folds (-5 dB, 400 one-turn sub-ints, random true
     offsets): pulls (TOA - truth)/sigma mean 0 and std 1 within 4 sigma; validateTOA reads
     the truth struct (all TOAs matched). H0 (no pulsar, 3000 sub-ints): detectPulsar's T0 mean 0, std 1.
  5. Speed: 20 h of a 100 Hz pulsar in 60 s sub-ints (1200 sub-ints of 6000 turns); pass
     below 30 s.
  Channels (C4a, 20 Oct), with channel-path time bins (0.96 us, Bnoise 3.0596 MHz) and
  illustrative noise coefficients [V X(1..7)] (sum rule V + 2*sum X ~ 1):
  6. Asymptotic weights with correlated time bins against foldProfile with the same
     NoiseCoeffs on a long synthetic series (30 turns per sub-int, so the grid pattern
     averages out): mean weight2/weight and weightX(d)/weight within 2e-3 relative, same D.
  7. Channel noise, 16 channels, 300 sub-ints without pulsar: variance 1, lag-d correlation
     = weightX(d)/weight2 (d = 1..D), neighbouring channels uncorrelated; within 4 sigma.
  8. 128 channels at -5 dB (3 x 60 one-turn sub-ints, ChannelGain 1), foldToTOA 'optimal':
     pulls mean 0, std 1 within 4 sigma; median sigma within 2 % of the one-channel
     full-band fold (same B_n in total and the sum rule: the channel path loses nothing
     without the forward taper).
  RFI profile (C4b, 20 Oct), a synthetic profile of 8 channels (no data files): channel 2
  level 4, channel 3 noise coefficients x1.3 (a narrowed spectrum), channel 4 valid
  fraction 0.9, channel 5 a strong carrier (level 5, coefficients x0.36):
  9. 600 sub-ints without pulsar: per-channel mean = level (4 sigma); the fold's weights =
     the clean model x valid fraction (exact); measured variance / the observer's model
     variance = the profile's noise ratio per channel (4 sigma).
  Schedule (C3a, 20 Oct):
  10. One pass covering turns 0..nTurns-1 with 'SubintTime' = SubintPeriods/f0 gives the
      same fold, bit for bit, as the nTurns interface; two passes [0 30; 100 130.004] s
      with 10 s sub-ints at 100 Hz: 3 + 3 sub-ints, every sub-int's turns inside its pass,
      the expected turnRef / tRef, truth.pass; foldToTOA on that fold at -5 dB gives 6 valid
      TOAs and validateTOA reads the truth (absolute times; true offset in pulseCenters).
  Scintillation (C3b, 20 Oct):
  11. FluxScale 0.25 in half the sub-ints (-5 dB, 400 one-turn sub-ints): fitted amplitude
      ratio 0.25 within 4 sigma. Scintillation field (64 channels of 6.25 MHz, scintle
      bandwidth 50 MHz, scintle time 60 s, 20000 sub-ints of 6 s): F mean 1, variance 1,
      P(F < 0.1) = 1 - exp(-0.1) (exponential statistics), normalized intensity ACF 1/e at
      the scintle time and 1/2 at the scintle bandwidth (within 0.06); one channel over
      400 MHz (the band average) keeps mean 1 with variance below 0.3. 'Scintillation', [] (the
      switch off, pipelineParams scintOn = false) = no option, bit for bit.
  Adaptive sub-ints (C3c, 20 Oct):
  12. mergeSubints: the one-turn foldProfile fold of test 1 merged in groups of 5 = the
      5-turn foldProfile fold (sum, weight, weight2, weightX, prof, turnRef, tRef; 1e-12).
      adaptiveSubints on 4000 one-turn sub-ints at -25 dB (SNR ~1.2 per turn): groups close
      at T0 >= 7; TOA pulls of the complete groups mean 0, std 1 (4 sigma) - the stopping
      rule does not bias the arrival times; with scintillation (2 s scintles) the same, and
      bright stretches need fewer sub-ints than faint ones; two passes with a gap: no group
      spans the gap.
  Dual polarization and groups over gaps (20 Oct):
  13. 'NPol', 1 = no option, bit for bit (with and without the RFI profile of test 9);
      'NPol', 2 = 'Bnoise' x 2 (1e-12); with the RFI profile and NPol 2, the variance against
      the observer's model with Bnoise x 2 = the profile's ratios (as test 9, 4 sigma).
      adaptiveSubints 'AcrossGaps': false = the default of test 12, unchanged; true on the
      two-pass fold of test 12: some groups span the gap, at most one incomplete group (the
      last), their TOAs within 5 sigma of the truth; 'TargetSNR' Inf with 'MaxSubints' 10:
      fixed groups of 10.
Errors at the end if any check fails.
%}

root = fullfile(fileparts(mfilename('fullpath')), '..');
addpath(fullfile(root, 'functions'));
oldDir = cd(root); restoreDir = onCleanup(@() cd(oldDir)); %#ok<NASGU>
tmp = fullfile(tempdir, 'testSimulateFold');
if ~isfolder(tmp), mkdir(tmp); end
cleanTmp = onCleanup(@() delete(fullfile(tmp, '*'))); %#ok<NASGU>
fails = {};

f0 = 100; T = 1/f0; dt = 1e-6; N = 2048; fwhm = 0.05; dphi0 = 0.0123;
TRef = 0.5 * T;

% ---------------------------------------------------------------------------------
% 1. Conventions against foldProfile (noise-free)
% ---------------------------------------------------------------------------------
nT = 20;                                            % turns of synthetic power
nb = round(nT * T / dt);
t  = (0.5 + (0:nb-1)) * dt;                         % bin centroids
sTime = fwhm / (2*sqrt(2*log(2))) * T;              % [s]
jj = -ceil(10*sTime/T):ceil(10*sTime/T);
peakNorm = sum(exp(-0.5 * (jj*T / sTime).^2));
P = zeros(size(t));
for j = -1:nT                                       % pulses near the data, box-averaged
    c = TRef + j*T + dphi0*T;
    P = P + sTime*sqrt(pi/2)/dt * (erf((t + dt/2 - c)/(sqrt(2)*sTime)) - ...
                                   erf((t - dt/2 - c)/(sqrt(2)*sTime)));
end
P = 1 + P / peakNorm;                               % baseline 1, pulse peak ~1 above it
pf = fullfile(tmp, 'power.dat');
fid = fopen(pf, 'w', 'ieee-le'); fwrite(fid, single(P), 'single'); fclose(fid);
info_det = struct('file', pf, 'nChan', 1, 'N', nb, 'binTime0', t(1), 'binDt', dt, ...
    'fullySupportedBins', [1 nb], 'byteOrder', 'ieee-le', 'chanFreqs', []);
[ifp, fp] = foldProfile(info_det, fullfile(tmp, 'fold.mat'), f0, 'TRef', TRef, 'NBin', N, ...
    'SubintPeriods', 1, 'SaveFile', false, 'Verbose', false);
[ifs, fs, ts] = simulateFold(f0, 0, 1, 'PhaseOffset', dphi0, 'Bnoise', Inf, 'Verbose', false);
full = find(fp.subint.nTimeBins == round(T/dt));    % complete sub-ints of the real fold
dMax = max(abs(fp.prof(:, full) - fs.prof), [], 'all');
[~, iPk] = max(fs.prof);
namesOk = all(ismember(fieldnames(fp), fieldnames(fs))) && ...
    all(ismember(fieldnames(fp.subint), fieldnames(fs.subint))) && ...
    all(ismember(fieldnames(ifp), fieldnames(ifs)));
pass = dMax < 5e-4 * ts.peak && namesOk;
fprintf(['1. conventions vs foldProfile (noise-free, offset %.4f): %d complete sub-ints, max |diff| ' ...
         '%.2e of the peak (limit 5e-4), peak at bin %d (phase %.5f); field names %d: %s\n'], ...
    dphi0, numel(full), dMax / ts.peak, iPk, (iPk-1)/N, namesOk, passStr(pass));
if ~pass, fails{end+1} = 'conventions'; end

% ---------------------------------------------------------------------------------
% 2. Noise statistics, asymptotic weights
% ---------------------------------------------------------------------------------
Bn = 391.6e6; nS = 4000;
[~, f2] = simulateFold(f0, -Inf, nS, 'Seed', 2, 'Verbose', false);
s2 = sqrt(f2.weight2 / (Bn*dt)) ./ f2.weight;       % s_j with m = 1
u  = (f2.prof - 1) ./ s2;
[vU, r1, r2] = deal(mean(u(:).^2), mean(u .* circshift(u, -1), 'all'), mean(u .* circshift(u, -2), 'all'));
nU = numel(u);
eV = 4*sqrt(2*1.125/nU); eR = 4*sqrt(1.5/nU);       % 4 sigma (correlated pairs, generous)
okH0 = abs(vU - 1) < eV && abs(r1 - 0.25) < eR && abs(r2) < eR;
[~, f3, t3] = simulateFold(f0, 0, nS, 'Seed', 3, 'Verbose', false);
[~, f3m] = simulateFold(f0, 0, 1, 'Bnoise', Inf, 'Verbose', false);   % mean profile
m3 = f3m.prof;
s3 = m3 .* sqrt(f3.weight2 / (Bn*dt)) ./ f3.weight;
u3 = (f3.prof - m3) ./ s3;
top = m3(:, 1) - 1 > 0.5 * t3.peak;                 % brightest bins (m up to 2)
vTop = mean(u3(top, :).^2, 'all'); nTop = nnz(top) * nS;
okPulse = abs(vTop - 1) < 4*sqrt(2*1.125/nTop);
pass = okH0 && okPulse;
fprintf(['2. noise, asymptotic weights (%d sub-ints): H0 var %.4f, lag-1 %.4f (0.25), lag-2 %+.4f ' ...
         '(limits +-%.4f / %.4f); 0 dB, bins above half peak (m %.2f-%.2f): var %.4f: %s\n'], ...
    nS, vU, r1, r2, eV, eR, min(m3(top)), max(m3(top)), vTop, passStr(pass));
if ~pass, fails{end+1} = 'noise asymptotic'; end

% ---------------------------------------------------------------------------------
% 3. Copied weights (non-uniform, from the foldProfile fold of test 1)
% ---------------------------------------------------------------------------------
Wc = struct('weight', fp.weight(:, full), 'weight2', fp.weight2(:, full), 'weightX', fp.weightX(:, full));
nC = numel(full); nD = 400;
U = zeros(N, nC, nD);
for d = 1:nD
    [~, fc] = simulateFold(f0, -Inf, nC, 'Weights', Wc, 'Seed', 100 + d, 'Verbose', false);
    U(:, :, d) = (fc.prof - 1) ./ (sqrt(fc.weight2 / (Bn*dt)) ./ fc.weight);
end
rModel = Wc.weightX ./ sqrt(Wc.weight2 .* circshift(Wc.weight2, -1));
vC = mean(U.^2, 'all');
dR = mean(mean(U .* circshift(U, -1), 3) - rModel, 'all');
nCU = numel(U);
pass = abs(vC - 1) < 4*sqrt(2*1.125/nCU) && abs(dR) < 4*sqrt(1.5/nCU);
fprintf(['3. copied weights (W %.2f-%.2f, model correlation %.3f-%.3f), %d draws: var %.4f, ' ...
         'mean (measured - model) correlation %+.5f: %s\n'], min(Wc.weight(:)), max(Wc.weight(:)), ...
    min(rModel(:)), max(rModel(:)), nD, vC, dR, passStr(pass));
if ~pass, fails{end+1} = 'noise copied weights'; end

% ---------------------------------------------------------------------------------
% 4. Estimator on simulated folds
% ---------------------------------------------------------------------------------
nE = 400;
rs = RandStream('mt19937ar', 'Seed', 4);
offs = 0.004 * (rand(rs, nE, 1) - 0.5);
[i4, f4, t4] = simulateFold(f0, -5, nE, 'PhaseOffset', offs, 'Seed', 5, 'Verbose', false);
toa = foldToTOA(f4, i4, fwhm, Bn, 'Verbose', false);
v = toa.valid;
z = (toa.toa(v) - t4.pulseCenters(v)) ./ toa.toaErr(v);
val = validateTOA(toa, t4, 'Plot', false);
okT = nnz(v) == nE && abs(mean(z)) < 4/sqrt(nnz(v)) && abs(std(z) - 1) < 4/sqrt(2*nnz(v)) && ...
    val.n == nnz(v);
[i5, f5] = simulateFold(f0, -Inf, 3000, 'Seed', 6, 'Verbose', false);
d5 = detectPulsar(f5, i5, gaussianTemplate(N, fwhm), 'Bnoise', Bn, 'Verbose', false);
T0 = d5.T0(d5.tested); n5 = numel(T0);
okH = abs(mean(T0)) < 4/sqrt(n5) && abs(std(T0) - 1) < 4/sqrt(2*n5);
pass = okT && okH;
fprintf(['4. estimator: -5 dB, %d TOAs, median SNR %.1f, median sigma %.3f us, pulls mean %+.3f ' ...
         'std %.3f (+- %.3f), validateTOA ratio %.3f; H0 %d sub-ints: T0 mean %+.3f std %.3f ' ...
         '(+- %.3f), false alarms %d known / %d unknown (P_FA 1e-3: %.1f): %s\n'], nnz(v), ...
    median(toa.snr(v)), median(toa.toaErr(v))*1e6, mean(z), std(z), 1/sqrt(2*nnz(v)), val.ratio, ...
    n5, mean(T0), std(T0), 1/sqrt(2*n5), nnz(d5.detectedKnown), nnz(d5.detectedUnknown), ...
    1e-3*n5, passStr(pass));
if ~pass, fails{end+1} = 'estimator'; end

% ---------------------------------------------------------------------------------
% 5. Speed: 20 h in 60 s sub-ints
% ---------------------------------------------------------------------------------
tic;
[i6, ~] = simulateFold(f0, -50, 20*3600*f0, 'SubintPeriods', 60*f0, 'Seed', 7, 'Verbose', false);
el = toc;
pass = el < 30;
fprintf('5. speed: %d sub-ints of %d turns (20 h) in %.2f s (limit 30 s): %s\n', ...
    i6.nSub, i6.subintPeriods, el, passStr(pass));
if ~pass, fails{end+1} = 'speed'; end

% ---------------------------------------------------------------------------------
% 6. Asymptotic weights with correlated time bins vs foldProfile (NoiseCoeffs)
% ---------------------------------------------------------------------------------
dtC = 0.96e-6; BnC = 3.0596e6;
ncC = [0.8877 0.0458 0.0051 0.0020 0.0012 0.0009 0.0007 0.0004];
nT6 = 90; nb6 = round(nT6 * T / dtC);
t6 = (0.5 + (0:nb6-1)) * dtC;
pf6 = fullfile(tmp, 'power6.dat');
fid = fopen(pf6, 'w', 'ieee-le'); fwrite(fid, ones(1, nb6, 'single'), 'single'); fclose(fid);
det6 = struct('file', pf6, 'nChan', 1, 'N', nb6, 'binTime0', t6(1), 'binDt', dtC, ...
    'fullySupportedBins', [1 nb6], 'byteOrder', 'ieee-le', 'chanFreqs', []);
[~, fp6] = foldProfile(det6, fullfile(tmp, 'fold6.mat'), f0, 'TRef', TRef, 'NBin', N, ...
    'SubintPeriods', 30, 'NoiseCoeffs', ncC, 'SaveFile', false, 'Verbose', false);
c6 = find(fp6.subint.nTimeBins == round(30 * T / dtC) | ...
          fp6.subint.nTimeBins == round(30 * T / dtC) + 1 | fp6.subint.nTimeBins == round(30 * T / dtC) - 1);
[i6s, fs6] = simulateFold(f0, -Inf, 30, 'SubintPeriods', 30, 'BinDt', dtC, 'Bnoise', BnC, ...
    'NoiseCoeffs', ncC, 'Verbose', false);
D6 = size(fp6.weightX, 4);
refW2 = mean(fp6.weight2(:, c6) ./ fp6.weight(:, c6), 'all');
simW2 = fs6.weight2(1) / fs6.weight(1);
refX = zeros(1, D6); simX = zeros(1, D6);
for d = 1:D6
    refX(d) = mean(fp6.weightX(:, c6, 1, d) ./ fp6.weight(:, c6), 'all');
    simX(d) = fs6.weightX(1, 1, 1, d) / fs6.weight(1);
end
relErr = max(abs([simW2 - refW2, simX - refX]) ./ abs([refW2, refX(1:min(2, D6)), refW2*ones(1, D6 - min(2, D6))]));
pass = i6s.covLags == D6 && relErr < 2e-3;
fprintf(['6. correlated time bins (delta %.3f phase bins, D %d): weight2/weight %.5f (foldProfile %.5f), ' ...
         'weightX/weight %s (%s), max rel. error %.1e (limit 2e-3; lags > 2 relative to weight2): %s\n'], ...
    f0*dtC*N, D6, simW2, refW2, mat2str(round(simX, 5)), mat2str(round(refX, 5)), relErr, passStr(pass));
if ~pass, fails{end+1} = 'correlated weights'; end

% ---------------------------------------------------------------------------------
% 7. Channel noise statistics (16 channels, no pulsar)
% ---------------------------------------------------------------------------------
nS7 = 300; nC7 = 16;
[~, f7] = simulateFold(f0, -Inf, nS7, 'NChan', nC7, 'BinDt', dtC, 'Bnoise', BnC, ...
    'NoiseCoeffs', ncC, 'Seed', 8, 'Verbose', false);
s7 = sqrt(f7.weight2 / (BnC*dtC)) ./ f7.weight;     % shared by the channels
u7 = (f7.prof - 1) ./ s7;
n7 = numel(u7);
v7 = mean(u7.^2, 'all');
D7 = size(f7.weightX, 4);
rM = squeeze(f7.weightX(1, 1, 1, :)).' / f7.weight2(1);
rD = arrayfun(@(d) mean(u7 .* circshift(u7, -d, 1), 'all'), 1:D7);
cc = mean(u7(:, :, 1:end-1) .* u7(:, :, 2:end), 'all');   % neighbouring channels
ok7 = abs(v7 - 1) < 4*sqrt(2*1.2/n7) && all(abs(rD - rM) < 4*sqrt(1.5/n7)) && abs(cc) < 4*sqrt(1.5/n7);
fprintf(['7. channel noise (%d channels x %d sub-ints): var %.4f, lag correlation %s (model %s), ' ...
         'neighbouring channels %+.4f: %s\n'], nC7, nS7, v7, mat2str(round(rD, 4)), mat2str(round(rM, 4)), ...
    cc, passStr(ok7));
if ~ok7, fails{end+1} = 'channel noise'; end

% ---------------------------------------------------------------------------------
% 8. 128 channels at -5 dB, optimal weighting, vs the one-channel full-band fold
% ---------------------------------------------------------------------------------
z8 = []; sg8 = []; sn8 = [];
for rep8 = 1:3
    offs8 = 0.004 * (rand(rs, 60, 1) - 0.5);
    [i8, f8, t8] = simulateFold(f0, -5, 60, 'NChan', 128, 'BinDt', dtC, 'Bnoise', BnC, ...
        'NoiseCoeffs', ncC, 'PhaseOffset', offs8, 'Seed', 20 + rep8, 'Verbose', false);
    toa8 = foldToTOA(f8, i8, fwhm, BnC, 'Weighting', 'optimal', 'Verbose', false);
    v8 = toa8.valid;
    z8  = [z8; (toa8.toa(v8) - t8.pulseCenters(v8)) ./ toa8.toaErr(v8)]; %#ok<AGROW>
    sg8 = [sg8; toa8.toaErr(v8)]; sn8 = [sn8; toa8.snr(v8)];              %#ok<AGROW>
    clear f8
end
sgFB = median(toa.toaErr(toa.valid));               % test 4: one channel, full band, -5 dB
ratio8 = median(sg8) / sgFB;
ok8 = numel(z8) == 180 && abs(mean(z8)) < 4/sqrt(180) && abs(std(z8) - 1) < 4/sqrt(360) && ...
    abs(ratio8 - 1) < 0.02;
fprintf(['8. 128 channels, -5 dB, optimal: %d TOAs, median SNR %.1f, median sigma %.3f us (full band ' ...
         '%.3f us, ratio %.4f, limit 2 %%), pulls mean %+.3f std %.3f (+- %.3f): %s\n'], numel(z8), ...
    median(sn8), median(sg8)*1e6, sgFB*1e6, ratio8, mean(z8), std(z8), 1/sqrt(360), passStr(ok8));
if ~ok8, fails{end+1} = 'channels optimal'; end

% ---------------------------------------------------------------------------------
% 9. RFI profile (synthetic): level, weights, true vs model noise per channel
% ---------------------------------------------------------------------------------
nC9 = 8;
P9 = struct('label', 'synthetic', 'nChan', nC9, 'level', ones(1, nC9), ...
    'noiseCoeffs', repmat(ncC, nC9, 1), 'modelCoeffs', ncC, 'validFraction', ones(1, nC9), ...
    'binDt', dtC, 'BnoiseChan', BnC);
P9.level(2) = 4;
P9.noiseCoeffs(3, :) = 1.3 * ncC;
P9.validFraction(4) = 0.9;
P9.level(5) = 5; P9.noiseCoeffs(5, :) = 0.36 * ncC;
nS9 = 600;
[i9, f9, t9] = simulateFold(f0, -Inf, nS9, 'NChan', nC9, 'RFI', P9, 'Seed', 9, 'Verbose', false);
[~, f9c] = simulateFold(f0, -Inf, 1, 'NChan', nC9, 'BinDt', dtC, 'Bnoise', BnC, 'NoiseCoeffs', ncC, ...
    'Verbose', false);                              % clean weights
okW = isequal(size(f9.weight), [N nS9 nC9]) && i9.nWeightChan == nC9 && ...
    max(abs(squeeze(f9.weight(1, 1, :)).' - f9c.weight(1) * P9.validFraction)) < 1e-12 && ...
    max(abs(squeeze(f9.weight2(1, 1, :)).' - f9c.weight2(1) * P9.validFraction)) < 1e-12;
mu9 = squeeze(mean(f9.prof, [1 2])).';
nPer9 = N * nS9;
rLR = zeros(1, nC9); okM = true; okV = true;
for c = 1:nC9
    sM = P9.level(c) * sqrt(f9.weight2(:, :, c) / (BnC*dtC)) ./ f9.weight(:, :, c);   % observer's model
    u = (f9.prof(:, :, c) - P9.level(c)) ./ sM;
    rLR(c) = mean(u(:).^2);                         % per-bin variance / model
    okM = okM && abs(mu9(c) - P9.level(c)) < 4 * P9.level(c) * sqrt(f9c.weight2(1) / (BnC*dtC)) / ...
        f9c.weight(1) / sqrt(nPer9 / 2);
end
% expected per-bin variance ratio: true weight2 / model weight2 (same kernel, scaled coefficients)
rExp9 = [1, 1, 1.3, 1, 0.36, 1, 1, 1];
okV = all(abs(rLR - rExp9) < 4 * rExp9 * sqrt(2*1.2 / nPer9));
pass = okW && okM && okV;
fprintf(['9. RFI profile (8 synthetic channels, %d sub-ints): weights = model x valid fraction %d; ' ...
         'means %s (levels %s); variance / model %s (expected %s); truth noiseRatio %s: %s\n'], nS9, okW, ...
    mat2str(round(mu9, 3)), mat2str(P9.level), mat2str(round(rLR, 3)), mat2str(rExp9), ...
    mat2str(round(t9.noiseRatio, 3)), passStr(pass));
if ~pass, fails{end+1} = 'RFI profile'; end

% ---------------------------------------------------------------------------------
% 10. Schedule: passes and sub-int time
% ---------------------------------------------------------------------------------
nT10 = 2000; nP10 = 50;
[~, fa, ta] = simulateFold(f0, -5, nT10, 'SubintPeriods', nP10, 'Seed', 11, 'Verbose', false);
pass1 = [TRef - 0.5/f0 - 1e-7, TRef + (nT10 - 0.5)/f0 + 1e-7];
[~, fb, tb] = simulateFold(f0, -5, [], 'Passes', pass1, 'SubintTime', nP10/f0, 'Seed', 11, ...
    'Verbose', false);
same10 = isequaln(fa, fb) && isequal(ta.pulseCenters, tb.pulseCenters);
P10 = [0 30; 100 130.004];
[ic, fc, tc] = simulateFold(f0, -5, [], 'Passes', P10, 'SubintTime', 10, 'PhaseOffset', 0.001, ...
    'Seed', 12, 'Verbose', false);
su = fc.subint;
tStartT = TRef + (su.turnFirst - 0.5)/f0; tEndT = TRef + (su.turnLast + 0.5)/f0;   % turn edges
inside = all(tStartT >= P10(tc.pass, 1) - 1e-9 & tEndT <= P10(tc.pass, 2) + 1e-9);
expRef = [500 1500 2500 10500 11500 12500].';
okS = ic.nSub == 6 && isequal(tc.pass, [1 1 1 2 2 2].') && inside && isequal(su.turnRef, expRef) && ...
    max(abs(su.tRef - (TRef + expRef/f0))) < 1e-12 && ic.subintPeriods == 1000;
toa10 = foldToTOA(fc, ic, fwhm, Bn, 'Verbose', false);
val10 = validateTOA(toa10, tc, 'Plot', false);
z10 = (toa10.toa - tc.pulseCenters) ./ toa10.toaErr;
okE = all(toa10.valid) && abs(mean(tc.pulseCenters - su.tRef) - 0.001/f0) < 1e-12 && val10.n == 6;
pass = same10 && okS && okE;
fprintf(['10. schedule: one pass = nTurns interface bit for bit %d; two passes -> %d sub-ints (passes %s), ' ...
         'turns inside the passes %d, turnRef %s, tRef %s s; TOAs at -5 dB: %d valid, pulls %s: %s\n'], ...
    same10, ic.nSub, mat2str(tc.pass.'), inside, mat2str(su.turnRef.'), mat2str(round(su.tRef.', 3)), ...
    nnz(toa10.valid), mat2str(round(z10.', 2)), passStr(pass));
if ~pass, fails{end+1} = 'schedule'; end

% ---------------------------------------------------------------------------------
% 11. Scintillation: FluxScale and the random field
% ---------------------------------------------------------------------------------
Fh = [0.25 * ones(200, 1); ones(200, 1)];
[i11, f11] = simulateFold(f0, -5, 400, 'FluxScale', Fh, 'Seed', 13, 'Verbose', false);
t11 = foldToTOA(f11, i11, fwhm, Bn, 'Verbose', false);
aLo = t11.amp(1:200); aHi = t11.amp(201:400);
rA = mean(aLo) / mean(aHi);
eA = rA * sqrt(var(aLo)/200/mean(aLo)^2 + var(aHi)/200/mean(aHi)^2);
okF = abs(rA - 0.25) < 4 * eA;
S11 = struct('bandwidth', 50e6, 'timescale', 60);
[~, ~, ts11] = simulateFold(f0, -Inf, 20000 * 600, 'SubintPeriods', 600, 'NBin', 64, 'NChan', 64, ...
    'Bnoise', Inf, 'Scintillation', S11, 'Seed', 14, 'Verbose', false);
F11 = ts11.fluxScale;                               % [20000 x 64]
mF = mean(F11(:)); vF = var(F11(:)); p01 = mean(F11(:) < 0.1);
acf = @(a, b) (mean(a .* b, 'all') - mean(a, 'all') * mean(b, 'all')) / (mean(a.^2, 'all') - mean(a, 'all')^2);
rT = acf(F11(1:end-10, :), F11(11:end, :));         % 10 sub-ints = 60 s
rF = acf(F11(:, 1:end-8), F11(:, 9:end));           % 8 channels = 50 MHz
[~, ~, tb11] = simulateFold(f0, -Inf, 20000 * 600, 'SubintPeriods', 600, 'NBin', 64, 'NChan', 1, ...
    'Bnoise', Inf, 'Scintillation', S11, 'Seed', 15, 'Verbose', false);
mB = mean(tb11.fluxScale); vB = var(tb11.fluxScale);
[~, fOff] = simulateFold(f0, -5, 20, 'Seed', 16, 'Verbose', false);           % switch off:
[~, fNo]  = simulateFold(f0, -5, 20, 'Scintillation', [], 'Seed', 16, 'Verbose', false);   % [] = off
okOff = isequaln(fOff, fNo);
okS = abs(mF - 1) < 0.05 && abs(vF - 1) < 0.1 && abs(p01 - (1 - exp(-0.1))) < 0.01 && ...
    abs(rT - exp(-1)) < 0.06 && abs(rF - 0.5) < 0.06 && abs(mB - 1) < 0.05 && vB < 0.3 && okOff;
pass = okF && okS;
fprintf(['11. scintillation: FluxScale 0.25 -> fitted amplitude ratio %.3f (+- %.3f); field: mean %.3f, ' ...
         'var %.3f, P(F<0.1) %.4f (exp: %.4f), ACF at the scintle time %.3f (1/e %.3f), at the ' ...
         'scintle bandwidth %.3f (0.5); band average: mean %.3f, var %.3f; ''Scintillation'', [] = ' ...
         'off, bit for bit %d: %s\n'], rA, eA, mF, vF, p01, 1 - exp(-0.1), rT, exp(-1), rF, mB, vB, ...
    okOff, passStr(pass));
if ~pass, fails{end+1} = 'scintillation'; end

% ---------------------------------------------------------------------------------
% 12. Adaptive sub-ints: mergeSubints exactness, adaptiveSubints statistics
% ---------------------------------------------------------------------------------
[~, fp5] = foldProfile(info_det, fullfile(tmp, 'fold5.mat'), f0, 'TRef', TRef, 'NBin', N, ...
    'SubintPeriods', 5, 'SaveFile', false, 'Verbose', false);
grp = floor((fp.subint.turnFirst - fp.subint.turnFirst(1)) / 5) + 1;
[fm, im] = mergeSubints(fp, ifp, grp);
rel = @(a, b) max(abs(a(:) - b(:))) / max(abs(b(:)));
okM = isequal(size(fm.sum), size(fp5.sum)) && rel(fm.sum, fp5.sum) < 1e-12 && ...
    rel(fm.weight, fp5.weight) < 1e-12 && rel(fm.weight2, fp5.weight2) < 1e-12 && ...
    rel(fm.weightX, fp5.weightX) < 1e-12 && max(abs(fm.prof(:) - fp5.prof(:)), [], 'omitnan') < 1e-12 && ...
    isequal(fm.subint.turnRef, fp5.subint.turnRef) && max(abs(fm.subint.tRef - fp5.subint.tRef)) < 1e-12 && ...
    im.nSub == size(fp5.sum, 2);
tmpl = gaussianTemplate(N, fwhm);
pullsA = cell(1, 2); nMed = zeros(1, 2); okA = true;
for kase = 1:2
    if kase == 1, sArg = {}; else, sArg = {'Scintillation', struct('bandwidth', 400e6, 'timescale', 2)}; end
    [ia, fa, ta] = simulateFold(f0, -25, 4000, 'PhaseOffset', 0.002, 'Seed', 30 + kase, sArg{:}, 'Verbose', false);
    [fA, iA, sA] = adaptiveSubints(fa, ia, tmpl, Bn, 'Verbose', false);
    tA = foldToTOA(fA, iA, fwhm, Bn, 'Verbose', false);
    c  = sA.complete & tA.valid;
    zA = (tA.toa(c) - (fA.subint.tRef(c) + 0.002/f0)) ./ tA.toaErr(c);
    pullsA{kase} = zA; nMed(kase) = median(sA.nMerged(sA.complete));
    okA = okA && all(sA.T0(sA.complete) >= 7) && abs(mean(zA)) < 4/sqrt(numel(zA)) && ...
        abs(std(zA) - 1) < 4/sqrt(2*numel(zA));
    if kase == 2                                    % bright stretches close faster
        Fg = accumarray(sA.group, ta.fluxScale(:, 1), [], @mean);
        Fg = Fg(sA.complete); nG2 = sA.nMerged(sA.complete);
        okA = okA && median(nG2(Fg > median(Fg))) < median(nG2(Fg <= median(Fg)));
    end
end
[ig, fg, tg] = simulateFold(f0, -25, [], 'Passes', [0 20; 30 50], 'Seed', 40, 'Verbose', false);
[~, ~, sG] = adaptiveSubints(fg, ig, tmpl, Bn, 'Verbose', false);
noSpan = all(accumarray(sG.group, tg.pass, [], @(p) numel(unique(p))) == 1);
pass = okM && okA && noSpan;
fprintf(['12. adaptive sub-ints: mergeSubints = 5-turn foldProfile fold %d; adaptiveSubints at -25 dB: ' ...
         '%d complete groups, median %g sub-ints, pulls mean %+.3f std %.3f; with scintillation: %d ' ...
         'groups, median %g, pulls mean %+.3f std %.3f; no group across the gap %d: %s\n'], okM, ...
    numel(pullsA{1}), nMed(1), mean(pullsA{1}), std(pullsA{1}), numel(pullsA{2}), nMed(2), ...
    mean(pullsA{2}), std(pullsA{2}), noSpan, passStr(pass));
if ~pass, fails{end+1} = 'adaptive sub-ints'; end

% ---------------------------------------------------------------------------------
% 13. NPol (dual polarization) and AcrossGaps (groups over the gap between passes)
% ---------------------------------------------------------------------------------
[~, fn0] = simulateFold(f0, -5, 50, 'Seed', 13, 'Verbose', false);
[~, fn1] = simulateFold(f0, -5, 50, 'Seed', 13, 'NPol', 1, 'Verbose', false);
[~, fr0] = simulateFold(f0, -Inf, 20, 'NChan', nC9, 'RFI', P9, 'Seed', 13, 'Verbose', false);
[~, fr1] = simulateFold(f0, -Inf, 20, 'NChan', nC9, 'RFI', P9, 'Seed', 13, 'NPol', 1, 'Verbose', false);
okP1 = isequal(fn0, fn1) && isequal(fr0, fr1);      % NPol 1 = no option, bit for bit
[~, fn2] = simulateFold(f0, -5, 50, 'Seed', 13, 'NPol', 2, 'Verbose', false);
[~, fnB] = simulateFold(f0, -5, 50, 'Seed', 13, 'Bnoise', 2 * Bn, 'Verbose', false);
dP2 = max(abs(fn2.prof(:) - fnB.prof(:))) / max(abs(fnB.prof(:)));   % NPol 2 = Bnoise x 2
[~, f13] = simulateFold(f0, -Inf, nS9, 'NChan', nC9, 'RFI', P9, 'Seed', 9, 'NPol', 2, 'Verbose', false);
rL2 = zeros(1, nC9);
for c = 1:nC9                                       % observer's model with Bnoise x 2
    sM = P9.level(c) * sqrt(f13.weight2(:, :, c) / (2*BnC*dtC)) ./ f13.weight(:, :, c);
    u = (f13.prof(:, :, c) - P9.level(c)) ./ sM;
    rL2(c) = mean(u(:).^2);
end
okP2 = dP2 < 1e-12 && all(abs(rL2 - rExp9) < 4 * rExp9 * sqrt(2*1.2 / nPer9));
[~, ~, sD] = adaptiveSubints(fg, ig, tmpl, Bn, 'AcrossGaps', false, 'Verbose', false);
[fX, iX, sX] = adaptiveSubints(fg, ig, tmpl, Bn, 'AcrossGaps', true, 'Verbose', false);
spanX = accumarray(sX.group, tg.pass, [], @(p) numel(unique(p))) > 1;
tX = foldToTOA(fX, iX, fwhm, Bn, 'Verbose', false);
zX = (tX.toa(spanX & tX.valid) - fX.subint.tRef(spanX & tX.valid)) ./ tX.toaErr(spanX & tX.valid);
[~, ~, sF] = adaptiveSubints(fg, ig, tmpl, Bn, 'TargetSNR', Inf, 'MaxSubints', 10, 'AcrossGaps', true, ...
    'Verbose', false);
okG = isequal(sD, sG) && any(spanX) && nnz(~sX.complete) <= 1 && nnz(~sD.complete) >= 1 && ...
    ~isempty(zX) && all(abs(zX) < 5) && all(sF.nMerged(1:end-1) == 10);
pass = okP1 && okP2 && okG;
fprintf(['13. NPol: NPol 1 = no option %d; NPol 2 = Bnoise x 2 (max rel diff %.1e); with the RFI profile, ' ...
         'variance / model (Bnoise x 2) %s (expected %s). AcrossGaps: default = before %d; %d of %d groups ' ...
         'span the gap, incomplete groups %d (default %d), spanning TOAs z %s; fixed groups of 10 %d: %s\n'], ...
    okP1, dP2, mat2str(round(rL2, 3)), mat2str(rExp9), isequal(sD, sG), nnz(spanX), numel(spanX), ...
    nnz(~sX.complete), nnz(~sD.complete), mat2str(round(zX.', 2)), all(sF.nMerged(1:end-1) == 10), passStr(pass));
if ~pass, fails{end+1} = 'NPol / AcrossGaps'; end

if isempty(fails)
    fprintf('testSimulateFold: ALL PASSED\n');
else
    error('testSimulateFold:failed', 'FAILED: %s', strjoin(fails, ', '));
end
end


function s = passStr(p)
if p, s = 'PASS'; else, s = 'FAIL'; end
end
