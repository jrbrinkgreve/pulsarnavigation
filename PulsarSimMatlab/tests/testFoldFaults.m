function testFoldFaults()
%TESTFOLDFAULTS  TOA fit and detection under faults of a real receiver: each case puts one fault into the fold.
%{
Run from the PulsarSimMatlab folder: run('tests/testFoldFaults.m'), or runAllTests("testFoldFaults").
~3 min, ~2 GB of memory; no files. The fast simulator (simulateFold) draws 128-channel folds as the channel path
makes them (the noise coefficients and noise bandwidth of the frozen seed-43 power file,
data/test/); each case changes them the way a real receiver would (truth side, in this test),
and foldToTOA / estimateTOA / detectPulsar run on them unchanged. Risk register:
logs/2026-11-02_risk-register.md; these cases: logs/2026-11-03_real-hardware-faults.md.

  1. FM-41, the bandpass. A real receiver's gain G_c varies across the band (antenna match, LNA
     and cable slopes, standing-wave ripple, the SDR's filter): it scales channel c's baseline
     and pulse alike. The 'optimal' fit models the pulse as b * s_c * T with s_c = ChannelGain,
     which foldToTOA (and so processIQ's chain) leaves at 1: the same absolute pulse in every
     channel. It weights each channel by 1/variance ~ 1/G_c^2, so the low-gain channels count
     most, where the pulse is smallest. For a weak pulsar the SNR keeps
       eta_opt = sum(1/G) / sqrt(n sum(1/G^2)),
     and the 'equal' sum keeps eta_eq = sum(G) / sqrt(n sum(G^2)) (both 1 for a flat band).
     The bandpass here: a 6 dB tilt across the band (+3 dB at 1.2 GHz, -3 dB at 1.6 GHz) and a
     +-1 dB ripple of 20 MHz period. Checked, on the same noise draw with and without it: the
     total fold's SNR ratio within 0.005 of eta_opt and of eta_eq; with ChannelGain = G_c the
     optimal fit is the flat fit (data and weights scale together): TOAs within 1e-6 sigma, SNR
     within 1e-6; with ChannelGain = each channel's mean power (what an observer measures) the
     SNR within 0.005 of the flat fold's; with ChannelGain 1 the error bars stay honest: pulls
     against the truth with rms within 4/sqrt(2n) of 1. Shown: eta_opt and eta_eq for tilts of
     3, 6 and 10 dB, and for the pulsar's own spectrum (index -1.6, flat bandpass).
  2. FM-42, common-mode gain noise. The LNA's and the receiver's gain fluctuate (1/f noise,
     temperature, supply ripple), and the fluctuation multiplies every channel's power at once.
     The noise model holds the radiometer noise only, independent between channels. A common
     fluctuation adds the same amount R (in units of the combined profile's radiometer variance
     per bin) to the combined profile, but only R/nChan to each channel: the per-channel red.
     chi^2 hardly sees it, while every statistic of the combined profile takes all of it. White
     gain noise with R = 1 is folded into 300 sub-ints without a pulsar and into case 1's 300 with
     it, against the same draws without the gain noise. Checked: detectPulsar's noise ratio (the
     combined profile, off-pulse) +R within 0.03; the variance of the known-phase statistic T0
     +R sum|C_k|^2 / sum|C_k|^2 lam_k, and of the TOAs (in sigma) +R sum k^2|S_k|^2 / sum
     k^2|S_k|^2 lam_k, both within 4 sqrt(2/300) (lam_k: the spectrum of the radiometer part, which
     the fold's linear assignment correlates between neighbouring bins; C_k, S_k: the template's
     harmonics); the per-channel red. chi^2 +R/128 within 0.001. Shown: unknown-phase false
     alarms against Rice's formula; for B0329+54, a 1/f gain spectrum with knee f_k (where it
     equals the full band's radiometer noise: R_k = f_k / (k f0) at harmonic k, analytic): the
     detection statistic's variance and false alarms, the TOA error, the noise ratio and the
     per-channel chi^2, for knees of 0.1, 1 and 10 Hz.
Errors at the end if any check fails. (3 Nov 2026, risk register FM-41, FM-42.)
%}

root = fullfile(fileparts(mfilename('fullpath')), '..');
addpath(fullfile(root, 'functions'));
oldDir = cd(root); restoreDir = onCleanup(@() cd(oldDir));
run(fullfile(fileparts(mfilename('fullpath')), 'testParams.m'));   % frozen test settings, data/test/ (not pipelineParams)
fails = {};
detR = loadInfo(fullfile(dataDir, 'chan', 'test_dedisp_chan_power.dat'));
nc = [detR.noise.V, detR.noise.X]; Bc = detR.noise.Bnoise; dtB = detR.binDt;
nChan = detR.nChan;
fC = fLow + ((1:nChan) - 0.5) * chanWidth;                % [Hz] channel centres
f0 = 100; fw = 0.05; nb = 512;
sim = {'NChan', nChan, 'BinDt', dtB, 'Bnoise', Bc, 'NoiseCoeffs', nc, 'NBin', nb, 'FWHM', fw, ...
       'Verbose', false};
tmpl = gaussianTemplate(nb, fw);

% ---------------------------------------------------------------------------------
% 1. FM-41: the bandpass against the fit's equal channel amplitudes
% ---------------------------------------------------------------------------------
tiltDB = @(dB) dB/2 - dB * (fC - fLow) / (fHigh - fLow);  % +dB/2 at the bottom, -dB/2 at the top
G = 10.^((tiltDB(6) + sin(2*pi * fC / 20e6)) / 10);      % 6 dB tilt, +-1 dB ripple of 20 MHz
etaOpt = @(g) sum(1 ./ g) / sqrt(numel(g) * sum(1 ./ g.^2));
etaEq  = @(g) sum(g) / sqrt(numel(g) * sum(g.^2));
nSub = 300; spt = 20;                                     % sub-ints of 20 turns, SNR ~20 each
[info, fold, truth] = simulateFold(f0, -20, nSub * spt, 'SubintPeriods', spt, 'Seed', 11, sim{:});
foldB = scaleChannels(fold, G);
q = {'Bnoise', Bc, 'Verbose', false};
t0  = estimateTOA(fold,  info, tmpl, q{:}, 'Weighting', 'optimal');
tB  = estimateTOA(foldB, info, tmpl, q{:}, 'Weighting', 'optimal');
t0e = estimateTOA(fold,  info, tmpl, q{:}, 'Weighting', 'equal');
tBe = estimateTOA(foldB, info, tmpl, q{:}, 'Weighting', 'equal');
tBG = estimateTOA(foldB, info, tmpl, q{:}, 'Weighting', 'optimal', 'ChannelGain', G);
gM  = mean(foldB.profTotal, 1);                           % each channel's mean power (observer)
tBM = estimateTOA(foldB, info, tmpl, q{:}, 'Weighting', 'optimal', 'ChannelGain', gM / mean(gM));
rOpt = tB.total.snr / t0.total.snr;
rEq  = tBe.total.snr / t0e.total.snr;
v = t0.valid & tB.valid & tBG.valid;
dG = max(abs(tBG.toa(v) - t0.toa(v)) ./ t0.toaErr(v));
sG = abs(tBG.total.snr / t0.total.snr - 1);
rM = tBM.total.snr / t0.total.snr;
pull = (tB.toa(tB.valid) - truth.pulseCenters(tB.valid)) ./ tB.toaErr(tB.valid);
pull0 = (t0.toa(t0.valid) - truth.pulseCenters(t0.valid)) ./ t0.toaErr(t0.valid);
nP = numel(pull);
ok = abs(rOpt - etaOpt(G)) < 0.005 && abs(rEq - etaEq(G)) < 0.005 && dG < 1e-6 && sG < 1e-6 && ...
     abs(rM - 1) < 0.005 && abs(rms(pull) - 1) < 4 / sqrt(2 * nP);
fprintf(['1. FM-41 bandpass (6 dB tilt, +-1 dB ripple; gains %.2f-%.2f): SNR x %.4f optimal (predicted ' ...
         '%.4f), x %.4f equal (predicted %.4f); ChannelGain = G: TOAs within %.1e sigma, SNR %.1e of the ' ...
         'flat fit; ChannelGain = measured channel power: SNR x %.4f; pulls with ChannelGain 1: rms %.3f ' ...
         '(flat %.3f; %d sub-ints, within %.3f of 1), median TOA error x %.3f: %s\n'], min(G), max(G), ...
    rOpt, etaOpt(G), rEq, etaEq(G), dG, sG, rM, rms(pull), rms(pull0), nP, 4 / sqrt(2 * nP), ...
    median(tB.toaErr(v) ./ t0.toaErr(v)), passStr(ok));
if ~ok, fails{end+1} = 'FM-41 bandpass'; end
clear foldB
tilts = [3 6 10];
spec = (fC / 1.4e9).^-1.6;                                % pulse fraction ~ the pulsar's spectrum
fprintf(['   shown: tilt only, %s dB: optimal x %s, equal x %s; the pulsar''s spectrum (index -1.6, flat ' ...
         'band; the pulse fraction varies, the noise does not): x %.4f for both\n'], mat2str(tilts), ...
    mat2str(arrayfun(@(d) etaOpt(10.^(tiltDB(d)/10)), tilts), 4), ...
    mat2str(arrayfun(@(d) etaEq(10.^(tiltDB(d)/10)), tilts), 4), sum(spec) / sqrt(nChan * sum(spec.^2)));

% ---------------------------------------------------------------------------------
% 2. FM-42: common-mode gain noise
% ---------------------------------------------------------------------------------
R = 1;                                                    % gain noise / combined radiometer noise, per bin
[iH, foldH] = simulateFold(f0, -Inf, nSub * spt, 'SubintPeriods', spt, 'Seed', 21, sim{:});   % no pulsar
vC = mean(fold.weight2(:, 1) ./ fold.weight(:, 1).^2) / (Bc * dtB);   % one channel's variance per bin
rsG = RandStream('mt19937ar', 'Seed', 31);
sG2 = sqrt(R * vC / nChan);                               % rms of the folded gain per bin
foldHg = gainFold(foldH, sG2 * randn(rsG, nb, nSub));
foldg  = gainFold(fold,  sG2 * randn(rsG, nb, nSub));
dq = {'Bnoise', Bc, 'Weighting', 'optimal', 'Verbose', false};
dH  = detectPulsar(foldH,  iH, tmpl, dq{:});
dHg = detectPulsar(foldHg, iH, tmpl, dq{:});
tg  = estimateTOA(foldg, info, tmpl, q{:}, 'Weighting', 'optimal');
% predictions: the radiometer part of the combined profile is correlated between neighbouring
% bins (the fold's linear assignment: eigenvalues lam_k), the gain part is white
rD = reshape(fold.weightX(1, 1, 1, :), 1, []) / fold.weight2(1, 1);
kk = (0:nb-1).'; ks = min(kk, nb - kk);                   % |signed harmonic|
lam = 1 + 2 * sum(rD .* cos(2*pi * kk * (1:numel(rD)) / nb), 2);
Cs = abs(fft(tmpl - mean(tmpl))).^2;
Sd = ks.^2 .* abs(fft(tmpl)).^2;                          % the TOA's derivative template
fT0  = R * sum(Cs) / sum(Cs .* lam);
fTOA = R * sum(Sd) / sum(Sd .* lam);
v = t0.valid & tg.valid;
nrH = median(dHg.noiseRatio - dH.noiseRatio);
vT0 = var(dHg.T0 - dH.T0);
vTOA = var((tg.toa(v) - t0.toa(v)) ./ t0.toaErr(v));
dChi = mean(tg.redChi2(v) - t0.redChi2(v));
tolV = 4 * sqrt(2 / nSub);
ok = abs(nrH - R) < 0.03 && abs(vT0 / fT0 - 1) < tolV && abs(vTOA / fTOA - 1) < tolV && ...
     abs(dChi - R / nChan) < 0.001;
fprintf(['2. FM-42 common-mode gain noise, white, R = %g x the combined radiometer variance per bin: ' ...
         'noise ratio +%.3f (predicted +%g), known-phase T0 variance +%.3f (predicted +%.3f), TOA ' ...
         'variance +%.3f sigma^2 (predicted +%.3f; both within %.0f %%), per-channel red. chi2 +%.4f ' ...
         '(predicted R/nChan = %.4f): %s\n'], R, nrH, R, vT0, fT0, vTOA, fTOA, 100 * tolV, dChi, ...
    R / nChan, passStr(ok));
if ~ok, fails{end+1} = 'FM-42 gain noise'; end
% shown: false alarms at unknown phase, Rice's formula with the inflated, whiter process
eta = dH.etaUnknown;
Qf = @(x) 0.5 * erfc(x / sqrt(2));
w2 = (2*pi * ks).^2;
pfa = @(lamR, s2) Qf(eta / sqrt(s2)) + sqrt(sum(w2 .* Cs .* lamR) / sum(Cs .* lamR)) / (2*pi) * ...
    exp(-eta^2 / (2 * s2));
fprintf(['   shown: unknown-phase false alarms in %d sub-ints without a pulsar: %d without gain noise ' ...
         '(Rice %.1f), %d with it (Rice %.1f)\n'], nSub, nnz(dH.detectedUnknown), nSub * pfa(lam, 1), ...
    nnz(dHg.detectedUnknown), nSub * pfa(lam + R, 1 + fT0));
% shown: B0329+54 with a 1/f gain spectrum, knee f_k where it equals the full band's radiometer
% noise: R_k = f_k / (k f0) at harmonic k (analytic, white radiometer noise)
comps = [0.159 -11.427 1.479; 0.439 -1.123 3.281; 0.585 0 1.055; 0.152 7.719 2.751; 0.170 9.482 0.847];
comps = [comps(:, 1), comps(:, 2) / 360, comps(:, 3) * 2*sqrt(2*log(2)) / 360];
nP = 1024; f0P = 1.399541538720;
TP = reshape(gaussianTemplate(nP, comps), [], 1);
kP = (1:nP/2-1).';
CP = abs(fft(TP - mean(TP))).^2; CP = CP(kP + 1);
SP = abs(fft(TP)).^2;            SP = kP.^2 .* SP(kP + 1);
for fk = [0.1 1 10]
    Rk = fk ./ (kP * f0P);
    sT = 1 + sum(CP .* Rk) / sum(CP);
    etaP = fzero(@(e) Qf(e) + sqrt(sum((2*pi*kP).^2 .* CP) / sum(CP)) / (2*pi) * exp(-e^2/2) - 1e-3, [0 40]);
    pP = Qf(etaP / sqrt(sT)) + sqrt(sum((2*pi*kP).^2 .* CP .* (1 + Rk)) / sum(CP .* (1 + Rk))) / (2*pi) * ...
        exp(-etaP^2 / (2 * sT));
    fprintf(['   shown, B0329+54 (5 components, NBin %d), 1/f gain noise with knee %4.1f Hz: detection ' ...
             'statistic variance x %.2f (false alarms per sub-int %.1e instead of 1e-3), TOA error x %.3f, ' ...
             'noise ratio %.3f, per-channel red. chi2 +%.4f\n'], nP, fk, sT, pP, ...
        sqrt(1 + sum(SP .* Rk) / sum(SP)), 1 + mean(Rk), mean(Rk) / nChan);
end

if isempty(fails)
    fprintf('testFoldFaults: ALL PASSED\n');
else
    error('testFoldFaults:failed', 'FAILED: %s', strjoin(fails, ', '));
end
end


% =================================================================================
function f = scaleChannels(f, g)
%SCALECHANNELS  A receiver gain g(c) on channel c: the fold's power (profiles, sums, total) x g(c);
% the weights count time bins and stay.
g3 = reshape(g, 1, 1, []);
f.prof = f.prof .* g3;
f.sum  = f.sum .* g3;
f.profTotal = f.profTotal .* reshape(g, 1, []);
end

function f = gainFold(f, dbar)
%GAINFOLD  A gain 1 + dbar(j, s), common to all channels, on the power of phase bin j in sub-int s
% (the receiver's gain fluctuation, folded); the weights count time bins and stay.
f.prof = f.prof .* (1 + dbar);
f.sum  = f.sum .* (1 + dbar);
f.profTotal = reshape(sum(f.sum, 2), size(f.sum, 1), []) ./ reshape(sum(f.weight, 2), size(f.sum, 1), []);
end

function s = passStr(p)
if p, s = 'PASS'; else, s = 'FAIL'; end
end
