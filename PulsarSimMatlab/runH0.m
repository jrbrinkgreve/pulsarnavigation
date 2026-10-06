%{
RUNH0 - long noise-only (H0) run: detection statistics and noise correlations

Tests the H0 side of detectPulsar with many independent noise-only profiles:
  - T0 (matched filter at the known phase) must be N(0,1): mean 0, std 1;
  - exceedances of the thresholds at several P_FA must match the nominal
    rates (known phase: Phi^-1; unknown phase: Rice);
  - the normalized profile (normProfile, in units of the model sigma) must
    be correlated only between neighbouring phase bins (linear assignment,
    weightX). Correlation at larger lags would make the sigma of
    template-weighted sums (T0, and the TOA error bars) too small.
Motivation: the phase E sweep gave T0 std 1.102 +- 0.075 under H0, and the
phase D error-bar ratios were 1.03-1.06 (all from the same 10 noise
realizations).

Per pass: zero sky (generator with A = 0, made once) + receiver noise
('NoiseStd' 1, own seed seed + 10000 + i, no overlap with the sweeps), the
same processing as main.m, then detectPulsar. With L = 0.1 s and 1 turn per
sub-int: 9 complete sub-ints per pass, about 8 s per pass.

Parameters from pipelineParams.m; files in data/mc with an h0_ prefix.
%}

nPass  = 100;                                    % noise-only passes (~9 profiles each)
pfaList = [1e-1 1e-2 1e-3];                      % false-alarm levels to test
maxLag = 200;                                    % [bins] lags shown in the figure

% Paths and parameters
scriptDir = fileparts(mfilename('fullpath'));
addpath(fullfile(scriptDir, 'functions'));
pipelineParams;

mcDir      = fullfile(dataDir, "mc");
h0Zero     = fullfile(mcDir, "h0_zero.dat");
h0Rx       = fullfile(mcDir, "h0_rx.dat");
h0IQ       = fullfile(mcDir, "h0_rx_IQ.dat");
h0Dedisp   = fullfile(mcDir, "h0_dedisp.dat");
h0Envelope = fullfile(mcDir, "h0_envelope.dat");
h0Fold     = fullfile(mcDir, "h0_fold.mat");     % not saved (SaveFile false)
h0Result   = fullfile(mcDir, "h0Run.mat");

template = gaussianTemplate(nBin, ephem.profileFWHM);
c = template / max(template); c = c - mean(c);  % zero-mean template, as in detectPulsar

fprintf('runH0: %d noise-only passes, L = %.3g s, %d turn(s) per sub-int\n', nPass, L, subintPeriods);
tAll = tic;
info_zero = generatePulsarSignal(h0Zero, T, f_in, 0, L, dutycycle, ...
    'Seed', seed, 'EnvelopeMode', genEnvMode, 'Verbose', false);

T0all = []; Tmall = []; nrAll = [];
acSum = zeros(nBin, 1);                          % sum of circular autocorrelations of normProfile
rho1Sum = 0; nProf = 0;                          % model lag-1 correlation (weightX)
for i = 1:nPass
    tRun = tic;
    ns = seed + 10000 + i;
    info_rx = addNoiseAndRFI(info_zero.file, h0Rx, info_zero.actualFsOut, ...
        'Band', [fLow fHigh], 'NoiseStd', 1, 'RFI', rfi, 'Seed', ns, 'Verbose', false);
    info_IQ = applyIQmodulation(info_rx.file, h0IQ, info_rx.actualFsOut, fs, fLO, ...
        'FilterOrder', filterOrder, 'Band', [fLow fHigh], 'Verbose', false);
    info_dedisp = applyInverseDispersion(info_IQ.file, h0Dedisp, info_IQ.actualFsOut, ...
        info_IQ.fLO, ephem.DM, fLow, fHigh, 'RefFreq', refFreq, 'Verbose', false);
    info_det = detectPower(info_dedisp.file, h0Envelope, info_dedisp.fs, info_dedisp.fLO, ...
        f_out, 'FullySupported', info_dedisp.fullySupported, 'Verbose', false);
    [info_fold, fold] = foldProfile(info_det, h0Fold, ephem.f0, ...
        'F1', ephem.F1, 'TRef', ephem.TRef, 'NBin', nBin, 'SubintPeriods', subintPeriods, ...
        'SaveFile', false, 'Verbose', false);
    Bnoise = noiseBandwidth(fLow, fHigh, info_dedisp.edgeWidth);
    [detection, info_detect] = detectPulsar(fold, info_fold, template, 'Bnoise', Bnoise, ...
        'Verbose', false);

    t = detection.tested;
    T0all = [T0all; detection.T0(t)];            %#ok<AGROW>
    Tmall = [Tmall; detection.Tmax(t)];          %#ok<AGROW>
    nrAll = [nrAll; detection.noiseRatio(t)];    %#ok<AGROW>
    for s = find(t).'
        z = detection.normProfile(:, s);
        acSum = acSum + real(ifft(abs(fft(z)).^2)) / nBin;   % circular autocorrelation
        W  = fold.weight(:, s);  W2 = fold.weight2(:, s);  WX = fold.weightX(:, s);
        rho1Sum = rho1Sum + mean(WX ./ sqrt(W2 .* circshift(W2, -1)));
        nProf = nProf + 1;
    end
    fprintf('  pass %3d/%d (seed %d): %d profiles, T0 mean %+.2f, max Tmax %.2f, %.0f s\n', ...
        i, nPass, ns, nnz(t), mean(detection.T0(t)), max(detection.Tmax(t)), toc(tRun));
end
fprintf('runH0: done in %.1f min\n', toc(tAll)/60);

% ---- Statistics -------------------------------------------------------------------------
n  = numel(T0all);
ac = acSum / nProf;                              % measured autocorrelation, lag 0..nBin-1
rho1 = rho1Sum / nProf;                          % model lag-1 correlation
acModel = zeros(nBin, 1); acModel(1) = 1; acModel(2) = rho1; acModel(end) = rho1;
Rc = real(ifft(abs(fft(c)).^2));                 % template autocorrelation (circular)
predRatio = sqrt(sum(Rc .* ac) / sum(Rc .* acModel));   % T0 std implied by the measured correlation

Q   = @(x) 0.5 * erfc(x / sqrt(2));              % 1 - Phi(x)
rf  = info_detect.riceFactor;
riceP = @(e) Q(e) + rf * exp(-e.^2 / 2);
etaK = sqrt(2) * erfcinv(2 * pfaList);
etaU = arrayfun(@(p) fzero(@(e) riceP(e) - p, [0 40]), pfaList);

fprintf('\nrunH0 summary (%d noise-only profiles)\n', n);
fprintf('  T0: mean %+.4f (expect 0 +- %.3f), std %.4f (expect 1 +- %.3f)\n', ...
    mean(T0all), 1/sqrt(n), std(T0all), 1/sqrt(2*n));
fprintf('  noise ratio: median %.4f, mean %.4f\n', median(nrAll), mean(nrAll));
fprintf('  normProfile correlation: lag 0 %.4f, lag 1 %.4f (model %.4f), lags 2-%d mean %+.5f, max |r| %.4f (noise ~%.4f)\n', ...
    ac(1), ac(2), rho1, maxLag, mean(ac(3:maxLag+1)), max(abs(ac(3:maxLag+1))), 1/sqrt(nBin*nProf));
fprintf('  T0 std implied by the measured correlation: %.4f (measured %.4f)\n', predRatio, std(T0all));
fprintf('  exceedances      P_FA   expected   T0 > eta_known (eta)   Tmax > eta_unknown (eta)\n');
for k = 1:numel(pfaList)
    fprintf('               %8.3g   %6.1f +- %4.1f   %5d (%.2f)            %5d (%.2f)\n', ...
        pfaList(k), n*pfaList(k), sqrt(n*pfaList(k)), nnz(T0all > etaK(k)), etaK(k), ...
        nnz(Tmall > etaU(k)), etaU(k));
end

save(h0Result, 'T0all', 'Tmall', 'nrAll', 'ac', 'acModel', 'rho1', 'Rc', 'predRatio', ...
    'pfaList', 'etaK', 'etaU', 'rf', 'nPass', 'L', 'subintPeriods');
fprintf('runH0: results saved to %s\n', h0Result);

% ---- Figure -------------------------------------------------------------------------------
figure('Name', 'H0 run', 'Color', 'w');
tl = tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, sprintf('Noise-only (H0) run: %d profiles', n));

ax1 = nexttile(tl);
histogram(ax1, T0all, 40, 'Normalization', 'pdf', 'FaceColor', [0.6 0.6 0.6]); hold(ax1, 'on');
zz = linspace(-5, 5, 400);
plot(ax1, zz, exp(-zz.^2/2)/sqrt(2*pi), 'r', 'LineWidth', 1.5); hold(ax1, 'off');
grid(ax1, 'on'); xlabel(ax1, 'T0 (known phase)'); ylabel(ax1, 'pdf');
title(ax1, sprintf('T0 vs N(0,1): mean %+.3f, std %.3f', mean(T0all), std(T0all)));

ax2 = nexttile(tl);
ee = linspace(0, 6, 300);
sm = sort(Tmall); s0 = sort(T0all);
semilogy(ax2, sm, (n:-1:1)/n, 'b.', 'DisplayName', 'T_{max} measured'); hold(ax2, 'on');
semilogy(ax2, ee, min(riceP(ee), 1), 'b-', 'DisplayName', 'Rice (unknown phase)');
semilogy(ax2, s0, (n:-1:1)/n, 'k.', 'DisplayName', 'T0 measured');
semilogy(ax2, ee, Q(ee), 'k-', 'DisplayName', 'Q(\eta) (known phase)');
hold(ax2, 'off'); grid(ax2, 'on'); ylim(ax2, [0.5/n 1]); xlim(ax2, [0 6]);
xlabel(ax2, 'threshold \eta'); ylabel(ax2, 'P(statistic > \eta)');
legend(ax2, 'Location', 'southwest'); title(ax2, 'False-alarm probability');

ax3 = nexttile(tl);
lag = (1:maxLag).';
plot(ax3, lag, ac(lag + 1), 'k.-', 'DisplayName', 'measured'); hold(ax3, 'on');
plot(ax3, lag, acModel(lag + 1), 'r-', 'LineWidth', 1.5, 'DisplayName', 'model (weightX)');
sigAc = 1/sqrt(nBin*nProf);
yline(ax3, 2*sigAc, 'b--', 'HandleVisibility', 'off'); yline(ax3, -2*sigAc, 'b--', 'HandleVisibility', 'off');
hold(ax3, 'off'); grid(ax3, 'on');
xlabel(ax3, 'lag [phase bins]'); ylabel(ax3, 'correlation of normProfile');
legend(ax3, 'Location', 'northeast');
title(ax3, sprintf('Bin-bin correlation (lag 0: %.3f); implied T0 std %.3f', ac(1), predRatio));
