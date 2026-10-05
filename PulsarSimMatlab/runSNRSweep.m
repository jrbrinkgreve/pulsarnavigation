%{
RUNSNRSWEEP - TOA precision and error bars versus SNR (phase D)

For each SNR in snrList, nReal independent realizations go through the
same stages as main.m; the pooled TOAs are compared with the ground truth
(validateTOA) and with the best achievable TOA error (expectedPowerModel,
exact band tapers). Two questions:
  - low SNR: below which SNR per sub-int do the error bars stop being
    valid (wrong correlation peak, outliers)?
  - high SNR: how far does FFTFIT fall behind the optimum once the
    pulsar's own noise (self-noise) dominates?

Unlike runMonteCarlo.m, the PULSAR realization changes too: at high SNR the
self-noise would otherwise be identical in every realization and missing
from the scatter. Per realization j the sky (generator + dispersion) is
made once and reused for all SNR values (it does not depend on SNR):
  pulsar seed = seed + 100*(j-1),  noise seed = pulsar seed + 1
so j = 1 equals main.m (42 / 43): at -5 dB its TOAs must match main.
The same noise seed is used for every SNR of a realization (the noise is
only rescaled), which makes the curves smoother; TOAs of one SNR point
stay independent across realizations.

Parameters come from pipelineParams.m (snrDB and seed are overridden
here). Files go to data/mc with a sweep_ prefix and are overwritten.
Cost at L = 0.1 s: about nReal * (15 s + numel(snrList) * 17 s).
%}

nReal   = 10;                                    % pulsar + noise realizations
snrList = [-25 -20 -15 -10 -5 0 10 20];          % [dB] S_peak/SEFD
replotOnly = false;    % true: only redraw the figure from data/mc/snrSweep.mat


% Paths and parameters
scriptDir = fileparts(mfilename('fullpath'));
addpath(fullfile(scriptDir, 'functions'));
pipelineParams;

mcDir      = fullfile(dataDir, "mc");
swRaw      = fullfile(mcDir, "sweep_raw.dat");
swDisp     = fullfile(mcDir, "sweep_dispersed.dat");
swRx       = fullfile(mcDir, "sweep_rx.dat");
swIQ       = fullfile(mcDir, "sweep_rx_IQ.dat");
swDedisp   = fullfile(mcDir, "sweep_dedisp.dat");
swEnvelope = fullfile(mcDir, "sweep_envelope.dat");
swFold     = fullfile(mcDir, "sweep_fold.mat");  % not saved (SaveFile false)
swResult   = fullfile(mcDir, "snrSweep.mat");
if replotOnly
    S = load(swResult, 'res');
    plotSweep(S.res);
    return
end

template = gaussianTemplate(nBin, ephem.profileFWHM);
nSNR = numel(snrList);
toaC = cell(nReal, nSNR);                        % toa struct per (realization, SNR)
genC = cell(nReal, 1);                           % info_gen per realization
optSigma = nan(1, nSNR);                         % best TOA error per sub-int [s]
optSNR   = nan(1, nSNR);                         % best SNR per sub-int

fprintf('runSNRSweep: %d realizations x %d SNR values (%s dB), L = %.3g s, %d turn(s) per sub-int\n', ...
    nReal, nSNR, num2str(snrList), L, subintPeriods);
tAll = tic;
for j = 1:nReal
    pulsarSeed = seed + 100*(j - 1);             % j = 1 is main's seed
    noiseSeedJ = pulsarSeed + 1;                 % j = 1 is main's noiseSeed

    % Sky (synthetic): new pulsar realization, dispersion
    tSky = tic;
    info_gen = generatePulsarSignal(swRaw, T, f_in, A, L, dutycycle, ...
        'Seed', pulsarSeed, 'EnvelopeMode', genEnvMode, 'Verbose', false);
    info_disp = applyDispersionStream(info_gen.file, swDisp, DM, ...
        info_gen.actualFsOut, fLow, fHigh, 'RefFreq', refFreq, 'Verbose', false);
    genC{j} = info_gen;
    fprintf('  realization %d/%d: pulsar seed %d, noise seed %d, sky %.0f s\n', ...
        j, nReal, pulsarSeed, noiseSeedJ, toc(tSky));

    for s = 1:nSNR
        tRun = tic;
        % Receiver (synthetic): noise at RF, IQ conversion
        info_rx = addNoiseAndRFI(info_disp.file, swRx, info_disp.actualFsOut, ...
            'Band', [fLow fHigh], 'SNRdB', snrList(s), 'SignalInfo', info_gen, ...
            'RFI', rfi, 'Seed', noiseSeedJ, 'Verbose', false);
        info_IQ = applyIQmodulation(info_rx.file, swIQ, info_rx.actualFsOut, fs, fLO, ...
            'FilterOrder', filterOrder, 'Band', [fLow fHigh], 'Verbose', false);

        % Processing (ephemeris + receiver settings only), same calls as main.m
        info_dedisp = applyInverseDispersion(info_IQ.file, swDedisp, info_IQ.actualFsOut, ...
            info_IQ.fLO, ephem.DM, fLow, fHigh, 'RefFreq', refFreq, 'Verbose', false);
        info_det = detectPower(info_dedisp.file, swEnvelope, info_dedisp.fs, info_dedisp.fLO, ...
            f_out, 'FullySupported', info_dedisp.fullySupported, 'Verbose', false);
        [info_fold, fold] = foldProfile(info_det, swFold, ephem.f0, ...
            'F1', ephem.F1, 'TRef', ephem.TRef, 'NBin', nBin, 'SubintPeriods', subintPeriods, ...
            'SaveFile', false, 'Verbose', false);
        Bnoise = noiseBandwidth(fLow, fHigh, info_dedisp.edgeWidth);
        toa = estimateTOA(fold, info_fold, template, 'Bnoise', Bnoise, 'Verbose', false);
        toaC{j, s} = toa;

        % Best achievable (ground truth, exact band tapers); same for every j
        if j == 1
            M = expectedPowerModel(info_gen, info_disp, info_IQ, info_dedisp, info_rx);
            optSigma(s) = M.toaErrPulse / sqrt(subintPeriods);
            optSNR(s)   = M.snrPulse * sqrt(subintPeriods);
        end

        v = toa.valid;
        fprintf('    %+4g dB: %d TOAs, median SNR %.1f, median TOA error %.4g us, %.0f s\n', ...
            snrList(s), nnz(v), median(toa.snr(v)), median(toa.toaErr(v))*1e6, toc(tRun));
    end
end
fprintf('runSNRSweep: done in %.1f min\n', toc(tAll)/60);

% Pooled validation per SNR point
res = struct('snrDB', snrList, 'optSNR', optSNR, 'optSigma', optSigma, ...
    'n', zeros(1, nSNR), 'rmsErr', nan(1, nSNR), 'robustErr', nan(1, nSNR), ...
    'rmsPred', nan(1, nSNR), 'ratio', nan(1, nSNR), 'redChi2', nan(1, nSNR), ...
    'frac1', nan(1, nSNR), 'fracOut', nan(1, nSNR), 'meanErr', nan(1, nSNR));
genAll = [genC{:}];
for s = 1:nSNR
    toaS = [toaC{:, s}];
    if ~any(vertcat(toaS.valid)), continue; end
    fprintf('\n--- %+g dB ---\n', snrList(s));
    val = validateTOA(toaS, genAll, 'Plot', false);
    res.n(s)         = val.n;
    res.meanErr(s)   = val.meanErr;
    res.rmsErr(s)    = val.rmsErr;
    res.robustErr(s) = 1.4826 * median(abs(val.err));   % rms of the Gaussian core
    res.rmsPred(s)   = val.rmsPred;
    res.ratio(s)     = val.ratio;
    res.redChi2(s)   = val.chi2 / val.dof;
    res.frac1(s)     = val.frac1;
    res.fracOut(s)   = mean(abs(val.errNorm) > 5);
end

% Summary table
fprintf('\nrunSNRSweep summary (%d realizations, %d turn(s) per sub-int; times in us)\n', ...
    nReal, subintPeriods);
fprintf('  snrDB  SNRopt    n      rms   robust     pred      opt  rms/pred  pred/opt  chi2red  |z|<1  |z|>5\n');
for s = 1:nSNR
    fprintf('  %+5g %7.1f %4d %8.4g %8.4g %8.4g %8.4g %9.3f %9.3f %8.2f %5.0f%% %5.1f%%\n', ...
        snrList(s), optSNR(s), res.n(s), res.rmsErr(s)*1e6, res.robustErr(s)*1e6, ...
        res.rmsPred(s)*1e6, optSigma(s)*1e6, res.ratio(s), res.rmsPred(s)/optSigma(s), ...
        res.redChi2(s), 100*res.frac1(s), 100*res.fracOut(s));
end

save(swResult, 'res', 'toaC', 'genAll', 'snrList', 'nReal', 'subintPeriods', 'L');
fprintf('runSNRSweep: results saved to %s\n', swResult);
plotSweep(res);


% ======================================================================
%  Local functions
%  ======================================================================
function plotSweep(res)
%PLOTSWEEP  TOA error vs SNR, error-bar honesty, FFTFIT efficiency, threshold.
ok  = res.n > 0;
x   = res.snrDB;
dr  = 1 ./ sqrt(2 * max(res.n, 1));                    % 1-sigma of an rms ratio
fig = figure('Name', 'SNR sweep', 'Color', 'w'); %#ok<NASGU>
tl  = tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, 'Phase D: TOA error versus SNR');

ax1 = nexttile(tl);
errorbar(ax1, x(ok), res.rmsErr(ok)*1e6, dr(ok).*res.rmsErr(ok)*1e6, 'ko-', ...
    'MarkerFaceColor', 'k', 'DisplayName', 'rms error'); hold(ax1, 'on');
plot(ax1, x(ok), res.robustErr(ok)*1e6, 'bs--', 'DisplayName', 'robust (1.48 MAD)');
plot(ax1, x(ok), res.rmsPred(ok)*1e6, 'm^-', 'DisplayName', 'rms predicted \sigma');
plot(ax1, x, res.optSigma*1e6, 'r-', 'LineWidth', 1.5, 'DisplayName', 'best achievable');
hold(ax1, 'off'); set(ax1, 'YScale', 'log'); grid(ax1, 'on');
xlabel(ax1, 'SNR [dB] (S_{peak}/SEFD)'); ylabel(ax1, '\sigma_{TOA} per sub-int [\mus]');
legend(ax1, 'Location', 'northeast'); title(ax1, 'TOA error');

ax2 = nexttile(tl);
yLim = [0.8 3];                                        % off-scale points: on the edge, labelled
rB = res.ratio;                   rB(~ok) = NaN;       % error-bar honesty
rM = res.rmsPred ./ res.optSigma; rM(~ok) = NaN;       % FFTFIT efficiency
[cB, offB] = clipToAxis(rB, yLim);
[cM, offM] = clipToAxis(rM, yLim);
eB = dr .* rB; eB(offB) = 0;
errorbar(ax2, x, cB, eB, 'ko-', 'MarkerFaceColor', 'k', ...
    'DisplayName', 'rms error / predicted'); hold(ax2, 'on');
plot(ax2, x, cM, 'm^-', 'DisplayName', 'predicted / best achievable');
yline(ax2, 1, 'r-', 'HandleVisibility', 'off');
labelOffScale(ax2, x(offB), rB(offB), cB(offB), 'k');
labelOffScale(ax2, x(offM), rM(offM), cM(offM), 'm');
hold(ax2, 'off'); ylim(ax2, yLim); grid(ax2, 'on');
xlabel(ax2, 'SNR [dB]'); ylabel(ax2, 'ratio');
legend(ax2, 'Location', 'northeast'); title(ax2, 'Ratios');

ax3 = nexttile(tl);
semilogx(ax3, res.optSNR(ok), res.frac1(ok), 'ko-', 'MarkerFaceColor', 'k', ...
    'DisplayName', '|z| < 1'); hold(ax3, 'on');
semilogx(ax3, res.optSNR(ok), res.fracOut(ok), 'rs-', 'DisplayName', '|z| > 5 (outliers)');
yline(ax3, 0.683, 'k--', 'HandleVisibility', 'off'); hold(ax3, 'off');
ylim(ax3, [0 1]); grid(ax3, 'on');
xlabel(ax3, 'best SNR per sub-int'); ylabel(ax3, 'fraction of TOAs');
legend(ax3, 'Location', 'east'); title(ax3, 'Threshold: where error bars fail');
end


function [yc, off] = clipToAxis(y, yLim)
%CLIPTOAXIS  Clip y to [yLim(1), yLim(2)]; off marks clipped points (NaN stays NaN).
yc  = min(max(y, yLim(1)), yLim(2));
yc(isnan(y)) = NaN;
off = y < yLim(1) | y > yLim(2);
end


function labelOffScale(ax, x, y, yc, color)
%LABELOFFSCALE  Write the true value next to points drawn on the axis edge.
for i = 1:numel(x)
    if y(i) > yc(i), va = 'top'; else, va = 'bottom'; end
    text(ax, x(i), yc(i), sprintf(' %.3g', y(i)), 'Color', color, ...
        'VerticalAlignment', va, 'HorizontalAlignment', 'left', 'FontSize', 8);
end
end
