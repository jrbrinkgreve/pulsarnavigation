%{
RUNSNRSWEEP - TOA precision, error bars and detection versus SNR (phases D, E)

For each SNR in snrList, nReal independent realizations go through the
same stages as main.m; the pooled TOAs are compared with the ground truth
(validateTOA) and with the best achievable TOA error (expectedPowerModel,
exact band tapers). Questions:
  - low SNR: below which SNR per sub-int do the error bars stop being
    valid (wrong correlation peak, outliers)?
  - high SNR: how far does FFTFIT fall behind the optimum once the
    pulsar's own noise (self-noise) dominates?
  - detection (detectPulsar): false-alarm rate under H0, detection
    probability vs SNR against theory, and whether the good TOAs
    (detected, chi^2 not flagged; as in main.m) are free of outliers.

H0 (no pulsar at all): per realization one extra pass on a zero sky (the
generator with A = 0, made once) plus receiver noise ('NoiseStd' 1; the
level does not matter, the detection statistics are normalized), with the
realization's noise seed. Only detectPulsar runs there (no TOAs to
validate). A very low SNR such as -50 dB is NOT used as H0: it is a real
operating point (small antennas), only needing long integrations.

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
Cost at L = 0.1 s: about nReal * (7 s + (numel(snrList) + 1) * 8 s).
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
swZero     = fullfile(mcDir, "sweep_zero.dat");  % H0 sky: no pulsar
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
detC = cell(nReal, nSNR);                        % detectPulsar result per (realization, SNR)
detH0C = cell(nReal, 1);                         % detectPulsar result of the H0 pass
genC = cell(nReal, 1);                           % info_gen per realization
optSigma = nan(1, nSNR);                         % best TOA error per sub-int [s]
optSNR   = nan(1, nSNR);                         % best SNR per sub-int

fprintf('runSNRSweep: %d realizations x %d SNR values (%s dB), L = %.3g s, %d turn(s) per sub-int\n', ...
    nReal, nSNR, num2str(snrList), L, subintPeriods);
tAll = tic;
% H0 sky: the generator with A = 0 (zeros, same length; dispersion of zeros is zeros)
info_zero = generatePulsarSignal(swZero, T, f_in, 0, L, dutycycle, ...
    'Seed', seed, 'EnvelopeMode', genEnvMode, 'Verbose', false);
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
        detection = detectPulsar(fold, info_fold, template, 'Bnoise', Bnoise, 'Verbose', false);
        toaC{j, s} = toa;
        detC{j, s} = detection;

        % Best achievable (ground truth, exact band tapers); same for every j
        if j == 1
            M = expectedPowerModel(info_gen, info_disp, info_IQ, info_dedisp, info_rx);
            optSigma(s) = M.toaErrPulse / sqrt(subintPeriods);
            optSNR(s)   = M.snrPulse * sqrt(subintPeriods);
        end

        v = toa.valid;
        fprintf(['    %+4g dB: %d TOAs, median SNR %.1f, median TOA error %.4g us, ' ...
                 'detected %d (unknown phase), chi2-flagged %d, %.0f s\n'], ...
            snrList(s), nnz(v), median(toa.snr(v)), median(toa.toaErr(v))*1e6, ...
            nnz(detection.detectedUnknown), nnz(toa.flagChi2), toc(tRun));
    end

    % H0 pass: receiver noise only (no pulsar), same noise seed, same processing
    tRun = tic;
    info_rx = addNoiseAndRFI(info_zero.file, swRx, info_zero.actualFsOut, ...
        'Band', [fLow fHigh], 'NoiseStd', 1, 'RFI', rfi, 'Seed', noiseSeedJ, 'Verbose', false);
    info_IQ = applyIQmodulation(info_rx.file, swIQ, info_rx.actualFsOut, fs, fLO, ...
        'FilterOrder', filterOrder, 'Band', [fLow fHigh], 'Verbose', false);
    info_dedisp = applyInverseDispersion(info_IQ.file, swDedisp, info_IQ.actualFsOut, ...
        info_IQ.fLO, ephem.DM, fLow, fHigh, 'RefFreq', refFreq, 'Verbose', false);
    info_det = detectPower(info_dedisp.file, swEnvelope, info_dedisp.fs, info_dedisp.fLO, ...
        f_out, 'FullySupported', info_dedisp.fullySupported, 'Verbose', false);
    [info_fold, fold] = foldProfile(info_det, swFold, ephem.f0, ...
        'F1', ephem.F1, 'TRef', ephem.TRef, 'NBin', nBin, 'SubintPeriods', subintPeriods, ...
        'SaveFile', false, 'Verbose', false);
    Bnoise = noiseBandwidth(fLow, fHigh, info_dedisp.edgeWidth);
    detection = detectPulsar(fold, info_fold, template, 'Bnoise', Bnoise, 'Verbose', false);
    detH0C{j} = detection;
    fprintf('       H0: %d sub-ints, T0 mean %+.2f, false alarms %d known / %d unknown phase, %.0f s\n', ...
        nnz(detection.tested), mean(detection.T0(detection.tested)), ...
        nnz(detection.detectedKnown), nnz(detection.detectedUnknown), toc(tRun));
end
fprintf('runSNRSweep: done in %.1f min\n', toc(tAll)/60);

% Pooled validation per SNR point
nanS = nan(1, nSNR);
res = struct('snrDB', snrList, 'optSNR', optSNR, 'optSigma', optSigma, ...
    'n', zeros(1, nSNR), 'rmsErr', nanS, 'robustErr', nanS, ...
    'rmsPred', nanS, 'ratio', nanS, 'redChi2', nanS, ...
    'frac1', nanS, 'fracOut', nanS, 'meanErr', nanS, ...
    'nTested', zeros(1, nSNR), 'T0mean', nanS, 'T0std', nanS, ...
    'pdKnown', nanS, 'pdUnknown', nanS, 'nFlagChi2', zeros(1, nSNR), ...
    'nGood', zeros(1, nSNR), 'fracOutGood', nanS, 'ratioGood', nanS, ...
    'etaKnown', detC{1, 1}.etaKnown, 'etaUnknown', detC{1, 1}.etaUnknown, ...
    'pfa', 1e-3);                                        % detectPulsar default
genAll = [genC{:}];
for s = 1:nSNR
    toaS = [toaC{:, s}];
    detS = [detC{:, s}];

    % Detection: pooled over realizations (tested = complete turns)
    tested = vertcat(detS.tested);
    T0     = vertcat(detS.T0);
    dK     = vertcat(detS.detectedKnown);
    dU     = vertcat(detS.detectedUnknown);
    res.nTested(s)   = nnz(tested);
    res.T0mean(s)    = mean(T0(tested));
    res.T0std(s)     = std(T0(tested));
    res.pdKnown(s)   = mean(dK(tested));
    res.pdUnknown(s) = mean(dU(tested));
    res.nFlagChi2(s) = nnz(vertcat(toaS.valid) & vertcat(toaS.flagChi2));

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

    % Good TOAs only (as in main.m): detected (unknown phase) and chi^2 not flagged
    toaG = toaS;
    for i = 1:nReal
        toaG(i).valid = toaS(i).valid & detS(i).detectedUnknown & ~toaS(i).flagChi2;
    end
    res.nGood(s) = nnz(vertcat(toaG.valid));
    if res.nGood(s) > 0
        fprintf('  good TOAs only:\n');
        valG = validateTOA(toaG, genAll, 'Plot', false);
        res.fracOutGood(s) = mean(abs(valG.errNorm) > 5);
        res.ratioGood(s)   = valG.ratio;
    end
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

Q = @(x) 0.5 * erfc(x / sqrt(2));                       % 1 - Phi(x)
fprintf(['\nDetection (P_FA %.3g per sub-int; thresholds %.2f known phase, %.2f unknown ' ...
         'phase); theory P_D = 1 - Phi(eta - SNRopt), unknown phase approximate\n'], ...
    res.pfa, res.etaKnown, res.etaUnknown);
fprintf('  snrDB  SNRopt  tested  T0 mean  T0 std   P_D known (th.)   P_D unknown (th.)  flagged  good  |z|>5 good  rms/pred good\n');
for s = 1:nSNR
    fprintf('  %+5g %7.2f %7d %8.2f %7.2f   %6.3f (%6.3f)   %6.3f (%6.3f)   %6d %5d %9.1f%% %12.3f\n', ...
        snrList(s), optSNR(s), res.nTested(s), res.T0mean(s), res.T0std(s), ...
        res.pdKnown(s), Q(res.etaKnown - optSNR(s)), res.pdUnknown(s), ...
        Q(res.etaUnknown - optSNR(s)), res.nFlagChi2(s), res.nGood(s), ...
        100*res.fracOutGood(s), res.ratioGood(s));
end
% H0 (no pulsar): pooled over the H0 passes
detH0 = [detH0C{:}];
tH0 = vertcat(detH0.tested);
T0h = vertcat(detH0.T0);  T0h = T0h(tH0);
Tmh = vertcat(detH0.Tmax); Tmh = Tmh(tH0);
nH0 = nnz(tH0);
res.h0 = struct('n', nH0, 'T0mean', mean(T0h), 'T0std', std(T0h), ...
    'faKnown', nnz(T0h > res.etaKnown), 'faUnknown', nnz(Tmh > res.etaUnknown), ...
    'Tmax', Tmh);
fprintf(['  H0 (no pulsar, %d sub-ints): T0 mean %+.3f, std %.3f (expect 0 +- %.3f, ' ...
         '1 +- %.3f); false alarms %d known, %d unknown phase (expect %.2g each)\n'], ...
    nH0, res.h0.T0mean, res.h0.T0std, 1/sqrt(nH0), 1/sqrt(2*nH0), ...
    res.h0.faKnown, res.h0.faUnknown, res.pfa*nH0);

save(swResult, 'res', 'toaC', 'detC', 'detH0C', 'genAll', 'snrList', 'nReal', 'subintPeriods', 'L');
fprintf('runSNRSweep: results saved to %s\n', swResult);
plotSweep(res);


% ======================================================================
%  Local functions
%  ======================================================================
function plotSweep(res)
%PLOTSWEEP  TOA error vs SNR, error-bar honesty, FFTFIT efficiency, threshold,
% detection (with the H0 false alarms in its title).
hasDet = isfield(res, 'pdKnown');                      % older results: no detection
ok  = res.n > 0;
x   = res.snrDB;
dr  = 1 ./ sqrt(2 * max(res.n, 1));                    % 1-sigma of an rms ratio
fig = figure('Name', 'SNR sweep', 'Color', 'w'); %#ok<NASGU>
tl  = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, 'Phases D/E: TOA error and detection versus SNR');

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
if hasDet
    semilogx(ax3, res.optSNR(ok), res.fracOutGood(ok), 'bd--', ...
        'DisplayName', '|z| > 5 among good TOAs');
end
yline(ax3, 0.683, 'k--', 'HandleVisibility', 'off'); hold(ax3, 'off');
ylim(ax3, [0 1]); grid(ax3, 'on');
xlabel(ax3, 'best SNR per sub-int'); ylabel(ax3, 'fraction of TOAs');
legend(ax3, 'Location', 'east'); title(ax3, 'Threshold: where error bars fail');

if ~hasDet, return; end
ax4 = nexttile(tl);
Q   = @(z) 0.5 * erfc(z / sqrt(2));                    % 1 - Phi(z)
sg  = logspace(-1, log10(max(res.optSNR)), 300);       % theory curves
nT  = max(res.nTested, 1);
okD = res.nTested > 0;
eK  = sqrt(res.pdKnown .* (1 - res.pdKnown) ./ nT);    % binomial 1-sigma
eU  = sqrt(res.pdUnknown .* (1 - res.pdUnknown) ./ nT);
semilogx(ax4, sg, Q(res.etaKnown - sg), 'k-', 'DisplayName', 'theory, known phase'); hold(ax4, 'on');
semilogx(ax4, sg, Q(res.etaUnknown - sg), 'b-', 'DisplayName', 'theory, unknown phase (approx.)');
errorbar(ax4, res.optSNR(okD), res.pdKnown(okD), eK(okD), 'ko', 'MarkerFaceColor', 'k', ...
    'DisplayName', 'measured, known phase');
errorbar(ax4, res.optSNR(okD), res.pdUnknown(okD), eU(okD), 'bs', 'MarkerFaceColor', 'b', ...
    'DisplayName', 'measured, unknown phase');
hold(ax4, 'off'); set(ax4, 'XScale', 'log'); ylim(ax4, [0 1.02]); grid(ax4, 'on');
xlabel(ax4, 'best SNR per sub-int'); ylabel(ax4, 'detection probability P_D');
legend(ax4, 'Location', 'southeast');
tit = sprintf('Detection (P_{FA} %.3g, \\eta = %.2f / %.2f)', res.pfa, res.etaKnown, res.etaUnknown);
if isfield(res, 'h0')                                  % false alarms without a pulsar
    tit = sprintf('%s; H0: %d/%d, %d/%d false alarms', tit, ...
        res.h0.faKnown, res.h0.n, res.h0.faUnknown, res.h0.n);
end
title(ax4, tit);
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
