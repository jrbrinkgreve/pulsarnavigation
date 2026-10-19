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
same processing as main.m (processIQ, IQ -> fold), then detectPulsar. Front
end from pipelineParams (default: channels with RFI excision, 'optimal'
weighting; 'fullband' the reference path). With L = 0.1 s and 1 turn per
sub-int: 9 complete sub-ints per pass on the full band (about 8 s per pass),
8 on the channel path (about 40 s per pass).

Model correlation of normProfile (the noise-only profile in units of its
model sigma) from the fold weights: full band WX/sqrt(W2*W2') at lag 1; with
channels the combination detectPulsar uses ('optimal': weights 1/variance per
channel; 'equal': the sum of the channel profiles), lags 1..3, assuming equal
noise levels in all channels (white receiver noise). With excision also the
blanked fraction per pass: on pure noise, the false-flag rate of detectRFI.

Parameters from pipelineParams.m; files in data/mc with an h0_ prefix (channel
files in data/mc/chan). (19 Oct 2026: on processIQ, both front ends; before:
full band only.)
%}

nPass  = 100;                                    % noise-only passes (~9 profiles each)
pfaList = [1e-1 1e-2 1e-3];                      % false-alarm levels to test
maxLag = 200;                                    % [bins] lags shown in the figure

% Paths and parameters
scriptDir = fileparts(mfilename('fullpath'));
addpath(fullfile(scriptDir, 'functions'));
pipelineParams;

mcDir      = fullfile(dataDir, "mc");
mcChan     = fullfile(mcDir, "chan");
h0Zero     = fullfile(mcDir, "h0_zero.dat");
h0Rx       = fullfile(mcDir, "h0_rx.dat");
h0IQ       = fullfile(mcDir, "h0_rx_IQ.dat");
h0Result   = fullfile(mcDir, "h0Run.mat");
% processIQ's output files (own names, so main's files in data/ stay)
h0Files = struct('dedisp', fullfile(mcDir, "h0_dedisp.dat"), ...
    'envelope', fullfile(mcDir, "h0_envelope.dat"), 'fold', fullfile(mcDir, "h0_fold.mat"), ...
    'chanIQ', fullfile(mcChan, "h0_rx_IQ_chan"), 'chanBlank', fullfile(mcChan, "h0_rx_IQ_chan_blanked"), ...
    'chanDedisp', fullfile(mcChan, "h0_dedisp_chan"), 'chanPower', fullfile(mcChan, "h0_dedisp_chan_power.dat"), ...
    'chanDedispB', fullfile(mcChan, "h0_dedispB_chan"), 'chanPowerB', fullfile(mcChan, "h0_dedispB_chan_power.dat"), ...
    'chanWeight', fullfile(mcChan, "h0_dedispB_chan_weights.dat"), 'foldChan', fullfile(mcDir, "h0_fold_chan.mat"));
if ~isfolder(mcChan), mkdir(mcChan); end

fullBand   = strcmp(frontEnd, 'fullband');
excisionOn = ~fullBand && excision;
if fullBand && ~isempty(rfi)
    warning('runH0:noExcision', ['RFI is on (rfiOn) but the full-band path has no excision: ' ...
        'the noise-only passes contain the RFI unremoved.']);
end

template = gaussianTemplate(nBin, ephem.profileFWHM);
c = template / max(template); c = c - mean(c);  % zero-mean template, as in detectPulsar

pathStr = 'full band';
if ~fullBand, pathStr = sprintf('channels, excision %d, weighting ''%s''', excisionOn, weighting); end
fprintf('runH0: %d noise-only passes, L = %.3g s, %d turn(s) per sub-int, %s\n', nPass, L, ...
    subintPeriods, pathStr);
tAll = tic;
info_zero = generatePulsarSignal(h0Zero, T, f_in, 0, L, dutycycle, ...
    'Seed', seed, 'EnvelopeMode', genEnvMode, 'Verbose', false);

T0all = []; Tmall = []; nrAll = [];
acSum = zeros(nBin, 1);                          % sum of circular autocorrelations of normProfile
rhoSum = 0; nProf = 0;                           % model correlation at lags 1..nLag (weightX)
blanked = nan(nPass, 1);                         % excision: blanked fraction per pass
for i = 1:nPass
    tRun = tic;
    ns = seed + 10000 + i;
    info_rx = addNoiseAndRFI(info_zero.file, h0Rx, info_zero.actualFsOut, ...
        'Band', [fLow fHigh], 'NoiseStd', 1, 'RFI', rfi, 'Seed', ns, 'Verbose', false);
    info_IQ = applyIQmodulation(info_rx.file, h0IQ, info_rx.actualFsOut, fs, fLO, ...
        'FilterOrder', filterOrder, 'Band', [fLow fHigh], 'Verbose', false);
    % Processing (ephemeris + receiver settings only), as in main.m; detection only
    out  = processIQ(info_IQ, rx, ephem, proc, h0Files, 'Verbose', false);
    fold = out.fold;
    [detection, info_detect] = detectPulsar(fold, out.info_fold, template, 'Bnoise', out.Bnoise, ...
        'Weighting', out.weighting, 'Verbose', false);
    if excisionOn, blanked(i) = mean(out.info_w.blankedFraction); end

    t = detection.tested;
    T0all = [T0all; detection.T0(t)];            %#ok<AGROW>
    Tmall = [Tmall; detection.Tmax(t)];          %#ok<AGROW>
    nrAll = [nrAll; detection.noiseRatio(t)];    %#ok<AGROW>
    for s = find(t).'
        z = detection.normProfile(:, s);
        acSum = acSum + real(ifft(abs(fft(z)).^2)) / nBin;   % circular autocorrelation
        if fullBand                              % one channel, lag 1 (as before 19 Oct)
            W  = fold.weight(:, s);  W2 = fold.weight2(:, s);  WX = fold.weightX(:, s);
            rhoSum = rhoSum + mean(WX ./ sqrt(W2 .* circshift(W2, -1)));
        else
            rhoSum = rhoSum + modelCorr(fold, s, out.weighting);
        end
        nProf = nProf + 1;
    end
    blankStr = '';
    if excisionOn, blankStr = sprintf(', blanked %.4f %%', 100*blanked(i)); end
    fprintf('  pass %3d/%d (seed %d): %d profiles, T0 mean %+.2f, max Tmax %.2f%s, %.0f s\n', ...
        i, nPass, ns, nnz(t), mean(detection.T0(t)), max(detection.Tmax(t)), blankStr, toc(tRun));
end
fprintf('runH0: done in %.1f min\n', toc(tAll)/60);

% ---- Statistics -------------------------------------------------------------------------
n  = numel(T0all);
ac = acSum / nProf;                              % measured autocorrelation, lag 0..nBin-1
rho  = rhoSum / nProf;                           % model correlation, lags 1..numel(rho)
rho1 = rho(1);
acModel = zeros(nBin, 1); acModel(1) = 1;
acModel(1 + (1:numel(rho))) = rho; acModel(end + 1 - (1:numel(rho))) = rho;
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
if numel(rho) > 1
    fprintf('  normProfile correlation lags 2-%d: measured %s, model %s\n', numel(rho), ...
        mat2str(round(ac(3:numel(rho)+1).', 4)), mat2str(round(rho(2:end), 4)));
end
if excisionOn
    fprintf('  blanked fraction (detectRFI false flags on pure noise): mean %.4f %% (per pass %.4f-%.4f %%)\n', ...
        100*mean(blanked), 100*min(blanked), 100*max(blanked));
end
fprintf('  T0 std implied by the measured correlation: %.4f (measured %.4f)\n', predRatio, std(T0all));
fprintf('  exceedances      P_FA   expected   T0 > eta_known (eta)   Tmax > eta_unknown (eta)\n');
for k = 1:numel(pfaList)
    fprintf('               %8.3g   %6.1f +- %4.1f   %5d (%.2f)            %5d (%.2f)\n', ...
        pfaList(k), n*pfaList(k), sqrt(n*pfaList(k)), nnz(T0all > etaK(k)), etaK(k), ...
        nnz(Tmall > etaU(k)), etaU(k));
end

save(h0Result, 'T0all', 'Tmall', 'nrAll', 'ac', 'acModel', 'rho1', 'rho', 'Rc', 'predRatio', ...
    'pfaList', 'etaK', 'etaU', 'rf', 'nPass', 'L', 'subintPeriods', 'frontEnd', 'blanked');
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


% ======================================================================
%  Local functions
% ======================================================================
function rho = modelCorr(fold, s, weighting)
%MODELCORR  Model correlation of detectPulsar's normProfile at lags 1..nLag for sub-int s of
% a channel fold, from the fold weights W, W2, WX (per channel, or shared), assuming equal
% noise levels in all channels (white receiver noise; the level cancels).
%   'optimal': z_j = sum_c w_cj (p_cj - pbar_c) / sqrt(sum_c w_cj), w = 1/var ~ W^2/W2:
%     rho_d(j) = sum_c W_cj W_cj+d WX_cjd / (W2_cj W2_cj+d) / sqrt(sum_c W_cj^2/W2_cj * (same at j+d))
%   'equal': z_j = (sum_c p_cj - a) / sqrt(v_j), var_cj ~ W2/W^2, cov_cjd ~ WX/(W W'):
%     rho_d(j) = sum_c WX_cjd/(W_cj W_cj+d) / sqrt(sum_c W2_cj/W_cj^2 * (same at j+d))
% One channel: both reduce to WX/sqrt(W2*W2'). Averaged over the bins j.
N    = size(fold.weight, 1);
W    = reshape(fold.weight(:, s, :), N, []);
W2   = reshape(fold.weight2(:, s, :), N, []);
nW   = size(W, 2);
nLag = size(fold.weightX, 4);
WX   = reshape(fold.weightX(:, s, :, :), N, nW, nLag);
rho  = zeros(1, nLag);
for d = 1:nLag
    Wd = circshift(W, -d, 1); W2d = circshift(W2, -d, 1);
    if strcmp(weighting, 'optimal')
        num = W .* Wd .* WX(:, :, d) ./ (W2 .* W2d);
        a = W.^2 ./ W2; b = Wd.^2 ./ W2d;
    else
        num = WX(:, :, d) ./ (W .* Wd);
        a = W2 ./ W.^2; b = W2d ./ Wd.^2;
    end
    num(~isfinite(num)) = 0; a(~isfinite(a)) = 0; b(~isfinite(b)) = 0;
    r = sum(num, 2) ./ sqrt(sum(a, 2) .* sum(b, 2));
    rho(d) = mean(r(isfinite(r)));
end
end
