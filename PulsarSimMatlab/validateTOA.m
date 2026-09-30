function [val, fig] = validateTOA(toa, info_gen, opts)
%VALIDATETOA  Compare estimated TOAs with the generator ground truth.
%{
For each valid sub-integration the true arrival time is the generator
pulse centre closest to the TOA (the dedispersion keeps the reference
frequency at its original time, so the ground truth is directly
info_gen.pulseCenters). Errors, normalized errors (error / predicted
uncertainty) and summary statistics are computed; for a correct estimator
with correct uncertainties the normalized errors are N(0,1).

  [val, fig] = validateTOA(toa, info_gen)
  [val, fig] = validateTOA([toa1 toa2 ...], [info_gen1 info_gen2 ...])   % Monte Carlo

Multiple runs: pass struct arrays of equal length (e.g. collected over
seeds); all TOAs are pooled.

Name-value options:
  'Plot'   show the figure (default true)
  'Label'  text for the figure title (default '')

Output val:
  n, err [s], errNorm, predErr [s], t (truth) [s], run index,
  meanErr, meanErrUnc, rmsErr, rmsPred (rms of predicted errors),
  ratio (rmsErr / rmsPred), chi2, dof, pLow / pHigh (chi^2 tail
  probabilities), frac1 / frac2 (fraction with |errNorm| < 1 / 2;
  expect 0.683 / 0.954), and the total-fold offset(s) with uncertainty.
%}

arguments
    toa       struct
    info_gen  struct
    opts.Plot  (1,1) logical = true
    opts.Label {mustBeTextScalar} = ''
end

if numel(toa) ~= numel(info_gen)
    error('validateTOA:runs', 'toa and info_gen must have the same number of runs.');
end

err = []; pred = []; tru = []; runIdx = []; tTOA = []; amp = []; snr = [];
totOff = nan(numel(toa), 1); totErr = nan(numel(toa), 1);
for r = 1:numel(toa)
    v  = toa(r).valid & isfinite(toa(r).toa);
    tt = toa(r).toa(v);
    pc = info_gen(r).pulseCenters(:);
    [~, iN] = min(abs(tt.' - pc), [], 1);                % nearest true centre
    truth = pc(iN(:));
    if any(abs(tt - truth) > 0.25 * info_gen(r).T)
        warning('validateTOA:far', ['Run %d: some TOAs are more than T/4 from a ' ...
            'true pulse centre; check the phase model.'], r);
    end
    err    = [err;    tt - truth];                        %#ok<AGROW>
    pred   = [pred;   toa(r).toaErr(v)];                  %#ok<AGROW>
    tru    = [tru;    truth];                             %#ok<AGROW>
    tTOA   = [tTOA;   tt];                                %#ok<AGROW>
    runIdx = [runIdx; r * ones(nnz(v), 1)];               %#ok<AGROW>
    amp    = [amp;    toa(r).amp(v)];                     %#ok<AGROW>
    snr    = [snr;    toa(r).snr(v)];                     %#ok<AGROW>
    if isfield(toa(r), 'total')
        totOff(r) = toa(r).total.timeOffset;
        totErr(r) = toa(r).total.timeOffsetErr;
    end
end

n = numel(err);
if n == 0
    error('validateTOA:none', 'No valid TOAs to validate.');
end
z = err ./ pred;
chi2 = sum(z.^2);
val = struct();
val.n        = n;
val.t        = tru;
val.run      = runIdx;
val.err      = err;
val.predErr  = pred;
val.errNorm  = z;
val.meanErr  = mean(err);
val.meanErrUnc = sqrt(mean(pred.^2) / n);
val.rmsErr   = sqrt(mean(err.^2));
val.rmsPred  = sqrt(mean(pred.^2));
val.ratio    = val.rmsErr / val.rmsPred;
val.chi2     = chi2;
val.dof      = n;
val.pLow     = gammainc(chi2/2, n/2);                     % P(chi2 <= observed)
val.pHigh    = 1 - val.pLow;                              % P(chi2 >= observed)
val.frac1    = mean(abs(z) < 1);
val.frac2    = mean(abs(z) < 2);
val.totalOffset    = totOff;
val.totalOffsetErr = totErr;

fprintf('validateTOA: %d TOAs from %d run(s)\n', n, numel(toa));
fprintf('  mean error %+.4f us (+- %.4f us expected), rms %.4f us, predicted rms %.4f us, ratio %.3f\n', ...
    val.meanErr*1e6, val.meanErrUnc*1e6, val.rmsErr*1e6, val.rmsPred*1e6, val.ratio);
fprintf('  chi2 = %.1f for %d dof (red. %.3f), P(<=) = %.3f, P(>=) = %.3f\n', ...
    chi2, n, chi2/n, val.pLow, val.pHigh);
fprintf('  |z| < 1: %.1f%% (expect 68.3%%), |z| < 2: %.1f%% (expect 95.4%%)\n', ...
    100*val.frac1, 100*val.frac2);
for r = 1:numel(toa)
    if isfinite(totOff(r))
        fprintf('  run %d total fold offset %+.4f +- %.4f us (%.2f sigma)\n', ...
            r, totOff(r)*1e6, totErr(r)*1e6, totOff(r)/totErr(r));
    end
end

fig = [];
if ~opts.Plot, return; end

% ---- Figure --------------------------------------------------------------------------------
fig = figure('Name', 'TOA validation', 'Color', 'w');
tl  = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, sprintf('%s  %d TOAs: mean %+.3f us, rms %.3f us (pred. %.3f), \\chi^2_{red} = %.2f', ...
    char(opts.Label), n, val.meanErr*1e6, val.rmsErr*1e6, val.rmsPred*1e6, chi2/n));

ax1 = nexttile(tl);
errorbar(ax1, tTOA*1e3, err*1e6, pred*1e6, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 4);
yline(ax1, 0, 'r-'); grid(ax1, 'on');
xlabel(ax1, 'TOA [ms]'); ylabel(ax1, 'TOA - truth [\mus]');
title(ax1, 'TOA errors (\pm1\sigma predicted)');

ax2 = nexttile(tl);
nb = max(8, min(40, round(sqrt(n)*2)));
histogram(ax2, z, nb, 'Normalization', 'pdf', 'FaceColor', [0.6 0.6 0.6]); hold(ax2, 'on');
zz = linspace(-4, 4, 400);
plot(ax2, zz, exp(-zz.^2/2)/sqrt(2*pi), 'r', 'LineWidth', 1.5); hold(ax2, 'off');
grid(ax2, 'on'); xlabel(ax2, 'error / predicted \sigma'); ylabel(ax2, 'pdf');
title(ax2, sprintf('Normalized errors vs N(0,1): %.0f%% within 1\\sigma', 100*val.frac1));

ax3 = nexttile(tl);
zs = sort(z);
q  = sqrt(2) * erfinv(2*((1:n).' - 0.5)/n - 1);
plot(ax3, q, zs, 'ko', 'MarkerSize', 4); hold(ax3, 'on');
lim = max(3, max(abs([q; zs])));
plot(ax3, [-lim lim], [-lim lim], 'r-'); hold(ax3, 'off');
axis(ax3, 'equal'); xlim(ax3, [-lim lim]); ylim(ax3, [-lim lim]); grid(ax3, 'on');
xlabel(ax3, 'N(0,1) quantile'); ylabel(ax3, 'normalized error');
title(ax3, 'Q-Q plot (points on the line = correct error model)');

ax4 = nexttile(tl);
yyaxis(ax4, 'left');  plot(ax4, tTOA*1e3, amp, 'o-', 'MarkerSize', 4); ylabel(ax4, 'amplitude b');
yyaxis(ax4, 'right'); plot(ax4, tTOA*1e3, snr, 's-', 'MarkerSize', 4); ylabel(ax4, 'SNR (b / \sigma_b)');
grid(ax4, 'on'); xlabel(ax4, 'TOA [ms]');
title(ax4, 'Fitted amplitude and SNR per sub-integration');
end