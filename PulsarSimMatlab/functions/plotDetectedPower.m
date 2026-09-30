function [fig, check] = plotDetectedPower(info_det, info_gen, tStart, tSpan, opts)
%PLOTDETECTEDPOWER  Check the square-law detector output against ground truth.
%{
Reads the whole detectPower output (it is small) and shows:

  1. overview of the full file (channel sum), with true pulse centres and
     the not-fully-supported region (grey)
  2. (only if nChan > 1) channel-vs-time waterfall of the zoom window;
     dedispersed pulses must be vertical in every channel
  3. zoom window: measured power vs the expected profile
  4. normalized residual (measured - expected) / sigma_expected on the pulses,
     which should be unit-variance noise without structure. The pulsar
     signal is itself noise, so sigma scales with the expected power
     ("self-noise"); dividing by it makes the check visible at any level.

The expected profile needs the upstream info structs ('InfoDisp',
'InfoIQ', 'InfoDedisp'). Without them, panels 3-4 show the measurement
only and the per-pulse check compares against the true centres.

It prints, for every fully supported pulse in the FILE (not just the zoom
window), the power-centroid offset from the expected centroid, its
noise-limited uncertainty, and the measured/expected pulse energy, plus a
summary. This is effectively a first TOA-vs-ground-truth test.

  plotDetectedPower(info_det, info_gen)                               % 0-20 ms
  [fig, check] = plotDetectedPower(info_det, info_gen, 0, 20e-3, ...
      'InfoDisp', info_disp, 'InfoIQ', info_IQ, 'InfoDedisp', info_dedisp);

Expected power per bin (no receiver noise):
  gain^2 * A^2 * <p(t)>_bin * Beff / fsIn, with p = G ('power' envelope)
  or G^2 ('amplitude'), Beff = integral of (W_fwd*W_inv)^2 over the band,
  and <.>_bin the average over the bin (boxcar), evaluated on a sub-grid.
%}

arguments
    info_det  struct
    info_gen  struct
    tStart    (1,1) double {mustBeNonnegative} = 0
    tSpan     (1,1) double {mustBePositive}    = 20e-3
    opts.InfoDisp   struct = struct([])
    opts.InfoIQ     struct = struct([])
    opts.InfoDedisp struct = struct([])
    opts.SubSamples (1,1) double {mustBeInteger, mustBePositive} = 16
    opts.OnPulseFrac (1,1) double {mustBePositive} = 0.01
end

% ---- Read ----------------------------------------------------------------------------
[fid, msg] = fopen(info_det.file, 'r', info_det.byteOrder);
if fid == -1
    error('plotDetectedPower:open', 'Could not open "%s": %s', info_det.file, msg);
end
X = fread(fid, [info_det.nChan, info_det.N], 'single=>double');
fclose(fid);
nChan = info_det.nChan;
N  = size(X, 2);
P  = sum(X, 1);                                   % channel sum
dt = info_det.binDt;
tb = info_det.binTime0 + (0:N-1) * dt;            % bin centroid times
sup = info_det.fullySupportedBins;
tSup = [tb(sup(1)) - dt/2, tb(sup(2)) + dt/2];

% ---- Expected profile --------------------------------------------------------------------
haveExp = ~isempty(opts.InfoDisp) && ~isempty(opts.InfoIQ) && ~isempty(opts.InfoDedisp);
if haveExp
    [Pexp, Beff, Bnoise] = expectedBinPower(tb, dt, opts.SubSamples, info_gen, ...
        opts.InfoDisp, opts.InfoIQ, opts.InfoDedisp);
    relStd = 1 / sqrt(dt * Bnoise);               % per-bin relative noise
else
    Pexp = []; Beff = NaN; Bnoise = NaN; relStd = NaN;
end

% ---- Per-pulse check over the whole file --------------------------------------------------------
if strcmpi(info_gen.envelopeMode, 'power')
    sigP = info_gen.sigma;
else
    sigP = info_gen.sigma / sqrt(2);
end
halfW = 4 * sigP;
tcAll = info_gen.pulseCenters;
check = struct('tc', {}, 'offset', {}, 'offsetNoise', {}, 'energyRatio', {});
for tc = tcAll
    if tc - halfW < tSup(1) || tc + halfW > tSup(2)
        continue
    end
    sel = tb >= tc - halfW & tb <= tc + halfW;
    cm = sum(tb(sel) .* P(sel)) / sum(P(sel));
    if haveExp
        w   = Pexp(sel);
        ce  = sum(tb(sel) .* w) / sum(w);
        er  = sum(P(sel)) / sum(w);
        sOf = sqrt(sum(((tb(sel) - ce) .* w * relStd).^2)) / sum(w);
    else
        ce = tc; er = NaN; sOf = NaN;
    end
    check(end+1) = struct('tc', tc, 'offset', cm - ce, ...
        'offsetNoise', sOf, 'energyRatio', er); %#ok<AGROW>
end

fprintf('plotDetectedPower: %d bins x %d chan @ %.6g Hz (%.3g us bins)\n', ...
    N, nChan, 1/dt, dt*1e6);
if isempty(check)
    fprintf('  no fully supported pulse in the file.\n');
else
    for c = check
        fprintf('  pulse %9.4f ms: offset %+9.3f us (noise ~%.3f us), energy ratio %.4f\n', ...
            c.tc*1e3, c.offset*1e6, c.offsetNoise*1e6, c.energyRatio);
    end
    off = [check.offset];
    fprintf('  summary over %d pulses: mean offset %+.3f us, rms %.3f us', ...
        numel(off), mean(off)*1e6, sqrt(mean(off.^2))*1e6);
    if haveExp
        fprintf(' (expected rms ~%.3f us), mean energy ratio %.4f\n', ...
            sqrt(mean([check.offsetNoise].^2))*1e6, mean([check.energyRatio]));
    else
        fprintf('\n');
    end
end

% ---- Normalized residual statistics (whole file, supported on-pulse bins) -----
if haveExp
    onP  = Pexp > opts.OnPulseFrac * max(Pexp);
    zAll = nan(size(P));
    zAll(onP) = (P(onP) - Pexp(onP)) ./ (relStd * Pexp(onP));
    use  = onP & tb >= tSup(1) & tb <= tSup(2);
    zU   = zAll(use);
    fprintf(['  normalized residual over %d on-pulse bins: mean %+.4f, std %.4f ' ...
             '(expect 0 +- %.4f, 1 +- %.4f)\n'], numel(zU), mean(zU), std(zU), ...
        1/sqrt(numel(zU)), 1/sqrt(2*numel(zU)));
else
    zAll = [];
end

% ---- Plot ---------------------------------------------------------------------------------------------
nRows = 3 + (nChan > 1);
fig = figure('Name', 'Detected power check', 'Color', 'w');
tl  = tiledlayout(fig, nRows, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, sprintf('Detected power: %d channel(s), %.3g \\mus bins, f_{out} = %.4g Hz', ...
    nChan, dt*1e6, 1/dt));
tMsAll = tb * 1e3;
winSel = tb >= tStart & tb <= tStart + tSpan;
tW = tMsAll(winSel);
yTop = 1.05 * max([P, Pexp]);

% 1. Overview
ax1 = nexttile(tl);
hold(ax1, 'on');
shadeOutside(ax1, [tMsAll(1) tMsAll(end)], tSup*1e3, yTop);
plot(ax1, tMsAll, P, 'k');
plot(ax1, tcAll*1e3, yTop*0.98*ones(size(tcAll)), 'rv', 'MarkerFaceColor', 'r', ...
    'MarkerSize', 4);
xline(ax1, [tStart, tStart + tSpan]*1e3, 'b-', 'HandleVisibility', 'off');
hold(ax1, 'off'); box(ax1, 'on'); grid(ax1, 'on');
xlim(ax1, [tMsAll(1) tMsAll(end)]); ylim(ax1, [0 yTop]);
ylabel(ax1, 'power'); xlabel(ax1, 'Time [ms]');
legend(ax1, {'not fully supported', 'measured', 'true centres'}, 'Location', 'northeast');
title(ax1, 'Whole file (blue lines: zoom window)');

% 2. Waterfall (sub-bands only)
if nChan > 1
    ax2 = nexttile(tl);
    imagesc(ax2, tW, info_det.chanFreqs/1e9, X(:, winSel)); axis(ax2, 'xy');
    hold(ax2, 'on');
    tcW = tcAll(tcAll >= tStart & tcAll <= tStart + tSpan);
    for tc = tcW
        xline(ax2, tc*1e3, 'w--');
    end
    hold(ax2, 'off');
    ylabel(ax2, 'Channel centre [GHz]');
    cb = colorbar(ax2); cb.Label.String = 'power';
    title(ax2, 'Channels vs time (dedispersed pulses should be vertical)');
end

% 3. Zoom: measured vs expected
ax3 = nexttile(tl);
hold(ax3, 'on');
shadeOutside(ax3, [tW(1) tW(end)], tSup*1e3, yTop);
plot(ax3, tW, P(winSel), 'k.-', 'MarkerSize', 6);
leg = {'not fully supported', 'measured'};
if haveExp
    plot(ax3, tW, Pexp(winSel), 'r', 'LineWidth', 1.2);
    leg{end+1} = 'expected';
end
tcW = tcAll(tcAll >= tStart & tcAll <= tStart + tSpan);
for tc = tcW
    xline(ax3, tc*1e3, 'r:', 'HandleVisibility', 'off');
end
hold(ax3, 'off'); box(ax3, 'on'); grid(ax3, 'on'); ylim(ax3, [0 yTop]);
ylabel(ax3, 'power'); legend(ax3, leg, 'Location', 'northeast');
title(ax3, 'Zoom: measured vs ground-truth expectation');

% 4. Normalized residual
ax4 = nexttile(tl);
if haveExp
    zW = zAll(winSel);
    hold(ax4, 'on');
    hR = plot(ax4, tW, zW, 'k');
    yline(ax4, 0, 'r-', 'HandleVisibility', 'off');
    h1 = yline(ax4, [-1 1], 'b-', 'LineWidth', 1);
    h3 = yline(ax4, [-3 3], 'b--', 'LineWidth', 1);
    hold(ax4, 'off'); box(ax4, 'on'); grid(ax4, 'on');
    ylim(ax4, [-5 5]);
    legend(ax4, [hR, h1(1), h3(1)], ...
        {'(measured - expected) / \sigma', '\pm1\sigma', '\pm3\sigma'}, ...
        'Location', 'northeast');
    title(ax4, sprintf(['Normalized residual on the pulses (bins with expected power ' ...
        '> %g%% of peak): should be unit-variance noise'], 100*opts.OnPulseFrac));
    ylabel(ax4, '\sigma');
else
    text(ax4, 0.5, 0.5, 'Pass InfoDisp, InfoIQ, InfoDedisp for the expected profile', ...
        'HorizontalAlignment', 'center', 'Units', 'normalized');
    axis(ax4, 'off');
end
xlabel(ax4, 'Time [ms]');

zoomAxes = [ax3, ax4];
if nChan > 1, zoomAxes = [ax2, zoomAxes]; end
linkaxes(zoomAxes, 'x');
xlim(ax3, [tW(1) tW(end)]);
end


% =========================================================================================
function [Pexp, Beff, Bnoise] = expectedBinPower(tb, dt, nSub, info_gen, info_disp, info_IQ, info_dedisp)
fLow = info_dedisp.fLow; fHigh = info_dedisp.fHigh;
f = linspace(fLow, fHigh, 200001);
W = taperW(f, info_disp.fLow, info_disp.fHigh, info_disp.edgeWidth) .* ...
    taperW(f, fLow, fHigh, info_dedisp.edgeWidth);
Beff = trapz(f, W.^2);                            % sets the mean power
Bnoise = Beff^2 / trapz(f, W.^4);                 % sets the per-bin variance

% Boxcar average of the envelope over each bin on a sub-grid
off = ((1:nSub) - 0.5) / nSub - 0.5;              % symmetric offsets in bins
p = zeros(size(tb));
for s = off
    p = p + envelopePower(tb + s*dt, info_gen);
end
p = p / nSub;
Pexp = info_IQ.gainFactor^2 * info_gen.A^2 * p * Beff / info_IQ.fsIn;
end


function p = envelopePower(t, info_gen)
T = info_gen.T; sig = info_gen.sigma;
nN = max(1, ceil(6*sig/T + 0.5));
peakNorm = sum(exp(-0.5*((-nN:nN)*T/sig).^2));
k0 = round(t/T - 0.5);
G = zeros(size(t));
for j = -nN:nN
    G = G + exp(-0.5*((t - (k0 + j + 0.5)*T)/sig).^2);
end
G = G / peakNorm;
if strcmpi(info_gen.envelopeMode, 'power')
    p = G;
else
    p = G.^2;
end
end


function W = taperW(f, fLow, fHigh, e)
W = zeros(size(f));
ib = f >= fLow & f <= fHigh;
W(ib) = 1;
lo = ib & f < fLow + e;   W(lo) = sin(pi/2 * (f(lo) - fLow) / e).^2;
hi = ib & f > fHigh - e;  W(hi) = sin(pi/2 * (fHigh - f(hi)) / e).^2;
end


function shadeOutside(ax, xr, sup, yTop)
c = [0.88 0.88 0.88];
if sup(1) > xr(1)
    patch(ax, [xr(1) sup(1) sup(1) xr(1)], [0 0 yTop yTop], c, 'EdgeColor', 'none');
else
    patch(ax, nan(1,4), nan(1,4), c, 'EdgeColor', 'none');      % legend entry
end
if sup(2) < xr(2)
    patch(ax, [sup(2) xr(2) xr(2) sup(2)], [0 0 yTop yTop], c, ...
        'EdgeColor', 'none', 'HandleVisibility', 'off');
end
end