function [fig, check] = plotDetectedPower(info_det, info_gen, tStart, tSpan, opts)
%PLOTDETECTEDPOWER  Check the square-law detector output against ground truth.
%{
Reads the whole detectPower output (it is small) and shows:

  1. overview of the full file (channel sum, display-averaged), with the true
     pulse centres and the not-fully-supported region (grey)
  2. (only if nChan > 1) channel-vs-time waterfall of the zoom window;
     dedispersed pulses must be vertical in every channel
  3. zoom window: measured power vs the expected power (pulsar + receiver
     noise baseline)
  4. normalized residual (measured - expected) / sigma_expected. With
     receiver noise this covers every bin; without it only the on-pulse bins
     (off-pulse the expected variance is zero). It should be unit-variance
     noise without structure. RFI is not in the model, so it shows up here.

Expected power and its variance come from expectedPowerModel (ground truth:
generator, both band tapers, IQ gain, receiver-noise level). Pass the
upstream info structs 'InfoDisp', 'InfoIQ', 'InfoDedisp' and, with noise,
'InfoRx' (from addNoiseAndRFI).

Printed, per fully supported pulse in the whole file: power-centroid offset
from the expected centroid (after subtracting the off-pulse baseline), its
predicted uncertainty, and measured/expected pulse energy; then summaries
and the normalized-residual statistics on-pulse and off-pulse. The centroid
is a quick check; it becomes very noisy below a per-pulse SNR of ~10, where
the fold and TOA checks are the meaningful ones.

  [fig, check] = plotDetectedPower(info_det, info_gen, 0, 20e-3, ...
      'InfoDisp', info_disp, 'InfoIQ', info_IQ, 'InfoDedisp', info_dedisp, ...
      'InfoRx', info_rx);

Option 'DisplayAverage': time bins averaged for DISPLAY only (default []:
automatic, so the pulse stands out of the noise where possible). All
statistics use the full resolution.
%}

arguments
    info_det  struct
    info_gen  struct
    tStart    (1,1) double {mustBeNonnegative} = 0
    tSpan     (1,1) double {mustBePositive}    = 20e-3
    opts.InfoDisp   struct = struct([])
    opts.InfoIQ     struct = struct([])
    opts.InfoDedisp struct = struct([])
    opts.InfoRx     struct = struct([])
    opts.SubSamples (1,1) double {mustBeInteger, mustBePositive} = 16
    opts.OnPulseFrac (1,1) double {mustBePositive} = 0.01
    opts.DisplayAverage double = []
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
inSup = tb >= tSup(1) & tb <= tSup(2);

if strcmpi(info_gen.envelopeMode, 'power')
    sigP = info_gen.sigma;
else
    sigP = info_gen.sigma / sqrt(2);
end
tcAll = info_gen.pulseCenters;

% ---- Expected power ---------------------------------------------------------------------
haveExp = ~isempty(opts.InfoDisp) && ~isempty(opts.InfoIQ) && ~isempty(opts.InfoDedisp);
if haveExp
    M  = expectedPowerModel(info_gen, opts.InfoDisp, opts.InfoIQ, opts.InfoDedisp, opts.InfoRx);
    off = ((1:opts.SubSamples) - 0.5)/opts.SubSamples - 0.5;   % boxcar sub-grid
    pS = zeros(size(tb));
    for s = off
        pS = pS + M.envelope(tb + s*dt);
    end
    pS   = pS / numel(off);
    Ps   = M.sigScale * pS;                       % pulsar part
    Pexp = Ps + M.Pn;                             % + receiver-noise baseline
    sdB  = sqrt(M.var(Ps, dt));                   % per-bin std
    offP = pS < 1e-6;                             % truly off-pulse bins
else
    M = struct('Pn', 0, 'hasNoise', false, 'hasRFI', false);
    Ps = []; Pexp = []; sdB = [];
    Tp   = info_gen.T;                            % distance to the nearest pulse
    dist = abs(mod(tb - tcAll(1) + Tp/2, Tp) - Tp/2);
    offP = dist > 6*sigP;
end
use0 = offP & inSup;
if any(use0), base = median(P(use0)); else, base = 0; end   % measured baseline

% ---- Per-pulse check over the whole file ------------------------------------------------
halfW = 4 * sigP;
check = struct('tc', {}, 'offset', {}, 'offsetNoise', {}, 'energyRatio', {});
for tc = tcAll
    if tc - halfW < tSup(1) || tc + halfW > tSup(2), continue; end
    sel = tb >= tc - halfW & tb <= tc + halfW;
    Q   = P(sel) - base;
    cm  = sum(tb(sel) .* Q) / sum(Q);
    if haveExp
        w   = Ps(sel);
        ce  = sum(tb(sel) .* w) / sum(w);
        er  = sum(Q) / sum(w);
        sOf = sqrt(sum(((tb(sel) - ce) .* sdB(sel)).^2)) / sum(w);
    else
        ce = tc; er = NaN; sOf = NaN;
    end
    check(end+1) = struct('tc', tc, 'offset', cm - ce, ...
        'offsetNoise', sOf, 'energyRatio', er); %#ok<AGROW>
end

fprintf('plotDetectedPower: %d bins x %d chan @ %.6g Hz (%.3g us bins)\n', ...
    N, nChan, 1/dt, dt*1e6);
if haveExp
    fprintf('  expected: pulse peak %.4g, noise baseline %.4g, measured off-pulse baseline %.4g', ...
        M.sigScale, M.Pn, base);
    if M.hasRFI, fprintf('  (RFI present: not in the model)'); end
    fprintf('\n');
end
if isempty(check)
    fprintf('  no fully supported pulse in the file.\n');
else
    if numel(check) <= 20
        for c = check
            fprintf('  pulse %9.4f ms: offset %+9.3f us (noise ~%.3f us), energy ratio %.4f\n', ...
                c.tc*1e3, c.offset*1e6, c.offsetNoise*1e6, c.energyRatio);
        end
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

% ---- Normalized residuals ------------------------------------------------------------------
zAll = [];
if haveExp
    zAll = (P - Pexp) ./ sdB;
    zAll(sdB == 0) = NaN;
    onP = Ps > opts.OnPulseFrac * M.sigScale;
    statLine('on-pulse ', zAll(onP & inSup));
    if M.hasNoise
        statLine('off-pulse', zAll(offP & inSup));
    end
end

% ---- Display averaging ---------------------------------------------------------------------
nd = opts.DisplayAverage;
if isempty(nd)
    nd = 1;
    if haveExp && M.hasNoise
        need = (3 * sqrt(M.cnn/dt) * M.Pn / M.sigScale)^2;   % bins for 3:1 per display bin
        nd = max(1, min(round(sigP/(8*dt)), ceil(need)));
    end
end
avg = @(v) mean(reshape(v(1:floor(numel(v)/nd)*nd), nd, []), 1);
tD  = avg(tb);
PD  = avg(P);
if haveExp
    PexpD = avg(Pexp);
    zD = (PD - PexpD) ./ (avg(sdB) / sqrt(nd));
    zD(avg(sdB) == 0) = NaN;
    if ~M.hasNoise
        zD(avg(Ps) <= opts.OnPulseFrac * M.sigScale) = NaN;
    end
end

% ---- Plot ---------------------------------------------------------------------------------------
nRows = 3 + (nChan > 1);
fig = figure('Name', 'Detected power check', 'Color', 'w');
tl  = tiledlayout(fig, nRows, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, sprintf('Detected power: %d channel(s), %.3g \\mus bins (display: %.3g \\mus)', ...
    nChan, dt*1e6, nd*dt*1e6));
tMsD = tD * 1e3;
winD = tD >= tStart & tD <= tStart + tSpan;
tW   = tMsD(winD);
yLo  = min(PD); yHi = max(PD);
if haveExp, yLo = min(yLo, min(PexpD)); yHi = max(yHi, max(PexpD)); end
pad  = 0.05 * (yHi - yLo + eps);
yl   = [yLo - pad, yHi + pad];

% 1. Overview
ax1 = nexttile(tl);
hold(ax1, 'on');
shadeOutside(ax1, [tMsD(1) tMsD(end)], tSup*1e3, yl);
plot(ax1, tMsD, PD, 'k');
if haveExp, plot(ax1, tMsD, PexpD, 'r'); end
plot(ax1, tcAll*1e3, yl(2)*ones(size(tcAll)) - pad, 'rv', 'MarkerFaceColor', 'r', 'MarkerSize', 4);
xline(ax1, [tStart, tStart + tSpan]*1e3, 'b-', 'HandleVisibility', 'off');
hold(ax1, 'off'); box(ax1, 'on'); grid(ax1, 'on');
xlim(ax1, [tMsD(1) tMsD(end)]); ylim(ax1, yl);
ylabel(ax1, 'power'); xlabel(ax1, 'Time [ms]');
leg = {'not fully supported', 'measured'};
if haveExp, leg{end+1} = 'expected'; end
leg{end+1} = 'true centres';
legend(ax1, leg, 'Location', 'northeast');
title(ax1, 'Whole file (blue lines: zoom window)');

% 2. Waterfall (sub-bands only)
if nChan > 1
    ax2 = nexttile(tl);
    XD = zeros(nChan, numel(tD));
    for c = 1:nChan, XD(c, :) = avg(X(c, :)); end
    imagesc(ax2, tW, info_det.chanFreqs/1e9, XD(:, winD)); axis(ax2, 'xy');
    hold(ax2, 'on');
    for tc = tcAll(tcAll >= tStart & tcAll <= tStart + tSpan)
        xline(ax2, tc*1e3, 'w--');
    end
    hold(ax2, 'off');
    ylabel(ax2, 'Channel centre [GHz]');
    cb = colorbar(ax2); cb.Label.String = 'power';
    title(ax2, 'Channels vs time (dedispersed pulses should be vertical)');
end

% 3. Zoom
ax3 = nexttile(tl);
hold(ax3, 'on');
shadeOutside(ax3, [tW(1) tW(end)], tSup*1e3, yl);
plot(ax3, tW, PD(winD), 'k.-', 'MarkerSize', 6);
leg = {'not fully supported', 'measured'};
if haveExp
    plot(ax3, tW, PexpD(winD), 'r', 'LineWidth', 1.2);
    leg{end+1} = 'expected';
end
for tc = tcAll(tcAll >= tStart & tcAll <= tStart + tSpan)
    xline(ax3, tc*1e3, 'r:', 'HandleVisibility', 'off');
end
hold(ax3, 'off'); box(ax3, 'on'); grid(ax3, 'on'); ylim(ax3, yl);
ylabel(ax3, 'power'); legend(ax3, leg, 'Location', 'northeast');
title(ax3, 'Zoom: measured vs ground-truth expectation');

% 4. Normalized residual
ax4 = nexttile(tl);
if haveExp
    hR = plot(ax4, tW, zD(winD), 'k'); hold(ax4, 'on');
    yline(ax4, 0, 'r-', 'HandleVisibility', 'off');
    h1 = yline(ax4, [-1 1], 'b-'); h3 = yline(ax4, [-3 3], 'b--');
    hold(ax4, 'off'); ylim(ax4, [-5 5]); grid(ax4, 'on');
    legend(ax4, [hR h1(1) h3(1)], {'(measured - expected)/\sigma', '\pm1\sigma', '\pm3\sigma'}, ...
        'Location', 'northeast');
    if M.hasNoise
        title(ax4, 'Normalized residual, all bins (unit-variance noise expected; RFI shows up here)');
    else
        title(ax4, 'Normalized residual on the pulses (unit-variance noise expected)');
    end
    ylabel(ax4, '\sigma');
else
    text(ax4, 0.5, 0.5, 'Pass InfoDisp, InfoIQ, InfoDedisp (and InfoRx) for the expected power', ...
        'HorizontalAlignment', 'center', 'Units', 'normalized');
    axis(ax4, 'off');
end
xlabel(ax4, 'Time [ms]');

zoomAxes = [ax3, ax4];
if nChan > 1, zoomAxes = [ax2, zoomAxes]; end
linkaxes(zoomAxes, 'x');
xlim(ax3, [tW(1) tW(end)]);

check = struct('pulses', check, 'baseline', base, 'displayAverage', nd, 'zAll', zAll);
end


% =========================================================================================
function statLine(name, z)
z = z(isfinite(z));
n = numel(z);
if n == 0, return; end
fprintf(['  normalized residual %s over %d bins: mean %+.4f, std %.4f ' ...
         '(expect 0 +- %.4f, 1 +- %.4f)\n'], name, n, mean(z), std(z), ...
    1/sqrt(n), 1/sqrt(2*n));
end


function shadeOutside(ax, xr, sup, yl)
c = [0.88 0.88 0.88];
if sup(1) > xr(1)
    patch(ax, [xr(1) sup(1) sup(1) xr(1)], [yl(1) yl(1) yl(2) yl(2)], c, 'EdgeColor', 'none');
else
    patch(ax, nan(1,4), nan(1,4), c, 'EdgeColor', 'none');      % legend entry
end
if sup(2) < xr(2)
    patch(ax, [sup(2) xr(2) xr(2) sup(2)], [yl(1) yl(1) yl(2) yl(2)], c, ...
        'EdgeColor', 'none', 'HandleVisibility', 'off');
end
end