function [fig, check] = plotFoldCheck(fold, info_fold, info_gen, opts)
%PLOTFOLDCHECK  Check folded profiles against the generator ground truth.
%{
Panels:
  1. total folded profile (channel sum) over the full turn, phase -0.5..0.5,
     with the expected profile
  2. zoom around phase 0: measured vs expected
  3. normalized residual (measured - expected) / sigma on the pulse
     (sigma from self-noise: per time bin 1/sqrt(binDt*Bnoise), combined
     with the fold weights) - should be unit-variance noise
  4. sub-integration stack (zoom): the pulse must sit at phase 0 in every row
  5. per-sub-integration power-centroid offset [us] with +-1 sigma; a
     quick preview of the TOAs (the proper estimator is the next module)

The expected profile assumes the fold phase model matches the generator
(phase 0 at a true pulse centre). It includes the averaging of the
detector time bin (boxcar, binDt) and of the phase assignment (uniform
over one bin for 'nearest', triangular over +-1 bin for 'linear').
Absolute level needs 'InfoDisp', 'InfoIQ', 'InfoDedisp'; without them the
expected shape is scaled to the measured profile by least squares and no
noise-normalized statistics are available.

  [fig, check] = plotFoldCheck(fold, info_fold, info_gen, ...
      'InfoDisp', info_disp, 'InfoIQ', info_IQ, 'InfoDedisp', info_dedisp);
%}

arguments
    fold      struct
    info_fold struct
    info_gen  struct
    opts.InfoDisp    struct = struct([])
    opts.InfoIQ      struct = struct([])
    opts.InfoDedisp  struct = struct([])
    opts.ZoomPhase   (1,1) double {mustBePositive} = 0.1
    opts.OnPulseFrac (1,1) double {mustBePositive} = 0.01
    opts.SubGrid     (1,1) double {mustBeInteger, mustBePositive} = 16
end

Nbin = info_fold.NBin;
f0   = info_fold.f0;
P    = 1 / f0;
nSub = info_fold.nSub;

% Phase axis centred on 0: bins reordered to phases -0.5 .. 0.5
sh    = floor(Nbin/2);
phC   = ((0:Nbin-1) - sh) / Nbin;                        % centred phases
reord = @(v) circshift(v, sh, 1);

profTot = reord(sum(fold.profTotal, 2));                  % channel sum, [Nbin x 1]
profSub = reord(squeeze(sum(fold.prof, 3)));              % [Nbin x nSub]
if nSub == 1, profSub = profSub(:); end
wSub  = reord(fold.weight);
w2Sub = reord(fold.weight2);
if isfield(fold, 'weightX')
    wxSub = reord(fold.weightX);
else
    wxSub = zeros(size(wSub));                        % older fold files
end
wTot  = sum(wSub, 2);
w2Tot = sum(w2Sub, 2);
wxTot = sum(wxSub, 2);

% ---- Expected profile ------------------------------------------------------------------
haveAbs = ~isempty(opts.InfoDisp) && ~isempty(opts.InfoIQ) && ~isempty(opts.InfoDedisp);
shape = expectedShape(phC, info_fold, info_gen, opts.SubGrid);   % unit-scale p(phase)
if haveAbs
    [scale, Bnoise] = powerScale(info_gen, opts.InfoDisp, opts.InfoIQ, opts.InfoDedisp);
    Pexp = scale * shape;
    relStd = 1 / sqrt(info_fold.binDt * Bnoise);          % per time bin
    sigTot = relStd * Pexp .* sqrt(w2Tot) ./ wTot;        % per phase bin, total fold
else
    ok = ~isnan(profTot);
    Pexp = shape * (shape(ok)' * profTot(ok)) / (shape(ok)' * shape(ok));
    relStd = NaN; sigTot = nan(size(Pexp));
end

% ---- Checks -------------------------------------------------------------------------------------
if strcmpi(info_gen.envelopeMode, 'power')
    sigPh = info_gen.sigma * f0;
else
    sigPh = info_gen.sigma * f0 / sqrt(2);
end
win = abs(phC(:)) <= 4*sigPh;

check = struct();
cm = sum(phC(win)' .* profTot(win)) / sum(profTot(win));
ce = sum(phC(win)' .* Pexp(win))    / sum(Pexp(win));
check.totalOffset = (cm - ce) * P;
check.totalEnergyRatio = sum(profTot(win)) / sum(Pexp(win));
if haveAbs
    aC = zeros(Nbin, 1);
    aC(win) = (phC(win)' - ce) / sum(Pexp(win));
    check.totalOffsetNoise = combNoise(aC, relStd * Pexp, wTot, w2Tot, wxTot) * P;
    onP = Pexp > opts.OnPulseFrac * max(Pexp);
    zr  = (profTot - Pexp) ./ sigTot;
    check.normResidMean = mean(zr(onP), 'omitnan');
    check.normResidStd  = std(zr(onP), 'omitnan');
    check.nOnPulse      = nnz(onP);
else
    check.totalOffsetNoise = NaN; zr = nan(size(Pexp)); onP = false(size(Pexp));
    check.normResidMean = NaN; check.normResidStd = NaN; check.nOnPulse = 0;
end

subOff = nan(nSub, 1); subNoise = nan(nSub, 1);
for s = 1:nSub
    p = profSub(:, s);
    if any(isnan(p(win))), continue; end
    subOff(s) = (sum(phC(win)' .* p(win)) / sum(p(win)) - ce) * P;
    if haveAbs
        aC = zeros(Nbin, 1);
        aC(win) = (phC(win)' - ce) / sum(Pexp(win));
        subNoise(s) = combNoise(aC, relStd * Pexp, wSub(:, s), w2Sub(:, s), wxSub(:, s)) * P;
    end
end
check.subOffset = subOff;
check.subOffsetNoise = subNoise;
check.subTime = fold.subint.tMean;

fprintf('plotFoldCheck: NBin = %d (%.3g us), %d sub-int(s), %s assignment\n', ...
    Nbin, P/Nbin*1e6, nSub, info_fold.assign);
fprintf('  total profile: centroid offset %+.4f us (noise ~%.4f us), energy ratio %.5f\n', ...
    check.totalOffset*1e6, check.totalOffsetNoise*1e6, check.totalEnergyRatio);
if haveAbs
    fprintf(['  normalized residual over %d on-pulse bins: mean %+.3f, std %.3f ' ...
             '(expect 0 +- %.3f, 1 +- %.3f)\n'], check.nOnPulse, check.normResidMean, ...
        check.normResidStd, 1/sqrt(check.nOnPulse), 1/sqrt(2*check.nOnPulse));
end
okS = ~isnan(subOff);
fprintf('  sub-int centroid offsets: mean %+.3f us, rms %.3f us', ...
    mean(subOff(okS))*1e6, sqrt(mean(subOff(okS).^2))*1e6);
if haveAbs
    fprintf(' (expected rms ~%.3f us)\n', sqrt(mean(subNoise(okS).^2))*1e6);
else
    fprintf('\n');
end

% ---- Plot -------------------------------------------------------------------------------------
fig = figure('Name', 'Fold check', 'Color', 'w');
tl  = tiledlayout(fig, 3, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, sprintf('Folded profile: P = %.6g ms, %d bins (%.3g \\mus), %d sub-int(s)', ...
    P*1e3, Nbin, P/Nbin*1e6, nSub));
zp = opts.ZoomPhase;
zSel = abs(phC) <= zp;

ax1 = nexttile(tl, 1);
plot(ax1, phC, profTot, 'k', phC, Pexp, 'r'); grid(ax1, 'on');
xlim(ax1, [-0.5 0.5]); xlabel(ax1, 'Phase [turns]'); ylabel(ax1, 'power');
legend(ax1, {'folded', 'expected'}, 'Location', 'northeast');
title(ax1, 'Total folded profile');

ax2 = nexttile(tl, 2);
plot(ax2, phC(zSel), profTot(zSel), 'k.-', phC(zSel), Pexp(zSel), 'r', 'MarkerSize', 6);
xline(ax2, 0, 'r:'); grid(ax2, 'on'); xlim(ax2, [-zp zp]);
xlabel(ax2, 'Phase [turns]'); ylabel(ax2, 'power');
title(ax2, 'Zoom around phase 0');

ax3 = nexttile(tl, 3);
if haveAbs
    zz = zr; zz(~onP) = NaN;
    hR = plot(ax3, phC, zz, 'k.-', 'MarkerSize', 6); hold(ax3, 'on');
    yline(ax3, 0, 'r-', 'HandleVisibility', 'off');
    h1 = yline(ax3, [-1 1], 'b-'); h3 = yline(ax3, [-3 3], 'b--');
    hold(ax3, 'off'); ylim(ax3, [-5 5]);
    legend(ax3, [hR h1(1) h3(1)], {'(folded - expected)/\sigma', '\pm1\sigma', '\pm3\sigma'}, ...
        'Location', 'northeast');
    title(ax3, 'Normalized residual on the pulse');
else
    plot(ax3, phC, profTot - Pexp, 'k');
    title(ax3, 'Residual (shape only: pass upstream info structs for \sigma)');
end
grid(ax3, 'on'); xlim(ax3, [-zp zp]); xlabel(ax3, 'Phase [turns]');

ax4 = nexttile(tl, 4);
imagesc(ax4, phC(zSel), 1:nSub, profSub(zSel, :).'); axis(ax4, 'xy');
xline(ax4, 0, 'w--'); xlabel(ax4, 'Phase [turns]'); ylabel(ax4, 'Sub-integration');
cb = colorbar(ax4); cb.Label.String = 'power';
title(ax4, 'Sub-integrations (pulse at phase 0 in every row)');

ax5 = nexttile(tl, 5, [1 2]);
tS = fold.subint.tMean * 1e3;
if haveAbs
    errorbar(ax5, tS(okS), subOff(okS)*1e6, subNoise(okS)*1e6, 'ko', 'MarkerFaceColor', 'k');
else
    plot(ax5, tS(okS), subOff(okS)*1e6, 'ko', 'MarkerFaceColor', 'k');
end
yline(ax5, 0, 'r-'); grid(ax5, 'on');
xlabel(ax5, 'Sub-integration mean time [ms]'); ylabel(ax5, 'offset [\mus]');
title(ax5, 'Per-sub-integration power-centroid offset (TOA preview, \pm1\sigma self-noise)');
end


% ========================================================================================
function shape = expectedShape(phC, info_fold, info_gen, nS)
% Unit-amplitude expected profile at the centred phases phC, averaged over
% the detector time bin and the phase-assignment kernel.
f0 = info_fold.f0; Nbin = info_fold.NBin; dt = info_fold.binDt;
if strcmp(info_fold.assign, 'nearest')
    u = ((1:nS) - 0.5)/nS - 0.5;  wu = ones(size(u));
else
    u = ((1:2*nS) - 0.5)/nS - 1;  wu = 1 - abs(u);
end
wu = wu / sum(wu);
v = ((1:8) - 0.5)/8 - 0.5;                                % time-bin boxcar
tc0 = info_fold.TRef - info_fold.Phi0 / f0;               % time of phase 0
shape = zeros(numel(phC), 1);
for iu = 1:numel(u)
    for iv = 1:numel(v)
        t = tc0 + (phC(:) + u(iu)/Nbin) / f0 + v(iv)*dt;
        shape = shape + wu(iu) * envelopePower(t, info_gen) / numel(v);
    end
end
end


function sd = combNoise(a, sTB, W, W2, WX)
%COMBNOISE  Std of sum_j a_j * prof_j, including the covariance of
% neighbouring phase bins that share time bins (linear assignment).
% sTB = per-time-bin noise std at each phase bin.
v  = sTB.^2 .* W2 ./ W.^2;
cv = sTB .* circshift(sTB, -1) .* WX ./ (W .* circshift(W, -1));
v(~isfinite(v)) = 0; cv(~isfinite(cv)) = 0;
sd = sqrt(sum(a.^2 .* v) + 2*sum(a .* circshift(a, -1) .* cv));
end


function [scale, Bnoise] = powerScale(info_gen, info_disp, info_IQ, info_dedisp)
fLow = info_dedisp.fLow; fHigh = info_dedisp.fHigh;
f = linspace(fLow, fHigh, 200001);
W = taperW(f, info_disp.fLow, info_disp.fHigh, info_disp.edgeWidth) .* ...
    taperW(f, fLow, fHigh, info_dedisp.edgeWidth);
Beff   = trapz(f, W.^2);
Bnoise = Beff^2 / trapz(f, W.^4);
scale  = info_IQ.gainFactor^2 * info_gen.A^2 * Beff / info_IQ.fsIn;
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