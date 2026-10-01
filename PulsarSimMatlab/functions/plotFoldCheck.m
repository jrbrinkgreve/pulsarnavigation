function [fig, check] = plotFoldCheck(fold, info_fold, info_gen, opts)
%PLOTFOLDCHECK  Check folded profiles against the generator ground truth.
%{
Panels:
  1. total folded profile (channel sum) over the full turn, phase -0.5..0.5,
     with the expected profile (pulsar + receiver-noise baseline)
  2. zoom around phase 0: measured vs expected
  3. normalized residual (measured - expected) / sigma; all bins when
     receiver noise is present, on-pulse only without it. Should be
     unit-variance noise; RFI (not in the model) shows up here
  4. sub-integration stack (baseline subtracted, zoom): pulse at phase 0
  5. per-sub-integration power-centroid offset [us] with +-1 sigma (a quick
     TOA preview; meaningless below ~10 SNR per sub-integration)

The expected profile assumes the fold phase model matches the generator
(phase 0 at a true pulse centre). It includes the averaging of the detector
time bin (boxcar) and of the phase assignment ('nearest': uniform over one
bin, 'linear': triangular over +-1 bin). Level and noise come from
expectedPowerModel: pass 'InfoDisp', 'InfoIQ', 'InfoDedisp' and, with
noise, 'InfoRx'. The noise of each phase bin combines the per-time-bin
variance with the fold weights, including the neighbour covariance of
linear assignment. Centroids use the profile minus its measured off-pulse
baseline (median of |phase| > 0.25).

  [fig, check] = plotFoldCheck(fold, info_fold, info_gen, ...
      'InfoDisp', info_disp, 'InfoIQ', info_IQ, 'InfoDedisp', info_dedisp, 'InfoRx', info_rx);
%}

arguments
    fold      struct
    info_fold struct
    info_gen  struct
    opts.InfoDisp    struct = struct([])
    opts.InfoIQ      struct = struct([])
    opts.InfoDedisp  struct = struct([])
    opts.InfoRx      struct = struct([])
    opts.ZoomPhase   (1,1) double {mustBePositive} = 0.1
    opts.OnPulseFrac (1,1) double {mustBePositive} = 0.01
    opts.SubGrid     (1,1) double {mustBeInteger, mustBePositive} = 16
end

Nbin = info_fold.NBin;
f0   = info_fold.f0;
Per  = 1 / f0;
nSub = info_fold.nSub;
dt   = info_fold.binDt;

% Phase axis centred on 0: bins reordered to phases -0.5 .. 0.5
sh    = floor(Nbin/2);
phC   = ((0:Nbin-1) - sh).' / Nbin;                       % column
reord = @(v) circshift(v, sh, 1);

profTot = reord(sum(fold.profTotal, 2));                   % [Nbin x 1]
profSub = reshape(reord(sum(fold.prof, 3)), Nbin, nSub);   % [Nbin x nSub]
wSub  = reshape(reord(fold.weight),  Nbin, nSub);
w2Sub = reshape(reord(fold.weight2), Nbin, nSub);
if isfield(fold, 'weightX')
    wxSub = reshape(reord(fold.weightX), Nbin, nSub);
else
    wxSub = zeros(Nbin, nSub);                             % older fold files
end
wTot = sum(wSub, 2); w2Tot = sum(w2Sub, 2); wxTot = sum(wxSub, 2);

offMask = abs(phC) > 0.25;
base    = median(profTot(offMask), 'omitnan');
baseSub = median(profSub(offMask, :), 1, 'omitnan');

% ---- Expected profile ------------------------------------------------------------------
shape   = expectedShape(phC, info_fold, info_gen, opts.SubGrid);   % unit-peak pulsar profile
haveAbs = ~isempty(opts.InfoDisp) && ~isempty(opts.InfoIQ) && ~isempty(opts.InfoDedisp);
if haveAbs
    M    = expectedPowerModel(info_gen, opts.InfoDisp, opts.InfoIQ, opts.InfoDedisp, opts.InfoRx);
    Ps   = M.sigScale * shape;
    Pexp = Ps + M.Pn;
    sTB  = sqrt(M.var(Ps, dt));                            % per time bin, at each phase
    sigTot = sTB .* sqrt(w2Tot) ./ wTot;                   % per phase bin, total fold
else
    M  = struct('Pn', 0, 'hasNoise', false, 'hasRFI', false, 'sigScale', NaN);
    ok = ~isnan(profTot);
    q  = profTot - base;
    Ps = shape * (shape(ok)' * q(ok)) / (shape(ok)' * shape(ok));   % shape-only fit
    Pexp = Ps + base;
    sTB = nan(Nbin, 1); sigTot = nan(Nbin, 1);
end

% ---- Checks -----------------------------------------------------------------------------
if strcmpi(info_gen.envelopeMode, 'power')
    sigPh = info_gen.sigma * f0;
else
    sigPh = info_gen.sigma * f0 / sqrt(2);
end
win = abs(phC) <= 4*sigPh;
ce  = sum(phC(win) .* Ps(win)) / sum(Ps(win));
aC  = zeros(Nbin, 1);
aC(win) = (phC(win) - ce) / sum(Ps(win));

check = struct();
q = profTot - base;
check.totalOffset      = (sum(phC(win) .* q(win)) / sum(q(win)) - ce) * Per;
check.totalEnergyRatio = sum(q(win)) / sum(Ps(win));
check.baseline         = base;
if haveAbs
    check.totalOffsetNoise = combNoise(aC, sTB, wTot, w2Tot, wxTot) * Per;
    zr  = (profTot - Pexp) ./ sigTot;
    onP = Ps > opts.OnPulseFrac * M.sigScale;
else
    check.totalOffsetNoise = NaN;
    zr  = nan(Nbin, 1); onP = false(Nbin, 1);
end
check.zr = zr;

subOff = nan(nSub, 1); subNoise = nan(nSub, 1);
for s = 1:nSub
    p = profSub(:, s) - baseSub(s);
    if any(isnan(p(win))) || sum(p(win)) <= 0, continue; end
    subOff(s) = (sum(phC(win) .* p(win)) / sum(p(win)) - ce) * Per;
    if haveAbs
        subNoise(s) = combNoise(aC, sTB, wSub(:, s), w2Sub(:, s), wxSub(:, s)) * Per;
    end
end
check.subOffset = subOff;
check.subOffsetNoise = subNoise;
check.subTime = fold.subint.tMean;

fprintf('plotFoldCheck: NBin = %d (%.3g us), %d sub-int(s), %s assignment\n', ...
    Nbin, Per/Nbin*1e6, nSub, info_fold.assign);
if haveAbs
    fprintf('  expected pulse peak %.4g, noise baseline %.4g; measured baseline %.4g', ...
        M.sigScale, M.Pn, base);
    if M.hasRFI, fprintf('  (RFI present: not in the model)'); end
    fprintf('\n');
end
fprintf('  total profile: centroid offset %+.4f us (noise ~%.4f us), energy ratio %.5f\n', ...
    check.totalOffset*1e6, check.totalOffsetNoise*1e6, check.totalEnergyRatio);
if haveAbs
    statLine('on-pulse ', zr(onP));
    if M.hasNoise, statLine('off-pulse', zr(~onP)); end
end
okS = ~isnan(subOff);
if any(okS)
    fprintf('  sub-int centroid offsets (%d usable): mean %+.3f us, rms %.3f us', ...
        nnz(okS), mean(subOff(okS))*1e6, sqrt(mean(subOff(okS).^2))*1e6);
    if haveAbs
        fprintf(' (expected rms ~%.3f us)\n', sqrt(mean(subNoise(okS).^2))*1e6);
    else
        fprintf('\n');
    end
end

% ---- Plot ---------------------------------------------------------------------------------
fig = figure('Name', 'Fold check', 'Color', 'w');
tl  = tiledlayout(fig, 3, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, sprintf('Folded profile: P = %.6g ms, %d bins (%.3g \\mus), %d sub-int(s)', ...
    Per*1e3, Nbin, Per/Nbin*1e6, nSub));
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
    zz = zr;
    if ~M.hasNoise, zz(~onP) = NaN; end
    hR = plot(ax3, phC, zz, 'k.-', 'MarkerSize', 5); hold(ax3, 'on');
    yline(ax3, 0, 'r-', 'HandleVisibility', 'off');
    h1 = yline(ax3, [-1 1], 'b-'); h3 = yline(ax3, [-3 3], 'b--');
    hold(ax3, 'off'); ylim(ax3, [-5 5]);
    legend(ax3, [hR h1(1) h3(1)], {'(folded - expected)/\sigma', '\pm1\sigma', '\pm3\sigma'}, ...
        'Location', 'northeast');
    if M.hasNoise
        xlim(ax3, [-0.5 0.5]); title(ax3, 'Normalized residual, all phase bins');
    else
        xlim(ax3, [-zp zp]);   title(ax3, 'Normalized residual on the pulse');
    end
else
    plot(ax3, phC, profTot - Pexp, 'k'); xlim(ax3, [-0.5 0.5]);
    title(ax3, 'Residual (shape only: pass upstream info structs for \sigma)');
end
grid(ax3, 'on'); xlabel(ax3, 'Phase [turns]');

ax4 = nexttile(tl, 4);
imagesc(ax4, phC(zSel), 1:nSub, (profSub(zSel, :) - baseSub).'); axis(ax4, 'xy');
xline(ax4, 0, 'w--'); xlabel(ax4, 'Phase [turns]'); ylabel(ax4, 'Sub-integration');
cb = colorbar(ax4); cb.Label.String = 'power - baseline';
title(ax4, 'Sub-integrations (pulse at phase 0 in every row)');

ax5 = nexttile(tl, 5, [1 2]);
tS = fold.subint.tMean * 1e3;
if haveAbs && any(okS)
    errorbar(ax5, tS(okS), subOff(okS)*1e6, subNoise(okS)*1e6, 'ko', 'MarkerFaceColor', 'k');
elseif any(okS)
    plot(ax5, tS(okS), subOff(okS)*1e6, 'ko', 'MarkerFaceColor', 'k');
end
yline(ax5, 0, 'r-'); grid(ax5, 'on');
xlabel(ax5, 'Sub-integration mean time [ms]'); ylabel(ax5, 'offset [\mus]');
title(ax5, 'Per-sub-integration power-centroid offset (TOA preview, \pm1\sigma)');
end


% ========================================================================================
function statLine(name, z)
z = z(isfinite(z));
n = numel(z);
if n == 0, return; end
fprintf(['  normalized residual %s over %d bins: mean %+.3f, std %.3f ' ...
         '(expect 0 +- %.3f, 1 +- %.3f)\n'], name, n, mean(z), std(z), ...
    1/sqrt(n), 1/sqrt(2*n));
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


function shape = expectedShape(phC, info_fold, info_gen, nS)
% Unit-peak expected pulsar profile at the centred phases phC, averaged over
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