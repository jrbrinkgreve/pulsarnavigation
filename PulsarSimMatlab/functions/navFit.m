function fit = navFit(toa, psrs, site, opts)
%NAVFIT  The position of a fixed site from TOA residuals: a robust batch fit (block D4a; observer side).
%{
The observer predicts every pulse for the site it ASSUMES (makePredictor). If the true site is
dx away (Earth-fixed, ECEF), each TOA arrives early or late by the offset projected onto the
pulsar's direction. This fits dx to the residuals of many TOAs.

  fit = navFit(toa, psrs, site)
  fit = navFit(toa, psrs, site, 'Start', "calibrated")
  fit = navFit(struct('t', t, 'sigma', s), psrs, site, 'Predict', true)   % precision only

A TOA merges sub-ints, often over hours, sometimes over two passes a day apart: its residual is
then the average over its sub-ints, not the value at one time. With toa.tSub (a cell, per TOA the
times of its sub-ints) each row is that average: a_i = mean_j a(tau_ij) (the sub-ints weigh
equally: equal length, constant flux), and the clock rate's column the mean time. Without it, a
TOA of the daily sine spread over 9 h of hour angle is damped to sinc(omega 4.5 h) = 0.79, one
across the daily gap to ~0.8 (found 22 Oct 2026, block D4d).

The model. A TOA of pulsar k at time t_i has the residual (observed - predicted)
  r_i = a_i . dx + c_k + noise,   a_i = -R_z(theta(t_i))' n_k / c,
with theta the Earth rotation angle (as delayModel), n_k the unit vector to the pulsar and c_k a
constant per pulsar: the clock offset plus the pulsar's own offset (profile reference point vs
timing model). The Earth's rotation turns the two equatorial components of dx (X, Y) into a
daily sine; the polar component (Z) gives a constant -dZ sin(dec_k)/c, the same for every TOA of
that pulsar - indistinguishable from c_k. Hence the switch 'Start' (Jasper's decision, 21 Oct):
  "self"        (default) the constants are unknown and fitted: dX, dY and one c_k per pulsar.
                The polar component stays unobservable for a fixed site.
  "calibrated"  the observer started from a known position and time, so the constants are known
                and already removed from r (c_k = 0): dX, dY, dZ - the full 3-D position.
Left out: the factor 1/(1 + Doppler) of the residual (~1e-4: 10 m on 100 km), and the Einstein
site term (31 ns per 100 km, 2e-4).

The clock ('Clock', block D4b). A real receiver clock adds its unestimated error db(t) to every
residual. It wanders: the crystal of makeClock by ~76 us a day and ~13 ms a month. TOAs close in
time share almost the same error, so it is not extra independent noise per TOA but correlated
noise: the fit becomes generalized least squares with the covariance
  Cov(r_i, r_j) = sigma_i^2 delta_ij + K_ij,
  K_ij = q1 min(tau_i, tau_j) + q2 (m^2 M / 2 - m^3 / 6),  m = min, M = max of tau_i, tau_j,
tau = t - t0 the time since the clock's start ('ClockStart'), q1 = h0/2 (white FM: the time error a
random walk), q2 = 2 pi^2 h_2 (random-walk FM) - the paper's two-state model (makeClock) for any two
times. Flicker FM (h_1) is left out of K: at daily scales it adds ~3 us against ~76 us. With "self"
the clock's offset merges into the constants c_k and its rate at t0 is a free parameter (one column
tau); with "calibrated" both were known at t0 and only the wander since then is unknown. As a
by-product the fit estimates the clock wander at every TOA (the GLS posterior mean K S^-1 (r - A p))
and predicts it at any later time with an error bar (clockAt, clockAtErr): what the observer feeds
back into its predictor (makePredictor 'Clock'). The position hardly suffers for a fixed site: the
clock is red, the position signal a daily sine.

Robustness. A few TOAs near the detection threshold land one pulse component off (the 100-day
scenario runs: 1.9 % beyond 3 sigma, 0.3-0.5 % beyond 4 sigma, an honest core). The noise model
is therefore a mixture: (1 - eps) N(0, sigma^2) + eps N(0, (k sigma)^2). The fit iterates: each
TOA gets the probability gamma that it belongs to the core, its weight becomes gamma / sigma^2 +
(1 - gamma) / (k sigma)^2, and the weighted least squares is solved again, until the estimate
changes by less than 1e-9 of its error (EM for the mixture with eps and k fixed). An outlier
then counts about k^2 = 16 times less, instead of pulling the fit. 'Robust' false: plain
weighted least squares.

Inputs:
  toa    struct of [1 x N] vectors: t (time of the TOA, s since J2000, TT), r (residual, s;
         not needed with 'Predict'), sigma (s), k (pulsar index into psrs; default all 1); optional
         tSub (a cell, per TOA the times of its sub-ints: the row averaged over them, above).
  psrs   struct array: raDeg, decDeg per pulsar.
  site   makeTrajectory('fixed', ...) struct: the ASSUMED site (ecef, enu, params.Lon).

Name-value options:
  'Start'     "self" (default) or "calibrated" (above).
  'Robust'    true (default): the mixture; false: weighted least squares.
  'OutlierFraction'  eps (default 0.02); 'OutlierScale'  k (default 4).
  'MaxIter'   iterations of the robust fit (default 50).
  'Clock'     [] (default: a perfect clock) or the clock's power-law coefficients [h0 h_1 h_2], or a
              makeClock struct (its field h): the clock wander becomes correlated noise (above).
  'ClockStart' [s since J2000] the clock's start t0 (default: the first TOA).
  'Predict'   true: no residuals; only the covariance the schedule would give (Gaussian
              errors), for planning which pulsars and how long.

Output fit: names, p (estimates: dX dY [dZ] [m], c_k [s]), cov; east, radial, polar [m] and their
errors (radial = in the equatorial plane, away from the Earth's axis; polar along it; NaN where
unobservable) and covEq (east/radial 2 x 2); enu and covENU (3-D, "calibrated" only; NaN
otherwise); resid (post-fit residuals r - A p, s; with a clock they contain its wander), z ((resid -
clock) / sigma), gamma (core probability), outlier (gamma < 0.5), chi2Core (mean z^2 of the core
TOAs), nIter, nTOA, start, robust, eps, k, made; with 'Clock' also clock (the clock's change since t0
at each TOA, s: rate tau + wander with "self", the wander with "calibrated"; its offset at t0 is in
the constants c_k), clockErr (including the uncertainty of the fitted parameters: universal
kriging), clockRate (the fitted rate, "self"; NaN otherwise), clockAt / clockAtErr (@(t) the same
predicted at any time t >= t0, e.g. the next day), clockStart, h.
%}

arguments
    toa struct
    psrs struct
    site struct
    opts.Start (1,1) string {mustBeMember(opts.Start, ["self" "calibrated"])} = "self"
    opts.Robust (1,1) logical = true
    opts.OutlierFraction (1,1) double {mustBeInRange(opts.OutlierFraction, 0, 0.5)} = 0.02
    opts.OutlierScale (1,1) double {mustBeGreaterThan(opts.OutlierScale, 1)} = 4
    opts.MaxIter (1,1) double {mustBeInteger, mustBePositive} = 50
    opts.Predict (1,1) logical = false
    opts.Clock = []
    opts.ClockStart (1,1) double = NaN
end

c = 299792458;
t = toa.t(:); sg = toa.sigma(:); N = numel(t);
kp = ones(N, 1); if isfield(toa, 'k') && ~isempty(toa.k), kp = toa.k(:); end
if ~isfield(site, 'ecef') || ~isfield(site, 'enu')
    error('navFit:site', 'site must be a makeTrajectory(''fixed'', ...) struct (fields ecef, enu).');
end
if numel(sg) ~= N || numel(kp) ~= N || any(kp < 1 | kp > numel(psrs) | kp ~= round(kp))
    error('navFit:toa', 'toa.t, toa.sigma and toa.k must have one value per TOA, k indexing psrs.');
end
if ~opts.Predict && (~isfield(toa, 'r') || numel(toa.r) ~= N)
    error('navFit:toa', 'toa.r (the residuals) is needed unless ''Predict'' is true.');
end

% design matrix: a_i = -R_z(theta_i)' n_k / c, columns X Y (Z), then one constant per pulsar
th = 2*pi * mod(0.7790572732640 + 1.00273781191135448 * (t - 69.184) / 86400, 1);   % as delayModel
n = zeros(3, numel(psrs));
for j = 1:numel(psrs)
    n(:, j) = [cosd(psrs(j).decDeg) * cosd(psrs(j).raDeg); cosd(psrs(j).decDeg) * sind(psrs(j).raDeg); sind(psrs(j).decDeg)];
end
nk = n(:, kp);                                         % [3 x N] direction of each TOA's pulsar
[cT, sT] = deal(cos(th), sin(th));                     % per TOA: the rotation, or its mean over the sub-ints
tm = t;                                                 % per TOA: its time, or its sub-ints' mean time
if isfield(toa, 'tSub') && ~isempty(toa.tSub)
    if numel(toa.tSub) ~= N, error('navFit:toa', 'toa.tSub must have one entry (sub-int times) per TOA.'); end
    for i = 1:N
        ts = toa.tSub{i};
        if isempty(ts), continue, end
        ths = 2*pi * mod(0.7790572732640 + 1.00273781191135448 * (ts(:) - 69.184) / 86400, 1);
        cT(i) = mean(cos(ths)); sT(i) = mean(sin(ths)); tm(i) = mean(ts);
    end
end
ax =  cT .* nk(1, :).' + sT .* nk(2, :).';             % R_z(theta)' n: x and y components (averaged)
ay = -sT .* nk(1, :).' + cT .* nk(2, :).';
A = -[ax, ay, nk(3, :).'] / c;                          % [N x 3]: d r / d (dX dY dZ)
used = unique(kp).';
if opts.Start == "self"
    A = A(:, 1:2);
    C = double(kp == used);                             % [N x nPsr] constants
    A = [A, C];
    names = [{'dX', 'dY'}, arrayfun(@(j) sprintf('c%d', j), used, 'UniformOutput', false)];
else
    names = {'dX', 'dY', 'dZ'};
end
nPos = 2 + (opts.Start == "calibrated");
K = []; h = []; t0 = NaN; q = [0 0];
if ~isempty(opts.Clock)                                 % the clock's wander as correlated noise
    if isstruct(opts.Clock), h = opts.Clock.h; else, h = opts.Clock; end
    if numel(h) ~= 3 || any(h < 0)
        error('navFit:clock', 'Clock must be [h0 h_1 h_2] >= 0 or a makeClock struct.');
    end
    t0 = opts.ClockStart; if isnan(t0), t0 = min(t); end
    tau = t - t0;
    if any(tau < 0), error('navFit:clock', 'ClockStart is after the first TOA.'); end
    q = [h(1) / 2, 2 * pi^2 * h(3)];
    K = clockCov(tau, tau.', q);
    if opts.Start == "self"                             % the clock rate at t0 is free
        A = [A, tm - t0]; names{end+1} = 'rate';
    end
end
sc = 1 ./ sqrt(sum(A.^2, 1));                           % column scales: metres and seconds differ by 1e9
A = A .* sc;                                            % solved in scaled units, scaled back below

% the fit (with a clock: generalized least squares, the wander estimated alongside)
w = 1 ./ sg.^2; gamma = ones(N, 1); nIter = 0; xh = zeros(N, 1);
if opts.Predict
    p = NaN(size(A, 2), 1); res = NaN(N, 1);
    [~, cov] = gls(A, [], w, K);
else
    r = toa.r(:);
    [p, cov, xh] = gls(A, r, w, K);
    if opts.Robust
        e = opts.OutlierFraction; k = opts.OutlierScale;
        for nIter = 1:opts.MaxIter
            u = (r - A * p - xh) ./ sg;                 % the white part of the residual
            core = (1 - e) * exp(-u.^2 / 2);
            wide = e / k * exp(-u.^2 / (2 * k^2));
            gamma = core ./ (core + wide);
            w = gamma ./ sg.^2 + (1 - gamma) ./ (k * sg).^2;
            [pNew, cov, xh] = gls(A, r, w, K);
            step = max(abs(pNew - p) ./ sqrt(diag(cov)));
            p = pNew;
            if step < 1e-9, break, end
        end
        if step >= 1e-9
            warning('navFit:converge', 'The robust fit did not converge in %d iterations (last step %.1e sigma).', ...
                opts.MaxIter, step);
        end
    end
    res = r - A * p;
end
pS = p;                                                 % scaled, for the clock prediction
p = sc.' .* p; cov = (sc.' * sc) .* cov;               % back to metres and seconds

% reported in the site's frame: east and radial (equatorial plane), polar (the Earth's axis)
lon = site.params.Lon;
Q = [-sind(lon), cosd(lon); cosd(lon), sind(lon)];     % rows: east, radial in terms of X, Y
eq = Q * p(1:2); covEq = Q * cov(1:2, 1:2) * Q';
fit = struct();
fit.names   = names;
fit.p       = p;
fit.cov     = cov;
fit.east    = eq(1);  fit.eastErr   = sqrt(covEq(1, 1));
fit.radial  = eq(2);  fit.radialErr = sqrt(covEq(2, 2));
fit.covEq   = covEq;
fit.polar   = NaN;    fit.polarErr  = NaN;
fit.enu     = NaN(3, 1); fit.covENU = NaN(3);
if nPos == 3
    fit.polar = p(3); fit.polarErr = sqrt(cov(3, 3));
    E = site.enu;                                       % columns east, north, up (ECEF)
    fit.enu = E' * p(1:3); fit.covENU = E' * cov(1:3, 1:3) * E;
end
fit.resid    = res.';
fit.z        = ((res - xh) ./ sg).';
fit.gamma    = gamma.';
fit.outlier  = (gamma < 0.5).';
fit.chi2Core = mean(fit.z(~fit.outlier).^2);
fit.nIter    = nIter;
fit.nTOA     = N;
fit.start    = char(opts.Start);
fit.robust   = opts.Robust;
fit.eps      = opts.OutlierFraction;
fit.k        = opts.OutlierScale;
if ~isempty(K)
    % the clock's change since t0 (rate tau + wander; with "calibrated" the wander alone) and its error:
    % universal kriging, var = K_qq - k' S^-1 k + g' cov g, g = a - A' S^-1 k (a: the rate column at q)
    S = K + diag(1 ./ w); L = chol(S, 'lower');
    covS = (sc.' * sc) .\ cov;                         % the parameter covariance in scaled units
    rateCol = opts.Start == "self";
    if opts.Predict
        alpha = zeros(N, 1); rate = 0;
    else
        alpha = L' \ (L \ (toa.r(:) - A * pS));          % S^-1 (r - A p)
        rate = 0; if rateCol, rate = p(end); end
    end
    aOf = @(tq) [zeros(size(A, 2) - rateCol, numel(tq)); repmat(sc(end) * (tq(:).' - t0), rateCol, 1)];
    errOf = @(kq, kqq, tq) sqrt(max(kqq - sum((L \ kq).^2, 1).' + ...
        sum((aOf(tq) - A' * (L' \ (L \ kq))) .* (covS * (aOf(tq) - A' * (L' \ (L \ kq)))), 1).', 0));
    fit.clock      = (xh + rate * tau).';
    fit.clockErr   = errOf(K, diag(K), t).';
    fit.clockRate  = NaN;
    if rateCol, fit.clockRate = p(end); end
    fit.clockAt    = @(tq) reshape(clockCov(tq(:) - t0, tau.', q) * alpha + rate * (tq(:) - t0), size(tq));
    fit.clockAtErr = @(tq) reshape(errOf(clockCov(tau, tq(:).' - t0, q), clockCov(tq(:) - t0, tq(:) - t0, q, true), tq), ...
        size(tq));
    fit.clockStart = t0;
    fit.h          = h;
end
fit.made     = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm'));
end


% =====================================================================================
function [p, cov, xh] = gls(A, r, w, K)
%GLS  Weighted (K empty) or generalized least squares with S = K + diag(1/w); xh = K S^-1 (r - A p),
% the clock wander's posterior mean (zeros without a clock). r = [] gives the covariance only.
N = size(A, 1); xh = zeros(N, 1); p = [];
if isempty(K)
    if ~isempty(r), p = lsq(A, r, w); end
    cov = inv(A' * (w .* A));
else
    L = chol(K + diag(1 ./ w), 'lower');
    Aw = L \ A;
    if ~isempty(r)
        p = Aw \ (L \ r);
        xh = K * (L' \ (L \ (r - A * p)));
    end
    cov = inv(Aw' * Aw);
end
end

function C = clockCov(ta, tb, q, diagOnly)
%CLOCKCOV  Covariance of the clock's time error at times ta [n x 1] and tb [1 x m] since its start:
% white FM q1 min + random-walk FM q2 (m^2 M / 2 - m^3 / 6) (m, M the smaller and larger time).
% diagOnly: ta = tb, only the variances [n x 1].
if nargin > 3 && diagOnly
    C = q(1) * ta + q(2) * ta.^3 / 3;
    return
end
m = min(ta, tb); M = max(ta, tb);
C = q(1) * m + q(2) * (m.^2 .* M / 2 - m.^3 / 6);
end

function p = lsq(A, r, w)
%LSQ  Weighted least squares: minimize sum w (r - A p)^2.
sw = sqrt(w);
p = (sw .* A) \ (sw .* r);
end
