function fit = navFit(toa, psrs, site, opts)
%NAVFIT  The position of a fixed site from TOA residuals: a robust batch fit (block D4a; observer side).
%{
The observer predicts every pulse for the site it ASSUMES (makePredictor). If the true site is
dx away (Earth-fixed, ECEF), each TOA arrives early or late by the offset projected onto the
pulsar's direction. This fits dx to the residuals of many TOAs.

  fit = navFit(toa, psrs, site)
  fit = navFit(toa, psrs, site, 'Start', "calibrated")
  fit = navFit(struct('t', t, 'sigma', s), psrs, site, 'Predict', true)   % precision only

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
The clock is perfect here (block D4b adds it as correlated noise). Left out: the factor
1/(1 + Doppler) of the residual (~1e-4: 10 m on 100 km), and the Einstein site term (31 ns per
100 km, 2e-4).

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
         not needed with 'Predict'), sigma (s), k (pulsar index into psrs; default all 1).
  psrs   struct array: raDeg, decDeg per pulsar.
  site   makeTrajectory('fixed', ...) struct: the ASSUMED site (ecef, enu, params.Lon).

Name-value options:
  'Start'     "self" (default) or "calibrated" (above).
  'Robust'    true (default): the mixture; false: weighted least squares.
  'OutlierFraction'  eps (default 0.02); 'OutlierScale'  k (default 4).
  'MaxIter'   iterations of the robust fit (default 50).
  'Predict'   true: no residuals; only the covariance the schedule would give (Gaussian
              errors), for planning which pulsars and how long.

Output fit: names, p (estimates: dX dY [dZ] [m], c_k [s]), cov; east, radial, polar [m] and their
errors (radial = in the equatorial plane, away from the Earth's axis; polar along it; NaN where
unobservable) and covEq (east/radial 2 x 2); enu and covENU (3-D, "calibrated" only; NaN
otherwise); resid (post-fit residuals, s), z (resid / sigma), gamma (core probability), outlier
(gamma < 0.5), chi2Core (mean z^2 of the core TOAs), nIter, nTOA, start, robust, eps, k, made.
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
ax =  cos(th) .* nk(1, :).' + sin(th) .* nk(2, :).';   % R_z(theta)' n: x and y components
ay = -sin(th) .* nk(1, :).' + cos(th) .* nk(2, :).';
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
sc = 1 ./ sqrt(sum(A.^2, 1));                           % column scales: metres and seconds differ by 1e9
A = A .* sc;                                            % solved in scaled units, scaled back below

% the fit
w = 1 ./ sg.^2; gamma = ones(N, 1); nIter = 0;
if opts.Predict
    p = NaN(size(A, 2), 1); res = NaN(N, 1);
else
    r = toa.r(:);
    p = lsq(A, r, w);
    if opts.Robust
        e = opts.OutlierFraction; k = opts.OutlierScale;
        for nIter = 1:opts.MaxIter
            u = (r - A * p) ./ sg;
            core = (1 - e) * exp(-u.^2 / 2);
            wide = e / k * exp(-u.^2 / (2 * k^2));
            gamma = core ./ (core + wide);
            w = gamma ./ sg.^2 + (1 - gamma) ./ (k * sg).^2;
            pNew = lsq(A, r, w);
            step = max(abs(pNew - p) ./ sqrt(diag(inv(A' * (w .* A)))));
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
cov = inv(A' * (w .* A));
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
fit.z        = (res ./ sg).';
fit.gamma    = gamma.';
fit.outlier  = (gamma < 0.5).';
fit.chi2Core = mean(fit.z(~fit.outlier).^2);
fit.nIter    = nIter;
fit.nTOA     = N;
fit.start    = char(opts.Start);
fit.robust   = opts.Robust;
fit.eps      = opts.OutlierFraction;
fit.k        = opts.OutlierScale;
fit.made     = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm'));
end


% =====================================================================================
function p = lsq(A, r, w)
%LSQ  Weighted least squares: minimize sum w (r - A p)^2.
sw = sqrt(w);
p = (sw .* A) \ (sw .* r);
end
