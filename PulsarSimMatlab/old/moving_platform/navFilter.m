function fit = navFilter(toa, psrs, traj, opts)
%NAVFILTER  Position, velocity and clock of a moving receiver from TOA residuals: a Kalman filter and
%smoother around its dead-reckoning track (block D4d; observer side).
%{
On a ship the observer knows its track only from dead reckoning (DR): speed log and compass. The DR
error grows - an unknown current of 0.3 m/s is 26 km a day. The observer predicts every pulse for the
DR track (makePredictor), so each TOA residual measures the DR error projected onto the pulsar's
direction, plus the clock. This filter estimates the corrections to the DR track over time:

  fit = navFilter(toa, psrs, trajDR)
  fit = navFilter(toa, psrs, trajDR, 'Clock', makeClock(...), 'Start', "calibrated")

The state, all in metres (the clock as c times its offset, as in GNSS receivers):
  dp    position error, true - DR, in the LOCAL frame of the DR track: east, north (and up if the
        height is not known) [m]
  dv    velocity error in the same local frame [m/s]: dp' = dv; dv a random walk ('VelocityNoise'),
        its start unknown to 'VelocityPrior' (an unknown current, a log bias)
The local frame rotates with the ship: a current pushes it along its own east and north, so a
constant current is a constant dv here. (A first version kept dv fixed in Earth-fixed coordinates;
over 35 deg of longitude a constant east current then gets a 37 % vertical part, which the height
constraint forbade: the east error was pinned and the filter overconfident by a factor 30.)
  b, d       the clock's offset and rate, c x [s] and c x [s/s] ('Clock': the two-state model of
             makeClock, q1 = h0/2, q2 = 2 pi^2 h_2, as navFit)
  J_k        one constant per pulsar ("self": unknown; "calibrated": known, 0)
The measurement of a TOA of pulsar k at time t (residual r, s): c r = a(t) . M(t)' dp + b + J_k + noise,
with a(t) = -R_z(theta(t))' n_k (the Earth rotation angle as delayModel) and M(t)' = [east north up]
the local frame of the DR track at t (WGS84 normal). This is linear in the state, so
the filter is an ordinary linear Kalman filter around the DR track (no extended filter needed); the
Rauch-Tung-Striebel smoother then gives the best estimate at every TOA from all of them.

The switch 'Start' (Jasper, 21 Oct), as navFit:
  "self"        the clock offset at t0 is part of the constants J_k, its rate is unknown;
  "calibrated"  a known position and time at the start (port): J_k = 0, and the clock's offset and
                rate at t0 known to 'ClockPrior' (default exactly: b = d = 0). A port calibration is
                never exact: 1 ms of offset is 300 km in c b, more than a TOA's error. An exact prior
                then states the clock too precisely (testNavFilter check 6: clock pulls rms 2.6; the
                position is hardly hurt, as the clock's own wander takes the offset within a day).
Height ('Height', default true): a ship stays at sea level, so dp and dv have no up component. That
is what lets one pulsar find north: the radial (equatorial) component then fixes it. The horizontal
error is taken in the tangent plane of the DR point; on the real surface an error d lies d^2/2R lower
(7 km at 300 km, ~14 us in a residual). In the closed loop (runShip) the filter re-centres the track
every day, keeping d small. 'Height' false adds the up components (a site of unknown height).

Robustness: a TOA's innovation nu (its residual against the prediction from all earlier TOAs) with
variance S is weighed with the mixture of navFit: gamma = core probability under (1 - eps) N(0, S) +
eps N(0, S - R + k^2 R), and the update uses the noise variance R / (gamma + (1 - gamma) / k^2). An
outlier then counts ~k^2 = 16 times less. One pass (the filter does not iterate, unlike navFit).

Inputs:
  toa    struct of [1 x N]: t (s since J2000, TT), r (residual against the prediction for traj, s; not
         needed with 'Predict'), sigma (s), k (pulsar index, default 1). Any order; sorted inside.
         Optional tSub: a cell, per TOA the times of its sub-ints. A TOA merges sub-ints over hours,
         sometimes over two passes, so its residual is their average: its row is then the mean over
         the sub-ints (the sub-ints weigh equally), with the position at the TOA's time t moved by the
         velocity: mean a M' (dp + dv (tau - t)), and the clock b + d mean(tau - t). Without it a TOA
         spread over 9 h of hour angle is damped to 0.79 of the daily signal (22 Oct 2026, D4d).
  psrs   struct array: raDeg, decDeg.
  traj   makeTrajectory struct ('ECEF' frame): the DR track (or a fixed site) the residuals refer to.

Name-value options:
  'Start'          "self" (default) | "calibrated".
  'Clock'          [] (a perfect clock) | [h0 h_1 h_2] | a makeClock struct.
  'ClockPrior'     [s, s/s] the sigmas of the clock's offset and rate at the first TOA with
                   "calibrated" (default [0 0]: known exactly); needs 'Clock'. Not used with "self"
                   (the offset is in the J_k, the rate free).
  'Height'         true (default) | false (a site or probe whose height is unknown).
  'PositionPrior'  [m] initial sigma of dp per axis (default 1e7: unknown).
  'VelocityPrior'  [m/s] initial sigma of dv per axis (default 0.3; 0 for a fixed site).
  'VelocityNoise'  [m^2/s^3] spectral density of dv's random walk per axis (default 1.2e-7: ~0.1 m/s
                   per day; 0 for a fixed site).
  'Robust', 'OutlierFraction', 'OutlierScale'  as navFit (true, 0.02, 4).
  'Predict'        true: covariances only (no residuals needed), for planning.

Output fit: t (sorted TOA times), order (the sort), names, x (smoothed states, [nState x N]), P
(smoothed covariances, [nState x nState x N]), xF, PF (filtered), enu ([3 x N] position error in the
local east-north-up frame of the DR track, m; up 0 with 'Height') and enuErr, clock (b / c, s) and clockErr, gamma,
outlier, nu (innovations / sqrt(S)), predict (@(t) state and covariance extrapolated from the last
filtered state, for feedback to the predictor), opts, made.
%}

arguments
    toa struct
    psrs struct
    traj struct
    opts.Start (1,1) string {mustBeMember(opts.Start, ["self" "calibrated"])} = "self"
    opts.Clock = []
    opts.ClockPrior (1,2) double {mustBeNonnegative} = [0 0]
    opts.Height (1,1) logical = true
    opts.PositionPrior (1,1) double {mustBePositive} = 1e7
    opts.VelocityPrior (1,1) double {mustBeNonnegative} = 0.3
    opts.VelocityNoise (1,1) double {mustBeNonnegative} = 1.2e-7
    opts.Robust (1,1) logical = true
    opts.OutlierFraction (1,1) double {mustBeInRange(opts.OutlierFraction, 0, 0.5)} = 0.02
    opts.OutlierScale (1,1) double {mustBeGreaterThan(opts.OutlierScale, 1)} = 4
    opts.Predict (1,1) logical = false
end

c = 299792458;
if ~isfield(traj, 'frame') || ~strcmpi(traj.frame, 'ECEF')
    error('navFilter:traj', 'traj must be a makeTrajectory struct in the ''ECEF'' frame.');
end
t = toa.t(:).'; N = numel(t); sg = toa.sigma(:).';
kp = ones(1, N); if isfield(toa, 'k') && ~isempty(toa.k), kp = toa.k(:).'; end
if numel(sg) ~= N || numel(kp) ~= N || any(kp < 1 | kp > numel(psrs) | kp ~= round(kp))
    error('navFilter:toa', 'toa.t, toa.sigma and toa.k must have one value per TOA, k indexing psrs.');
end
haveR = ~opts.Predict;
if haveR && (~isfield(toa, 'r') || numel(toa.r) ~= N)
    error('navFilter:toa', 'toa.r (the residuals) is needed unless ''Predict'' is true.');
end
[t, order] = sort(t); sg = sg(order); kp = kp(order);
if haveR, z = c * toa.r(order); else, z = zeros(1, N); end

% the state layout
haveClk = ~isempty(opts.Clock);
nPsr = numel(psrs);
nD = 3 - opts.Height;                                   % local axes: east, north (, up)
iP = 1:nD; iV = nD + (1:nD); nS = 2 * nD;
iB = []; iD = []; iJ = [];
if haveClk, iB = nS + 1; iD = nS + 2; nS = nS + 2; end
if opts.Start == "self", iJ = nS + (1:nPsr); nS = nS + nPsr; end
ax = {'E', 'N', 'U'};
names = [strcat('d', ax(1:nD)), strcat('v', ax(1:nD)), repmat({'b', 'd'}, 1, haveClk), ...
    arrayfun(@(j) sprintf('J%d', j), 1:numel(iJ), 'UniformOutput', false)];
q = [0 0];
if haveClk
    if isstruct(opts.Clock), h = opts.Clock.h; else, h = opts.Clock; end
    if numel(h) ~= 3 || any(h < 0), error('navFilter:clock', 'Clock must be [h0 h_1 h_2] >= 0 or a makeClock struct.'); end
    q = c^2 * [h(1) / 2, 2 * pi^2 * h(3)];             % [m^2/s] and [m^2/s^3]
elseif any(opts.ClockPrior > 0)
    error('navFilter:clock', '''ClockPrior'' needs ''Clock'' (the clock''s states).');
end

% the prior at t0 (the first TOA)
x = zeros(nS, 1);
P0 = zeros(nS);
P0(iP, iP) = opts.PositionPrior^2 * eye(nD);
P0(iV, iV) = opts.VelocityPrior^2 * eye(nD);
if haveClk && opts.Start == "self", P0(iD, iD) = (c * 1e-5)^2; end   % the rate unknown (well beyond 2 ppm)
if haveClk && opts.Start == "calibrated"                % the port calibration's own error
    P0(iB, iB) = (c * opts.ClockPrior(1))^2; P0(iD, iD) = (c * opts.ClockPrior(2))^2;
end
if ~isempty(iJ), P0(iJ, iJ) = 1e7^2 * eye(nPsr); end
P = P0;

% the measurement rows and the local frame along the track
n = zeros(3, nPsr);
for j = 1:nPsr
    n(:, j) = [cosd(psrs(j).decDeg) * cosd(psrs(j).raDeg); cosd(psrs(j).decDeg) * sind(psrs(j).raDeg); sind(psrs(j).decDeg)];
end
th = 2*pi * mod(0.7790572732640 + 1.00273781191135448 * (t - 69.184) / 86400, 1);   % as delayModel
[E, Nn, U] = localFrame(traj.pos(t));
% the rows of the position, velocity and clock-rate states: at each TOA's time, or averaged over its sub-ints
rowP = zeros(N, nD); rowV = zeros(N, nD); rowD = zeros(N, 1);
for i = 1:N
    rowP(i, :) = posRow(th(i), n(:, kp(i)), [E(:, i), Nn(:, i), U(:, i)], nD);
end
if isfield(toa, 'tSub') && ~isempty(toa.tSub)
    if numel(toa.tSub) ~= N, error('navFilter:toa', 'toa.tSub must have one entry (sub-int times) per TOA.'); end
    tS = toa.tSub(order);
    for i = find(~cellfun(@isempty, tS(:).'))
        ts = tS{i}(:).';
        thS = 2*pi * mod(0.7790572732640 + 1.00273781191135448 * (ts - 69.184) / 86400, 1);
        [Es, Ns, Us] = localFrame(traj.pos(ts));
        Rs = zeros(numel(ts), nD);
        for j = 1:numel(ts), Rs(j, :) = posRow(thS(j), n(:, kp(i)), [Es(:, j), Ns(:, j), Us(:, j)], nD); end
        rowP(i, :) = mean(Rs, 1); rowV(i, :) = mean(Rs .* (ts.' - t(i)), 1); rowD(i) = mean(ts - t(i));
    end
end

% forward pass
[xF, xPr] = deal(zeros(nS, N)); [PF, PPr] = deal(zeros(nS, nS, N));
[gamma, nuN] = deal(ones(1, N), zeros(1, N));
e = opts.OutlierFraction; kO = opts.OutlierScale;
tPrev = t(1);
for i = 1:N
    [F, Q] = transition(t(i) - tPrev, nS, iP, iV, iB, iD, opts.VelocityNoise, q);
    x = F * x; P = F * P * F' + Q;
    xPr(:, i) = x; PPr(:, :, i) = P;
    H = zeros(1, nS);
    H(iP) = rowP(i, :); H(iV) = rowV(i, :);
    if haveClk, H(iB) = 1; H(iD) = rowD(i); end
    if ~isempty(iJ), H(iJ(kp(i))) = 1; end
    R = (c * sg(i))^2;
    nu = z(i) - H * x; S = H * P * H' + R;
    Reff = R;
    if haveR && opts.Robust
        Sw = S - R + kO^2 * R;
        core = (1 - e) * exp(-nu^2 / (2 * S)) / sqrt(S);
        wide = e * exp(-nu^2 / (2 * Sw)) / sqrt(Sw);
        gamma(i) = core / (core + wide);
        Reff = R / (gamma(i) + (1 - gamma(i)) / kO^2);
    end
    nuN(i) = nu / sqrt(S);
    if haveR, [x, P] = update(x, P, H, z(i), Reff); else, [~, P] = update(x, P, H, 0, Reff); end
    xF(:, i) = x; PF(:, :, i) = P;
    tPrev = t(i);
end

% backward pass (Rauch-Tung-Striebel)
xS = xF; PS = PF;
for i = N-1:-1:1
    F = transition(t(i+1) - t(i), nS, iP, iV, iB, iD, opts.VelocityNoise, q);
    C = smootherGain(PF(:, :, i) * F', PPr(:, :, i+1));
    xS(:, i) = xF(:, i) + C * (xS(:, i+1) - xPr(:, i+1));
    PS(:, :, i) = PF(:, :, i) + C * (PS(:, :, i+1) - PPr(:, :, i+1)) * C';
    PS(:, :, i) = (PS(:, :, i) + PS(:, :, i)') / 2;
end

% the position error in the local frame of the track (the states themselves), and the clock
[enu, enuErr] = deal(zeros(3, N));
enu(1:nD, :) = xS(iP, :);
for i = 1:N
    enuErr(1:nD, i) = sqrt(max(diag(PS(iP, iP, i)), 0));
end
fit = struct();
fit.t = t; fit.order = order; fit.names = names;
fit.x = xS; fit.P = PS; fit.xF = xF; fit.PF = PF;
fit.enu = enu; fit.enuErr = enuErr;
fit.clock = NaN(1, N); fit.clockErr = NaN(1, N);
if haveClk
    fit.clock = xS(iB, :) / c; fit.clockErr = sqrt(squeeze(PS(iB, iB, :)).') / c;
end
fit.gamma = gamma; fit.outlier = gamma < 0.5; fit.nu = nuN;
xL = xF(:, end); PL = PF(:, :, end); tL = t(end); vn = opts.VelocityNoise;
fit.predict = @(tq) extrapolate(tq, tL, xL, PL, nS, iP, iV, iB, iD, vn, q);
fit.index = struct('P', iP, 'V', iV, 'B', iB, 'D', iD, 'J', iJ);
fit.opts = opts;
fit.made = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm'));
end


% =====================================================================================
function h = posRow(th, n, M, nD)
%POSROW  A TOA's row for the position error in the local frame: the Earth-fixed row -R_z(theta)' n turned
% into the local frame M = [east north up].
aE = -[cos(th) * n(1) + sin(th) * n(2), -sin(th) * n(1) + cos(th) * n(2), n(3)];
h = aE * M(:, 1:nD);
end


function [F, Q] = transition(dt, nS, iP, iV, iB, iD, qv, q)
%TRANSITION  State transition over dt: position integrates velocity (a random walk), the clock is the
% two-state model (offset integrates rate; white FM q1, random-walk FM q2), the constants stay.
F = eye(nS); Q = zeros(nS);
nD = numel(iP);
F(iP, iV) = dt * eye(nD);
Q(iP, iP) = qv * dt^3 / 3 * eye(nD); Q(iP, iV) = qv * dt^2 / 2 * eye(nD);
Q(iV, iP) = Q(iP, iV); Q(iV, iV) = qv * dt * eye(nD);
if ~isempty(iB)
    F(iB, iD) = dt;
    Q([iB iD], [iB iD]) = [q(1) * dt + q(2) * dt^3 / 3, q(2) * dt^2 / 2; q(2) * dt^2 / 2, q(2) * dt];
end
end

function [x, P] = update(x, P, H, z, R)
%UPDATE  One scalar Kalman update in Joseph form (keeps P symmetric and positive).
S = H * P * H' + R;
K = P * H' / S;
x = x + K * (z - H * x);
I = eye(numel(x)) - K * H;
P = I * P * I' + K * R * K';
end

function C = smootherGain(A, Pp)
%SMOOTHERGAIN  A / Pp on the states with variance (a fixed site has none in its velocity), with Pp scaled
% to a correlation matrix first: its entries span ~15 orders of magnitude (1e14 m^2 down to 1e-2 m^2/s^2).
C = zeros(size(A));
k = find(diag(Pp) > 0);
if isempty(k), return, end
D = sqrt(diag(Pp(k, k))).';
C(:, k) = ((A(:, k) ./ D) / (Pp(k, k) ./ (D.' * D))) ./ D;
end

function [E, N, U] = localFrame(r)
%LOCALFRAME  East, north, up unit vectors [3 x M] at Earth-fixed positions r [3 x M] (WGS84 normal,
% Bowring's latitude).
a = 6378137; f = 1 / 298.257223563; e2 = f * (2 - f); b = a * (1 - f); ep2 = (a^2 - b^2) / b^2;
p = hypot(r(1, :), r(2, :)); lon = atan2(r(2, :), r(1, :));
th = atan2(r(3, :) * a, p * b);
lat = atan2(r(3, :) + ep2 * b * sin(th).^3, p - e2 * a * cos(th).^3);
E = [-sin(lon); cos(lon); zeros(size(lon))];
N = [-sin(lat) .* cos(lon); -sin(lat) .* sin(lon); cos(lat)];
U = [cos(lat) .* cos(lon); cos(lat) .* sin(lon); sin(lat)];
end

function [x, P] = extrapolate(tq, tL, xL, PL, nS, iP, iV, iB, iD, qv, q)
%EXTRAPOLATE  The state and covariance at time tq >= the last TOA, from the last filtered state.
[F, Q] = transition(tq - tL, nS, iP, iV, iB, iD, qv, q);
x = F * xL; P = F * PL * F' + Q;
end
