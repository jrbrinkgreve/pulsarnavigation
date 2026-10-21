function traj = makeTrajectory(type, opts)
%MAKETRAJECTORY  Where the receiver is: a trajectory for the delay model (block D1a).
%{
The phase prediction and the navigation fit need the receiver's position at any time. This
builds it as one struct with a position and a velocity function, whatever the receiver is: a
fixed site, a ship (waypoints), a deep-space probe (barycentric positions), or anything else
that can give positions - e.g. a Kalman filter's propagated state ('custom'): the entry point
for sensor fusion (Jasper, 21 Oct 2026). delayModel turns it into delays; the simulator's truth
and the observer's prediction are two trajectories.

  traj = makeTrajectory('fixed', 'Lat', 52, 'Lon', 5)
  traj = makeTrajectory('fixed', 'Lat', 52, 'Lon', 5, 'OffsetENU', [100e3 0 0])   % 100 km east
  traj = makeTrajectory('track', 'Time', t, 'Lat', lat, 'Lon', lon)               % ship waypoints
  traj = makeTrajectory('barycentric', 'Time', t, 'Pos', P)                         % probe, [3 x K] m
  traj = makeTrajectory('custom', 'Frame', 'ECEF', 'PosFcn', @(t) ..., 'VelFcn', @(t) ...)

Conventions:
  time   t [s] since J2000.0 (2000-01-01 12:00 TT), the time the Earth model uses; row vectors.
  frame  'ECEF': Earth-fixed (rotates with the Earth; x to lat 0 / lon 0, z to the north pole);
         'SSB': barycentric, equatorial axes (as ICRS). Positions [m], [3 x N] for N times.
  vel    relative to the trajectory's own frame (a fixed site: 0); delayModel adds the Earth's
         rotation for 'ECEF'.
  Geodetic coordinates on the WGS84 ellipsoid (a = 6378137 m, 1/f = 298.257223563): latitude
  and longitude [deg] (east positive), height [m] above the ellipsoid.

Types and their options:
  'fixed'        'Lat', 'Lon', 'Height' (0): one point, at all times. 'OffsetENU' [east north
                 up] [m] (default [0 0 0]) moves it in the local horizon frame - e.g. the true
                 site 100 km east of the assumed one.
  'track'        'Time' [1 x K] (increasing), 'Lat', 'Lon', 'Height' (0 or [1 x K]): waypoints,
                 converted to ECEF and interpolated per coordinate ('Interp': 'linear' (default;
                 constant velocity per leg), 'pchip' (smooth, no overshoot at turns) or 'spline'
                 (smoothest; for orbits)). Longitudes may wrap; the
                 interpolation runs in ECEF, so 180 deg is no problem (chord vs arc between
                 waypoints 10 km apart: 2 m). Outside [Time(1), Time(end)]: error.
  'barycentric'  'Time' [1 x K], 'Pos' [3 x K] [m] in the 'SSB' frame; 'Interp' as 'track'.
  'custom'       'Frame' ('ECEF' or 'SSB'), 'PosFcn' @(t) -> [3 x N], 'VelFcn' @(t) -> [3 x N]
                 (default: central difference of PosFcn over 'VelStep' s, default 1).

Output traj: type, frame, pos (@(t) -> [3 x N] m), vel (@(t) -> [3 x N] m/s), tRange [first
last] time it is valid for ([-Inf Inf] for 'fixed' and 'custom'), params (the inputs), and for
'fixed' also ecef (the point, [3 x 1]) and enu (the local east, north, up unit vectors, columns).
%}

arguments
    type {mustBeTextScalar}
    opts.Lat       double = []
    opts.Lon       double = []
    opts.Height    double = 0
    opts.OffsetENU (1,3) double = [0 0 0]
    opts.Time      double = []
    opts.Pos       double = []
    opts.Interp    {mustBeTextScalar} = 'linear'
    opts.Frame     {mustBeTextScalar} = 'ECEF'
    opts.PosFcn    = []
    opts.VelFcn    = []
    opts.VelStep   (1,1) double {mustBePositive} = 1
end

type = lower(char(type));
traj = struct('type', type, 'frame', '', 'pos', [], 'vel', [], 'tRange', [-Inf Inf], 'params', opts);
switch type
    case 'fixed'
        need(opts, {'Lat', 'Lon'}, type);
        [r, E] = geodeticToEcef(opts.Lat, opts.Lon, opts.Height);
        r = r + E * opts.OffsetENU(:);
        traj.frame = 'ECEF';
        traj.pos = @(t) repmat(r, 1, numel(t));
        traj.vel = @(t) zeros(3, numel(t));
        traj.ecef = r;
        traj.enu = E;
    case 'track'
        need(opts, {'Time', 'Lat', 'Lon'}, type);
        tk = checkTimes(opts.Time);
        K = numel(tk);
        if numel(opts.Lat) ~= K || numel(opts.Lon) ~= K || ~any(numel(opts.Height) == [1 K])
            error('makeTrajectory:size', 'Lat, Lon (and Height) need one value per waypoint.');
        end
        h = opts.Height(:).' .* ones(1, K);
        P = zeros(3, K);
        for k = 1:K
            P(:, k) = geodeticToEcef(opts.Lat(k), opts.Lon(k), h(k));
        end
        [traj.pos, traj.vel] = interpolant(tk, P, opts.Interp);
        traj.frame = 'ECEF';
        traj.tRange = tk([1 end]);
    case 'barycentric'
        need(opts, {'Time', 'Pos'}, type);
        tk = checkTimes(opts.Time);
        if ~isequal(size(opts.Pos), [3 numel(tk)])
            error('makeTrajectory:size', 'Pos must be [3 x K], one column per time.');
        end
        [traj.pos, traj.vel] = interpolant(tk, opts.Pos, opts.Interp);
        traj.frame = 'SSB';
        traj.tRange = tk([1 end]);
    case 'custom'
        fr = upper(char(opts.Frame));
        if ~any(strcmp(fr, {'ECEF', 'SSB'}))
            error('makeTrajectory:frame', 'Frame must be ''ECEF'' or ''SSB''.');
        end
        if ~isa(opts.PosFcn, 'function_handle')
            error('makeTrajectory:fcn', 'PosFcn must be a function handle @(t) -> [3 x N].');
        end
        traj.frame = fr;
        traj.pos = opts.PosFcn;
        if isa(opts.VelFcn, 'function_handle')
            traj.vel = opts.VelFcn;
        else
            f = opts.PosFcn; hs = opts.VelStep;
            traj.vel = @(t) (f(t + hs) - f(t - hs)) / (2 * hs);
        end
    otherwise
        error('makeTrajectory:type', 'Type must be ''fixed'', ''track'', ''barycentric'' or ''custom''.');
end
end


% =====================================================================================
function [r, E] = geodeticToEcef(lat, lon, h)
%GEODETICTOECEF  WGS84 geodetic [deg, deg, m] -> ECEF [m] (column), and the local east, north,
% up unit vectors as the columns of E.
a = 6378137; f = 1 / 298.257223563; e2 = f * (2 - f);
p = deg2rad(lat); l = deg2rad(lon);
N = a / sqrt(1 - e2 * sin(p)^2);
r = [(N + h) * cos(p) * cos(l); (N + h) * cos(p) * sin(l); (N * (1 - e2) + h) * sin(p)];
E = [-sin(l), -sin(p) * cos(l), cos(p) * cos(l);
      cos(l), -sin(p) * sin(l), cos(p) * sin(l);
      0,       cos(p),          sin(p)];
end

function [posF, velF] = interpolant(tk, P, method)
%INTERPOLANT  Position and velocity functions through the points P [3 x K] at times tk:
% piecewise linear (velocity constant per leg), pchip or cubic spline (velocity: the derivative
% of the cubic pieces, exact). Times outside [tk(1), tk(end)]: error.
method = lower(char(method));
switch method
    case 'linear'
        posF = @(t) linPos(tk, P, inRange(t, tk));
        velF = @(t) linVel(tk, P, inRange(t, tk));
    case {'pchip', 'spline'}
        if strcmp(method, 'pchip'), pp = pchip(tk, P); else, pp = spline(tk, P); end   % 3 at once
        [br, cf, ~, ord, dim] = unmkpp(pp);
        dpp = mkpp(br, cf(:, 1:ord-1) .* (ord-1:-1:1), dim);
        posF = @(t) ppval(pp, inRange(t, tk));
        velF = @(t) ppval(dpp, inRange(t, tk));
    otherwise
        error('makeTrajectory:interp', 'Interp must be ''linear'', ''pchip'' or ''spline''.');
end
end

function r = linPos(tk, P, t)
r = zeros(3, numel(t));
for c = 1:3, r(c, :) = interp1(tk, P(c, :), t, 'linear'); end
end

function v = linVel(tk, P, t)
% velocity of the leg each time falls in (the last leg for t = tk(end))
k = discretize(t, tk);
k(isnan(k) & t == tk(end)) = numel(tk) - 1;
dt = diff(tk);
V = diff(P, 1, 2) ./ dt;                                % [3 x K-1]
v = V(:, k);
end

function t = inRange(t, tk)
if any(t < tk(1) | t > tk(end))
    error('makeTrajectory:range', 'Time outside the trajectory [%g, %g] s.', tk(1), tk(end));
end
end

function tk = checkTimes(t)
tk = t(:).';
if numel(tk) < 2 || any(diff(tk) <= 0)
    error('makeTrajectory:time', 'Time needs at least two strictly increasing values.');
end
end

function need(opts, names, type)
for i = 1:numel(names)
    if isempty(opts.(names{i}))
        error('makeTrajectory:missing', 'Type ''%s'' needs ''%s''.', type, names{i});
    end
end
end
