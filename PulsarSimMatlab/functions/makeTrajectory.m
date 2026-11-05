function traj = makeTrajectory(type, opts)
%MAKETRAJECTORY  Where the receiver is: a trajectory for the delay model (block D1a).
%{
The phase prediction and the navigation fit need the receiver's position at any time. This builds
it as one struct with a position and a velocity function. The array stands still on the Earth, so
in the pipeline this is one point ('fixed'). delayModel turns it into delays; the simulator's truth
and the observer's prediction are two trajectories. A moving platform can be added later without
changing the delay model or the predictor (old/moving_platform keeps the earlier ship and probe types).

  traj = makeTrajectory('fixed', 'Lat', 52, 'Lon', 5)
  traj = makeTrajectory('fixed', 'Lat', 52, 'Lon', 5, 'OffsetENU', [100e3 0 0])   % 100 km east
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
  'custom'       'Frame' ('ECEF' or 'SSB'), 'PosFcn' @(t) -> [3 x N], 'VelFcn' @(t) -> [3 x N]
                 (default: central difference of PosFcn over 'VelStep' s, default 1). For the
                 tests: the geocentre, and probes with a closed-form truth.

Output traj: type, frame, pos (@(t) -> [3 x N] m), vel (@(t) -> [3 x N] m/s), params (the inputs),
and for 'fixed' also ecef (the point, [3 x 1]) and enu (the local east, north, up unit vectors,
columns).
%}

arguments
    type {mustBeTextScalar}
    opts.Lat       double = []
    opts.Lon       double = []
    opts.Height    double = 0
    opts.OffsetENU (1,3) double = [0 0 0]
    opts.Frame     {mustBeTextScalar} = 'ECEF'
    opts.PosFcn    = []
    opts.VelFcn    = []
    opts.VelStep   (1,1) double {mustBePositive} = 1
end

type = lower(char(type));
traj = struct('type', type, 'frame', '', 'pos', [], 'vel', [], 'params', opts);
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
        error('makeTrajectory:type', 'Type must be ''fixed'' or ''custom''.');
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

function need(opts, names, type)
for i = 1:numel(names)
    if isempty(opts.(names{i}))
        error('makeTrajectory:missing', 'Type ''%s'' needs ''%s''.', type, names{i});
    end
end
end
