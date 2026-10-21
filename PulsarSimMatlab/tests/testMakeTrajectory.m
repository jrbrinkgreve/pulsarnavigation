function testMakeTrajectory()
%TESTMAKETRAJECTORY  Unit tests for makeTrajectory (block D1a: where the receiver is).
%{
Run from the PulsarSimMatlab folder: run('tests/testMakeTrajectory.m'), or
runAllTests("testMakeTrajectory"). ~1 s; no files. Truth from closed forms and an independent
inverse (ECEF -> geodetic by Bowring's iteration, not the code's forward formula).

  1. Fixed site: lat 0 / lon 0 -> [a 0 0], the north pole -> [0 0 b] (b = a(1-f)), 1e-6 m;
     52 N 5 E 30 m: the inverse returns lat, lon, h (1e-9 deg, 1e-6 m); the east / north / up
     vectors orthonormal (1e-12), up = d(position)/d(height) (1e-9); 'OffsetENU' [100 km 0 0]
     moves the point 100 km, perpendicular to north and up (1e-9); velocity 0.
  2. Ship track, linear: an eastward then a northward leg at 52 N (10 km per 1000 s): the
     waypoints themselves (= 'fixed' at those points, 1e-6 m), a leg's midpoint = the mean of
     its ends, velocity = (end - start) / duration per leg, ~10 m/s; a time outside the track
     is refused.
  3. Ship track, pchip: waypoints exact; velocity = central difference of the position
     (0.01 s, 1e-6 m/s); continuous at an interior waypoint (the linear track jumps there).
  4. Across 180 deg longitude (179.9 E -> 179.9 W at 20 N): the midpoint lies on 180 deg (1e-6
     deg), at 20 N (1e-4 deg: the chord midpoint is ~3e-5 deg higher) and below the ellipsoid
     by the chord sag (0-20 m; ~9 m for this 21 km leg).
  5. Barycentric: a circular orbit of 1 AU, one year, sampled hourly, 'spline': position at
     half-hours within 1 m of the circle, velocity within 1e-5 of 2 pi AU / yr, frame 'SSB'
     (pchip, second-order slopes, would be off by metres here).
  6. Custom (the entry point for a filter's state): uniform motion as PosFcn; the default
     velocity (central difference) equals the true one (1e-9); frame checked.
  7. Input checks: unknown type, missing Lat, non-increasing Time, bad Frame -> refused with
     their ids.
Errors at the end if any check fails. (21 Oct 2026, block D1a.)
%}

root = fullfile(fileparts(mfilename('fullpath')), '..');
addpath(fullfile(root, 'functions'));
fails = {};
a = 6378137; f = 1 / 298.257223563; b = a * (1 - f);

% ---------------------------------------------------------------------------------
% 1. Fixed site
% ---------------------------------------------------------------------------------
t3 = [0 1e6 -5e8];
r0 = makeTrajectory('fixed', 'Lat', 0, 'Lon', 0).pos(t3);
rp = makeTrajectory('fixed', 'Lat', 90, 'Lon', 0).pos(0);
s1 = makeTrajectory('fixed', 'Lat', 52, 'Lon', 5, 'Height', 30);
[la, lo, hh] = ecefToGeodetic(s1.ecef);
E = s1.enu;
up = makeTrajectory('fixed', 'Lat', 52, 'Lon', 5, 'Height', 31).ecef - s1.ecef;
s1e = makeTrajectory('fixed', 'Lat', 52, 'Lon', 5, 'Height', 30, 'OffsetENU', [100e3 0 0]);
d1 = s1e.ecef - s1.ecef;
ok1 = max(abs(r0 - [a; 0; 0]), [], 'all') < 1e-6 && max(abs(rp - [0; 0; b])) < 1e-6 && ...
      abs(la - 52) < 1e-9 && abs(lo - 5) < 1e-9 && abs(hh - 30) < 1e-6 && ...
      max(abs(E.' * E - eye(3)), [], 'all') < 1e-12 && norm(up - E(:, 3)) < 1e-9 && ...
      abs(norm(d1) - 100e3) < 1e-9 * 100e3 && abs(d1.' * E(:, 2)) < 1e-9 * 100e3 && ...
      abs(d1.' * E(:, 3)) < 1e-9 * 100e3 && isequal(s1.vel(t3), zeros(3, 3)) && strcmp(s1.frame, 'ECEF');
fprintf(['1. fixed: [a 0 0] %.1e m, pole %.1e m; 52 N 5 E 30 m back as %.10f / %.10f / %.6f m; ENU ' ...
         'orthonormal %.1e, up = d r / d h %.1e; 100 km east: %.6f km, . north %.1e, . up %.1e m: %s\n'], ...
    max(abs(r0(:, 1) - [a; 0; 0])), max(abs(rp - [0; 0; b])), la, lo, hh, ...
    max(abs(E.' * E - eye(3)), [], 'all'), norm(up - E(:, 3)), norm(d1) / 1e3, d1.' * E(:, 2), ...
    d1.' * E(:, 3), passStr(ok1));
if ~ok1, fails{end+1} = 'fixed site'; end

% ---------------------------------------------------------------------------------
% 2. Ship track, linear
% ---------------------------------------------------------------------------------
latK = 52 + [0 0 10e3 / 111.2e3]; lonK = 5 + [0 10e3 / (111.32e3 * cosd(52)) 10e3 / (111.32e3 * cosd(52))];
tK = [0 1000 2000];
tr = makeTrajectory('track', 'Time', tK, 'Lat', latK, 'Lon', lonK);
Pk = zeros(3, 3);
for k = 1:3, Pk(:, k) = makeTrajectory('fixed', 'Lat', latK(k), 'Lon', lonK(k)).ecef; end
mid = tr.pos(500);
v1 = tr.vel([100 900]); v2 = tr.vel(1500);
ok2 = max(abs(tr.pos(tK) - Pk), [], 'all') < 1e-6 && norm(mid - (Pk(:, 1) + Pk(:, 2)) / 2) < 1e-6 && ...
      max(abs(v1 - (Pk(:, 2) - Pk(:, 1)) / 1000), [], 'all') < 1e-9 && ...
      norm(v2 - (Pk(:, 3) - Pk(:, 2)) / 1000) < 1e-9 && abs(norm(v1(:, 1)) - 10) < 0.05 && ...
      strcmp(errId(@() tr.pos(2500)), 'makeTrajectory:range');
fprintf(['2. track, linear: waypoints %.1e m, leg midpoint %.1e m, leg velocities %.4f / %.4f m/s ' ...
         '(per leg exact), outside the track refused %d: %s\n'], max(abs(tr.pos(tK) - Pk), [], 'all'), ...
    norm(mid - (Pk(:, 1) + Pk(:, 2)) / 2), norm(v1(:, 1)), norm(v2), ...
    strcmp(errId(@() tr.pos(2500)), 'makeTrajectory:range'), passStr(ok2));
if ~ok2, fails{end+1} = 'track linear'; end

% ---------------------------------------------------------------------------------
% 3. Ship track, pchip
% ---------------------------------------------------------------------------------
tp = makeTrajectory('track', 'Time', tK, 'Lat', latK, 'Lon', lonK, 'Interp', 'pchip');
tq = [10 400 999 1001 1700];
vNum = (tp.pos(tq + 0.01) - tp.pos(tq - 0.01)) / 0.02;
jumpP = norm(tp.vel(1000 - 1e-6) - tp.vel(1000 + 1e-6));
jumpL = norm(tr.vel(1000 - 1e-6) - tr.vel(1000 + 1e-6));
ok3 = max(abs(tp.pos(tK) - Pk), [], 'all') < 1e-6 && max(abs(tp.vel(tq) - vNum), [], 'all') < 1e-6 && ...
      jumpP < 1e-6 && jumpL > 1;
fprintf(['3. track, pchip: waypoints %.1e m, velocity vs central difference %.1e m/s; velocity jump at ' ...
         'the turn %.1e m/s (linear %.2f m/s): %s\n'], max(abs(tp.pos(tK) - Pk), [], 'all'), ...
    max(abs(tp.vel(tq) - vNum), [], 'all'), jumpP, jumpL, passStr(ok3));
if ~ok3, fails{end+1} = 'track pchip'; end

% ---------------------------------------------------------------------------------
% 4. Across 180 deg longitude
% ---------------------------------------------------------------------------------
tw = makeTrajectory('track', 'Time', [0 1000], 'Lat', [20 20], 'Lon', [179.9 -179.9]);
[lw, ow, hw] = ecefToGeodetic(tw.pos(500));
ok4 = abs(abs(ow) - 180) < 1e-6 && abs(lw - 20) < 1e-4 && hw < 0 && hw > -20;
fprintf('4. across 180 deg: midpoint lon %.8f deg, lat %.8f deg, height %.2f m (chord): %s\n', ow, lw, hw, ...
    passStr(ok4));
if ~ok4, fails{end+1} = 'date line'; end

% ---------------------------------------------------------------------------------
% 5. Barycentric: circular orbit, 1 AU, 1 year, hourly samples, pchip
% ---------------------------------------------------------------------------------
AU = 149597870700; yr = 365.25 * 86400; w = 2*pi / yr;
tb = 0:3600:30 * 86400;
pb = makeTrajectory('barycentric', 'Time', tb, 'Pos', AU * [cos(w * tb); sin(w * tb); zeros(size(tb))], ...
    'Interp', 'spline');
th = tb(1:end-1) + 1800;
eP = max(vecnorm(pb.pos(th) - AU * [cos(w * th); sin(w * th); zeros(size(th))]));
eV = max(vecnorm(pb.vel(th) - AU * w * [-sin(w * th); cos(w * th); zeros(size(th))])) / (AU * w);
ok5 = eP < 1 && eV < 1e-5 && strcmp(pb.frame, 'SSB');
fprintf('5. barycentric orbit (spline, hourly): position %.3f m from the circle, velocity %.1e relative: %s\n', ...
    eP, eV, passStr(ok5));
if ~ok5, fails{end+1} = 'barycentric'; end

% ---------------------------------------------------------------------------------
% 6. Custom: the filter's entry point
% ---------------------------------------------------------------------------------
r0c = [1e6; 2e6; 3e6]; vc = [7; -3; 0.5];
tc = makeTrajectory('custom', 'Frame', 'SSB', 'PosFcn', @(t) r0c + vc * t);
ok6 = max(abs(tc.vel([0 1e3 -2e5]) - vc), [], 'all') < 1e-9 * norm(vc) * 1e3 && strcmp(tc.frame, 'SSB');
fprintf('6. custom (uniform motion): default velocity vs true %.1e m/s, frame %s: %s\n', ...
    max(abs(tc.vel([0 1e3 -2e5]) - vc), [], 'all'), tc.frame, passStr(ok6));
if ~ok6, fails{end+1} = 'custom'; end

% ---------------------------------------------------------------------------------
% 7. Input checks
% ---------------------------------------------------------------------------------
ids = {errId(@() makeTrajectory('boat', 'Lat', 1, 'Lon', 1)), 'makeTrajectory:type'; ...
       errId(@() makeTrajectory('fixed', 'Lon', 1)), 'makeTrajectory:missing'; ...
       errId(@() makeTrajectory('track', 'Time', [0 0], 'Lat', [1 1], 'Lon', [1 1])), 'makeTrajectory:time'; ...
       errId(@() makeTrajectory('custom', 'Frame', 'XYZ', 'PosFcn', @(t) t)), 'makeTrajectory:frame'};
ok7 = all(cellfun(@strcmp, ids(:, 1), ids(:, 2)));
fprintf('7. input checks: %d of %d refused with the right id: %s\n', nnz(cellfun(@strcmp, ids(:, 1), ids(:, 2))), ...
    size(ids, 1), passStr(ok7));
if ~ok7, fails{end+1} = 'input checks'; end

if isempty(fails)
    fprintf('testMakeTrajectory: ALL PASSED\n');
else
    error('testMakeTrajectory:failed', 'FAILED: %s', strjoin(fails, ', '));
end
end


% =================================================================================
function [lat, lon, h] = ecefToGeodetic(r)
% ECEF [m] -> WGS84 geodetic [deg, deg, m] by Bowring's iteration (independent of the code).
a = 6378137; f = 1 / 298.257223563; e2 = f * (2 - f);
x = r(1); y = r(2); z = r(3);
lon = atan2d(y, x);
p = hypot(x, y);
lat = atan2(z, p * (1 - e2));
for it = 1:10
    N = a / sqrt(1 - e2 * sin(lat)^2);
    h = p / cos(lat) - N;
    lat = atan2(z, p * (1 - e2 * N / (N + h)));
end
N = a / sqrt(1 - e2 * sin(lat)^2);
h = p / cos(lat) - N;
lat = rad2deg(lat);
end

function id = errId(f)
% The identifier of the error f() raises ('' if none).
id = '';
try
    f();
catch err
    id = err.identifier;
end
end

function s = passStr(p)
if p, s = 'PASS'; else, s = 'FAIL'; end
end
