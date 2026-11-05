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
  2. Custom: uniform motion as PosFcn; the default velocity (central difference) equals the true
     one (1e-9); frame checked.
  3. Input checks: unknown type, missing Lat, bad Frame -> refused with their ids.
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
% 2. Custom
% ---------------------------------------------------------------------------------
r0c = [1e6; 2e6; 3e6]; vc = [7; -3; 0.5];
tc = makeTrajectory('custom', 'Frame', 'SSB', 'PosFcn', @(t) r0c + vc * t);
ok2 = max(abs(tc.vel([0 1e3 -2e5]) - vc), [], 'all') < 1e-9 * norm(vc) * 1e3 && strcmp(tc.frame, 'SSB');
fprintf('2. custom (uniform motion): default velocity vs true %.1e m/s, frame %s: %s\n', ...
    max(abs(tc.vel([0 1e3 -2e5]) - vc), [], 'all'), tc.frame, passStr(ok2));
if ~ok2, fails{end+1} = 'custom'; end

% ---------------------------------------------------------------------------------
% 3. Input checks
% ---------------------------------------------------------------------------------
ids = {errId(@() makeTrajectory('boat', 'Lat', 1, 'Lon', 1)), 'makeTrajectory:type'; ...
       errId(@() makeTrajectory('fixed', 'Lon', 1)), 'makeTrajectory:missing'; ...
       errId(@() makeTrajectory('custom', 'Frame', 'XYZ', 'PosFcn', @(t) t)), 'makeTrajectory:frame'};
ok3 = all(cellfun(@strcmp, ids(:, 1), ids(:, 2)));
fprintf('3. input checks: %d of %d refused with the right id: %s\n', nnz(cellfun(@strcmp, ids(:, 1), ids(:, 2))), ...
    size(ids, 1), passStr(ok3));
if ~ok3, fails{end+1} = 'input checks'; end

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
