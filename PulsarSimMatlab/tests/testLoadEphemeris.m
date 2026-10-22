function testLoadEphemeris()
%TESTLOADEPHEMERIS  Unit tests for loadEphemeris (block D5a: JPL's DE440 from its binary SPK file).
%{
Run from the PulsarSimMatlab folder: run('tests/testLoadEphemeris.m'), or runAllTests("testLoadEphemeris").
~2 s. Needs data/ephem/de440s.bsp (naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/).

  1. JPL Horizons (independent: JPL's own service, DE441, queried 22 Oct 2026): the geocentre
     relative to the SSB, ICRF, at 2000-01-01 12:00, 2026-03-01 00:00 and 2045-01-19 18:00 TDB.
     Position within 10 m, velocity within 1 mm/s (DE440 and DE441 agree far better for the
     Earth today; Horizons prints 16 digits).
  2. The velocity is the derivative of the position: a central difference over +-60 s at 20 times
     agrees within 1e-5 m/s (Earth: 30 km/s).
  3. The Earth and the Moon about their barycentre: r_Earth = -r_Moon / EMRAT at every time (the
     two are separate fits in the file), with DE440's Earth/Moon mass ratio 81.30057 (1e-6).
  4. TDB - TT at the geocentre (block D5c), integrated from the file, against the USNO series
     (Kaplan 2005, USNO Circular 179 eq. 2.6, stated good to ~10 us; typed into this test
     independently of loadEphemeris' copy) over 2000-2050, daily: rms below 5 us, largest
     difference below 10 us (the constant is the series' mean by construction; the shape is the
     test). Shown: delayModel's analytic two terms against the same series.
  5. Bad inputs refused: a missing file ('loadEphemeris:file'), a file that is not SPK
     ('loadEphemeris:format'), a body not in the file ('loadEphemeris:body'), a time outside the
     file ('loadEphemeris:range').
Errors at the end if any check fails. (22 Oct 2026, block D5a; check 4 for D5c.)
%}

root = fullfile(fileparts(mfilename('fullpath')), '..');
addpath(fullfile(root, 'functions'));
fails = {};
file = fullfile(root, 'data', 'ephem', 'de440s.bsp');
eph = loadEphemeris(file);

% 1. JPL Horizons: JD (TDB), x y z [km], vx vy vz [km/s]
H = [2451545.0,  -2.756674048281145E+07,  1.323613811535491E+08,  5.741865328625385E+07, ...
                 -2.978494749851088E+01, -5.029753814928081E+00, -2.180645069035755E+00;
     2461100.5,  -1.397051888927500E+08,  4.564847304932671E+07,  1.980848586205852E+07, ...
                 -1.062918609124083E+01, -2.578446268581494E+01, -1.117620362141082E+01;
     2468000.25, -7.276672036836250E+07,  1.182115469657561E+08,  5.123886500423498E+07, ...
                 -2.644499451883116E+01, -1.355332237866367E+01, -5.875016466349400E+00];
t = (H(:, 1).' - 2451545.0) * 86400;
[r, v] = eph.earth(t);
dr = max(vecnorm(r - 1e3 * H(:, 2:4).')); dv = max(vecnorm(v - 1e3 * H(:, 5:7).'));
ok1 = dr < 10 && dv < 1e-3;
fprintf('1. the geocentre against JPL Horizons at three dates: largest difference %.3f m, %.2e m/s: %s\n', dr, dv, passStr(ok1));
if ~ok1, fails{end+1} = 'Horizons'; end

% 2. Velocity = derivative of the position
tq = 9555.5 * 86400 + linspace(0, 365.25 * 86400, 20) + 1234.5;
[~, vq] = eph.earth(tq);
[rP, ~] = eph.earth(tq + 60); [rM, ~] = eph.earth(tq - 60);
d2 = max(vecnorm(vq - (rP - rM) / 120));
ok2 = d2 < 1e-5;
fprintf('2. velocity = central difference of the position: largest difference %.1e m/s: %s\n', d2, passStr(ok2));
if ~ok2, fails{end+1} = 'velocity'; end

% 3. Earth and Moon about their barycentre
rE = eph.state(399, 3, tq); rMn = eph.state(301, 3, tq);
q = -dot(rE, rMn) ./ dot(rMn, rMn);                  % = 1 / EMRAT
perp = max(vecnorm(rE + q .* rMn) ./ vecnorm(rE));
ok3 = abs(mean(q) * 81.30057 - 1) < 1e-6 && std(q) / mean(q) < 1e-9 && perp < 1e-9;
fprintf('3. Earth = -Moon / EMRAT about the EMB: EMRAT %.6f (spread %.1e, off the line %.1e): %s\n', 1 / mean(q), ...
    std(q) / mean(q), perp, passStr(ok3));
if ~ok3, fails{end+1} = 'EMRAT'; end

% 4. TDB - TT against the USNO series
t5 = 0:86400:50 * 365.25 * 86400; T = t5 / (36525 * 86400);
usno = 0.001657 * sin(628.3076 * T + 6.2401) + 0.000022 * sin(575.3385 * T + 4.2970) + ...
       0.000014 * sin(1256.6152 * T + 6.1969) + 0.000005 * sin(606.9777 * T + 4.0212) + ...
       0.000005 * sin(52.9691 * T + 0.4444) + 0.000002 * sin(21.3299 * T + 5.5431) + 0.000010 * T .* sin(628.3076 * T + 4.2490);
e5 = eph.tdbMinusTT(t5) - usno;
g = deg2rad(357.53 + 0.98560028 * (t5 / 86400)); e2 = 1.657e-3 * sin(g) + 1.4e-5 * sin(2*g) - usno; e2 = e2 - mean(e2);
ok5 = rms(e5) < 5e-6 && max(abs(e5)) < 10e-6;
fprintf('4. TDB - TT from the file against the USNO series, 2000-2050: rms %.2f us, largest %.2f us (the analytic two terms: rms %.1f, largest %.1f us): %s\n', ...
    1e6 * rms(e5), 1e6 * max(abs(e5)), 1e6 * rms(e2), 1e6 * max(abs(e2)), passStr(ok5));
if ~ok5, fails{end+1} = 'TDB - TT'; end

% 5. Bad inputs
bad = [tempname, '.bsp']; fb = fopen(bad, 'w'); fwrite(fb, zeros(1, 2048, 'uint8')); fclose(fb);
ok4 = strcmp(errId(@() loadEphemeris(fullfile(root, 'data', 'ephem', 'none.bsp'))), 'loadEphemeris:file') && ...
      strcmp(errId(@() loadEphemeris(bad)), 'loadEphemeris:format') && ...
      strcmp(errId(@() eph.state(12345, 0, 0)), 'loadEphemeris:body') && ...
      strcmp(errId(@() eph.earth(eph.span(2) + 86400)), 'loadEphemeris:range');
delete(bad);
fprintf('5. missing file, not SPK, unknown body, time outside refused %d: %s\n', ok4, passStr(ok4));
if ~ok4, fails{end+1} = 'input checks'; end

if isempty(fails)
    fprintf('testLoadEphemeris: ALL PASSED\n');
else
    error('testLoadEphemeris:failed', 'FAILED: %s', strjoin(fails, ', '));
end
end


% =================================================================================
function id = errId(f)
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
