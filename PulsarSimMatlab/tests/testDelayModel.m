function testDelayModel()
%TESTDELAYMODEL  Unit tests for delayModel (block D1b: receiver -> solar-system barycentre).
%{
Run from the PulsarSimMatlab folder: run('tests/testDelayModel.m'), or
runAllTests("testDelayModel"). ~2 s; no files. Each term against a closed form computed here,
independently of the function's code paths (the Earth rotation angle and WGS84 by their
standard formulas). Times in 2026 (s since J2000 TT).

  1. Earth's orbit: without Moon and planets the EMB's distance from the Sun (= the SSB) over a
     year spans a(1 - e) ... a(1 + e) (1e-6 AU); with planets the Sun's offset from the SSB is
     0.002-0.008 AU (Jupiter ~0.005, Saturn ~0.003); with the Moon the Earth is 4424-4932 km from
     the EMB (0.01215 x 364,096-405,906 km).
  2. Annual Roemer, exact: no Moon / planets, the geocentre, a pulsar in the ecliptic towards the
     EMB's perihelion: r.n = r cos(true anomaly), so peak-to-peak = a(1 - e) + a(1 + e) = 2a/c
     (998.0 s), with e and a at the perihelion and aphelion dates (e falls 4.4e-5 per century:
     half a year apart that is 1.1e-4 s; extremes refined on a 10 s grid; 1e-5 s).
     Shown: B0329+54 (ecliptic latitude 34.3 deg: ~2a cos(beta)/c).
  3. Daily Roemer: a site at 52 N 5 E minus the geocentre = (p cos(dec) cos(theta + lon - ra)
     + z sin(dec)) / c, theta the Earth rotation angle, p and z the site's distance from the
     axis and from the equator (1e-12 s over two days); its maximum at the pulsar's transit
     (theta + lon = ra; 1 s). B0329+54: amplitude 7.6 ms.
  4. The deliverable's measurement equation: the same site 100 km east (OffsetENU) minus the
     assumed one = (R(theta) dr . n) / c (roemer, 1e-12 s): a daily sine of amplitude
     100 km cos(dec) / c (193 us for B0329+54; 1e-5, the 1-minute grid) with zero mean (the east vector has no polar
     part); the full model adds only the Einstein site term v_E . dr / c^2 (< 40 ns).
  5. Ships: a ship standing still (two equal waypoints) = the fixed site (1e-12 s, all terms); a
     ship sailing east at 10 m/s: d(roemer)/dt (numerical) = vObs . n / c (1e-6 relative).
  6. A probe ('SSB', no planets: the Sun at the SSB) at 1 AU on the x axis: roemer = AU n_x / c
     (1e-12 s), einstein 0, shapiro = 2 T_sun ln(1 + n_x) (1e-15 s) - 0 at n_x = 0, -45 us at
     n_x = -0.99 (the Sun nearly in front of the pulsar).
  7. Einstein: over a year peak-to-peak 2 x 1.657 ms (5e-5 s, the 14 us harmonic and the site
     term included); the site term at 52 N between 0.5 and 2.1 us (|v_E| |r_site| / c^2 = 2.1 us).
  8. Shapiro for B0329+54 over a year between -20 and +8 us (expected -17 ... +6 us).
  9. A trajectory with an unknown frame is refused.
 10. JPL's DE440 behind the same interface ('Ephemeris', block D5b; data/ephem/de440s.bsp): the
     geocentre is the file's at TDB = TT + eph.tdbMinusTT (1 cm: the time's last
     bit, 1.2e-7 s at 8e8 s, is 4 mm of the Earth's motion; velocity 1e-5 m/s); the analytic
     orbit against it over a year: the Earth relative to the Sun within 2e4 km, relative to the
     SSB within 4e5 km (the analytic Sun leaves out Uranus and Neptune); the navigation's
     measurement (the site 100 km east minus the site, all terms, two days) is the same with both
     to 1e-11 s; Shapiro for B0329+54 within -20 ... +8 us; the Einstein term at the geocentre is
     the file's TDB - TT (D5c, 1e-15 s); a bad 'Ephemeris' refused.
Errors at the end if any check fails. (21 Oct 2026, block D1b; check 10 22 Oct, D5b.)
%}

root = fullfile(fileparts(mfilename('fullpath')), '..');
addpath(fullfile(root, 'functions'));
fails = {};
c = 299792458; AU = 149597870700; Tsun = 1.32712440018e20 / c^3;
t0 = 9496.5 * 86400;                                  % 2026-01-01 00:00 TT [s since J2000]
ty = t0 + (0:3600:366 * 86400);                       % a year, hourly
td = t0 + 100 * 86400 + (0:60:2 * 86400);             % two days, per minute
geo = makeTrajectory('custom', 'Frame', 'ECEF', 'PosFcn', @(t) zeros(3, numel(t)), ...
    'VelFcn', @(t) zeros(3, numel(t)));
b0329 = struct('raDeg', 53.24738, 'decDeg', 54.57878);
site = makeTrajectory('fixed', 'Lat', 52, 'Lon', 5);

% ---------------------------------------------------------------------------------
% 1. Earth's orbit
% ---------------------------------------------------------------------------------
T = (t0 + 0.5 * 365.25 * 86400) / (36525 * 86400);
aT = 1.00000261 + 0.00000562 * T; eT = 0.01671123 - 0.00004392 * T;
d0 = delayModel(ty, geo, b0329, 'Moon', false, 'Planets', false);
rr = vecnorm(d0.rEarth) / AU;
dP = delayModel(ty, geo, b0329, 'Moon', false);
rs = vecnorm(dP.rSun) / AU;
dM = delayModel(ty, geo, b0329, 'Planets', false);
rm = vecnorm(dM.rEarth - d0.rEarth) / 1e3;
ok1 = abs(min(rr) - aT * (1 - eT)) < 1e-6 && abs(max(rr) - aT * (1 + eT)) < 1e-6 && ...
      all(rs > 0.002 & rs < 0.008) && all(rm > 4424 - 1 & rm < 4932 + 1);
fprintf(['1. orbit: EMB-Sun %.7f ... %.7f AU (a(1-e) %.7f, a(1+e) %.7f); Sun-SSB %.4f ... %.4f AU; ' ...
         'Earth-EMB %.0f ... %.0f km: %s\n'], min(rr), max(rr), aT * (1 - eT), aT * (1 + eT), min(rs), ...
    max(rs), min(rm), max(rm), passStr(ok1));
if ~ok1, fails{end+1} = 'orbit'; end

% ---------------------------------------------------------------------------------
% 2. Annual Roemer, exact
% ---------------------------------------------------------------------------------
wbar = 102.93768193 + 0.32327364 * T;                 % longitude of perihelion [deg]
psrP = eclipticToRaDec(wbar, 0);
dR = delayModel(ty, geo, psrP, 'Moon', false, 'Planets', false, 'Terms', "roemer");
[~, iP] = max(dR.roemer); [~, iA] = min(dR.roemer);
% refine the extremes on a 10 s grid (hourly samples miss them by up to 1 - cos(3.7e-4 rad): 34 us)
fP = ty(iP) + (-7200:10:7200); fA = ty(iA) + (-7200:10:7200);
[rMax, jP] = max(delayModel(fP, geo, psrP, 'Moon', false, 'Planets', false, 'Terms', "roemer").roemer);
[rMin, jA] = min(delayModel(fA, geo, psrP, 'Moon', false, 'Planets', false, 'Terms', "roemer").roemer);
pp = rMax - rMin;
Tp = fP(jP) / (36525 * 86400); Ta = fA(jA) / (36525 * 86400);     % e and a at perihelion, aphelion
ppT = ((1.00000261 + 0.00000562 * Tp) * (1 - (0.01671123 - 0.00004392 * Tp)) + ...
       (1.00000261 + 0.00000562 * Ta) * (1 + (0.01671123 - 0.00004392 * Ta))) * AU / c;
dB = delayModel(ty, geo, b0329, 'Moon', false, 'Planets', false, 'Terms', "roemer");
ok2 = abs(pp - ppT) < 1e-5;
fprintf(['2. annual Roemer towards perihelion: peak-to-peak %.6f s (a(1-e) + a(1+e) at the two dates %.6f s; ' ...
         '2a/c %.5f s); B0329+54 %.2f s (2a cos(beta)/c %.2f s): %s\n'], pp, ppT, 2 * aT * AU / c, ...
    max(dB.roemer) - min(dB.roemer), 2 * aT * AU / c * cosd(34.3), passStr(ok2));
if ~ok2, fails{end+1} = 'annual Roemer'; end

% ---------------------------------------------------------------------------------
% 3. Daily Roemer (Earth rotation) and transit
% ---------------------------------------------------------------------------------
[p52, z52] = siteAxes(52, 0);
ra = deg2rad(b0329.raDeg); de = deg2rad(b0329.decDeg);
ds = delayModel(td, site, b0329, 'Terms', "roemer");
dg = delayModel(td, geo, b0329, 'Terms', "roemer");
th = era(td);
expD = (p52 * cos(de) * cos(th + deg2rad(5) - ra) + z52 * sin(de)) / c;
e3 = max(abs(ds.roemer - dg.roemer - expD));
[~, im] = max(ds.roemer(1:1441) - dg.roemer(1:1441));      % first day
tTr = transitTime(td(im), deg2rad(5), ra);
ok3 = e3 < 1e-12 && abs(td(im) - tTr) <= 60;
fprintf(['3. daily Roemer, 52 N: site - geocentre vs closed form %.1e s; amplitude %.3f ms (p cos(dec)/c); ' ...
         'maximum %.0f s from the transit (grid 60 s): %s\n'], e3, p52 * cos(de) / c * 1e3, td(im) - tTr, ...
    passStr(ok3));
if ~ok3, fails{end+1} = 'daily Roemer'; end

% ---------------------------------------------------------------------------------
% 4. The deliverable's measurement equation: the site 100 km east
% ---------------------------------------------------------------------------------
siteE = makeTrajectory('fixed', 'Lat', 52, 'Lon', 5, 'OffsetENU', [100e3 0 0]);
dr = siteE.ecef - site.ecef;
dRe = delayModel(td, siteE, b0329, 'Terms', "roemer").roemer - ds.roemer;
expE = ((cos(th) * dr(1) - sin(th) * dr(2)) * cos(de) * cos(ra) + (sin(th) * dr(1) + cos(th) * dr(2)) ...
        * cos(de) * sin(ra) + dr(3) * sin(de)) / c;
dAll = delayModel(td, siteE, b0329).total - delayModel(td, site, b0329).total;
ampE = (max(dRe) - min(dRe)) / 2;
ok4 = max(abs(dRe - expE)) < 1e-12 && abs(ampE / (100e3 * cos(de) / c) - 1) < 1e-5 && ...   % 1-min grid
      abs(mean(dRe)) < 0.01 * ampE && max(abs(dAll - dRe)) < 4e-8;
fprintf(['4. site 100 km east: Roemer change vs (R dr . n)/c %.1e s; amplitude %.2f us (100 km cos(dec)/c ' ...
         '%.2f us), mean %.2f us; full model - Roemer %.1f ns (Einstein site term): %s\n'], ...
    max(abs(dRe - expE)), ampE * 1e6, 100e3 * cos(de) / c * 1e6, mean(dRe) * 1e6, ...
    max(abs(dAll - dRe)) * 1e9, passStr(ok4));
if ~ok4, fails{end+1} = 'site offset'; end

% ---------------------------------------------------------------------------------
% 5. Ships
% ---------------------------------------------------------------------------------
still = makeTrajectory('track', 'Time', [td(1) td(end)], 'Lat', [52 52], 'Lon', [5 5]);
eS = max(abs(delayModel(td, still, b0329).total - delayModel(td, site, b0329).total));
dlon = 10 * (td(end) - td(1)) / (111.32e3 * cosd(52));
ship = makeTrajectory('track', 'Time', [td(1) td(end)], 'Lat', [52 52], 'Lon', [5 5 + dlon]);
tm = td(2:end-1);
dsh = delayModel(tm, ship, b0329, 'Terms', "roemer");
rate = (delayModel(tm + 0.5, ship, b0329, 'Terms', "roemer").roemer - ...
        delayModel(tm - 0.5, ship, b0329, 'Terms', "roemer").roemer);
eV = max(abs(rate - dsh.nHat.' * dsh.vObs / c)) / max(abs(dsh.nHat.' * dsh.vObs / c));
ok5 = eS < 1e-12 && eV < 1e-6;
fprintf('5. ships: standing still = fixed site %.1e s; sailing east, d(roemer)/dt vs vObs.n/c %.1e relative: %s\n', ...
    eS, eV, passStr(ok5));
if ~ok5, fails{end+1} = 'ships'; end

% ---------------------------------------------------------------------------------
% 6. A probe in the barycentric frame (Sun at the SSB)
% ---------------------------------------------------------------------------------
probe = makeTrajectory('custom', 'Frame', 'SSB', 'PosFcn', @(t) repmat([AU; 0; 0], 1, numel(t)), ...
    'VelFcn', @(t) zeros(3, numel(t)));
nx = [0 0.5 -0.99];
e6 = 0; sh = zeros(size(nx));
for k = 1:numel(nx)
    ps = struct('raDeg', acosd(nx(k)), 'decDeg', 0);             % n = [nx, sqrt(1-nx^2), 0]
    dp = delayModel(t0, probe, ps, 'Planets', false);
    sh(k) = dp.shapiro;
    e6 = max([e6, abs(dp.roemer - AU * nx(k) / c) / 1e-12, abs(dp.shapiro - 2 * Tsun * log(1 + nx(k))) / 1e-15, ...
              abs(dp.einstein) * 1e20]);
end
ok6 = e6 < 1;
fprintf(['6. probe at 1 AU: roemer, shapiro (2 T ln(1 + n_x)) and einstein (0) to their tolerances %d; ' ...
         'shapiro at n_x = 0 / 0.5 / -0.99: %.2f / %.2f / %.2f us: %s\n'], e6 < 1, sh * 1e6, passStr(ok6));
if ~ok6, fails{end+1} = 'probe'; end

% ---------------------------------------------------------------------------------
% 7. Einstein
% ---------------------------------------------------------------------------------
dE = delayModel(ty, site, b0329, 'Terms', "einstein");
dE0 = delayModel(ty, geo, b0329, 'Terms', "einstein");
siteT = dE.einstein - dE0.einstein;
ok7 = abs((max(dE0.einstein) - min(dE0.einstein)) - 2 * 1.657e-3) < 5e-5 && ...
      max(abs(siteT)) > 0.5e-6 && max(abs(siteT)) < 2.1e-6;
fprintf('7. Einstein: peak-to-peak %.4f ms (2 x 1.657); site term up to %.2f us: %s\n', ...
    (max(dE0.einstein) - min(dE0.einstein)) * 1e3, max(abs(siteT)) * 1e6, passStr(ok7));
if ~ok7, fails{end+1} = 'Einstein'; end

% ---------------------------------------------------------------------------------
% 8. Shapiro over a year (B0329+54)
% ---------------------------------------------------------------------------------
dS = delayModel(ty, site, b0329, 'Terms', "shapiro");
ok8 = min(dS.shapiro) > -20e-6 && max(dS.shapiro) < 8e-6 && min(dS.shapiro) < -10e-6;
fprintf('8. Shapiro, B0329+54, a year: %.2f ... %.2f us (expected -17 ... +6): %s\n', ...
    min(dS.shapiro) * 1e6, max(dS.shapiro) * 1e6, passStr(ok8));
if ~ok8, fails{end+1} = 'Shapiro'; end

% ---------------------------------------------------------------------------------
% 9. Input check
% ---------------------------------------------------------------------------------
bad = site; bad.frame = 'GALACTIC';
ok9 = strcmp(errId(@() delayModel(t0, bad, b0329)), 'delayModel:frame');
fprintf('9. unknown frame refused %d: %s\n', ok9, passStr(ok9));
if ~ok9, fails{end+1} = 'input check'; end

% 10. JPL's DE440
% ---------------------------------------------------------------------------------
eph = loadEphemeris(fullfile(root, 'data', 'ephem', 'de440s.bsp'));
dJ = delayModel(ty, geo, b0329, 'Ephemeris', eph);
dA = delayModel(ty, geo, b0329);
[rF, vF] = eph.earth(ty + eph.tdbMinusTT(ty));
eEin = max(abs(dJ.einstein - eph.tdbMinusTT(ty)));   % geo: no site term
eF = max(vecnorm(dJ.rEarth - rF)); eV = max(vecnorm(dJ.vEarth - vF));
eHel = max(vecnorm((dA.rEarth - dA.rSun) - (dJ.rEarth - dJ.rSun)));
eBar = max(vecnorm(dA.rEarth - dJ.rEarth));
navJ = delayModel(td, siteE, b0329, 'Ephemeris', eph).total - delayModel(td, site, b0329, 'Ephemeris', eph).total;
eNav = max(abs(navJ - dAll));
shJ = delayModel(ty, site, b0329, 'Ephemeris', eph, 'Terms', "shapiro").shapiro;
ok10 = eF < 1e-2 && eV < 1e-5 && eEin < 1e-15 && eHel < 2e7 && eBar < 4e8 && eNav < 1e-11 && min(shJ) > -20e-6 && max(shJ) < 8e-6 && ...
       strcmp(errId(@() delayModel(t0, site, b0329, 'Ephemeris', struct('x', 1))), 'delayModel:ephemeris');
fprintf(['10. DE440: geocentre = the file''s to %.1e m, %.1e m/s; the analytic Earth off by %.0f km from the Sun, %.0f km from the SSB; ' ...
         'the 100 km measurement the same to %.1e s; Shapiro %.1f ... %.1f us; Einstein = the file''s TDB - TT to %.0e s: %s\n'], ...
    eF, eV, eHel / 1e3, eBar / 1e3, eNav, 1e6 * min(shJ), 1e6 * max(shJ), eEin, passStr(ok10));
if ~ok10, fails{end+1} = 'DE440'; end

if isempty(fails)
    fprintf('testDelayModel: ALL PASSED\n');
else
    error('testDelayModel:failed', 'FAILED: %s', strjoin(fails, ', '));
end
end


% =================================================================================
function psr = eclipticToRaDec(lamDeg, betDeg)
% Ecliptic J2000 longitude / latitude -> RA / Dec [deg] (obliquity 23.439281 deg).
ep = 23.439281;
x = cosd(betDeg) * cosd(lamDeg);
y = cosd(ep) * cosd(betDeg) * sind(lamDeg) - sind(ep) * sind(betDeg);
z = sind(ep) * cosd(betDeg) * sind(lamDeg) + cosd(ep) * sind(betDeg);
psr = struct('raDeg', mod(atan2d(y, x), 360), 'decDeg', asind(z));
end

function [p, z] = siteAxes(lat, h)
% WGS84: the site's distance from the rotation axis and from the equatorial plane [m].
a = 6378137; f = 1 / 298.257223563; e2 = f * (2 - f);
N = a / sqrt(1 - e2 * sind(lat)^2);
p = (N + h) * cosd(lat); z = (N * (1 - e2) + h) * sind(lat);
end

function th = era(t)
% Earth rotation angle (IERS 2003), UT1 = TT - 69.184 s.
th = 2*pi * mod(0.7790572732640 + 1.00273781191135448 * (t - 69.184) / 86400, 1);
end

function tT = transitTime(tGuess, lon, ra)
% The time near tGuess when era + lon = ra (mod 2 pi): the pulsar on the local meridian.
w = 2*pi * 1.00273781191135448 / 86400;
dphi = mod(ra - lon - era(tGuess) + pi, 2*pi) - pi;
tT = tGuess + dphi / w;
end

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
