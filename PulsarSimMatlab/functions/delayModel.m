function d = delayModel(t, traj, psr, opts)
%DELAYMODEL  Delay from the receiver to the solar-system barycentre (block D1b).
%{
A pulsar's rotation phase is defined at the solar-system barycentre (SSB). This gives, for a
receiver on a trajectory (makeTrajectory) and a pulsar direction, the time to ADD to the
receiver's arrival time to get the arrival time at the SSB:

  delta(t) = t_SSB - t_obs = roemer + curvature + einstein + shapiro

  d = delayModel(t, traj, psr)
  d = delayModel(t, makeTrajectory('fixed', 'Lat', 52, 'Lon', 5), struct('raDeg', 53.2474, 'decDeg', 54.5788))

Inputs:
  t     [1 x N] s since J2000.0 (2000-01-01 12:00 TT), as makeTrajectory.
  traj  makeTrajectory struct: 'ECEF' (on the Earth: fixed site, ship) or 'SSB' (a probe).
  psr   struct: raDeg, decDeg (J2000 / ICRS direction), optional distance [m] (wave-front
        curvature term; default: none).

Name-value options:
  'Terms'    which terms to include, any of "roemer" "curvature" "einstein" "shapiro"
             (default all; curvature only with psr.distance).
  'Moon'     true (default): the Earth's offset from the Earth-Moon barycentre (EMB).
  'Planets'  true (default): the Sun's offset from the SSB (Jupiter, Saturn). Moon and Planets
             false leave the EMB on its Keplerian orbit around a Sun at the SSB (for tests).
  'Ephemeris' [] (default): the analytic Earth below. A loadEphemeris struct (JPL's DE440): the
             geocentre, its velocity and the Sun from the file, read at TDB = TT + the Einstein
             term (block D5b), and that term (TDB - TT) integrated from the file (D5c, to ~2 us,
             against the analytic two terms' ~60 us); 'Moon' and 'Planets' then do not apply.

The terms (signs: contributions to t_SSB - t_obs):
  roemer     r_obs . n / c: the light-travel time across the observer's offset from the SSB
             (n: unit vector to the pulsar). Up to +-499 s (Earth's orbit), +-21 ms (rotation).
  curvature  -(|r_obs|^2 - (r_obs . n)^2) / (2 c D): the wave front is a sphere, not a plane
             (B0329+54 at ~1 kpc: ~1 us annual).
  einstein   TDB - TT for a clock on the Earth: 1.657 ms sin(g) + 14 us sin(2g) (g: the Earth's
             mean anomaly; the clock runs at a varying rate along the eccentric orbit) plus the
             site term v_E . r_site / c^2 (~1 us, daily). 0 for an 'SSB' trajectory (a probe's
             own clock: later).
  shapiro    +2 T_sun ln((|r| + r . n) / AU), r the observer relative to the Sun, T_sun = GM/c^3
             = 4.925 us: light from the pulsar passing the Sun is delayed by -2 T_sun ln(|r| + r.n)
             (+ a constant); the contribution is minus that delay - strongly negative when the
             Sun is near the line of sight (B0329+54: -17 ... +6 us over a year).

Earth model (analytic; decision 1 of block D, 21 Oct 2026: DE440 behind the same interface at
the end of D, now 'Ephemeris'): the EMB on a Keplerian orbit with JPL's approximate elements and their rates
(Standish, "Keplerian Elements for Approximate Positions of the Major Planets", 1800-2050,
~20" for the EMB); the Sun around the SSB from Jupiter and Saturn on the same kind of orbits
(the Sun moves ~0.007 AU, +-2.5 s and +-1.4 s); the Earth around the EMB from a low-precision
lunar series (+-4671 km, +-15.6 ms; ~0.3 deg). Ecliptic J2000 -> equator J2000 (obliquity
23.439281 deg). Earth rotation: the Earth rotation angle from UT1 ~ TT - 69.184 s, without
precession, nutation and polar motion (they matter against real data, not between the
simulator's truth and prediction, which share this model). Absolute accuracy (against DE440, 22
Oct 2026): the Earth relative to the Sun within ~1.5e4 km, but the Sun relative to the SSB is off
by up to ~3e5 km, as Uranus and Neptune are left out; so the Earth is up to 3.2e5 km from where
DE440 puts it, and B0329+54's Roemer delay is off by up to 0.66 s, 71 ms of it varying within two
years. The differences that the navigation uses are exact within the model.

Output d: total, roemer, curvature, einstein, shapiro [1 x N] (s); rObs, vObs [3 x N] (the
receiver relative to the SSB, equatorial J2000, m and m/s); rEarth, vEarth (geocentre), rSun
(Sun) [3 x N]; nHat [3 x 1]; frame of the trajectory; model (text).
%}

arguments
    t (1,:) double
    traj struct
    psr struct
    opts.Terms (1,:) string = ["roemer" "curvature" "einstein" "shapiro"]
    opts.Moon (1,1) logical = true
    opts.Planets (1,1) logical = true
    opts.Ephemeris = []
end

c = 299792458; AU = 149597870700; Tsun = 1.32712440018e20 / c^3;
use = @(name) any(opts.Terms == name);
a = deg2rad(psr.raDeg); dl = deg2rad(psr.decDeg);
n = [cos(dl) * cos(a); cos(dl) * sin(a); sin(dl)];
N = numel(t);

if isempty(opts.Ephemeris)
    tdbTT = tdbMinusTT(t);                                % TDB - TT at the geocentre
    [rE, rS] = earthAndSun(t, opts.Moon, opts.Planets);   % geocentre and Sun wrt the SSB
    h = 60;                                               % velocity by central difference
    vE = (earthAndSun(t + h, opts.Moon, opts.Planets) - earthAndSun(t - h, opts.Moon, opts.Planets)) / (2*h);
else
    eph = opts.Ephemeris;
    if ~isstruct(eph) || ~isfield(eph, 'earth') || ~isfield(eph, 'sun') || ~isfield(eph, 'tdbMinusTT')
        error('delayModel:ephemeris', '''Ephemeris'' must be a loadEphemeris struct, or [] for the analytic orbit.');
    end
    tdbTT = eph.tdbMinusTT(t);                            % TDB - TT from the ephemeris (D5c)
    tdb = t + tdbTT;                                      % the ephemeris' time scale (1.7 ms: 50 m of the Earth's motion)
    [rE, vE] = eph.earth(tdb); rS = eph.sun(tdb);
end

switch upper(traj.frame)
    case 'ECEF'
        rT = traj.pos(t); vT = traj.vel(t);
        th = eraAngle(t);                             % Earth rotation angle [rad]
        cs = cos(th); sn = sin(th);
        rSite = [cs .* rT(1, :) - sn .* rT(2, :); sn .* rT(1, :) + cs .* rT(2, :); rT(3, :)];
        w = 2*pi * 1.00273781191135448 / 86400;       % d(angle)/dt [rad/s]
        vRot = [cs .* vT(1, :) - sn .* vT(2, :); sn .* vT(1, :) + cs .* vT(2, :); vT(3, :)];
        vSite = vRot + w * [-rSite(2, :); rSite(1, :); zeros(1, N)];
        rObs = rE + rSite; vObs = vE + vSite;
    case 'SSB'
        rObs = traj.pos(t); vObs = traj.vel(t);
        rSite = zeros(3, N);
    otherwise
        error('delayModel:frame', 'Trajectory frame must be ''ECEF'' or ''SSB''.');
end

rn = n.' * rObs;                                       % r_obs . n [m]
z = zeros(1, N);
d = struct('total', z, 'roemer', z, 'curvature', z, 'einstein', z, 'shapiro', z);
if use("roemer"), d.roemer = rn / c; end
if use("curvature") && isfield(psr, 'distance') && ~isempty(psr.distance)
    d.curvature = -(sum(rObs.^2, 1) - rn.^2) / (2 * c * psr.distance);
end
if use("einstein") && strcmpi(traj.frame, 'ECEF')
    d.einstein = tdbTT + sum(vE .* rSite, 1) / c^2;
end
if use("shapiro")
    rs = rObs - rS;                                    % observer relative to the Sun
    d.shapiro = 2 * Tsun * log((sqrt(sum(rs.^2, 1)) + n.' * rs) / AU);
end
d.total  = d.roemer + d.curvature + d.einstein + d.shapiro;
d.rObs   = rObs; d.vObs = vObs; d.rEarth = rE; d.vEarth = vE; d.rSun = rS;
d.nHat   = n; d.frame = upper(traj.frame);
if isempty(opts.Ephemeris)
    d.model = sprintf(['analytic (Keplerian EMB, Moon %d, planets %d; Earth rotation angle, no ' ...
                       'precession/nutation)'], opts.Moon, opts.Planets);
else
    d.model = sprintf('JPL ephemeris %s; Earth rotation angle, no precession/nutation', opts.Ephemeris.file);
end
end


% =====================================================================================
function dt = tdbMinusTT(t)
%TDBMINUSTT  TDB - TT at the geocentre [s]: its two largest terms (the clock's varying rate along the
% Earth's eccentric orbit; g the Earth's mean anomaly).
g = deg2rad(357.53 + 0.98560028 * (t / 86400));
dt = 1.657e-3 * sin(g) + 1.4e-5 * sin(2*g);
end


% =====================================================================================
function [rE, rS] = earthAndSun(t, withMoon, withPlanets)
%EARTHANDSUN  Geocentre and Sun relative to the SSB [3 x N] m, equatorial J2000.
AU = 149597870700; T = t / (36525 * 86400);           % Julian centuries since J2000
% Standish's elements (a [AU], e, I, L, long. perihelion, long. node [deg]) and rates per century
emb = [1.00000261 0.01671123 -0.00001531 100.46457166 102.93768193 0;
       0.00000562 -0.00004392 -0.01294668 35999.37244981 0.32327364 0];
jup = [5.20288700 0.04838624 1.30439695 34.39644051 14.72847983 100.47390909;
       -0.00011607 -0.00013253 -0.00183714 3034.74612775 0.21252668 0.20469106];
sat = [9.53667594 0.05386179 2.48599187 49.95424423 92.59887831 113.66242448;
       -0.00125060 -0.00050991 0.00193609 1222.49362201 -0.41897216 -0.28867794];
rEMB = kepler(emb, T) * AU;                            % heliocentric, ecliptic
rSun = zeros(size(rEMB));
if withPlanets                                         % the Sun around the SSB (masses / M_sun)
    mE = 1/328900.56; mJ = 1/1047.348644; mS = 1/3497.9018;
    rSun = -(mE * rEMB + mJ * kepler(jup, T) * AU + mS * kepler(sat, T) * AU) / (1 + mE + mJ + mS);
end
rEc = rEMB + rSun;                                     % EMB wrt the SSB
if withMoon                                            % the Earth around the EMB
    rEc = rEc - 0.0121505856 * moonGeocentric(t);
end
rS = eclToEq(rSun);
rE = eclToEq(rEc);
end

function r = kepler(el, T)
%KEPLER  Heliocentric ecliptic-J2000 position [AU] from elements el (row 1) and rates per century
% (row 2) at T centuries: mean anomaly, Kepler's equation (Newton), orbit plane -> ecliptic.
e  = el(1, 2) + el(2, 2) * T;
a  = el(1, 1) + el(2, 1) * T;
I  = deg2rad(el(1, 3) + el(2, 3) * T);
L  = el(1, 4) + el(2, 4) * T;
wb = el(1, 5) + el(2, 5) * T;
Om = deg2rad(el(1, 6) + el(2, 6) * T);
M  = deg2rad(mod(L - wb + 180, 360) - 180);
w  = deg2rad(wb) - Om;                                 % argument of perihelion
E  = M + e .* sin(M);
for it = 1:8
    E = E - (E - e .* sin(E) - M) ./ (1 - e .* cos(E));
end
xp = a .* (cos(E) - e); yp = a .* sqrt(1 - e.^2) .* sin(E);
cw = cos(w); sw = sin(w); cO = cos(Om); sO = sin(Om); cI = cos(I); sI = sin(I);
r = [(cw.*cO - sw.*sO.*cI) .* xp + (-sw.*cO - cw.*sO.*cI) .* yp;
     (cw.*sO + sw.*cO.*cI) .* xp + (-sw.*sO + cw.*cO.*cI) .* yp;
     (sw.*sI) .* xp + (cw.*sI) .* yp];
end

function r = moonGeocentric(t)
%MOONGEOCENTRIC  Low-precision geocentric Moon [3 x N] m, ecliptic (Astronomical Almanac series:
% ~0.3 deg, ~1000 km).
dd = t / 86400;
Lp = 218.316 + 13.176396 * dd; Mp = 134.963 + 13.064993 * dd; F = 93.272 + 13.229350 * dd;
lam = deg2rad(Lp + 6.289 * sind(Mp)); bet = deg2rad(5.128 * sind(F));
D = (385001 - 20905 * cosd(Mp)) * 1e3;
r = [D .* cos(bet) .* cos(lam); D .* cos(bet) .* sin(lam); D .* sin(bet)];
end

function r = eclToEq(r)
%ECLTOEQ  Ecliptic J2000 -> equatorial J2000 (obliquity 23.439281 deg).
ep = deg2rad(23.439281);
r = [r(1, :); cos(ep) * r(2, :) - sin(ep) * r(3, :); sin(ep) * r(2, :) + cos(ep) * r(3, :)];
end

function th = eraAngle(t)
%ERAANGLE  Earth rotation angle [rad] (IERS 2003) with UT1 ~ TT - 69.184 s.
Du = (t - 69.184) / 86400;
th = 2*pi * mod(0.7790572732640 + 1.00273781191135448 * Du, 1);
end
