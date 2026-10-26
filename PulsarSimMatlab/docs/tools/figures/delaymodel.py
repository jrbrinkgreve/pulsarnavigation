import numpy as np
# A line-by-line port of functions/delayModel.m's analytic branch (no 'Ephemeris') and of makeTrajectory's
# 'fixed' and linear 'track' types, for the figures of localisation.html and fold.html#prediction. Same constants,
# same formulas; times in s since J2000.0 (TT), positions in m, equatorial J2000.
C = 299792458.0; AU = 149597870700.0; TSUN = 1.32712440018e20 / C**3
B0329 = dict(raDeg=53.2474, decDeg=54.5787)
EPOCH = 9555.5 * 86400                                   # 2026-03-01 00:00 TT, as runScenario and the tests

def geodetic_to_ecef(lat, lon, h=0.0):
    a = 6378137.0; f = 1 / 298.257223563; e2 = f * (2 - f)
    p, l = np.radians(lat), np.radians(lon)
    N = a / np.sqrt(1 - e2 * np.sin(p)**2)
    return np.array([(N + h) * np.cos(p) * np.cos(l), (N + h) * np.cos(p) * np.sin(l), (N * (1 - e2) + h) * np.sin(p)])

def fixed(lat, lon, h=0.0):
    r = geodetic_to_ecef(lat, lon, h)
    return dict(frame='ECEF', pos=lambda t: np.repeat(r[:, None], np.size(t), 1), vel=lambda t: np.zeros((3, np.size(t))))

def track(tk, lat, lon):
    """makeTrajectory('track', ...) with the default 'linear' interpolation (in ECEF)."""
    tk = np.asarray(tk, float); P = np.stack([geodetic_to_ecef(a, b) for a, b in zip(lat, lon)], 1)
    def pos(t):
        t = np.atleast_1d(t); return np.stack([np.interp(t, tk, P[i]) for i in range(3)])
    def vel(t):
        t = np.atleast_1d(t); k = np.clip(np.searchsorted(tk, t, side='right') - 1, 0, len(tk) - 2)
        return (P[:, k + 1] - P[:, k]) / (tk[k + 1] - tk[k])
    return dict(frame='ECEF', pos=pos, vel=vel)

def tdb_minus_tt(t):
    g = np.radians(357.53 + 0.98560028 * (t / 86400))
    return 1.657e-3 * np.sin(g) + 1.4e-5 * np.sin(2 * g)

def kepler(el, T):
    e = el[0][1] + el[1][1] * T; a = el[0][0] + el[1][0] * T
    I = np.radians(el[0][2] + el[1][2] * T); L = el[0][3] + el[1][3] * T
    wb = el[0][4] + el[1][4] * T; Om = np.radians(el[0][5] + el[1][5] * T)
    M = np.radians(np.mod(L - wb + 180, 360) - 180); w = np.radians(wb) - Om
    E = M + e * np.sin(M)
    for _ in range(8):
        E = E - (E - e * np.sin(E) - M) / (1 - e * np.cos(E))
    xp = a * (np.cos(E) - e); yp = a * np.sqrt(1 - e**2) * np.sin(E)
    cw, sw, cO, sO, cI, sI = np.cos(w), np.sin(w), np.cos(Om), np.sin(Om), np.cos(I), np.sin(I)
    return np.stack([(cw*cO - sw*sO*cI) * xp + (-sw*cO - cw*sO*cI) * yp,
                     (cw*sO + sw*cO*cI) * xp + (-sw*sO + cw*cO*cI) * yp,
                     (sw*sI) * xp + (cw*sI) * yp])

def moon_geocentric(t):
    dd = t / 86400
    Lp = 218.316 + 13.176396 * dd; Mp = 134.963 + 13.064993 * dd; F = 93.272 + 13.229350 * dd
    lam = np.radians(Lp + 6.289 * np.sin(np.radians(Mp))); bet = np.radians(5.128 * np.sin(np.radians(F)))
    D = (385001 - 20905 * np.cos(np.radians(Mp))) * 1e3
    return np.stack([D * np.cos(bet) * np.cos(lam), D * np.cos(bet) * np.sin(lam), D * np.sin(bet)])

def ecl_to_eq(r):
    ep = np.radians(23.439281)
    return np.stack([r[0], np.cos(ep) * r[1] - np.sin(ep) * r[2], np.sin(ep) * r[1] + np.cos(ep) * r[2]])

EMB = [[1.00000261, 0.01671123, -0.00001531, 100.46457166, 102.93768193, 0],
       [0.00000562, -0.00004392, -0.01294668, 35999.37244981, 0.32327364, 0]]
JUP = [[5.20288700, 0.04838624, 1.30439695, 34.39644051, 14.72847983, 100.47390909],
       [-0.00011607, -0.00013253, -0.00183714, 3034.74612775, 0.21252668, 0.20469106]]
SAT = [[9.53667594, 0.05386179, 2.48599187, 49.95424423, 92.59887831, 113.66242448],
       [-0.00125060, -0.00050991, 0.00193609, 1222.49362201, -0.41897216, -0.28867794]]

def earth_and_sun(t):
    T = t / (36525 * 86400)
    rEMB = kepler(EMB, T) * AU
    mE, mJ, mS = 1 / 328900.56, 1 / 1047.348644, 1 / 3497.9018
    rSun = -(mE * rEMB + mJ * kepler(JUP, T) * AU + mS * kepler(SAT, T) * AU) / (1 + mE + mJ + mS)
    rEc = rEMB + rSun - 0.0121505856 * moon_geocentric(t)
    return ecl_to_eq(rEc), ecl_to_eq(rSun)

def era_angle(t, dT=69.184):
    return 2 * np.pi * np.mod(0.7790572732640 + 1.00273781191135448 * (t - dT) / 86400, 1)

def n_hat(psr):
    a, d = np.radians(psr['raDeg']), np.radians(psr['decDeg'])
    return np.array([np.cos(d) * np.cos(a), np.cos(d) * np.sin(a), np.sin(d)])

def delay(t, traj, psr, distance=None):
    """delayModel(t, traj, psr): the terms [s] as a dict, plus rObs, rSite, vE."""
    t = np.atleast_1d(np.asarray(t, float)); n = n_hat(psr)
    rE, rS = earth_and_sun(t)
    h = 60.0; vE = (earth_and_sun(t + h)[0] - earth_and_sun(t - h)[0]) / (2 * h)
    rT = traj['pos'](t); th = era_angle(t); cs, sn = np.cos(th), np.sin(th)
    rSite = np.stack([cs * rT[0] - sn * rT[1], sn * rT[0] + cs * rT[1], rT[2]])
    rObs = rE + rSite; rn = n @ rObs
    d = dict(roemer=rn / C, curvature=np.zeros_like(t), einstein=tdb_minus_tt(t) + np.sum(vE * rSite, 0) / C**2)
    if distance:
        d['curvature'] = -(np.sum(rObs**2, 0) - rn**2) / (2 * C * distance)
    rs = rObs - rS
    d['shapiro'] = 2 * TSUN * np.log((np.sqrt(np.sum(rs**2, 0)) + n @ rs) / AU)
    d['total'] = d['roemer'] + d['curvature'] + d['einstein'] + d['shapiro']
    d.update(rObs=rObs, rSite=rSite, rEarth=rE, rSun=rS, vE=vE, n=n)
    return d
