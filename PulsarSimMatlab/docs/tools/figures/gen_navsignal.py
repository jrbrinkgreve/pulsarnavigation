import numpy as np
import scipy.io as sio
from pathlib import Path
from plotsvg import Ax
# The navigation signal (navigation.html #signal): B0329+54's simulated TOAs from 52 N 5 E, the true site 100 km east of the
# assumed one, perfect clock, nominal flux (runScenario motion run of 21 Oct 2026; local only:
# data/scenario/pre-22oct-tag/scenario_B0329_nominal_clean_30d_motion_perfect.mat, needs scipy).
# The curve: navFit's row r = -R_z(theta)' n . dx / c for dx = 100 km east, theta the Earth rotation angle.
root = Path(__file__).resolve().parents[3]
R = sio.loadmat(root / 'data/scenario/pre-22oct-tag/scenario_B0329_nominal_clean_30d_motion_perfect.mat',
                squeeze_me=True, struct_as_record=False)['R']
G = [g for g in np.atleast_1d(R.groups) if g.good]
t = np.array([g.tAbs for g in G]); r = np.array([g.resid for g in G]) * 1e6; s = np.array([g.sigma for g in G]) * 1e6
c = 299792458.0; lat, lon = np.radians(52), np.radians(5)
ra, dec = np.radians(R.motionTruth.psr.raDeg), np.radians(R.motionTruth.psr.decDeg)
n = np.array([np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra), np.sin(dec)])
dx = 100e3 * np.array([-np.sin(lon), np.cos(lon), 0])
era = lambda tt: 2 * np.pi * np.mod(0.7790572732640 + 1.00273781191135448 * (tt - 69.184) / 86400, 1)
def model(tt):
    th = era(tt)
    return -((np.cos(th) * n[0] + np.sin(th) * n[1]) * dx[0] + (-np.sin(th) * n[0] + np.cos(th) * n[1]) * dx[1]) / c * 1e6
w = 1 / s**2
c0 = np.sum(w * (r - model(t))) / np.sum(w)              # the per-pulsar constant (profile reference), as the fit's c_k
r = r - c0
truth = np.array([g.trueResid for g in G]) * 1e6
print(f'{len(t)} good TOAs; sigma median {np.median(s):.0f} us; constant {c0:.0f} us; model vs trueResid max {np.max(np.abs(model(t) - truth)):.2f} us')
o = []
# left: three days of TOAs
d0 = np.floor(t[0] / 86400) * 86400                     # the first TOA's day
td = (t - t[0]) / 3600; sel = td <= 72
L = Ax(80, 50, 400, 230, (0, 72), (-2600, 2600))
H = lambda tt: np.mod(era(tt) + lon - ra + np.pi, 2 * np.pi) - np.pi
tc = np.linspace(t[0], t[0] + 72 * 3600, 2000)
el = np.degrees(np.arcsin(np.sin(lat) * np.sin(dec) + np.cos(lat) * np.cos(dec) * np.cos(H(tc))))
below = el < 20; edges = np.flatnonzero(np.diff(below.astype(int)))
segs = []; start = 0 if below[0] else None
for e in edges:
    if below[e + 1]: start = e + 1
    elif start is not None: segs.append((start, e)); start = None
if start is not None: segs.append((start, len(tc) - 1))
for a, b in segs:
    xa, xb = float(L.X((tc[a] - t[0]) / 3600)), float(L.X((tc[b] - t[0]) / 3600))
    o.append(f'<rect class="fillfg" x="{xa:.1f}" y="{L.y0}" width="{xb - xa:.1f}" height="{L.h}"><title>below 20°: no data</title></rect>')
o.append(f'<path class="grid" d="M{L.x0} {float(L.Y(0)):.1f} H{L.x0 + L.w}"/>')
o.append(L.frame(xticks=[0, 12, 24, 36, 48, 60, 72], yticks=[-2000, -1000, 0, 1000, 2000], xlab='hours since the first TOA', ylab='residual [µs]'))
for x, y, e in zip(td[sel], r[sel], s[sel]):
    X, Y1, Y2 = float(L.X(x)), float(L.Y(y - e)), float(L.Y(y + e))
    o.append(f'<path class="err" d="M{X:.1f} {Y1:.1f} V{Y2:.1f}"/>')
for x, y, e in zip(td[sel], r[sel], s[sel]):
    o.append(f'<circle class="s1f" cx="{float(L.X(x)):.1f}" cy="{float(L.Y(y)):.1f}" r="4"><title>a TOA: {y:+.0f} ± {e:.0f} µs</title></circle>')
o.append(f'<path class="s2" d="{L.path((tc - t[0]) / 3600, model(tc))}"/>')
o.append('<text class="ttl" x="80" y="24">Three days of TOAs</text>')
# right: 30 days folded on the sidereal day, 2 h bins of hour angle
Hh = np.degrees(H(t)) / 15
Rp = Ax(560, 50, 270, 230, (-12, 12), (-500, 500))
o.append(f'<path class="grid" d="M{Rp.x0} {float(Rp.Y(0)):.1f} H{Rp.x0 + Rp.w}"/>')
o.append(Rp.frame(xticks=[-12, -6, 0, 6, 12], yticks=[-400, -200, 0, 200, 400], xlab='hour angle [h]', ylab='mean residual [µs]'))
Hc = np.linspace(-12, 12, 400)
# the curve against hour angle: the same model at times with that hour angle
tH = t[0] + (Hc - np.degrees(H(t[0])) / 15) * 3600 / 1.00273781191135448
o.append(f'<path class="s2" d="{Rp.path(Hc, model(tH))}"/>')
for b0 in range(-12, 12, 2):
    m = (Hh >= b0) & (Hh < b0 + 2)
    if m.sum() < 3: continue
    mu = np.sum(w[m] * r[m]) / np.sum(w[m]); err = 1 / np.sqrt(np.sum(w[m]))
    X = float(Rp.X(b0 + 1))
    o.append(f'<path class="err" d="M{X:.1f} {float(Rp.Y(mu - err)):.1f} V{float(Rp.Y(mu + err)):.1f}"/>')
    o.append(f'<circle class="s1f" cx="{X:.1f}" cy="{float(Rp.Y(mu)):.1f}" r="5"><title>{m.sum()} TOAs, hour angle {b0} to {b0 + 2} h: {mu:+.0f} ± {err:.0f} µs</title></circle>')
o.append(f'<text class="ttl" x="560" y="24">All {len(t)} TOAs of 30 days, by hour angle</text>')
lx, ly = 80, 326
o.append(f'<circle class="s1f" cx="{lx}" cy="{ly}" r="4"/><text class="sub" x="{lx + 10}" y="{ly + 4}">measured residuals, ±1σ</text>')
o.append(f'<path class="s2" d="M{lx + 210} {ly} h20"/><text class="sub" x="{lx + 236}" y="{ly + 4}">the signal of a site 100 km east: a daily sine of 193 µs</text>')
svg = '''<figure class="scroll">
<svg class="wide" viewBox="0 0 860 340" role="img" aria-label="Simulated TOA residuals of B0329+54 for a site 100 km east: three days of single TOAs, and all 30 days averaged by hour angle, against the daily sine of 193 microseconds">
  ''' + '\n  '.join(o) + f'''
</svg>
<figcaption><b>Figure 10.1.</b> The navigation signal in the simulated data. Left: three days of B0329+54's TOAs for a site 100 km east of the assumed one
(grey: the pulsar below 20°). Each TOA's error bar, about {np.median(s):.0f} µs, is four times the signal, a daily sine of 193 µs. Right: all {len(t)} TOAs of
30 days averaged in 2-hour bins of hour angle, where the sine stands out. Data: the 30-day run of <code>runScenario</code> with
<code>motion.on</code>; the constant offset of the profile reference is removed.</figcaption>
</figure>'''
open('navsignal.svg', 'w').write(svg)
