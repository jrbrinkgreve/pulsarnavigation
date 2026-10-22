import numpy as np
import scipy.io as sio
from pathlib import Path
from plotsvg import Ax
# navFit converging (navigation.html #navfit): the east component of a site 100 km east of the assumed one, fitted to the
# TOAs of the first N days, N = 1 ... 100 (runScenario motion run, B0329+54, nominal flux, perfect clock, 22 Oct 2026; local only:
# data/scenario/scenario_B0329_nominal_clean_100d_motion_perfect.mat, needs scipy). The fit is navFit's "self" model - dX, dY and
# one constant - by weighted least squares (navFit without its outlier model; the rows at each TOA's time: for these TOAs the
# averaged rows change the result by under 1 km).
root = Path(__file__).resolve().parents[3]
R = sio.loadmat(root / 'data/scenario/scenario_B0329_nominal_clean_100d_motion_perfect.mat', squeeze_me=True, struct_as_record=False)['R']
G = [g for g in np.atleast_1d(R.groups) if g.good]
t = np.array([g.tAbs for g in G]); r = np.array([g.resid for g in G]); s = np.array([g.sigma for g in G])
c = 299792458.0; lon = np.radians(5)
ra, dec = np.radians(R.motionTruth.psr.raDeg), np.radians(R.motionTruth.psr.decDeg)
n = np.array([np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra), np.sin(dec)])
th = 2 * np.pi * np.mod(0.7790572732640 + 1.00273781191135448 * (t - 69.184) / 86400, 1)
A = -np.column_stack([np.cos(th) * n[0] + np.sin(th) * n[1], -np.sin(th) * n[0] + np.cos(th) * n[1]]) / c
A = np.column_stack([A, np.ones_like(t)])
e = np.array([-np.sin(lon), np.cos(lon), 0.0])          # east in terms of dX, dY (and the constant)
day = (t - t[0]) / 86400
days = np.arange(1, 101); est = []; err = []
for N in days:
    m = day < N
    W = 1 / s[m]**2
    C = np.linalg.inv(A[m].T @ (A[m] * W[:, None]))
    p = C @ (A[m].T @ (W * r[m]))
    est.append(e @ p / 1e3); err.append(np.sqrt(e @ C @ e) / 1e3)
est, err = np.array(est), np.array(err)
print(f'{len(t)} TOAs; day 30: {est[29]:.1f} +- {err[29]:.1f} km; day 100: {est[99]:.1f} +- {err[99]:.1f} km ({est[99] / err[99]:.1f} sigma)')
o = []
P = Ax(90, 40, 640, 240, (0, 100), (-150, 250))
for v in (0, 100, 200):
    o.append(f'<path class="grid" d="M{P.x0} {float(P.Y(v)):.1f} H{P.x0 + P.w}"/>')
o.append(P.frame(xticks=[0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100], yticks=[-100, 0, 100, 200], xlab='days of TOAs',
                 ylab='east component of the offset [km]'))
band = np.concatenate([est + err, (est - err)[::-1]]); dd = np.concatenate([days, days[::-1]])
o.append(f'<path class="s1band" d="M' + ' L'.join(f'{float(P.X(a)):.1f} {float(P.Y(b)):.1f}' for a, b in zip(dd, band)) + ' Z"/>')
o.append(f'<path class="s1" d="{P.path(days, est)}"><title>the fit to the first N days</title></path>')
o.append(f'<path class="s2" d="M{P.x0} {float(P.Y(100)):.1f} H{P.x0 + P.w}" style="stroke-dasharray:6 4"/>')
o.append(f'<text class="sub" x="{P.x0 + P.w + 6}" y="{float(P.Y(100)) + 4:.1f}">truth</text>')
for N in (30, 100):
    x, y = float(P.X(N)), float(P.Y(est[N - 1]))
    o.append(f'<circle class="s1f" cx="{x:.1f}" cy="{y:.1f}" r="5"><title>{N} days: {est[N - 1]:.1f} ± {err[N - 1]:.1f} km</title></circle>')
    o.append(f'<text class="sub" x="{x - 8:.1f}" y="{float(P.Y(est[N - 1] - err[N - 1])) + 16:.1f}" text-anchor="{"middle" if N < 100 else "end"}">{N} days: {est[N - 1]:.0f} ± {err[N - 1]:.0f} km</text>')
lx, ly = 90, 340
o.append(f'<path class="s1" d="M{lx} {ly} h20"/><text class="sub" x="{lx + 26}" y="{ly + 4}">the fit to the first N days, ±1σ</text>')
o.append(f'<path class="s2" d="M{lx + 260} {ly} h20" style="stroke-dasharray:6 4"/><text class="sub" x="{lx + 286}" y="{ly + 4}">the true offset, 100 km east</text>')
svg = '''<figure class="scroll">
<svg class="wide" viewBox="0 0 860 352" role="img" aria-label="The east component of the site offset fitted to the first N days of simulated TOAs, converging on the true 100 km">
  ''' + '\n  '.join(o) + f'''
</svg>
<figcaption>How the fit converges: the east component of the site offset, fitted to the first <var>N</var> days of B0329+54's simulated
TOAs (nominal flux, perfect clock; the truth 100 km east). The error shrinks as 1/√<var>N</var>: ±{err[29]:.0f} km after 30 days and
±{err[99]:.0f} km after 100, where the offset stands out at {est[99] / err[99]:.1f}σ (<code>navFit</code> itself, with its outlier model: 80 ± 16 km, 4.9σ). After 30 days this noise draw happens to sit 2σ low; the
later days pull it back. Data: <code>runScenario</code>'s 100-day motion run; the fit is <code>navFit</code>'s model by weighted least
squares.</figcaption>
</figure>'''
open('navfit.svg', 'w').write(svg)
