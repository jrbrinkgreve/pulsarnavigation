import numpy as np
import scipy.io as sio
from pathlib import Path
from plotsvg import Ax
# The ship (navigation.html #ship): runShip's results of 22 Oct 2026 (local only: data/scenario/ship_14d_nominal[_cal][_seed<n>].mat,
# needs scipy). (A) one ship: the drift off the dead-reckoning track, east and north, the truth against the smoothed filter with
# its 1-sigma band. (B) all seeds: the rms track error of dead reckoning, the filter with a self-calibrated start and with a
# calibrated one.
root = Path(__file__).resolve().parents[3] / 'data/scenario'
SEED = 7                                               # the ship shown in (A)
def load(seed, cal):
    f = root / f"ship_14d_nominal{'_cal' if cal else ''}{'' if seed == 1 else f'_seed{seed}'}.mat"
    return sio.loadmat(f, squeeze_me=True, struct_as_record=False)['R'] if f.exists() else None
o = []
# (A)
R = load(SEED, False); F = R.final
d = (F.t - F.t[0]) / 86400
for j, (comp, x0) in enumerate((('east', 70), ('north', 500))):
    tr, est, err = F.truth[j] / 1e3, F.enu[j] / 1e3, F.enuErr[j] / 1e3
    lo = min(tr.min(), (est - err).min()); hi = max(tr.max(), (est + err).max())
    pad = 0.08 * (hi - lo); P = Ax(x0, 50, 330, 210, (0, 14), (lo - pad, hi + pad))
    step = 100 if hi - lo > 300 else 50
    yt = list(np.arange(np.ceil((lo - pad) / step) * step, hi + pad, step))
    o.append(P.frame(xticks=[0, 2, 4, 6, 8, 10, 12, 14], yticks=yt, xlab='days', ylab=f'{comp}: off the dead-reckoning track [km]', yfmt='{:.0f}'))
    if lo < 0 < hi:
        o.append(f'<path class="grid" d="M{P.x0} {float(P.Y(0)):.1f} H{P.x0 + P.w}"/>')
    band = np.concatenate([est + err, (est - err)[::-1]]); dd = np.concatenate([d, d[::-1]])
    pts = ' L'.join(f'{float(P.X(a)):.1f} {float(P.Y(b)):.1f}' for a, b in zip(dd, band))
    o.append(f'<path class="s1band" d="M{pts} Z"/>')
    o.append(f'<path class="s1" d="{P.path(d, est)}"><title>the filter (smoothed): {comp}</title></path>')
    o.append(f'<path class="s2" d="{P.path(d, tr)}"><title>the truth: {comp}</title></path>')
    yT, yE = float(P.Y(tr[-1])), float(P.Y(est[-1]))     # end labels inside the panel: the upper one above, the lower below
    up, dn = (yT, yE) if yT < yE else (yE, yT)
    lab = {yT: f'truth {tr[-1]:.0f} km', yE: f'filter {est[-1]:.0f} ± {err[-1]:.0f} km'}
    o.append(f'<text class="sub" x="{float(P.X(13.8)):.1f}" y="{up - 10:.1f}" text-anchor="end">{lab[up]}</text>')
    o.append(f'<text class="sub" x="{float(P.X(13.8)):.1f}" y="{dn + 18:.1f}" text-anchor="end">{lab[dn]}</text>')
o.append(f'<text class="ttl" x="70" y="24">One ship (seed {SEED}): where the current pushed it, and what the pulsars say</text>')
lx, ly = 70, 316
o.append(f'<path class="s2" d="M{lx} {ly} h20"/><text class="sub" x="{lx + 26}" y="{ly + 4}">truth</text>')
o.append(f'<path class="s1" d="M{lx + 90} {ly} h20"/><text class="sub" x="{lx + 116}" y="{ly + 4}">the filter, smoothed over all TOAs, ±1σ</text>')
svgA = '''<figure class="scroll">
<svg class="wide" viewBox="0 0 860 330" role="img" aria-label="One ship over 14 days: the true drift off its dead-reckoning track, east and north, against the navigation filter's estimate">
  ''' + '\n  '.join(o) + f'''
</svg>
<figcaption>One simulated ship (<code>runShip</code>, seed {SEED}, self-calibrated start): over 14 days an unknown current pushes it
{np.hypot(*F.truth[:, -1]) / 1e3:.0f} km off its dead-reckoning track. B0329+54 and B0950+08 alone, with a drifting SDR clock, keep the smoothed
track within its error band. The band widens towards the end, where no later TOAs help.</figcaption>
</figure>'''
open('navship.svg', 'w').write(svgA)
# (B) all seeds
rows = []
for seed in range(1, 9):
    a, b = load(seed, False), load(seed, True)
    if a is None: continue
    dr = np.sqrt(np.mean(np.sum(a.final.truth[:2] ** 2, 0))) / 1e3
    e = lambda Rx: np.sqrt(np.mean(np.sum((Rx.final.enu[:2] - Rx.final.truth[:2]) ** 2, 0))) / 1e3
    rows.append((seed, dr, e(a), e(b) if b is not None else np.nan))
o = []
P = Ax(120, 40, 600, 30 * len(rows), (5, 2000), (len(rows) + 0.5, 0.5), logy=False)
lg = lambda v: np.log10(v)
X = lambda v: P.x0 + (lg(v) - lg(5)) / (lg(2000) - lg(5)) * P.w
for v in (10, 30, 100, 300, 1000):
    o.append(f'<path class="grid" d="M{X(v):.1f} {P.y0} V{P.y0 + P.h}"/><text class="sub" x="{X(v):.1f}" y="{P.y0 + P.h + 18}" text-anchor="middle">{v}</text>')
o.append(f'<text class="sub" x="{P.x0 + P.w / 2}" y="{P.y0 + P.h + 36}" text-anchor="middle">rms error of the track over 14 days [km] (log scale)</text>')
for i, (seed, dr, es, ec) in enumerate(rows):
    y = P.y0 + 15 + 30 * i
    o.append(f'<text class="sub" x="{P.x0 - 10}" y="{y + 4}" text-anchor="end">seed {seed}</text>')
    lo, hi = np.nanmin([es, ec, dr]), np.nanmax([es, ec, dr])
    o.append(f'<path class="thin" d="M{X(lo):.1f} {y} H{X(hi):.1f}"/>')
    o.append(f'<circle class="fillfg" cx="{X(dr):.1f}" cy="{y}" r="6" style="opacity:0.7"><title>seed {seed}, dead reckoning: {dr:.0f} km</title></circle>')
    o.append(f'<circle class="s1f" cx="{X(es):.1f}" cy="{y}" r="6"><title>seed {seed}, filter, self-calibrated start: {es:.0f} km</title></circle>')
    if np.isfinite(ec):
        o.append(f'<circle class="s2f" cx="{X(ec):.1f}" cy="{y}" r="6"><title>seed {seed}, filter, calibrated start: {ec:.0f} km</title></circle>')
med = np.nanmedian(np.array(rows)[:, 1:], 0)
ly = P.y0 + P.h + 58
o.append(f'<circle class="fillfg" cx="{P.x0}" cy="{ly}" r="6" style="opacity:0.7"/><text class="sub" x="{P.x0 + 12}" y="{ly + 4}">dead reckoning (median {med[0]:.0f} km)</text>')
o.append(f'<circle class="s1f" cx="{P.x0 + 230}" cy="{ly}" r="6"/><text class="sub" x="{P.x0 + 242}" y="{ly + 4}">filter, self-calibrated ({med[1]:.0f} km)</text>')
o.append(f'<circle class="s2f" cx="{P.x0 + 440}" cy="{ly}" r="6"/><text class="sub" x="{P.x0 + 452}" y="{ly + 4}">filter, calibrated start ({med[2]:.0f} km)</text>')
H = int(ly + 20)
svgB = f'''<figure class="scroll">
<svg class="wide" viewBox="0 0 860 {H}" role="img" aria-label="Track errors of eight simulated ships: dead reckoning against the navigation filter with a self-calibrated and a calibrated start">
  ''' + '\n  '.join(o) + f'''
</svg>
<figcaption>Eight ships with different currents, clocks and noise (<code>runShip</code>, seeds 1–8, 14 days): the rms distance of the
track from the truth for dead reckoning alone, and for the smoothed filter with a self-calibrated and with a calibrated start. Medians
{med[0]:.0f}, {med[1]:.0f} and {med[2]:.0f} km.</figcaption>
</figure>'''
open('navships.svg', 'w').write(svgB)
print('rows', [(r[0], round(r[1]), round(r[2]), None if np.isnan(r[3]) else round(r[3])) for r in rows])
