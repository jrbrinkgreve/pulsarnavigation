import numpy as np
from plotsvg import Ax
# A TOA is a time-average (navigation.html #averaged): B0329+54 at 52 N, the site 100 km east of the assumed one.
# navFit's row a = -R_z(theta)' n / c gives the residual r = A sin(H), H the hour angle, A = 100 km cos(dec) / c.
# Left: a TOA from 10 h of one pass. Right: a TOA from two 5 h stretches a sidereal day apart at the same hour angle.
c = 299792458.0; dec = np.radians(54.5787)
A = 100e3 * np.cos(dec) / c * 1e6                      # [us]
r = lambda H: A * np.sin(2 * np.pi * H / 24)           # H: hour angle [sidereal hours]
gapHalf = (24 - 19.6) / 2                              # B0329+54 below 20 deg: 4.4 h around anti-transit
o = []

def panel(x0, xlim, title):
    P = Ax(x0, 58, 330, 220, xlim, (-230, 230))
    o.append(f'<text class="ttl" x="{x0}" y="22">{title}</text>')
    for gc in np.arange(12, xlim[1] + 1, 24):          # the daily gap (pulsar below 20 deg)
        a, b = max(gc - gapHalf, xlim[0]), min(gc + gapHalf, xlim[1])
        if b > a:
            o.append(f'<rect class="fillfg" x="{float(P.X(a)):.1f}" y="{P.y0}" width="{float(P.X(b) - P.X(a)):.1f}" height="{P.h}"><title>below 20°: no data</title></rect>')
    o.append(f'<path class="grid" d="M{P.x0} {float(P.Y(0)):.1f} H{P.x0 + P.w}"/>')
    H = np.linspace(*xlim, 600)
    o.append(f'<path class="thin" d="{P.path(H, r(H))}"/>')
    o.append(P.frame(xticks=list(range(int(xlim[0]), int(xlim[1]) + 1, 6)), yticks=[-200, -100, 0, 100, 200],
                     xlab='hour angle [h] (0 = transit)', ylab='residual [µs]'))
    return P

def dot(P, x, y, cls, tip):
    X, Y = float(P.X(x)), float(P.Y(y))
    o.append(f'<circle class="{cls}" cx="{X:.1f}" cy="{Y:.1f}" r="5"><title>{tip}</title></circle>')
    return X, Y

# left: one pass, 10 h of data
L = panel(70, (-12, 24), 'One pass: 10 hours of data')
H1 = np.linspace(-1, 9, 1201); avg1 = r(H1).mean()
o.append(f'<rect class="s1band" x="{float(L.X(-1)):.1f}" y="{L.y0}" width="{float(L.X(9) - L.X(-1)):.1f}" height="{L.h}"/>')
o.append(f'<path class="s1" d="M{float(L.X(-1)):.1f} {float(L.Y(avg1)):.1f} H{float(L.X(9)):.1f}"><title>what the TOA measures: the average, {avg1:.0f} µs</title></path>')
x, y = dot(L, 4, r(4), 's2f', f'one time per TOA (its mean hour angle, 4 h): {r(4):.0f} µs')
o.append(f'<text class="sub" x="{x - 10:.1f}" y="{y + 4:.1f}" text-anchor="end">one time: {r(4):.0f} µs</text>')
o.append(f'<text class="sub" x="{float(L.X(-0.6)):.1f}" y="{float(L.Y(avg1)) + 18:.1f}">average: {avg1:.0f} µs</text>')
o.append(f'<text class="sub" x="{float(L.X(-1)):.1f}" y="{L.y0 + L.h - 8}">data</text>')

# right: two passes a day apart, 5 h each at the same hour angle
R = panel(500, (0, 36), 'Two passes a day apart')
spans = [(3.5, 8.5), (27.5, 32.5)]
Hs = np.concatenate([np.linspace(a, b, 601) for a, b in spans]); avg2 = r(Hs).mean()
for a, b in spans:
    o.append(f'<rect class="s1band" x="{float(R.X(a)):.1f}" y="{R.y0}" width="{float(R.X(b) - R.X(a)):.1f}" height="{R.h}"/>')
    o.append(f'<path class="s1" d="M{float(R.X(a)):.1f} {float(R.Y(avg2)):.1f} H{float(R.X(b)):.1f}"><title>what the TOA measures: the average, {avg2:.0f} µs</title></path>')
o.append(f'<path class="dash" d="M{float(R.X(8.5)):.1f} {float(R.Y(avg2)):.1f} H{float(R.X(27.5)):.1f}"/>')
mid = (3.5 + 32.5) / 2
x, y = dot(R, mid, r(mid), 's2f', f'one time at the middle of the span ({mid:.0f} h): {r(mid):.0f} µs, the opposite sign')
o.append(f'<path class="thin" d="M{x:.1f} {float(R.Y(-80)):.1f} V{y - 7:.1f}"/>')
o.append(f'<text class="sub" x="{x:.1f}" y="{float(R.Y(-62)):.1f}" text-anchor="middle">middle of the span: {r(mid):.0f} µs</text>')
o.append(f'<text class="sub" x="{float(R.X(18)):.1f}" y="{float(R.Y(avg2)) - 8:.1f}" text-anchor="middle">average: +{avg2:.0f} µs</text>')
o.append(f'<text class="sub" x="{float(R.X(3.5)):.1f}" y="{R.y0 + R.h - 8}">data</text>')
o.append(f'<text class="sub" x="{float(R.X(27.5)):.1f}" y="{R.y0 + R.h - 8}">data</text>')

# legend (two series; direct labels as well)
lx, ly = 70, 326
o.append(f'<path class="s1" d="M{lx} {ly} h20"/><text class="sub" x="{lx + 26}" y="{ly + 4}">what the TOA measures: the average over its data</text>')
o.append(f'<circle class="s2f" cx="{lx + 370}" cy="{ly}" r="5"/><text class="sub" x="{lx + 382}" y="{ly + 4}">the model at one time per TOA</text>')
o.append(f'<path class="thin" d="M{lx + 600} {ly} h20"/><text class="sub" x="{lx + 626}" y="{ly + 4}">the daily signal</text>')
print(f'A = {A:.1f} us; one pass: average {avg1:.1f} vs one time {r(4):.1f} ({avg1 / r(4):.2f}); two passes: average {avg2:.1f}, middle {r(mid):.1f}')
svg = '''<figure class="scroll">
<svg class="wide" viewBox="0 0 860 340" role="img" aria-label="The daily position signal of B0329+54 for a site 100 km east, with the time spans that two TOAs average over, and the value one time per TOA would give">
  ''' + '\n  '.join(o) + f'''
</svg>
<figcaption><b>Figure 10.2.</b> Why a TOA's row must be averaged. Grey: the residual a site 100 km east gives B0329+54, a daily sine of {A:.0f} µs against hour
angle (grey bands: the pulsar below 20°). Left: a TOA from 10 hours of one pass measures the average over its data, {avg1:.0f} µs; the model at
the data's mean hour angle says {r(4):.0f} µs, so the fit sees the offset {avg1 / r(4):.2f} times too small. Right: a TOA from 5 hours on each of two days at
the same hour angle measures +{avg2:.0f} µs, but the middle of its span lies at the opposite hour angle, where the model says {r(mid):.0f} µs. With
<code>toa.tSub</code> the fit averages the model over the same data, and both cases are exact.</figcaption>
</figure>'''
open('navavg.svg', 'w').write(svg)
