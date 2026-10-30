import numpy as np
from plotsvg import Ax
from delaymodel import delay, fixed, C, B0329
# The delays, term by term (localisation.html #delays): delayModel's analytic terms for B0329+54 seen from 52 N 5 E,
# over 2026 (the annual terms) and over two days from 1 March 2026 (the rotation), computed by delaymodel.py, a port
# of functions/delayModel.m. c: the size of every term on one log scale, against the TOA noise of the two pulsars
# (B0329+54 734 us, J0437 13 us: localisation.html #fix) and the signal navigation wants (100 km along the line of
# sight, 333 us). Curvature at 1 kpc (the page's table).
DAY = 86400.0
site = fixed(52, 5)
t0 = 9496.5 * DAY                                         # 2026-01-01 00:00 TT
ty = t0 + np.arange(0, 365 * DAY + 1, 6 * 3600.0); doy = (ty - t0) / DAY
dy = delay(ty, site, B0329, distance=3.0857e19)
td = 9555.5 * DAY + np.arange(0, 2 * DAY + 1, 300.0); hr = (td - td[0]) / 3600
dd = delay(td, site, B0329)
geo = (dy['n'] @ dy['rEarth']) / C                      # the geocentre's Roemer delay
rot = (dd['n'] @ dd['rSite']) / C                       # the site's part
rot0 = rot.mean(); rot = rot - rot0                     # its constant part (the axis) merges with the clock
months = [(0, 'Jan'), (90, 'Apr'), (181, 'Jul'), (273, 'Oct'), (365, 'Jan')]
o = []

def xt(A, ticks):
    for v, lab in ticks:
        x = float(A.X(v)); o.append(f'<path class="axis" d="M{x:.1f} {A.y0 + A.h} v5"/>')
        o.append(f'<text class="sub" x="{x:.1f}" y="{A.y0 + A.h + 17}" text-anchor="middle">{lab}</text>')

def panel(A, x, y, yt, unit, ttl, xticks, note):
    for v in yt:
        o.append(f'<path class="grid" d="M{A.x0} {float(A.Y(v)):.1f} H{A.x0 + A.w}"/>')
    o.append(A.frame(yticks=yt, ylab=unit))
    xt(A, xticks)
    o.append(f'<path class="s1" d="{A.path(x, y)}"/>')
    o.append(f'<text class="ttl" x="{A.x0}" y="{A.y0 - 26}">{ttl}</text>')
    o.append(f'<text class="sub" x="{A.x0 + A.w}" y="{A.y0 - 8}" text-anchor="end">{note}</text>')

A = Ax(80, 60, 320, 120, (0, 365), (-480, 480))
panel(A, doy, geo, [-400, -200, 0, 200, 400], 's', 'a · Roemer, the Earth\'s orbit', months,
      f'{geo.min():.0f} to +{geo.max():.0f} s over 2026')
B = Ax(520, 60, 320, 120, (0, 48), (-9, 9))
panel(B, hr, rot * 1e3, [-8, -4, 0, 4, 8], 'ms', 'b · Roemer, the Earth\'s rotation', [(h, f'{h} h') for h in (0, 12, 24, 36, 48)],
      f'±{(rot.max() - rot.min()) / 2 * 1e3:.2f} ms, once per sidereal day')
Cc = Ax(80, 270, 320, 120, (0, 365), (-2, 2))
panel(Cc, doy, dy['einstein'] * 1e3, [-2, -1, 0, 1, 2], 'ms', 'c · Einstein: TDB − TT', months,
      f'±{(dy["einstein"].max() - dy["einstein"].min()) / 2 * 1e3:.2f} ms')
D = Ax(520, 270, 320, 120, (0, 365), (-20, 10))
panel(D, doy, dy['shapiro'] * 1e6, [-20, -10, 0, 10], 'µs', 'd · Shapiro: the Sun\'s gravity', months,
      f'{dy["shapiro"].min() * 1e6:.1f} to +{dy["shapiro"].max() * 1e6:.1f} µs')
k = int(np.argmin(dy['shapiro']))
o.append(f'<circle class="s1f" cx="{float(D.X(doy[k])):.1f}" cy="{float(D.Y(dy["shapiro"][k] * 1e6)):.1f}" r="4"/>')
o.append(f'<text class="sub" x="{float(D.X(doy[k])) + 9:.1f}" y="{float(D.Y(dy["shapiro"][k] * 1e6)) + 4:.1f}">Sun closest to the line of sight</text>')

# e: every term on one log scale, one row each
E = Ax(240, 470, 590, 120, (-7, 3), (0, 1))
yA = E.y0 + E.h
for e in range(-7, 4):
    x = float(E.X(e)); o.append(f'<path class="grid" d="M{x:.1f} {E.y0 - 6} V{yA}"/>')
    o.append(f'<path class="axis" d="M{x:.1f} {yA} v5"/>')
    lab = {-7: '100 ns', -6: '1 µs', -5: '10 µs', -4: '100 µs', -3: '1 ms', -2: '10 ms', -1: '100 ms', 0: '1 s', 1: '10 s', 2: '100 s', 3: '1000 s'}[e]
    o.append(f'<text class="sub" x="{x:.1f}" y="{yA + 17}" text-anchor="middle">{lab}</text>')
o.append(f'<path class="axis" d="M{E.x0} {yA} H{E.x0 + E.w}"/>')
o.append(f'<text class="ttl" x="80" y="452">e · every term on one scale: its largest size</text>')
curv = np.max(np.abs(dy['curvature']))
terms = [(-geo.min(), 'Roemer, orbit', f'{-geo.min():.0f} s', 's1f'), ((rot.max() - rot.min()) / 2, 'Roemer, rotation', '7.6 ms', 's1f'),
         (1.66e-3, 'Einstein', '1.66 ms', 's1f'), (333e-6, 'ship, 100 km', '333 µs', 's2f'),
         (-dy['shapiro'].min(), 'Shapiro', '17 µs', 's1f'), (curv, 'curvature, 1 kpc', f'{curv * 1e6:.1f} µs', 's1f')]
for i, (v, lab, val, cls) in enumerate(terms):
    y = E.y0 + 8 + i * 20; x = float(E.X(np.log10(v)))
    col = ' style="fill:var(--series-2)"' if cls == 's2f' else ''
    o.append(f'<text class="sub" x="80" y="{y + 4}"{col}>{lab}</text>')
    o.append(f'<path class="thin" d="M{E.x0} {y} H{x:.1f}"/>')
    o.append(f'<circle class="{cls}" cx="{x:.1f}" cy="{y}" r="5"><title>{lab}: {val}</title></circle>')
    xv, anc = (x - 9, 'end') if cls == 's2f' else (x + 9, 'start')    # the signal's value left of its dot, clear of the dashed line
    o.append(f'<text class="sub" x="{xv:.1f}" y="{y + 4}" text-anchor="{anc}"{col}>{val}</text>')
for v, lab, anc, dx in ((734e-6, 'B0329+54 TOA noise', 'start', 5), (13e-6, 'J0437 TOA noise', 'end', -5)):
    x = float(E.X(np.log10(v)))
    o.append(f'<path class="dash" d="M{x:.1f} {E.y0 - 6} V{yA + 26}"/>')
    o.append(f'<text class="sub" x="{x + dx:.1f}" y="{yA + 36}" text-anchor="{anc}">{lab}</text>')

svg = '''<figure class="scroll">
<svg class="wide" viewBox="0 0 860 650" role="img" aria-label="The terms of the delay for B0329+54 from 52 degrees north: the orbital Roemer delay over 2026, the rotational Roemer delay over two days, the Einstein and Shapiro delays over 2026, and every term on one log scale against the TOA noise of B0329+54 and J0437">
  ''' + '\n  '.join(o) + '''
</svg>
<figcaption><b>Figure L.3.</b> The terms of the delay for B0329+54 seen from 52° N, computed with <code>delayModel</code>'s analytic Earth. (a) The Earth's
orbit, over 2026; not symmetric about zero, because the Sun itself moves around the barycentre. (b) The Earth's rotation, over two days from
1 March 2026, about its mean; the mean, 13.6 ms, is the site's offset along the Earth's axis, which does not turn. (c) The Einstein term and
(d) the Shapiro term over 2026; the Shapiro term dips when the Sun passes closest to the line of sight. (e) The largest size of each term on
one scale, with the 100 km signal that navigation wants in orange, and the TOA noise of the two pulsars dashed. Generated by
<code>docs/tools/figures/gen_delays.py</code>, with a port of <code>delayModel</code> (<code>delaymodel.py</code>).</figcaption>
</figure>'''
open('delays.svg', 'w').write(svg)
print('orbit', geo.min(), geo.max(), 'rot', (rot.max() - rot.min()) / 2, 'rot0', rot0, 'einstein', dy['einstein'].min(), dy['einstein'].max(),
      'shapiro', dy['shapiro'].min(), dy['shapiro'].max(), 'day', doy[k], 'curv', curv)
