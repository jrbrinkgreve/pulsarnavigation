import numpy as np
from plotsvg import Ax, pow10
P = [('J0437−4715', 21, 12, 1, (8, -8), 'start'), ('Vela', 0.3, 146, 1, (8, -8), 'start'), ('B0329+54', 4.2, 566, 1, (8, -8), 'start'),
     ('B0950+08', 67, 738, 0, (-8, -9), 'end'), ('B1749−28', 96, 566, 0, (8, 14), 'start'), ('B1642−03', 278, 326, 0, (8, 4), 'start'),
     ('B1929+10', 571, 480, 0, (8, -6), 'start'), ('B1937+21', 3600, 5, 0, (-8, 4), 'end'), ('J1713+0747', 20000, 26, 0, (8, 14), 'start'),
     ('J2145−0750', 16000, 34, 0, (-8, -9), 'end')]
for n, t, s, *_ in P: print(n, round(s*np.sqrt(t)))
A = Ax(80, 30, 720, 250, (0.1, 1e5), (2, 3000), logy=True)
# log x: wrap X
import math
lx0, lx1 = math.log10(0.1), math.log10(1e5)
X = lambda t: A.x0 + (np.log10(t) - lx0)/(lx1 - lx0)*A.w
o = []
o.append(f'<path class="axis" d="M{A.x0} {A.y0+A.h} H{A.x0+A.w} M{A.x0} {A.y0} V{A.y0+A.h}"/>')
for t in [0.1, 1, 10, 100, 1000, 1e4, 1e5]:
    x = X(t); o.append(f'<path class="axis" d="M{x:.1f} {A.y0+A.h} v5"/><text class="sub" x="{x:.1f}" y="{A.y0+A.h+18}" text-anchor="middle">{t:g}</text>')
for s in [10, 100, 1000]:
    y = float(A.Y(s)); o.append(f'<path class="axis" d="M{A.x0} {y:.1f} h-5"/><text class="sub" x="{A.x0-8}" y="{y+4:.1f}" text-anchor="end">{s:g}</text>')
o.append(f'<text class="sub" x="{A.x0+A.w/2}" y="{A.y0+A.h+36}" text-anchor="middle">time per TOA, worst case [h] (SNR 7, dual polarization)</text>')
o.append(f'<text class="sub" x="{A.x0}" y="{A.y0-10}">σ_TOA [µs]</text>')
tt = np.logspace(-1, 5, 200)
for F, lab in [(55, 'σ·√t = 55 µs·√h (J0437)'), (1160, 'σ·√t = 1160 (B0329+54)')]:
    ss = F/np.sqrt(tt); m = (ss >= 2) & (ss <= 3000)
    pts = ' L'.join(f'{X(a):.1f} {float(A.Y(b)):.1f}' for a, b in zip(tt[m], ss[m]))
    o.append(f'<path class="dash" d="M{pts}"/>')
    tl = 2.0 if F == 55 else 600
    o.append(f'<text class="sub" x="{X(tl):.1f}" y="{float(A.Y(F/np.sqrt(tl)))-8:.1f}" transform="rotate(-18 {X(tl):.1f} {float(A.Y(F/np.sqrt(tl)))-8:.1f})">{lab}</text>')
xd = X(24); o.append(f'<path class="thin" stroke-dasharray="2 3" d="M{xd:.1f} {A.y0} V{A.y0+A.h}"/><text class="sub" x="{xd+4:.1f}" y="{A.y0+12}">1 day</text>')
for n, t, s, bold, (dx, dy), anc in P:
    x, y = X(t), float(A.Y(s))
    o.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{5.5 if bold else 4}" class="{"dot2" if bold else "dot"}"/>')
    o.append(f'<text class="{"lbl" if bold else "sub"}" x="{x+dx:.1f}" y="{y+dy:.1f}" text-anchor="{anc}"{" font-weight=\"600\"" if bold else ""}>{n}</text>')
svg = '''<figure class="scroll">
<svg class="wide" viewBox="0 0 860 330" role="img" aria-label="Time per TOA against TOA precision for the candidate pulsars at the worst-case design point">
  ''' + '\n  '.join(o) + '''
</svg>
<figcaption>The table above as a map (worst case, both axes logarithmic). Down is more precise, left is faster. Along a dashed line the
figure of merit σ<sub>TOA</sub>·√<var>t</var> is constant: averaging <var>N</var> TOAs trades time for precision along it. J0437 and Vela
lead; B0329+54 is the only northern source with TOAs within a pass; the faint millisecond pulsars would be precise but take months.</figcaption>
</figure>'''
open('scenario.svg', 'w').write(svg)
