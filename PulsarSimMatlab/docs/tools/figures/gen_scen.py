import numpy as np
from plotsvg import Ax, pow10
# Time per TOA [h] and sigma_TOA [us] of scenario.html's table R2: Eq. (R.1) at the array's design point of 23 Nov 2026
# (SEFD 6.44e4 Jy, 150 MHz, two polarizations, SNR 7, W_eq = 2 W50). Until 23 Nov 2026 this figure showed the worst case of
# the 0.61 m^2 receiver (7.4e5 Jy, 400 MHz): 21, 0.3, 4.2, 67, 96, 278, 571, 3600, 20000, 16000 h.
P = [('J0437−4715', 0.42, 12, 1, (8, -8), 'start'), ('Vela', 0.007, 146, 1, (8, -8), 'start'), ('B0329+54', 0.09, 566, 1, (8, -8), 'start'),
     ('B0950+08', 1.3, 738, 0, (-8, -9), 'end'), ('B1749−28', 1.9, 566, 0, (8, 14), 'start'), ('B1642−03', 5.6, 326, 0, (8, 4), 'start'),
     ('B1929+10', 11.5, 480, 0, (8, -6), 'start'), ('B1937+21', 72, 5, 0, (-8, 4), 'end'), ('J1713+0747', 400, 26, 0, (8, 14), 'start'),
     ('J2145−0750', 320, 34, 0, (-8, -9), 'end')]
for n, t, s, *_ in P: print(n, round(s*np.sqrt(t)))
A = Ax(80, 30, 720, 250, (0.001, 1e3), (2, 3000), logy=True)
# log x: wrap X
import math
lx0, lx1 = math.log10(0.001), math.log10(1e3)
X = lambda t: A.x0 + (np.log10(t) - lx0)/(lx1 - lx0)*A.w
o = []
o.append(f'<path class="axis" d="M{A.x0} {A.y0+A.h} H{A.x0+A.w} M{A.x0} {A.y0} V{A.y0+A.h}"/>')
for t in [0.001, 0.01, 0.1, 1, 10, 100, 1000]:
    x = X(t); o.append(f'<path class="axis" d="M{x:.1f} {A.y0+A.h} v5"/><text class="sub" x="{x:.1f}" y="{A.y0+A.h+18}" text-anchor="middle">{t:g}</text>')
for s in [10, 100, 1000]:
    y = float(A.Y(s)); o.append(f'<path class="axis" d="M{A.x0} {y:.1f} h-5"/><text class="sub" x="{A.x0-8}" y="{y+4:.1f}" text-anchor="end">{s:g}</text>')
o.append(f'<text class="sub" x="{A.x0+A.w/2}" y="{A.y0+A.h+36}" text-anchor="middle">time per TOA [h] (SNR 7, dual polarization, 150 MHz)</text>')
o.append(f'<text class="sub" x="{A.x0}" y="{A.y0-10}">σ_TOA [µs]</text>')
tt = np.logspace(-3, 3, 200)
for F, lab in [(8, 'σ·√t = 8 µs·√h (J0437)'), (170, 'σ·√t = 170 (B0329+54)')]:
    ss = F/np.sqrt(tt); m = (ss >= 2) & (ss <= 3000)
    pts = ' L'.join(f'{X(a):.1f} {float(A.Y(b)):.1f}' for a, b in zip(tt[m], ss[m]))
    o.append(f'<path class="dash" d="M{pts}"/>')
    tl = 0.04 if F == 8 else 0.4
    o.append(f'<text class="sub" x="{X(tl):.1f}" y="{float(A.Y(F/np.sqrt(tl)))-8:.1f}" transform="rotate(-18 {X(tl):.1f} {float(A.Y(F/np.sqrt(tl)))-8:.1f})">{lab}</text>')
xd = X(24); o.append(f'<path class="thin" stroke-dasharray="2 3" d="M{xd:.1f} {A.y0} V{A.y0+A.h}"/><text class="sub" x="{xd+4:.1f}" y="{A.y0+12}">1 day</text>')
for n, t, s, bold, (dx, dy), anc in P:
    x, y = X(t), float(A.Y(s))
    o.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{5.5 if bold else 4}" class="{"dot2" if bold else "dot"}"/>')
    o.append(f'<text class="{"lbl" if bold else "sub"}" x="{x+dx:.1f}" y="{y+dy:.1f}" text-anchor="{anc}"{" font-weight=\"600\"" if bold else ""}>{n}</text>')
svg = '''<figure class="scroll">
<svg class="wide" viewBox="0 0 860 330" role="img" aria-label="Time per TOA against TOA precision for the candidate pulsars at the array's design point">
  ''' + '\n  '.join(o) + '''
</svg>
<figcaption><b>Figure R.1.</b> The table above as a map (the design point, both axes logarithmic). Down is more precise, left is faster. Along a dashed line the
figure of merit σ<sub>TOA</sub>·√<var>t</var> is constant: averaging <var>N</var> TOAs trades time for precision along it. J0437 and Vela
lead; B0329+54 is by far the fastest northern source; the faint millisecond pulsars would be precise but need tens to hundreds of hours
of data.</figcaption>
</figure>'''
open('scenario.svg', 'w').write(svg)
