import numpy as np
from plotsvg import Ax
# Which extra pulsars help (navigation.html #plan): the planning study's numbers (tests/expNavPlan.m, logged 22 Oct 2026):
# the polar component's 1-sigma error with a calibrated start and the crystal clock, nominal flux, from 52 N.
sets = ['B0329+54 alone', '+ B0950+08', '+ B0950+08, B1642−03, B1929+10']
polar = {30: [63.7, 49.1, 39.8], 100: [63.6, 34.3, 23.7]}       # [km], crystal clock
perfect100 = 8.9                                                # [km], B0329+54 alone, perfect clock, 100 days
o = []
A = Ax(110, 40, 650, 230, (0, 3), (0, 70))
for v in (20, 40, 60):
    o.append(f'<path class="grid" d="M{A.x0} {float(A.Y(v)):.1f} H{A.x0 + A.w}"/>')
o.append(A.frame(xticks=[], yticks=[0, 20, 40, 60], ylab='polar error, 1σ [km]'))
bw = 0.26
for i, name in enumerate(sets):
    for j, (days, cls) in enumerate(((30, 's1f'), (100, 's2f'))):
        v = polar[days][i]
        x0 = i + 0.5 + (j - 1) * bw + 0.01 * (1 if j else -1)
        X0, X1, Y = float(A.X(x0)), float(A.X(x0 + bw)), float(A.Y(v))
        yb = float(A.Y(0)); rr = 4
        o.append(f'<path class="{cls}" d="M{X0:.1f} {yb:.1f} V{Y + rr:.1f} Q{X0:.1f} {Y:.1f} {X0 + rr:.1f} {Y:.1f} H{X1 - rr:.1f} '
                 f'Q{X1:.1f} {Y:.1f} {X1:.1f} {Y + rr:.1f} V{yb:.1f} Z"><title>{name}, {days} days: ±{v:.1f} km</title></path>')
        o.append(f'<text class="sub" x="{(X0 + X1) / 2:.1f}" y="{Y - 6:.1f}" text-anchor="middle">{v:.0f}</text>')
    o.append(f'<text class="sub" x="{float(A.X(i + 0.5)):.1f}" y="{A.y0 + A.h + 18}" text-anchor="middle">{name}</text>')
y = float(A.Y(perfect100))
o.append(f'<path class="dash" d="M{A.x0} {y:.1f} H{A.x0 + A.w}"/>')
o.append(f'<text class="sub" x="{A.x0 + A.w + 6}" y="{y + 4:.1f}">perfect clock:</text>')
o.append(f'<text class="sub" x="{A.x0 + A.w + 6}" y="{y + 18:.1f}">{perfect100:.0f} km</text>')
lx, ly = A.x0 + A.w - 190, A.y0 + 2
o.append(f'<rect class="s1f" x="{lx}" y="{ly - 6}" width="14" height="12" rx="3"/><text class="sub" x="{lx + 20}" y="{ly + 4}">30 days</text>')
o.append(f'<rect class="s2f" x="{lx + 90}" y="{ly - 6}" width="14" height="12" rx="3"/><text class="sub" x="{lx + 110}" y="{ly + 4}">100 days</text>')
svg = '''<figure class="scroll">
<svg class="wide" viewBox="0 0 860 310" role="img" aria-label="The polar error with a crystal clock for B0329+54 alone and with extra pulsars, after 30 and 100 days">
  ''' + '\n  '.join(o) + '''
</svg>
<figcaption>What extra pulsars buy: the error of the polar component (along the Earth's axis), with a calibrated start and the crystal
clock, from 52° N at nominal flux. B0329+54 alone cannot tell the polar offset from the clock's wander, and more days do not help
(64 km after 30 and after 100 days). Pulsars at other declinations separate the two: B0950+08 brings it to 34 km in 100 days, all
three extra pulsars to 24 km. With a perfect clock B0329+54 alone already reaches 9 km. Numbers: the planning study
(<code>tests/expNavPlan.m</code>), from <code>navFit</code>'s <code>'Predict'</code>.</figcaption>
</figure>'''
open('navplan.svg', 'w').write(svg)
