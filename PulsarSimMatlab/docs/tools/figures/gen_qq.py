import numpy as np
from scipy.stats import norm
from plotsvg import Ax
# Reading validateTOA (validation.html): a: Q-Q plots of normalized errors z for three simulated cases (illustrative):
# honest error bars, error bars 1.5x too small, and 5 % outliers; b: the uncertainty of the rms ratio, ~1/sqrt(2n).
rng = np.random.default_rng(3)
n = 90
q = norm.ppf((np.arange(1, n + 1) - 0.5)/n)
cases = [('honest: on the line', np.sort(rng.standard_normal(n)), 'dot'),
         ('error bars 1.5× too small: steeper', np.sort(1.5*rng.standard_normal(n)), 'dot2'),
         ('5 % outliers: tails bend away', None, 'bad')]
z3 = rng.standard_normal(n); idx = rng.choice(n, 5, replace=False); z3[idx] = rng.choice([-1, 1], 5)*rng.uniform(4, 7, 5)
cases[2] = (cases[2][0], np.sort(z3), 'bad')
o = []
A = Ax(70, 46, 330, 230, (-3, 3), (-6, 6))
o.append(A.frame(xticks=[-3, -2, -1, 0, 1, 2, 3], yticks=[-6, -3, 0, 3, 6], xlab='N(0,1) quantile', ylab='sorted z = error ÷ stated error'))
o.append(f'<path class="dash" d="M{A.X(-3):.1f} {A.Y(-3):.1f} L{A.X(3):.1f} {A.Y(3):.1f}"/>')
for i, (lab, z, cls) in enumerate(cases):
    for a, b in zip(q, z):
        b = np.clip(b, -6, 6)
        if cls == 'bad':
            o.append(f'<rect class="bad" x="{A.X(a) - 2:.1f}" y="{A.Y(b) - 2:.1f}" width="4" height="4"/>')
        else:
            o.append(f'<path class="{cls}" d="M{A.X(a) - 2.2:.1f} {A.Y(b):.1f} a2.2 2.2 0 1 0 4.4 0 a2.2 2.2 0 1 0 -4.4 0"/>')
    ly = A.y0 + 10 + 15*i
    sym = (f'<rect class="bad" x="{A.x0 + 8}" y="{ly - 6}" width="5" height="5"/>' if cls == 'bad'
           else f'<path class="{cls}" d="M{A.x0 + 8} {ly - 3} a2.5 2.5 0 1 0 5 0 a2.5 2.5 0 1 0 -5 0"/>')
    o.append(sym + f'<text class="sub" x="{A.x0 + 18}" y="{ly}">{lab}</text>')
o.append(f'<text class="ttl" x="{A.x0}" y="20">a · the Q-Q plot (90 simulated TOAs)</text>')
B = Ax(520, 46, 300, 230, (5, 500), (0, 0.45))
nn = np.logspace(np.log10(5), np.log10(500), 200)
lxb = lambda v: B.x0 + (np.log10(v) - np.log10(5))/(np.log10(500) - np.log10(5))*B.w
o.append(B.frame(yticks=[0, 0.1, 0.2, 0.3, 0.4], ylab='1σ uncertainty of the ratio', yfmt=lambda v: f'±{int(v*100)} %' if v else '0'))
for v in (5, 10, 20, 50, 100, 200, 500):
    o.append(f'<path class="axis" d="M{lxb(v):.1f} {B.y0 + B.h} v5"/><text class="sub" x="{lxb(v):.1f}" y="{B.y0 + B.h + 18}" text-anchor="middle">{v}</text>')
o.append(f'<text class="sub" x="{B.x0 + B.w/2:.1f}" y="{B.y0 + B.h + 34}" text-anchor="middle">number of TOAs</text>')
pts = ' L'.join(f'{lxb(a):.1f} {float(B.Y(1/np.sqrt(2*a))):.1f}' for a in nn)
o.append(f'<path class="curve" d="M{pts}"/>')
for a, lab in ((8, '8 TOAs (main.m): ±25 %'), (90, '90 TOAs (Monte Carlo): ±7.5 %')):
    x, y = lxb(a), float(B.Y(1/np.sqrt(2*a)))
    o.append(f'<path class="dot2" d="M{x - 4:.1f} {y:.1f} a4 4 0 1 0 8 0 a4 4 0 1 0 -8 0"/>')
    anc = 'end' if a > 50 else 'start'
    o.append(f'<text class="sub" x="{x + (-8 if anc == "end" else 8):.1f}" y="{y - 8:.1f}" text-anchor="{anc}">{lab}</text>')
o.append(f'<text class="ttl" x="{B.x0}" y="20">b · how well a few TOAs test the error bars</text>')
svg = '''<figure class="scroll">
<svg class="wide" viewBox="0 0 860 330" role="img" aria-label="Reading validateTOA: Q-Q plots for honest, too small and outlier-ridden error bars, and the uncertainty of the rms ratio against the number of TOAs">
  ''' + '\n  '.join(o) + '''
</svg>
<figcaption><b>Figure 9.2.</b> How to read <code>validateTOA</code> (illustrative, not pipeline output). (a) Q-Q plots of normalized errors: honest error
bars lie on the line, too small ones make it steeper, outliers bend the tails away. (b) The rms ratio is itself uncertain by about
1/√(2<var>n</var>), so 8 TOAs can only show a ratio far from 1, while a Monte Carlo of 90 TOAs tests it to ±7.5 %. The validation record
therefore relies on Monte Carlo runs and sweeps, not on one <code>main.m</code> run.</figcaption>
</figure>'''
open('qq.svg', 'w').write(svg)
