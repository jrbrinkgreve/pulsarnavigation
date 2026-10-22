import os
import numpy as np
import scipy.io as sio
from plotsvg import Ax
# Experiment 3a (tests/expBlankingVariance.m, 6 Oct; its saved results data/mc/expBlankingVariance.mat):
# variance of a detected bin after blanking + per-channel dedispersion, grouped by valid fraction W.
# For each group: the two simple rules (var0*W^2, var0*W) and the measurement, divided by the exact value.
m = sio.loadmat(os.path.join(os.path.dirname(os.path.abspath(__file__)), '../../../data/mc/expBlankingVariance.mat'), squeeze_me=True, struct_as_record=False)['results']
o = []
A = Ax(90, 40, 720, 250, (0, 1), (0.3, 3.2), logy=True)
o.append(A.frame(xticks=[0, 0.25, 0.5, 0.75, 1], yticks=[0.3, 0.5, 1, 2, 3], xlab='valid fraction W of the detected bin',
                 ylab='variance ÷ exact variance', yfmt=lambda v: f'{v:g}'))
o.append(f'<path class="dash" d="M{A.x0} {A.Y(1):.1f} H{A.x0 + A.w}"/>')
worst = {'w2': 1, 'w1': 1}
for r in m:
    lo, hi = r.Wrange
    if lo < 0.02 or hi > 0.9999:
        continue                                     # skip the empty-bin and unblanked groups
    x = float(A.X((lo + hi)/2))
    for key, val, cls in (('w2', r.varW2, 'dot2'), ('w1', r.varW1, 'bad'), ('meas', r.varMeas, 'dot')):
        q = val/r.varExact
        if key in worst and lo >= 0.25:
            worst[key] = min(worst[key], q) if key == 'w2' else max(worst[key], q)
        if cls == 'bad':
            o.append(f'<rect class="bad" x="{x - 3:.1f}" y="{A.Y(q) - 3:.1f}" width="6" height="6"/>')
        else:
            o.append(f'<path class="{cls}" d="M{x - 3.5:.1f} {A.Y(q):.1f} a3.5 3.5 0 1 0 7 0 a3.5 3.5 0 1 0 -7 0"/>')
print('for W >= 0.25: W^2 rule down to', worst['w2'], ', W rule up to', worst['w1'])
lx, ly = A.x0 + 470, A.y0 + 14
o.append(f'<rect class="bad" x="{lx}" y="{ly - 7}" width="6" height="6"/><text class="sub" x="{lx + 12}" y="{ly}">rule "var ∝ W" (whole samples lost)</text>')
o.append(f'<path class="dot" d="M{lx - 0.5} {ly + 15} a3.5 3.5 0 1 0 7 0 a3.5 3.5 0 1 0 -7 0"/><text class="sub" x="{lx + 12}" y="{ly + 19}">measured</text>')
o.append(f'<path class="dot2" d="M{lx - 0.5} {ly + 34} a3.5 3.5 0 1 0 7 0 a3.5 3.5 0 1 0 -7 0"/><text class="sub" x="{lx + 12}" y="{ly + 38}">rule "var ∝ W²" (amplitude scaled)</text>')
o.append(f'<text class="sub" x="{A.X(0.5):.1f}" y="{A.Y(1) - 6:.1f}" text-anchor="middle">exact (blankingWeights)</text>')
svg = '''<figure class="scroll">
<svg class="wide" viewBox="0 0 860 330" role="img" aria-label="Experiment 3a: the variance of blanked detected bins against the two simple rules and the exact formula">
  ''' + '\\n  '.join(o) + '''
</svg>
<figcaption>The experiment behind this choice (<code>tests/expBlankingVariance.m</code>, from its saved results). It shows the
variance of a detected bin after blanking and per-channel dedispersion, for three masks at DM 5 and DM 100, grouped by the bin's valid
fraction <var>W</var>. Divided by the exact value, the measurement sits on 1. The two simple rules do not, and which one comes closer
depends on the mask. Left out: groups with <var>W</var> &lt; 0.02, where the <var>W</var> rule is 30× too high, and unblanked bins.</figcaption>
</figure>'''
open('blankvar.svg', 'w').write(svg.replace('\\n', '\n'))
