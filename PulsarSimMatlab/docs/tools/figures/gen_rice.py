import numpy as np
from scipy.stats import norm
from scipy.optimize import brentq
from plotsvg import Ax, pow10
# detectPulsar's thresholds against the false-alarm probability (detection.html): known phase eta = Phi^-1(1 - PFA);
# unknown phase from Rice's formula P(Tmax > eta) = Q(eta) + (sqrt(lambda2)/2pi) exp(-eta^2/2), Rice factor 5.51 for the
# 5 % Gaussian template; and two naive alternatives: 5.5 or 2048 independent trials (Bonferroni).
rice = 5.51
pfa = np.logspace(-1, -9, 200)
known = norm.isf(pfa)
unk = np.array([brentq(lambda e: norm.sf(e) + rice*np.exp(-e**2/2) - p, 0.5, 12) for p in pfa])
ind55 = norm.isf(pfa/5.5)
ind2048 = norm.isf(pfa/2048)
print('PFA 1e-3: known', norm.isf(1e-3), 'unknown', brentq(lambda e: norm.sf(e) + rice*np.exp(-e**2/2) - 1e-3, 0.5, 12),
      '5.5 trials', norm.isf(1e-3/5.5), '2048 bins', norm.isf(1e-3/2048))
o = []
A = Ax(90, 40, 560, 250, (-9, -1), (1, 7))
o.append(A.frame(xticks=[-9, -8, -7, -6, -5, -4, -3, -2, -1], yticks=[1, 2, 3, 4, 5, 6, 7], xlab='false-alarm probability per profile (log₁₀)',
                 ylab='threshold η (in noise standard deviations)', xfmt='{:g}'))
lp = np.log10(pfa)
o.append(f'<rect class="fillwhy" x="{A.X(-3):.1f}" y="{A.y0}" width="{A.X(-2) - A.X(-3):.1f}" height="{A.h}"/>')
o.append(f'<text class="sub" x="{A.X(-2.5):.1f}" y="{A.y0 + 12}" text-anchor="middle">validated</text>')
o.append(f'<path class="thin" d="M{A.X(-6):.1f} {A.y0} V{A.y0 + A.h}"/>')
o.append(f'<text class="sub" x="{A.X(-6) + 4:.1f}" y="{A.y0 + 12}">10⁻⁶: navigation (extrapolated)</text>')
for y, cls, lab, dash in ((known, 'curve', 'known phase (T0)', ''), (unk, 'curve2', 'unknown phase (Tmax, Rice)', ''),
                          (ind55, 'thin', '5.5 independent trials (too low)', ' stroke-dasharray="4 3"'), (ind2048, 'thin', '2048 independent bins (too high)', ' stroke-dasharray="1 3"')):
    o.append(f'<path class="{cls}"{dash} d="{A.path(lp, y)}"/>')
    i = [known is y, unk is y, ind55 is y, ind2048 is y].index(True)
    ly = float(A.Y(3.4)) + 16*i
    o.append(f'<path class="{cls}"{dash} d="M{A.X(-8.8):.1f} {ly - 4:.1f} h22"/><text class="sub" x="{A.X(-8.8) + 28:.1f}" y="{ly:.1f}">{lab}</text>')
for p in (1e-3,):
    for e, cls in ((norm.isf(p), 'dot'), (brentq(lambda e: norm.sf(e) + rice*np.exp(-e**2/2) - p, 0.5, 12), 'dot2')):
        o.append(f'<path class="{cls}" d="M{A.X(-3) - 4:.1f} {A.Y(e):.1f} a4 4 0 1 0 8 0 a4 4 0 1 0 -8 0"/>')
        o.append(f'<text class="sub" x="{A.X(-3) - 8:.1f}" y="{A.Y(e) + 4:.1f}" text-anchor="end">{e:.2f}</text>')
nx = 672
notes = [('ttl', 'The price of not knowing'), ('ttl', 'the phase'), ('sub', 'Searching all 2048 shifts gives'), ('sub', 'noise many chances, so the'),
         ('sub', 'threshold must rise: 3.09 → 4.15'), ('sub', 'at 10⁻³. Neighbouring shifts are'), ('sub', 'almost the same test, so 2048'),
         ('sub', 'independent trials would be far'), ('sub', 'too strict; Rice counts how often'), ('sub', 'a smooth noise curve crosses η.')]
for i, (c, line) in enumerate(notes):
    o.append(f'<text class="{c}" x="{nx}" y="{58 + 17*i}">{line}</text>')
svg = '''<figure class="scroll">
<svg class="wide" viewBox="0 0 860 330" role="img" aria-label="Detection thresholds against false-alarm probability: known phase, unknown phase from Rice's formula, and two naive counts of independent trials">
  ''' + '\n  '.join(o) + '''
</svg>
<figcaption>The thresholds of <code>detectPulsar</code> against the false-alarm probability, for the 5 % Gaussian template in 2048 phase bins
(Rice factor 5.51). Known phase: one Gaussian test. Unknown phase: Rice's formula for the maximum over all shifts. The two dotted lines show
why neither naive count works: 5.5 independent trials gives too low a threshold (3.57 at 10⁻³: too many false alarms), 2048 gives too high
a one (4.90: detections lost). The false-alarm rate is validated at 10⁻²–10⁻³ (green); 10⁻⁶ relies on the Gaussian and Rice tails.</figcaption>
</figure>'''
open('rice.svg', 'w').write(svg)
