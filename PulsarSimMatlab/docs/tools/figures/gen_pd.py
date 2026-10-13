import numpy as np
from scipy.stats import norm
from plotsvg import Ax
# Detection probability against the matched-filter SNR of a sub-int (detection.html): P_D = Phi(SNR - eta) with the
# known-phase (3.09) and unknown-phase (4.15) thresholds at PFA 1e-3 (the unknown-phase curve is this Gaussian approximation),
# and the fraction of good TOAs measured in the threshold sweep of 6 Oct (results.html): 40 % at 3.8, 92 % at 6.0, 100 % at 7.5.
snr = np.linspace(0, 10, 401)
o = []
A = Ax(90, 40, 560, 230, (0, 10), (0, 1.05))
o.append(A.frame(xticks=[0, 2, 4, 6, 8, 10], yticks=[0, 0.25, 0.5, 0.75, 1], xlab='SNR of one sub-int (matched filter)',
                 ylab='probability'))
o.append(f'<rect class="fillwhy" x="{A.X(6):.1f}" y="{A.y0}" width="{A.X(7) - A.X(6):.1f}" height="{A.h}"/>')
o.append(f'<text class="sub" x="{A.X(6.5):.1f}" y="{A.y0 + A.h - 22}" text-anchor="middle">design rule</text>')
o.append(f'<text class="sub" x="{A.X(6.5):.1f}" y="{A.y0 + A.h - 8}" text-anchor="middle">SNR 6–7</text>')
for eta, cls, lab in ((3.09, 'curve', 'detected, known phase (η = 3.09)'), (4.15, 'curve2', 'detected, unknown phase (η = 4.15)')):
    o.append(f'<path class="{cls}" d="{A.path(snr, norm.cdf(snr - eta))}"/>')
meas = [(3.8, 0.40), (6.0, 0.92), (7.5, 1.00)]
for s, f in meas:
    o.append(f'<path class="dot" d="M{A.X(s) - 4.5:.1f} {A.Y(f):.1f} a4.5 4.5 0 1 0 9 0 a4.5 4.5 0 1 0 -9 0"/>')
    o.append(f'<text class="sub" x="{A.X(s) + 8:.1f}" y="{A.Y(f) + 14:.1f}">{int(f*100)} %</text>')
lx, ly = A.x0 + 10, A.y0 + 20
o.append(f'<path class="curve" d="M{lx} {ly} h18"/><text class="sub" x="{lx + 24}" y="{ly + 4}">detected, known phase (η = 3.09)</text>')
o.append(f'<path class="curve2" d="M{lx} {ly + 16} h18"/><text class="sub" x="{lx + 24}" y="{ly + 20}">detected, unknown phase (η = 4.15)</text>')
o.append(f'<path class="dot" d="M{lx + 5} {ly + 32} a4 4 0 1 0 8 0 a4 4 0 1 0 -8 0"/><text class="sub" x="{lx + 24}" y="{ly + 36}">good TOAs, measured (6 Oct)</text>')
nx = 672
notes = [('ttl', 'Where SNR ≥ 6–7 comes from'), ('sub', 'Below ~4 a sub-int is often not'), ('sub', 'detected: its fitted peak may be a'),
         ('sub', 'noise peak, so the TOA is dropped.'), ('sub', 'At 6–7 almost every sub-int is'), ('sub', 'detected and its TOA is good, with'),
         ('sub', 'no outliers among the good ones.'), ('sub', 'So a sub-int is made long enough'), ('sub', 'to reach SNR 6–7 (fold page).')]
for i, (c, line) in enumerate(notes):
    o.append(f'<text class="{c}" x="{nx}" y="{58 + 17*i}">{line}</text>')
svg = '''<figure class="scroll">
<svg class="wide" viewBox="0 0 860 310" role="img" aria-label="Detection probability against sub-int SNR for known and unknown phase, with the measured fraction of good TOAs">
  ''' + '\n  '.join(o) + '''
</svg>
<figcaption>Detection probability against the SNR of a sub-int at a false-alarm probability of 10⁻³, from <var>P</var><sub>D</sub> =
Φ(SNR − η) (for the unknown phase this is an approximation), and the fraction of good TOAs measured in the threshold sweep of 6 Oct
(<a href="results.html#threshold">results</a>). The good TOAs lie close to the unknown-phase curve, because a good TOA must pass the
unknown-phase test (and the χ² flag).</figcaption>
</figure>'''
open('pd.svg', 'w').write(svg)
