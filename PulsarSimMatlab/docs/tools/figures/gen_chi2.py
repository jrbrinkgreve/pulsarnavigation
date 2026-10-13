import numpy as np
from scipy.stats import chi2
from plotsvg import Ax
# The chi^2 quality flag of estimateTOA (toa.html section 4): the distribution of the reduced chi^2 for a clean profile
# (~260 on-pulse bins, 257 degrees of freedom), and where the measured cases sit, against the threshold MaxRedChi2 = 2.
dof = 257
x = np.logspace(np.log10(0.3), np.log10(3000), 2000)
pdf = chi2.pdf(x*dof, dof)*dof
o = []
A = Ax(60, 70, 780, 160, (0.3, 3000), (0, 1.15))
lx = lambda v: A.x0 + (np.log10(v) - np.log10(0.3))/(np.log10(3000) - np.log10(0.3))*A.w
o.append(f'<path class="axis" d="M{A.x0} {A.y0 + A.h} H{A.x0 + A.w}"/>')
for v in (0.3, 0.5, 1, 2, 5, 10, 100, 1000, 3000):
    o.append(f'<path class="axis" d="M{lx(v):.1f} {A.y0 + A.h} v5"/><text class="sub" x="{lx(v):.1f}" y="{A.y0 + A.h + 18}" text-anchor="middle">{v:g}</text>')
o.append(f'<text class="sub" x="{A.x0 + A.w/2:.1f}" y="{A.y0 + A.h + 34}" text-anchor="middle">reduced χ² of the fit (on-pulse bins, log scale)</text>')
pk = pdf.max()
pts = ' L'.join(f'{lx(a):.1f} {float(A.Y(0.45*b/pk)):.1f}' for a, b in zip(x, pdf) if b/pk > 1e-3)
o.append(f'<path class="fillsoft" d="M{lx(x[pdf/pk > 1e-3][0]):.1f} {A.y0 + A.h} L{pts} L{lx(x[pdf/pk > 1e-3][-1]):.1f} {A.y0 + A.h} Z"/>')
o.append(f'<path class="curve" d="M{pts}"/>')
o.append(f'<text class="sub" x="{A.x0 + 4}" y="{A.Y(0.3):.1f}">clean profile:</text>')
o.append(f'<text class="sub" x="{A.x0 + 4}" y="{A.Y(0.3) + 14:.1f}">1 ± 0.09</text>')
o.append(f'<rect class="fillwarn" x="{lx(2):.1f}" y="{A.y0}" width="{A.x0 + A.w - lx(2):.1f}" height="{A.h}"/>')
o.append(f'<path class="curve4" d="M{lx(2):.1f} {A.y0 - 20} V{A.y0 + A.h}"/>')
o.append(f'<text class="lbl" x="{lx(2) + 6:.1f}" y="{A.y0 - 8}">threshold 2: flagged above</text>')
marks = [((1.05, 1.1), 'high-SNR excess', 0), ((1.04, 1.04), 'lopsided LTE channels', 1), ((1.3, 1.3), 'with a carrier', 2),
         ((939, 1625), 'radar in the profile (full band, 5 Oct)', 0)]
for (a, b), lab, row in marks:
    y = A.y0 + 16 + 24*row
    if a == b:
        o.append(f'<path class="dot2" d="M{lx(a) - 4:.1f} {y} a4 4 0 1 0 8 0 a4 4 0 1 0 -8 0"/>')
    else:
        o.append(f'<path class="curve2" stroke-width="5" d="M{lx(a):.1f} {y} H{lx(b):.1f}"/>')
    anchor = 'end'
    xt = lx(a) - 8
    o.append(f'<text class="sub" x="{xt:.1f}" y="{y + 4}" text-anchor="{anchor}">{lab}</text>')
svg = '''<figure class="scroll">
<svg class="wide" viewBox="0 0 860 280" role="img" aria-label="Reduced chi-squared of the TOA fit: the clean distribution, the measured cases and the threshold 2">
  ''' + '\n  '.join(o) + '''
</svg>
<figcaption>Why the χ² flag sits at 2. Blue: the distribution of the reduced χ² for a clean profile (~260 on-pulse bins). Orange: values seen in
the validation runs. The threshold is far above clean profiles and the small known excesses, so it does not throw away good TOAs, and
far below a radar spike in the profile, which it must catch. A carrier (~1.3) is deliberately not flagged: those TOAs are still valid.</figcaption>
</figure>'''
open('chi2.svg', 'w').write(svg)
