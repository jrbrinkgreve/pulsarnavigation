import numpy as np
from plotsvg import Ax, pow10
db   = np.array([-25, -20, -15, -10, -5, 0, 10, 20])
rms  = np.array([2430, 995, 26.3, 8.60, 3.07, 1.30, 0.561, 0.489])
rob  = np.array([1470, 83.9, 27.6, 8.53, 2.80, 1.25, 0.537, 0.482])
pred = np.array([217, 82.2, 25.6, 8.28, 2.91, 1.22, 0.538, 0.472])
opt  = np.array([249, 79.0, 25.3, 8.28, 2.89, 1.17, 0.358, 0.198])
def circ(ax, x, y, cls, r=4.5, fill=True):
    X, Y = float(ax.X(x)), float(ax.Y(y))
    if fill: return f'<circle cx="{X:.1f}" cy="{Y:.1f}" r="{r}" class="{cls}"/>'
    return f'<circle cx="{X:.1f}" cy="{Y:.1f}" r="{r}" fill="none" stroke="var(--warn)" stroke-width="1.6"/>'
o = []
A = Ax(70, 50, 450, 220, (-27, 22), (0.1, 5000), logy=True)
o.append(A.frame(xticks=[-25, -20, -15, -10, -5, 0, 10, 20], yticks=[0.1, 1, 10, 100, 1000], yfmt=lambda v: f'{v:g}', xlab='snrDB (pulse peak / receiver noise)', ylab='TOA error per sub-int [µs]'))
o.append(f'<rect class="fillwarn" x="{float(A.X(-27)):.1f}" y="{A.y0}" width="{float(A.X(-17.5)-A.X(-27)):.1f}" height="{A.h}"/>')
o.append(f'<text class="sub" x="{float(A.X(-26.5)):.1f}" y="{float(A.Y(0.17)):.1f}">outliers</text>')
o.append(f'<path class="curve" d="{A.path(db, pred)}"/>')
o.append(f'<path class="curve3" stroke-dasharray="6 4" d="{A.path(db, opt)}"/>')
for x, y in zip(db, rob): o.append(circ(A, x, y, '', fill=False, r=6))
for x, y in zip(db, rms): o.append(circ(A, x, y, 'dot2'))
lx, ly = 345, 72
o.append(f'<path class="curve" d="M{lx} {ly} h24"/><text class="sub" x="{lx+30}" y="{ly+4}">predicted error bar</text>')
o.append(f'<path class="curve3" stroke-dasharray="6 4" d="M{lx} {ly+18} h24"/><text class="sub" x="{lx+30}" y="{ly+22}">optimum (radiometer)</text>')
o.append(f'<circle cx="{lx+12}" cy="{ly+36}" r="4.5" class="dot2"/><text class="sub" x="{lx+30}" y="{ly+40}">measured rms</text>')
o.append(f'<circle cx="{lx+12}" cy="{ly+54}" r="6" fill="none" stroke="var(--warn)" stroke-width="1.6"/><text class="sub" x="{lx+30}" y="{ly+58}">robust (core only)</text>')
o.append(f'<text class="sub" x="{float(A.X(19)):.1f}" y="{float(A.Y(0.13)):.1f}" text-anchor="end">self-noise: FFTFIT 2.4× the optimum</text>')
o.append('<text class="ttl" x="70" y="16">a · phase D sweep: 10 realizations × 9 sub-ints per point</text>')
# panel b: threshold sweep
snr = np.array([3.8, 4.8, 6.0, 7.5, 9.4, 11.8])
yld = np.array([40, 72, 92, 100, 100, 100]); outl = np.array([5.0, 1.7, 0, 0, 0, 0])
B = Ax(600, 50, 240, 220, (3, 12.5), (0, 105))
o.append(B.frame(xticks=[4, 6, 8, 10, 12], yticks=[0, 25, 50, 75, 100], xlab='SNR per sub-int', ylab='%'))
o.append(f'<rect class="fillwhy" x="{float(B.X(6)):.1f}" y="{B.y0}" width="{float(B.X(7)-B.X(6)):.1f}" height="{B.h}"/>')
o.append(f'<text class="sub" x="{float(B.X(6.5)):.1f}" y="{float(B.Y(30)):.1f}" text-anchor="middle">6–7</text>')
o.append(f'<path class="curve3" d="{B.path(snr, yld)}"/>')
for x, y in zip(snr, yld): o.append(circ(B, x, y, 'ok'))
o.append(f'<path class="curve4" d="{B.path(snr, outl)}"/>')
for x, y in zip(snr, outl): o.append(circ(B, x, y, 'bad', r=3.5))
o.append(f'<text class="sub" x="{float(B.X(8.2)):.1f}" y="{float(B.Y(88)):.1f}">good TOAs</text>')
o.append(f'<text class="sub" x="{float(B.X(5.2)):.1f}" y="{float(B.Y(9)):.1f}">outliers (all TOAs)</text>')
o.append('<text class="ttl" x="600" y="16">b · threshold sweep (6 Oct)</text>')
svg = '''<figure class="scroll">
<svg class="wide" viewBox="0 0 860 320" role="img" aria-label="TOA error versus SNR, and the yield of good TOAs versus SNR per sub-integration">
  ''' + '\n  '.join(o) + '''
</svg>
<figcaption>a: the numbers of the table above (full band, 1 turn per sub-int). From −15 dB up the measured scatter follows the predicted
error bar; below, outliers (fits on noise peaks) blow up the rms while the core stays on the prediction. At high SNR FFTFIT stops
improving (self-noise) while the optimum keeps falling. b: the threshold sweep of section 9: the detector and χ² remove every outlier,
and above SNR 6–7 per sub-int nearly every TOA is usable.</figcaption>
</figure>'''
open('results.svg', 'w').write(svg)
