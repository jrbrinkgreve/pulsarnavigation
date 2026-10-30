import numpy as np
from plotsvg import Ax
# powerCovariance for the pipeline's channels (3.125 MHz, 62.5 kHz sin^2 edges, 4.1667 MHz, 4 samples per bin),
# recomputed with its formula: a: how long the detected power stays correlated (|R(l)|^2 against lag);
# b: the variance and neighbour covariances of a 0.96 us bin, against the independent-bin radiometer model.
def pc(chanW, edgeW, fs, n, Nf=2**18, maxLag=40):
    f = (np.arange(Nf) - Nf/2)*fs/Nf; inn = np.abs(f) <= chanW/2
    W = np.zeros(Nf); W[inn] = 1
    lo = inn & (f < -chanW/2 + edgeW); W[lo] = np.sin(np.pi/2*(f[lo] + chanW/2)/edgeW)**2
    hi = inn & (f > chanW/2 - edgeW); W[hi] = np.sin(np.pi/2*(chanW/2 - f[hi])/edgeW)**2
    R = np.fft.fft(np.fft.ifftshift(W**2)); R2 = np.abs(R/R[0])**2
    rad = R2.sum()/n; d = np.arange(-(n - 1), n); cnt = n - np.abs(d)
    C = np.array([np.sum(cnt*R2[(d + L*n) % Nf])/n**2 for L in range(maxLag + 1)])
    return C[0]/rad, C[1:]/rad, R2
fs = 800e6/192
V, X, R2 = pc(3.125e6, 62.5e3, fs, 4)
print('V', V, 'X', X[:7], 'sum V+2X(1..7)', V + 2*X[:7].sum())
o = []
A = Ax(80, 46, 330, 200, (0, 4.0), (0, 1.05))
o.append(A.frame(xticks=[0, 1, 2, 3, 4], yticks=[0, 0.5, 1], xlab='time between two samples [µs]', ylab='correlation of their power, |R|²'))
for k in range(5):                                   # the 0.96 us bins
    if k % 2 == 0:
        xa = float(A.X(k*4/fs*1e6)); xb = min(float(A.X((k + 1)*4/fs*1e6)), A.x0 + A.w)
        o.append(f'<rect class="fillfg" x="{xa:.1f}" y="{A.y0}" width="{xb - xa:.1f}" height="{A.h}"/>')
lags = np.arange(0, 17)
for l in lags:
    x = float(A.X(l/fs*1e6)); y = float(A.Y(R2[l]))
    o.append(f'<path class="line" d="M{x:.1f} {A.Y(0):.1f} V{y:.1f}"/>')
    o.append(f'<path class="dot" d="M{x-3:.1f} {y:.1f} a3 3 0 1 0 6 0 a3 3 0 1 0 -6 0"/>')
o.append(f'<text class="sub" x="{A.X(1.1):.1f}" y="{A.Y(0.85):.1f}">a channel of 3.125 MHz forgets in ~1/B ≈ 0.3 µs:</text>')
o.append(f'<text class="sub" x="{A.X(1.1):.1f}" y="{A.Y(0.85) + 14:.1f}">about one sample (0.24 µs) apart.</text>')
o.append(f'<text class="sub" x="{A.X(1.1):.1f}" y="{A.Y(0.85) + 28:.1f}">Grey/white bands: the 0.96 µs bins (4 samples).</text>')
o.append(f'<text class="ttl" x="{A.x0}" y="20">a · how long the detected power remembers</text>')
B = Ax(520, 46, 310, 200, (-0.6, 7.6), (0, 1.05))
o.append(B.frame(yticks=[0, 0.5, 1], ylab='relative to m²/(B·Δt)'))
labels = ['V', 'X₁', 'X₂', 'X₃', 'X₄', 'X₅', 'X₆', 'X₇']
vals = [V] + list(X[:7])
for i, (lab, v) in enumerate(zip(labels, vals)):
    x0 = float(B.X(i - 0.3)); x1 = float(B.X(i + 0.3))
    o.append(f'<rect class="box sky" x="{x0:.1f}" y="{B.Y(v):.1f}" width="{x1 - x0:.1f}" height="{B.Y(0) - B.Y(v):.1f}"/>')
    o.append(f'<text class="sub" x="{(x0 + x1)/2:.1f}" y="{B.y0 + B.h + 16}" text-anchor="middle">{lab}</text>')
    if i < 3:
        o.append(f'<text class="sub" x="{(x0 + x1)/2:.1f}" y="{B.Y(v) - 5:.1f}" text-anchor="middle">{v:.4f}</text>')
o.append(f'<path class="curve2" stroke-dasharray="6 4" d="M{B.X(-0.45):.1f} {B.Y(1):.1f} H{B.X(0.45):.1f}"/>')
o.append(f'<text class="sub" x="{B.X(0.6):.1f}" y="{B.Y(1) + 4:.1f}">independent bins: V = 1, every X = 0</text>')
o.append(f'<text class="sub" x="{B.X(1.6):.1f}" y="{B.Y(0.55):.1f}">V + 2·ΣX = 1: the total is right,</text>')
o.append(f'<text class="sub" x="{B.X(1.6):.1f}" y="{B.Y(0.55) + 14:.1f}">but split over neighbouring bins</text>')
o.append(f'<text class="sub" x="{B.x0 + B.w/2:.1f}" y="{B.y0 + B.h + 34}" text-anchor="middle">one 0.96 µs bin: its variance and covariances</text>')
o.append(f'<text class="ttl" x="{B.x0}" y="20">b · so neighbouring bins share noise</text>')
svg = '''<figure class="scroll">
<svg class="wide" viewBox="0 0 860 300" role="img" aria-label="Noise of detected bins in a narrow channel: the power correlation against lag, and the variance and covariances of one bin">
  ''' + '\n  '.join(o) + '''
</svg>
<figcaption><b>Figure 6b.1.</b> <code>powerCovariance</code> for the pipeline's channels, recomputed with its formula. a: the squared correlation of the channel
voltage, which is how strongly the detected power at two times is related (Isserlis). It falls within about one sample but has a slow
tail from the sharp channel edges. b: so a 0.96 µs bin has less variance than the radiometer formula (0.888) and shares some noise
with its neighbours. Summed over all lags the total is exactly the radiometer value; the fold carries <var>V</var> and <var>X</var> so the
error bars are exact phase bin by phase bin.</figcaption>
</figure>'''
open('powcov.svg', 'w').write(svg)
