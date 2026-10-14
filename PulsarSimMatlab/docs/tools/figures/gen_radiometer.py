import numpy as np
from plotsvg import Ax, pow10
# Radiometer equation (physics.html, section 5). a: power |z|^2 of complex white noise (mean 1), raw and averaged over 16 and
# 256 samples. b: relative scatter of an n-sample average vs n, simulated (2^20 samples) for white noise (band = sample rate)
# and for noise filtered to half the band; lines: exact Isserlis sums (1/n^2) sum |R(a-b)|^2, with R(k) = 1 (k = 0) or 0 for
# white and R(k) = sinc(k/2) for the half band.
rng = np.random.default_rng(11)
N = 2**20
z = (rng.standard_normal(N) + 1j*rng.standard_normal(N))/np.sqrt(2)
Z = np.fft.fft(z)
f = np.fft.fftfreq(N)
zh = np.fft.ifft(np.where(np.abs(f) < 0.25, Z, 0))
zh /= np.sqrt(np.mean(np.abs(zh)**2))
p, ph = np.abs(z)**2, np.abs(zh)**2
o = []
# ---------- a ----------
A = Ax(60, 40, 330, 250, (0, 120), (0, 3))
o.append(f'<text class="ttl" x="60" y="20">a · detected power, averaged over n samples</text>')
rows = [(1, 'n = 1: scatter 100 %'), (16, 'n = 16: 25 %'), (256, 'n = 256: 6 %')]
hh = 250/3
for i, (n, lab) in enumerate(rows):
    S = Ax(60, 40 + i*hh, 330, hh - 22, (0, 120), (0, 4 if n == 1 else 2))
    m = p[:120*n].reshape(120, n).mean(1)
    o.append(f'<path class="axis" d="M{S.x0} {S.y0 + S.h} H{S.x0 + S.w} M{S.x0} {S.y0} V{S.y0 + S.h}"/>')
    o.append(f'<path class="dash" d="M{S.x0} {float(S.Y(1)):.1f} H{S.x0 + S.w}"/>')
    o.append(f'<path class="curve" d="{S.path(np.arange(120) + 0.5, m)}"/>')
    o.append(f'<text class="sub" x="{S.x0 + S.w + 8}" y="{float(S.Y(1)) + 4:.1f}">mean 1</text>')
    o.append(f'<text class="lbl" x="{S.x0 + 6}" y="{S.y0 + 12}">{lab}</text>')
o.append(f'<text class="sub" x="{60 + 165}" y="{40 + 250 + 6}" text-anchor="middle">consecutive averages (time →)</text>')
# ---------- b ----------
B = Ax(520, 40, 310, 230, (0, np.log10(1024)), (1e-2, 1.5), logy=True)
lx = lambda n: B.x0 + np.log10(n)/np.log10(1024)*B.w
o.append(f'<path class="axis" d="M{B.x0} {B.y0 + B.h} H{B.x0 + B.w} M{B.x0} {B.y0} V{B.y0 + B.h}"/>')
for n in (1, 4, 16, 64, 256, 1024):
    o.append(f'<path class="axis" d="M{lx(n):.1f} {B.y0 + B.h} v5"/><text class="sub" x="{lx(n):.1f}" y="{B.y0 + B.h + 18}" text-anchor="middle">{n}</text>')
for v, lab in ((1, '100 %'), (0.1, '10 %'), (0.01, '1 %')):
    o.append(f'<path class="axis" d="M{B.x0} {float(B.Y(v)):.1f} h-5"/><text class="sub" x="{B.x0 - 8}" y="{float(B.Y(v)) + 4:.1f}" text-anchor="end">{lab}</text>')
o.append(f'<text class="sub" x="{B.x0 + B.w/2:.1f}" y="{B.y0 + B.h + 34}" text-anchor="middle">samples per average n</text>')
o.append(f'<text class="sub" x="{B.x0}" y="{B.y0 - 8}">relative scatter σ/mean</text>')
ns = 2**np.arange(0, 11)
nn = np.unique(np.round(np.logspace(0, np.log10(1024), 120)).astype(int))
k = np.arange(-1100, 1101)
R2 = np.sinc(k/2)**2
exact_h = []
for n in nn:
    lags = np.arange(-(n - 1), n)
    w = n - np.abs(lags)
    exact_h.append(np.sqrt(np.sum(w*R2[lags + 1100]))/n)
o.append(f'<path class="curve" d="M{" L".join(f"{lx(n):.1f} {float(B.Y(1/np.sqrt(n))):.1f}" for n in nn)}"/>')
o.append(f'<path class="curve2" d="M{" L".join(f"{lx(n):.1f} {float(B.Y(e)):.1f}" for n, e in zip(nn, exact_h))}"/>')
for n in ns:
    for arr, cls in ((p, 'dot'), (ph, 'dot2')):
        m = arr[:(N//n)*n].reshape(-1, n).mean(1)
        v = m.std()/m.mean()
        x, y = lx(n), float(B.Y(v))
        o.append(f'<path class="{cls}" d="M{x - 3.5:.1f} {y:.1f} a3.5 3.5 0 1 0 7 0 a3.5 3.5 0 1 0 -7 0"/>')
o.append(f'<path class="curve" d="M{B.x0 + 150} {B.y0 + 12} h22"/><text class="sub" x="{B.x0 + 178}" y="{B.y0 + 16}">full band: 1/√n</text>')
o.append(f'<path class="curve2" d="M{B.x0 + 150} {B.y0 + 30} h22"/><text class="sub" x="{B.x0 + 178}" y="{B.y0 + 34}">half band: 1/√(n/2)</text>')
o.append(f'<text class="sub" x="{B.x0 + 150}" y="{B.y0 + 52}">dots: simulated; lines: Isserlis</text>')
o.append(f'<text class="ttl" x="520" y="20">b · the scatter falls as 1/√(B·Δt)</text>')
svg = '''<figure class="scroll">
<svg class="wide" viewBox="0 0 860 340" role="img" aria-label="Radiometer equation: averaged power scatters less with more samples, as one over the square root of bandwidth times bin length">
  ''' + '\n  '.join(o) + '''
</svg>
<figcaption>The radiometer equation, simulated with complex Gaussian noise. a: one power sample scatters as much as its mean
(100 %); an average over <var>n</var> samples scatters 1/√<var>n</var> as much. b: what counts is not the number of samples but how
many are independent: <var>B</var>·Δ<var>t</var>. Noise filtered to half the band has neighbouring samples that are correlated, so
an average of <var>n</var> samples holds only ~<var>n</var>/2 independent ones (orange). The lines are the exact sums of Isserlis'
theorem, the formula behind the whole noise model.</figcaption>
</figure>'''
open('radiometer.svg', 'w').write(svg)
