import numpy as np
from plotsvg import Ax
# How simulateFold draws a fold (simulator.html #how), with its formulas and illustrative numbers: 128 phase bins (the code uses
# 2048; fewer bins make the smoothing and the correlation visible), the generator's unit-peak periodic Gaussian (FWHM 5 %), five
# time bins per phase bin, rho = 0 dB (the pulse peak doubles the baseline), noise 8 % of the baseline per phase bin.
#   mean  m_j = 1 + rho p~(phi_j - dphi),  p~_k = p_k sinc(k f0 dt) sinc^2(k/N)   (time-bin average, linear assignment)
#   noise sd m_j s0, neighbour correlation 1/4 (evenly spread time bins: weightX / weight2 = (n/6) / (2n/3)),
#         drawn as z = IFFT(sqrt(lambda) FFT(e)), lambda_k = 1 + 2 (1/4) cos(2 pi k / N), e white.
N = 128; fwhm = 0.05; rho = 1.0; dphi = 0.35; s0 = 0.08; f0dt = 1 / (5 * N)
k = np.fft.fftfreq(N, 1 / N); sg = fwhm / (2 * np.sqrt(2 * np.log(2)))
pk = np.exp(-2 * np.pi**2 * sg**2 * k**2); pk = pk / np.sum(np.exp(-2 * np.pi**2 * sg**2 * np.arange(-N, N)**2))   # peak 1
pt = pk * np.sinc(k * f0dt) * np.sinc(k / N)**2
phi = np.arange(N) / N
m = 1 + rho * N * np.real(np.fft.ifft(pt * np.exp(-2j * np.pi * k * dphi)))
lam = 1 + 0.5 * np.cos(2 * np.pi * k / N)
rng = np.random.default_rng(3)
def draw():
    e = rng.standard_normal(N); z = np.real(np.fft.ifft(np.sqrt(lam) * np.fft.fft(e))); return e, z
o = []
# a: mean and its noise band
A = Ax(70, 50, 230, 220, (0, 1), (0.6, 2.4))
for v in (1.0, 1.5, 2.0):
    o.append(f'<path class="grid" d="M{A.x0} {float(A.Y(v)):.1f} H{A.x0 + A.w}"/>')
o.append(A.frame(xticks=[0, 0.25, 0.5, 0.75, 1], yticks=[1, 1.5, 2], xlab='pulse phase [turns]', ylab='power / receiver noise'))
ph = np.r_[phi, 1.0]; mm = np.r_[m, m[0]]
up = [f'{float(A.X(a)):.1f} {float(A.Y(v + s0 * v)):.1f}' for a, v in zip(ph, mm)]
dn = [f'{float(A.X(a)):.1f} {float(A.Y(v - s0 * v)):.1f}' for a, v in zip(ph[::-1], mm[::-1])]
o.append(f'<path class="s1band" d="M{" L".join(up + dn)} Z"/>')
o.append(f'<path class="s1" d="{A.path(ph, mm)}"/>')
o.append(f'<path class="dash" d="M{float(A.X(dphi)):.1f} {A.y0 + A.h} V{A.y0 + 10}"/>')
o.append(f'<text class="sub" x="{float(A.X(dphi)) + 5:.1f}" y="{A.y0 + 20}">true phase Δφ</text>')
o.append(f'<text class="sub" x="{float(A.X(0.62)):.1f}" y="{float(A.Y(1.0 + 2.2 * s0)):.1f}" style="fill:var(--series-1)">mean ± 1σ</text>')
o.append('<text class="ttl" x="70" y="22">a · what the model knows</text>')
# b: three drawn sub-ints
B = Ax(370, 50, 230, 220, (0, 1), (0.4, 6.6))
draws = []
for i in range(3):
    e, z = draw(); draws.append((e, z)); prof = m + s0 * m * z; off = 2.0 * i
    o.append(f'<path class="grid" d="M{B.x0} {float(B.Y(1 + off)):.1f} H{B.x0 + B.w}"/>')
    o.append(f'<path class="s1" style="stroke-width:1.3" d="{B.path(ph, np.r_[prof, prof[0]] + off)}"/>')
    o.append(f'<text class="sub" x="{B.x0 + B.w + 6}" y="{float(B.Y(1 + off)) + 4:.1f}">sub-int {i + 1}</text>')
o.append(B.frame(xticks=[0, 0.25, 0.5, 0.75, 1], xlab='pulse phase [turns]', left=False))
o.append('<text class="ttl" x="370" y="22">b · three drawn sub-ints</text>')
# c: the correlation of the drawn noise against lag: the model, and 4000 draws
zz = np.array([draw()[1] for _ in range(4000)])
lags = np.arange(-3, 4); meas = np.array([np.mean(zz * np.roll(zz, -L, 1)) for L in lags]) / np.mean(zz**2)
r1 = meas[lags == 1][0]
model = np.where(lags == 0, 1.0, np.where(np.abs(lags) == 1, 0.25, 0.0))
Cc = Ax(680, 50, 160, 220, (-3.6, 3.6), (-0.1, 1.1))
for v in (0.25, 0.5, 1.0):
    o.append(f'<path class="grid" d="M{Cc.x0} {float(Cc.Y(v)):.1f} H{Cc.x0 + Cc.w}"/>')
o.append(Cc.frame(xticks=[-3, -2, -1, 0, 1, 2, 3], yticks=[0, 0.25, 0.5, 1], xlab='lag [phase bins]', ylab='noise correlation'))
for L, mv, me in zip(lags, model, meas):
    x = float(Cc.X(L))
    o.append(f'<rect class="fillsoft" x="{x - 9:.1f}" y="{float(Cc.Y(max(mv, 0))):.1f}" width="18" height="{abs(float(Cc.Y(0)) - float(Cc.Y(mv))):.1f}"><title>model: {mv:g}</title></rect>')
    o.append(f'<circle class="s1f" cx="{x:.1f}" cy="{float(Cc.Y(me)):.1f}" r="4"><title>4000 draws: {me:.3f}</title></circle>')
o.append(f'<text class="sub" x="{float(Cc.X(1)) + 10:.1f}" y="{float(Cc.Y(0.25)) - 6:.1f}">¼</text>')
o.append('<text class="ttl" x="680" y="22">c · neighbours share noise</text>')

svg = '''<figure class="scroll">
<svg class="wide" viewBox="0 0 860 340" role="img" aria-label="How simulateFold draws a fold: the mean profile with its one-sigma band, which widens on the pulse; three drawn sub-ints; and the noise correlation against lag: a quarter between neighbouring bins, none beyond">
  ''' + '\n  '.join(o) + '''
</svg>
<figcaption><b>Figure F.1.</b> How <code>simulateFold</code> draws one sub-int, with its own formulas and illustrative numbers (128 phase bins, a pulse that
doubles the baseline, noise of 8 % of the baseline). (a) The mean, the smoothed pulse at its true phase on the baseline, with a noise band
that is wider on the pulse, because the pulsar adds its own noise. (b) Three sub-ints drawn from that model. (c) The noise correlation of
neighbouring bins: ¼, as in a real fold, and none beyond; bars are the model, dots 4000 draws (''' + f'{r1:.3f}' + ''' at lag 1). Generated by
<code>docs/tools/figures/gen_simdraw.py</code>.</figcaption>
</figure>'''
open('simdraw.svg', 'w').write(svg)
print('neighbour correlation', r1, 'mean peak', m.max(), 'at', phi[np.argmax(m)])
