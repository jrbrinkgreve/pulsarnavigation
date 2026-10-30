import numpy as np
from math import erfc, sqrt
from plotsvg import Ax
rng = np.random.default_rng(5)
N = 2048; ph = np.arange(N)/N; ph = np.where(ph >= 0.5, ph - 1, ph)
sg = 0.05/(2*np.sqrt(2*np.log(2)))
s = np.exp(-ph**2/(2*sg**2)); c = s - s.mean(); nc = np.linalg.norm(c)
Cf = np.conj(np.fft.fft(c))
k = np.arange(1, N//2); Sk = np.abs(np.fft.fft(s)[k])**2
rice = np.sqrt(np.sum((2*np.pi*k)**2*Sk)/np.sum(Sk))/(2*np.pi)
print('Rice factor', rice)
snr = 5.0; b = snr/nc
def tmax(nd, pulse):
    out0 = []; outm = []
    for _ in range(nd//1000):
        n = rng.standard_normal((1000, N)) + (b*s if pulse else 0)
        T = np.real(np.fft.ifft(np.fft.fft(n, axis=1)*Cf, axis=1))/nc
        out0.append(T[:, 0]); outm.append(T.max(axis=1))
    return np.concatenate(out0), np.concatenate(outm)
h0_0, h0_m = tmax(40000, False); h1_0, h1_m = tmax(40000, True)
eta0, etam = 3.09, 4.15
print('H0: mean T0 %.3f std %.3f; Tmax median %.2f; P(Tmax>4.15) %.4f; P(T0>3.09) %.4f' % (h0_0.mean(), h0_0.std(), np.median(h0_m), (h0_m > etam).mean(), (h0_0 > eta0).mean()))
pd0 = (h1_0 > eta0).mean(); pdm = (h1_m > etam).mean()
print('H1 SNR 5: PD T0 %.3f (theory %.3f), PD Tmax %.3f, median Tmax %.2f' % (pd0, 1 - 0.5*erfc((eta0 - snr)/sqrt(2))*0 - (1-0.5*erfc((eta0-snr)/sqrt(2)))+1 if False else 0.5*erfc((eta0-snr)/sqrt(2)), pdm, np.median(h1_m)))
def dens(x, lo, hi, nb=64):
    e = np.linspace(lo, hi, nb + 1); h, _ = np.histogram(x, e, density=True); return (e[:-1] + e[1:])/2, h
xs = np.linspace(-4, 9, 1000)
phi = lambda x, m: np.exp(-(x - m)**2/2)/np.sqrt(2*np.pi)
o = []
A = Ax(60, 40, 340, 150, (-4, 9), (0, 0.5))
o.append('<text class="ttl" x="60" y="22">known phase (tracking): T0</text>')
o.append(A.frame(xticks=[-4, -2, 0, 2, 4, 6, 8], xlab='T0 [noise σ]', yticks=[]))
m1 = xs >= eta0
o.append(f'<path class="fillwhy" d="{A.area(xs[m1], phi(xs[m1], snr))}"/>')
o.append(f'<path class="curve" d="{A.path(xs, phi(xs, 0))}"/>')
o.append(f'<path class="curve3" d="{A.path(xs, phi(xs, snr))}"/>')
o.append(f'<path class="curve4" d="M{float(A.X(eta0)):.1f} {A.y0+A.h} V{A.y0+8}"/>')
o.append(f'<text class="sub" x="{float(A.X(eta0))+5:.1f}" y="{A.y0+12}">η₀ = 3.09</text>')
o.append(f'<text class="sub" x="{float(A.X(-1.1)):.1f}" y="{float(A.Y(0.32)):.1f}" text-anchor="end">H0: N(0, 1)</text>')
o.append(f'<text class="sub" x="{float(A.X(6.1)):.1f}" y="{float(A.Y(0.33)):.1f}">H1: SNR 5</text>')
o.append(f'<text class="lbl" x="{float(A.X(5)):.1f}" y="{float(A.Y(0.1)):.1f}" text-anchor="middle">P_D = {0.5*erfc((eta0-snr)/sqrt(2)):.2f}</text>')
x0, d0 = dens(h0_m, -1, 7); x1, d1 = dens(h1_m, 0, 9)
B = Ax(500, 40, 340, 150, (-4, 9), (0, 1.08*max(d0.max(), d1.max())))
o.append('<text class="ttl" x="500" y="22">unknown phase (acquisition): Tmax</text>')
o.append(B.frame(xticks=[-4, -2, 0, 2, 4, 6, 8], xlab='Tmax [noise σ]', yticks=[]))
m1 = x1 >= etam
o.append(f'<path class="fillwhy" d="{B.area(x1[m1], d1[m1])}"/>')
o.append(f'<path class="curve" d="{B.path(x0, d0)}"/>')
o.append(f'<path class="curve3" d="{B.path(x1, d1)}"/>')
o.append(f'<path class="curve4" d="M{float(B.X(etam)):.1f} {B.y0+B.h} V{B.y0+8}"/>')
o.append(f'<text class="sub" x="{float(B.X(etam))+4:.1f}" y="{B.y0+16}">η = 4.15 (Rice)</text>')
o.append(f'<text class="sub" x="{float(B.X(np.median(h0_m)))-24:.1f}" y="{float(B.Y(d0.max()*0.7)):.1f}" text-anchor="end">H0: the maximum</text>')
o.append(f'<text class="sub" x="{float(B.X(np.median(h0_m)))-24:.1f}" y="{float(B.Y(d0.max()*0.7))+15:.1f}" text-anchor="end">of pure noise</text>')
o.append(f'<text class="sub" x="{float(B.X(6.6)):.1f}" y="{float(B.Y(d1.max()))-4:.1f}">H1: SNR 5</text>')
o.append(f'<text class="lbl" x="{float(B.X(5.4)):.1f}" y="{float(B.Y(0.14)):.1f}" text-anchor="middle">P_D = {pdm:.2f}</text>')
svg = '''<figure class="scroll">
<svg class="wide" viewBox="0 0 860 236" role="img" aria-label="Distributions of the detection statistics under H0 and H1, with the thresholds for a false-alarm probability of 1e-3">
  ''' + '\n  '.join(o) + '''
</svg>
<figcaption><b>Figure D.2.</b> The two tests at <var>P</var><sub>FA</sub> = 10⁻³, for a 5 % Gaussian template in 2048 phase bins and a pulse with SNR 5
(left: exact normal densities; right: 40,000 simulated profiles per curve, white noise). The maximum over all shifts of pure noise is already
~2σ, so η rises from 3.09 to 4.15 and the detection probability (green area) drops from ''' + f'{0.5*erfc((eta0-snr)/sqrt(2)):.2f} to {pdm:.2f}' + '''. Simulated false-alarm rate of the
4.15 threshold: ''' + f'{(h0_m > etam).mean()*1e3:.1f}' + '''·10⁻³.</figcaption>
</figure>'''
open('detection.svg', 'w').write(svg)
