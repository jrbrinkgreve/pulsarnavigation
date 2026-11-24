import numpy as np
from plotsvg import Ax
# How noise moves the peak of the FFTFIT correlation (toa.html, Section 8.3): the correlation and its slope without
# noise and for one noise draw, and the scatter of the fitted shift over many draws against Eq. (8.11).
# The same profile as gen_toa.py (2048 bins, 5 % Gaussian, white noise), with more noise so that the shift can be seen.
rng = np.random.default_rng(3)
N = 2048; ph = np.arange(N)/N
sg = 0.05/(2*np.sqrt(2*np.log(2)))
def tmpl(x):
    x = (x + 0.5) % 1 - 0.5; return np.exp(-x**2/(2*sg**2))
tau0, a0, b0, sn = 0.004, 1.0, 0.5, 0.3
s = tmpl(ph); S = np.fft.fft(s); k = np.arange(1, N//2); w1 = 2*np.pi*k; Sk = S[k]
model = a0 + b0*tmpl(ph - tau0)
def Cd(p, t, o=0):                       # C(tau) and its derivatives, Eq. (8.5): Re sum_k P_k S_k* e^{2 pi i k tau}
    X = np.fft.fft(p)[k]*np.conj(Sk)
    return np.real((((1j*w1)**o)*np.exp(1j*np.outer(np.atleast_1d(t), w1))) @ X)
def dvec(t):                             # C'(tau) = sum_j d_j p_j, as in estimateTOA's fitOne
    D = np.zeros(N, complex); D[k] = 1j*w1*np.conj(Sk)*np.exp(1j*w1*t); return np.real(np.fft.fft(D))
def fit(p):
    t = tau0
    for _ in range(30):
        t -= Cd(p, t, 1)[0]/Cd(p, t, 2)[0]
    return t
d0 = dvec(tau0); C2 = Cd(model, tau0, 2)[0]
sig = sn*np.sqrt(np.sum(d0**2))/abs(C2)               # Eq. (8.11) for white noise
M = 4000; dt = np.empty(M); keep = None
for i in range(M):
    n = sn*rng.standard_normal(N); th = fit(model + n); dt[i] = th - tau0
    if keep is None and 1.4 < dt[i]/sig < 1.6: keep = (n, th)
n1, th1 = keep; dC = float(np.sum(d0*n1))
T = 10e3                                               # µs per turn (the 10 ms pulsar)
snr = b0*np.sqrt(np.sum(s**2))/sn
print(f'SNR {snr:.1f}; Eq. (8.11) {sig*T:.2f} us; scatter of {M} draws {dt.std()*T:.2f} us; ratio {dt.std()/sig:.3f}')
print(f'draw: shift {(th1-tau0)*T:.2f} us; -dC\'/C\'\' {-dC/C2*T:.2f} us')
o = []
tg = tau0 + np.linspace(-0.02, 0.02, 401); x = (tg - tau0)*T
Cm = Cd(model, tg); Cn = Cd(model + n1, tg); sc = Cm.max()
A = Ax(60, 40, 230, 150, (-200, 200), (min(Cm.min(), Cn.min())/sc - 0.02, max(Cm.max(), Cn.max())/sc + 0.07))
o.append('<text class="ttl" x="60" y="20">a · the peak of C(τ)</text>')
o.append(A.frame(xticks=[-200, 0, 200], xlab='shift τ − τ₀ [µs]', yticks=[]))
o.append(f'<path class="curve" d="{A.path(x, Cm/sc)}"/>')
o.append(f'<path class="curve2" d="{A.path(x, Cn/sc)}"/>')
for tt, CC, cls in ((tau0, Cm, 'curve'), (th1, Cn, 'curve2')):
    xx = float(A.X((tt - tau0)*T)); yy = float(A.Y(np.interp((tt - tau0)*T, x, CC/sc)))
    o.append(f'<path class="{cls}" d="M{xx-4:.1f} {yy:.1f} a4 4 0 1 0 8 0 a4 4 0 1 0 -8 0"/>')
    o.append(f'<path class="thin" d="M{xx:.1f} {yy+5:.1f} V{A.y0+A.h}"/>')
o.append(f'<text class="sub" x="{float(A.X(-165)):.1f}" y="{A.y0+A.h-22}">blue: no noise</text>')
o.append(f'<text class="sub" x="{float(A.X(-165)):.1f}" y="{A.y0+A.h-8}">orange: one draw</text>')
mb = np.abs(x) <= 100                                   # panel b shows the middle half of panel a
c1m = Cd(model, tg, 1)[mb]; c1n = Cd(model + n1, tg, 1)[mb]; xb = x[mb]; s1 = np.max(np.abs(c1m))
B = Ax(340, 40, 230, 150, (-100, 100), (min(c1m.min(), c1n.min())/s1 - 0.25, max(c1m.max(), c1n.max())/s1 + 0.2))
o.append('<text class="ttl" x="340" y="20">b · its slope C′(τ)</text>')
o.append(B.frame(xticks=[-100, 0, 100], xlab='shift τ − τ₀ [µs]', yticks=[], left=False))
y0 = float(B.Y(0)); x0 = float(B.X(0)); xh = float(B.X((th1 - tau0)*T)); yd = float(B.Y(dC/s1))
o.append(f'<path class="axis" d="M{B.x0} {y0:.1f} H{B.x0+B.w}"/>')
o.append(f'<path class="curve" d="{B.path(xb, c1m/s1)}"/>')
o.append(f'<path class="curve2" d="{B.path(xb, c1n/s1)}"/>')
o.append(f'<path class="dash" d="M{x0:.1f} {y0:.1f} V{yd:.1f}"/>')                       # the lift dC' at the true shift
o.append(f'<text class="sub" x="{x0+7:.1f}" y="{yd-7:.1f}">lift δC′</text>')
for xx, cls in ((x0, 'curve'), (xh, 'curve2')):                                        # the two zero crossings
    o.append(f'<path class="{cls}" d="M{xx-4:.1f} {y0:.1f} a4 4 0 1 0 8 0 a4 4 0 1 0 -8 0"/>')
o.append(f'<text class="sub" x="{B.x0+B.w-4}" y="{B.y0+14}" text-anchor="end">gradient C″</text>')
o.append(f'<text class="sub" x="{B.x0+4}" y="{B.y0+B.h-22}">the zero crossing</text>')
o.append(f'<text class="sub" x="{B.x0+4}" y="{B.y0+B.h-8}">moves by −δC′/C″</text>')
Cc = Ax(620, 40, 230, 150, (-80, 80), (0, 1.15))
o.append(f'<text class="ttl" x="620" y="20">c · the shift over {M} draws</text>')
o.append(Cc.frame(xticks=[-80, -40, 0, 40, 80], xlab='fitted shift τ̂ − τ₀ [µs]', yticks=[], left=False))
edges = np.linspace(-80, 80, 33); h, _ = np.histogram(dt*T, edges); wbin = edges[1] - edges[0]
g = lambda u: np.exp(-u**2/(2*(sig*T)**2)); hn = h/(M*wbin)*np.sqrt(2*np.pi)*sig*T
for e0, v in zip(edges[:-1], hn):
    xa, xb = float(Cc.X(e0)), float(Cc.X(e0 + wbin)); yv = float(Cc.Y(v))
    o.append(f'<path class="fillsoft" d="M{xa+0.5:.1f} {Cc.y0+Cc.h} V{yv:.1f} H{xb-0.5:.1f} V{Cc.y0+Cc.h} Z"/>')
ug = np.linspace(-80, 80, 201)
o.append(f'<path class="curve" d="{Cc.path(ug, g(ug))}"/>')
o.append(f'<text class="sub" x="{Cc.x0+4}" y="{Cc.y0+14}">curve: Eq. (8.11), {sig*T:.1f} µs</text>')
o.append(f'<text class="sub" x="{Cc.x0+4}" y="{Cc.y0+28}">bars: scatter {dt.std()*T:.1f} µs</text>')
svg = '''<figure class="scroll">
<svg class="wide" viewBox="0 0 860 236" role="img" aria-label="How noise moves the peak of the correlation: the peak with and without noise, the slope whose zero crossing moves, and the scatter of the fitted shift against Eq. 8.11">
  ''' + '\n  '.join(o) + '''
</svg>
<figcaption><b>Figure 8.2.</b> How noise moves the peak, for the profile of Figure 8.1 with more noise. (a) The correlation near its
peak without noise and for one noise draw: the maximum has moved. (b) The slope of both. Without noise it crosses zero at the true shift.
Noise lifts the line by δ<var>C</var>′, and the crossing moves by −δ<var>C</var>′/<var>C</var>″, Eq. (8.10). (c) The fitted shift over
many noise draws, against a Gaussian of the width that Eq. (8.11) gives.</figcaption>
</figure>'''
open('toaerr.svg', 'w').write(svg)
