import numpy as np
from plotsvg import Ax
rng = np.random.default_rng(11)
N = 2048; ph = np.arange(N)/N; ph = np.where(ph >= 0.5, ph - 1, ph)    # bin j at (j-1)/N, wrapped
sg = 0.05/(2*np.sqrt(2*np.log(2)))
def tmpl(x):
    x = (x + 0.5) % 1 - 0.5; return np.exp(-x**2/(2*sg**2))
tau0, a0, b0, sn = 0.004, 1.0, 0.5, 0.12
p = a0 + b0*tmpl(ph - tau0) + sn*rng.standard_normal(N)
s = tmpl(ph)
P = np.fft.fft(p); S = np.fft.fft(s); k = np.arange(1, N//2)
X = P[k]*np.conj(S[k])
def C(t): return np.real(np.exp(2j*np.pi*np.outer(t, k)) @ X)
tg = np.linspace(-0.05, 0.05, 2001); Cg = C(tg)
th = tg[np.argmax(Cg)]
for _ in range(30):          # Newton on C'
    e = np.exp(2j*np.pi*k*th); c1 = np.real((2j*np.pi*k*e) @ X); c2 = np.real((-(2*np.pi*k)**2*e) @ X); th -= c1/c2
bh = C(np.array([th]))[0]/np.sum(np.abs(S[k])**2)
ah = (np.real(P[0]) - bh*np.real(S[0]))/N
# error: tau = where C'(tau)=0; C' = sum_j d_j p_j with d_j = d/dtau of the correlation weights
d = np.real(np.fft.ifft(np.r_[0, (2j*np.pi*k*np.exp(2j*np.pi*k*th))*np.conj(S[k]), np.zeros(N - N//2)]))*N   # C' = sum d_j p_j
d = np.real(np.fft.fft(np.r_[0, 0*k, 0*np.arange(N - N//2)]))  # placeholder, replaced below
# direct: C'(tau) = Re sum_k 2πik e^{2πikτ} P_k conj(S_k), P_k = sum_j p_j e^{-2πikj/N}
jj = np.arange(N)
D = np.exp(-2j*np.pi*np.outer(k, jj)/N)            # P_k = D @ p
w = (2j*np.pi*k*np.exp(2j*np.pi*k*th)*np.conj(S[k]))
d = np.real(w @ D)
c2 = np.real((-(2*np.pi*k)**2*np.exp(2j*np.pi*k*th)) @ X)
sig_tau = sn*np.sqrt(np.sum(d**2))/abs(c2)
T = 10e3   # µs
print(f'tau0 {tau0*T:.2f} us, tauhat {th*T:.2f} us, sigma {sig_tau*T:.2f} us, a {ah:.3f}, b {bh:.3f}')
o = []
A = Ax(60, 40, 230, 150, (-0.1, 0.1), (0.5, 1.9))
o.append('<text class="ttl" x="60" y="20">a · profile and fitted model</text>')
o.append(A.frame(xticks=[-0.1, 0, 0.1], xlab='phase [turns]', yticks=[]))
m = np.abs(ph) <= 0.1; ix = np.argsort(ph[m])
o.append(f'<path class="thin" d="{A.path(ph[m][ix], p[m][ix])}"/>')
o.append(f'<path class="curve2" d="{A.path(ph[m][ix], ah + bh*tmpl(ph[m][ix] - th))}"/>')
o.append(f'<path class="dash" d="M{float(A.X(0)):.1f} {A.y0} V{A.y0+A.h}"/>')
o.append(f'<text class="sub" x="{float(A.X(0))-4:.1f}" y="{A.y0+12}" text-anchor="end">predicted</text>')
o.append(f'<text class="sub" x="{float(A.X(0.02)):.1f}" y="{A.y0+12}">a + b·s(φ − τ̂)</text>')
B = Ax(340, 40, 230, 150, (-0.05, 0.05), (min(Cg.min(), 0) - 0.55*Cg.max(), Cg.max()*1.15))
o.append('<text class="ttl" x="340" y="20">b · cross-correlation C(τ)</text>')
o.append(B.frame(xticks=[-0.05, 0, 0.05], xlab='shift τ [turns]', yticks=[]))
o.append(f'<path class="curve" d="{B.path(tg, Cg)}"/>')
xh = float(B.X(th)); yh = float(B.Y(C(np.array([th]))[0]))
o.append(f'<path class="dot2" d="M{xh-4:.1f} {yh:.1f} a4 4 0 1 0 8 0 a4 4 0 1 0 -8 0"/>')
o.append(f'<path class="thin" d="M{xh:.1f} {yh+5:.1f} V{B.y0+B.h}"/>')
o.append(f'<text class="sub" x="{xh+8:.1f}" y="{yh+4:.1f}">peak: C′(τ̂) = 0</text>')
o.append(f'<text class="sub" x="{float(B.X(0)):.1f}" y="{B.y0+B.h-24}" text-anchor="middle">τ̂ = {th*T:.1f} µs (true {tau0*T:.0f} µs)</text>')
o.append(f'<text class="sub" x="{float(B.X(0)):.1f}" y="{B.y0+B.h-8}" text-anchor="middle">σ_τ = {sig_tau*T:.1f} µs from its sharpness</text>')
Cc = Ax(620, 40, 230, 150, (-0.1, 0.1), (-1.1, 1.1))
o.append('<text class="ttl" x="620" y="20">c · weights d in C′ = Σ d_j·p_j</text>')
o.append(Cc.frame(xticks=[-0.1, 0, 0.1], xlab='phase [turns]', yticks=[], left=False))
o.append(f'<path class="axis" d="M{Cc.x0} {float(Cc.Y(0)):.1f} H{Cc.x0+Cc.w}"/>')
dn = d/np.max(np.abs(d))
o.append(f'<path class="curve3" d="{Cc.path(ph[m][ix], dn[m][ix])}"/>')
o.append(f'<text class="sub" x="{Cc.x0+4}" y="{Cc.y0+14}">flanks: all the timing</text>')
o.append(f'<text class="sub" x="{Cc.x0+4}" y="{Cc.y0+28}">peak, baseline: none</text>')
svg = '''<figure class="scroll">
<svg class="wide" viewBox="0 0 860 236" role="img" aria-label="FFTFIT: the profile, the cross-correlation with the template and its peak, and the weights that set the TOA error">
  ''' + '\n  '.join(o) + '''
</svg>
<figcaption><b>Figure 8.1.</b> FFTFIT on a simulated profile (2048 bins, 5 % Gaussian, true shift 40 µs, white noise). (a) The data and the fitted
<var>a</var> + <var>b</var>·<var>s</var>(φ − τ̂). (b) The cross-correlation with the template, exact between bins from the Fourier phases;
τ̂ is its peak. (c) The weights <var>d</var> of Eq. (8.4): the template slope, so only the pulse flanks set the error.</figcaption>
</figure>'''
open('toa.svg', 'w').write(svg)
