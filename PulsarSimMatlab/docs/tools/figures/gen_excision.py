import numpy as np
from plotsvg import Ax
rng = np.random.default_rng(7)
fs = 800e6/192; dt = 1/fs*1e6        # µs per channel sample
n = 600; i = np.arange(n)
x = (rng.standard_normal(n) + 1j*rng.standard_normal(n))/np.sqrt(2)
p0 = np.abs(x)**2
b0, bl = 300, 8                        # a 2 µs radar pulse: 8 samples at ~+15 dB in this channel
s = np.zeros(n, complex); s[b0:b0+bl] = np.sqrt(30)*np.exp(1j*rng.uniform(0, 2*np.pi, bl))
v = x + s; p = np.abs(v)**2
eta = {1: 13.82, 2: 8.39, 4: 5.39, 8: 3.68, 16: 2.69}
flag = np.zeros(n, bool)
for w, e in eta.items():
    st = max(1, w//2)
    for a in range(0, n - w + 1, st):
        if p[a:a+w].mean() > e: flag[a:a+w] = True
g = 10
mask = flag.copy()
for k in np.flatnonzero(flag): mask[max(0, k-g):k+g+1] = True
fl = np.flatnonzero(flag); mk = np.flatnonzero(mask)
print('flagged', fl.min(), fl.max(), len(fl), 'mask', mk.min(), mk.max(), len(mk))
pb = np.where(mask, 0, p)
# dedispersion filter of a channel at 1.4 GHz (delay to 1.6 GHz left out)
K, DM = 4.148808e3, 5.0; fc = 1400.0; cw = 3.125
N = 4096; fbb = np.fft.fftfreq(N, d=1/fs)/1e6     # MHz
W = (np.abs(fbb) <= cw/2).astype(float); ed = 0.02*cw
lo = (fbb >= -cw/2) & (fbb < -cw/2 + ed); hi = (fbb > cw/2 - ed) & (fbb <= cw/2)
W[lo] = np.sin(np.pi/2*(fbb[lo] + cw/2)/ed)**2; W[hi] = np.sin(np.pi/2*(cw/2 - fbb[hi])/ed)**2
f = fc + fbb; ftop = fc + cw/2
H = W*np.exp(-1j*2*np.pi*K*DM*(f - ftop)**2/(f*ftop**2)*1e6)
h2 = np.abs(np.fft.ifft(H))**2; h2 /= h2.sum()
lag = np.arange(N); lag = np.where(lag > N/2, lag - N, lag)
c = np.round((lag*h2).sum()).astype(int)
sup = lag[h2 > 0.5*h2.max()]; print('filter: centroid lag', c, 'support (half max)', sup.min(), sup.max(), sup.max()-sup.min())
h2c = np.roll(h2, -c)                  # centred on lag 0
keep = (~mask).astype(float)
Wt = np.real(np.fft.ifft(np.fft.fft(np.r_[keep, np.ones(N - n)]) * np.fft.fft(h2c)))[:n]
print('W min', Wt.min())
o = []
t = i*dt; T = t[-1]
def panel(y, h, ylim, title, ylab, xt=False):
    A = Ax(80, y, 740, h, (0, T), ylim)
    o.append(f'<text class="ttl" x="80" y="{y-10}">{title}</text>')
    return A
A = panel(34, 80, (0, 40), 'a · detectRFI: power |x|² of one raw channel (baseline m = 1)', '')
o.append(A.frame(yticks=[0, 13.82, 40], yfmt=lambda v: {0: '0', 13.82: '13.8', 40: '40'}[v]))
o.append(f'<path class="thin" d="{A.path(t, p)}"/>')
o.append(f'<path class="dash" stroke="var(--todo)" d="M{A.x0} {A.Y(13.82):.1f} H{A.x0+A.w}"/>')
o.append(f'<text class="sub" x="{A.x0+6}" y="{A.Y(13.82)-5:.1f}">η₁ = 13.82·m (one sample, PFA 10⁻⁶)</text>')
o.append(f'<text class="sub" x="{float(A.X(t[b0+bl]))+8:.1f}" y="{A.Y(34):.1f}">radar pulse: 2 µs = 8 samples</text>')
B = panel(160, 30, (0, 1), 'b · mask: flagged windows (dark) widened by the guard of 10 samples', '')
o.append(f'<path class="axis" d="M{B.x0} {B.y0+B.h} H{B.x0+B.w}"/>')
def spans(m):
    idx = np.flatnonzero(m); out = []
    if len(idx) == 0: return out
    st = idx[0]; pr = idx[0]
    for k in idx[1:]:
        if k != pr + 1: out.append((st, pr)); st = k
        pr = k
    out.append((st, pr)); return out
for a, b in spans(mask):
    o.append(f'<rect class="fillwarn" x="{float(B.X(t[a]-dt/2)):.1f}" y="{B.y0}" width="{float(B.X(t[b]+dt/2)-B.X(t[a]-dt/2)):.1f}" height="{B.h}"/>')
for a, b in spans(flag):
    o.append(f'<rect class="fillfg" x="{float(B.X(t[a]-dt/2)):.1f}" y="{B.y0}" width="{float(B.X(t[b]+dt/2)-B.X(t[a]-dt/2)):.1f}" height="{B.h}" style="opacity:0.7"/>')
xa = float(B.X(t[mk.max()]+dt/2)) + 8
o.append(f'<text class="sub" x="{xa:.1f}" y="{B.y0+20}">{len(mk)} samples blanked ({len(mk)*dt:.1f} µs)</text>')
C = panel(230, 60, (0, 40), 'c · blankChannels: masked voltages set to exactly 0 (no RFI, no noise, no pulsar left)', '')
o.append(C.frame(yticks=[0, 40]))
o.append(f'<path class="thin" d="{C.path(t, pb)}"/>')
D = panel(336, 80, (0.6, 1.03), 'd · after dedispersion: each output mixes ~180 input samples, so the blank becomes a wide dip', '')
o.append(D.frame(yticks=[0.6, 0.8, 1.0], xticks=[0, 20, 40, 60, 80, 100, 120, 140], xlab='time [µs] (channel samples at 4.17 MHz; the channel\'s delay to 1.6 GHz is left out)', yfmt='{:.1f}'))
o.append(f'<path class="curve3" d="{D.path(t, Wt)}"/>')
o.append(f'<text class="sub" x="{float(D.X(t[300])):.1f}" y="{D.Y(Wt.min())+18:.1f}" text-anchor="middle">W = {Wt.min():.2f}: valid fraction (blankingWeights)</text>')
svg = '''<figure class="scroll">
<svg class="wide" viewBox="0 0 860 470" role="img" aria-label="Excision step by step: a radar pulse in one channel is flagged, masked, blanked, and becomes a dip in the valid fraction after dedispersion">
  ''' + '\n  '.join(o) + '''
</svg>
<figcaption>One channel through the excision chain, simulated with the code's own numbers. The input is white noise
with a 2 µs pulse at about 30× the baseline. The detector uses its exact thresholds η₁…η₁₆ = 13.82 … 2.69 for windows of 1–16
samples and a guard of 10 samples. The filter is the dedispersion filter of a 3.125 MHz channel at 1.4 GHz. Panel d shows
<var>W</var>: the fraction of each output sample's filter window that was valid. <code>blankingWeights</code> computes <var>W</var>
exactly for every detected bin, together with two numbers for its noise (<var>V</var> and <var>X</var>(<var>L</var>), section 4).
The fold weights each bin by them.</figcaption>
</figure>'''
open('excision.svg', 'w').write(svg)
