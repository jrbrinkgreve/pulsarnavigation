import numpy as np
from plotsvg import Ax
# A toy version of generatePulsarSignal: same formula x = A*sqrt(G)*n, but a low sample rate
# (300 kHz instead of 4 GHz) so that single samples can be drawn.
T, duty, A = 10.0, 5.0, 1.0                      # ms, %, signal scale
sig = duty/100*T/(2*np.sqrt(2*np.log(2)))
def G(t):                                         # nearest pulse only (exact at 5 % duty)
    k0 = np.round(t/T - 0.5)
    return np.exp(-(t - (k0 + 0.5)*T)**2/(2*sig**2))
rng = np.random.default_rng(1)
o = []
# panel a: three turns
fa = 50.0                                         # samples per ms for the overview
t = np.arange(0, 30, 1/fa); x = A*np.sqrt(G(t))*rng.standard_normal(t.size)
P = Ax(70, 40, 760, 90, (0, 30), (-3.2, 3.2))
o.append(P.frame(xticks=[0, 5, 10, 15, 20, 25, 30], yticks=[-3, 0, 3], xlab='time [ms]', ylab='voltage x(t)'))
o.append(f'<path class="thin" d="{P.path(t, x)}"/>')
for c in (5, 15, 25):
    o.append(f'<text class="sub" x="{P.X(c):.1f}" y="{P.y0-6:.1f}" text-anchor="middle">pulse centre {c:g} ms</text>')
o.append(f'<text class="ttl" x="{P.x0}" y="20">a · the voltage over three turns: silent between pulses, a burst of noise once per turn</text>')
# panel b: zoom on one pulse, the samples and the envelope
fb = 300.0
t = np.arange(4.0, 6.0, 1/fb); x = A*np.sqrt(G(t))*rng.standard_normal(t.size)
Q = Ax(70, 215, 330, 120, (4, 6), (-3.2, 3.2))
o.append(Q.frame(xticks=[4, 4.5, 5, 5.5, 6], yticks=[-3, 0, 3], xlab='time [ms]', ylab='x(t), one pulse'))
o.append(f'<path class="thin" d="{Q.path(t, x)}"/>')
te = np.linspace(4, 6, 401); env = A*np.sqrt(G(te))
o.append(f'<path class="curve2" stroke-dasharray="5 3" d="{Q.path(te, env)}"/>')
o.append(f'<path class="curve2" stroke-dasharray="5 3" d="{Q.path(te, -env)}"/>')
o.append(f'<text class="sub" x="{Q.X(5.28):.1f}" y="{Q.Y(2.6):.1f}">dashed: ± A·√G(t), the local std</text>')
o.append(f'<text class="ttl" x="{Q.x0}" y="{Q.y0-26}">b · zoom: random samples, no waveform</text>')
# panel c: power, averaged over many turns, follows G
nTurn = 400; bw = 0.04                            # 40 us bins
edges = np.arange(4.0, 6.0 + 1e-9, bw); cen = edges[:-1] + bw/2
acc = np.zeros(cen.size)
ts = np.arange(4.0, 6.0, 1/fb)
idx = np.minimum(((ts - 4.0)/bw).astype(int), cen.size - 1)
for _ in range(nTurn):
    x = A*np.sqrt(G(ts))*rng.standard_normal(ts.size)
    acc += np.bincount(idx, x**2, cen.size)
cnt = np.bincount(idx, minlength=cen.size)*nTurn
pw = acc/cnt
R = Ax(510, 215, 320, 120, (4, 6), (0, 1.25))
o.append(R.frame(xticks=[4, 4.5, 5, 5.5, 6], yticks=[0, 0.5, 1], xlab='time [ms]', ylab='power x²'))
o.append(f'<path class="fillsoft" d="{R.area(te, A**2*G(te))}"/>')
o.append(f'<path class="curve" d="{R.path(te, A**2*G(te))}"/>')
for c, p in zip(cen, pw):
    o.append(f'<path class="dot2" d="M{R.X(c)-2.5:.1f} {R.Y(p):.1f} a2.5 2.5 0 1 0 5 0 a2.5 2.5 0 1 0 -5 0"/>')
o.append(f'<text class="sub" x="{R.X(5.3):.1f}" y="{R.Y(1.0):.1f}">line: A²·G(t)</text>')
o.append(f'<text class="sub" x="{R.X(5.3):.1f}" y="{R.Y(0.86):.1f}">dots: mean x², {nTurn} turns</text>')
o.append(f'<text class="ttl" x="{R.x0}" y="{R.y0-26}">c · the pulse is in the power</text>')
svg = '''<figure class="scroll">
<svg class="wide" viewBox="0 0 860 382" role="img" aria-label="The generated voltage: noise whose variance rises once per turn; its power follows the profile G">
  ''' + '\n  '.join(o) + '''
</svg>
<figcaption><b>Figure 1.1.</b> What <code>generatePulsarSignal</code> writes, by Eq. (1.1) but at a toy sample rate (300 kHz instead of 4 GHz) so that
single samples can be seen. The voltage only gets louder once per turn (b). Squared and averaged, its power follows
<var>A</var>²·<var>G</var>(<var>t</var>) (c). Off-pulse the generated signal is zero; receiver noise is added in stage 3.</figcaption>
</figure>'''
open('voltage.svg', 'w').write(svg)
