import numpy as np
from plotsvg import Ax
rng = np.random.default_rng(3)
T = 10.0; dtb = 0.04                     # ms; drawing resolution 40 µs
t = np.arange(0, 40, dtb) + dtb/2
sig = 0.5/(2*np.sqrt(2*np.log(2)))
phi = (t - T/2)/T                        # ephemeris phase: TRef = T/2, f0 = 1/T
turn = np.round(phi); frac = phi - turn
G = np.exp(-(frac*T)**2/(2*sig**2))
P = 1 + 1.0*G + 0.16*rng.standard_normal(t.size)   # baseline 1, noisy detected power
o = []
A = Ax(80, 46, 740, 90, (0, 40), (0.4, 2.4))
o.append('<text class="ttl" x="80" y="18">a · detected power vs time: each bin gets its phase φ(t) = f0·(t − TRef); dashed: the cuts at φ = k + ½</text>')
o.append(A.frame(xticks=[0, 5, 10, 15, 20, 25, 30, 35, 40], xlab='time [ms]', yticks=[]))
for k in range(4):
    xb = float(A.X(T*(k + 1)))           # phase k + 1/2 at t = (k+1)·T
    if k < 3:
        o.append(f'<path class="dash" d="M{xb:.1f} {A.y0} V{A.y0+A.h}"/>')
    o.append(f'<text class="sub" x="{float(A.X(T*k + T/2)):.1f}" y="{A.y0-4}" text-anchor="middle">turn {k}</text>')
o.append(f'<path class="thin" d="{A.path(t, P)}"/>')
# panel b: the turns stacked
rows = []
o.append('<text class="ttl" x="80" y="196">b · sub-integrations (here 1 turn each)</text>')
for k in range(4):
    B = Ax(80, 206 + k*30, 330, 26, (-0.5, 0.5), (0.4, 2.4))
    m = turn == k
    o.append(f'<path class="thin" d="{B.path(frac[m], P[m])}"/>')
    o.append(f'<text class="sub" x="{B.x0-8}" y="{B.y0+17}" text-anchor="end">{k}</text>')
    rows.append(P[m])
Bx = Ax(80, 206, 330, 3*30+26, (-0.5, 0.5), (0, 1))
o.append(f'<path class="axis" d="M80 {206+3*30+26} H410"/>')
for v in (-0.5, 0, 0.5):
    o.append(f'<text class="sub" x="{float(Bx.X(v)):.1f}" y="{206+3*30+26+16}" text-anchor="middle">{v:g}</text>')
o.append(f'<text class="sub" x="245" y="{206+3*30+26+32}" text-anchor="middle">phase (frac = φ − round(φ)) [turns]</text>')
# panel c: the profile
prof = np.mean(rows, axis=0); fr = frac[turn == 0]
C = Ax(520, 206, 300, 116, (-0.5, 0.5), (0.6, 2.2))
o.append('<text class="ttl" x="520" y="196">c · profile = sum / weight per phase bin</text>')
o.append(C.frame(xticks=[-0.5, 0, 0.5], yticks=[], xlab='phase [turns]'))
o.append(f'<path class="curve" d="{C.path(fr, prof)}"/>')
o.append(f'<path class="curve2" stroke-width="1.6" stroke-dasharray="5 3" d="{C.path(fr, 1 + 1.0*np.exp(-(fr*T)**2/(2*sig**2)))}"/>')
o.append(f'<text class="sub" x="{float(C.X(0.08)):.1f}" y="{C.y0+22}">noise ÷ √(turns): ÷ 2 here</text>')
o.append(f'<text class="sub" x="{float(C.X(0.08)):.1f}" y="{C.y0+38}">dashed: the true pulse</text>')
o.append('<path class="line" marker-end="url(#ar-fo)" d="M430 266 H508"/>')
o.append('<text class="sub" x="469" y="258" text-anchor="middle">average</text>')
svg = '''<figure class="scroll">
<svg class="wide" viewBox="0 0 860 380" role="img" aria-label="Folding: the time series is cut into turns by the ephemeris phase, the turns are stacked by phase and averaged into a profile">
  <defs>
    <marker id="ar-fo" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse">
      <path class="ah" d="M0 0 L10 5 L0 10 z"/>
    </marker>
  </defs>
  ''' + '\n  '.join(o) + '''
</svg>
<figcaption>Folding in three steps (simulated, coarse 40 µs bins, a bright pulse for visibility). Every time bin gets its phase from the
ephemeris; the integer part is its turn, the rest its place in the profile. Turns are cut halfway between pulses, so each pulse stays
whole. In the stored fold, bin <var>j</var> is centred at phase (<var>j</var> − 1)/<code>NBin</code> and negative phases wrap to the end, so the
pulse sits at both ends of the array; it is drawn here from −½ to ½.</figcaption>
</figure>'''
open('fold.svg', 'w').write(svg)
