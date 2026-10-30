import numpy as np
from plotsvg import Ax
from pfb import design, response, K, dF
# One prototype, shifted to every channel: the 256 FFT slots across the IQ band (which rows are kept),
# and the real responses of six channels around the LO.
fLO, fLow = 1300.0, 1200.0
h, n, w, I = design(192)
beta = 0.5                                    # (fLow - fLO)/dF + 0.5 = -31.5 -> kFirst -32, beta 1/2

o = []
# (a) the 256 slots
A = Ax(80, 64, 750, 40, (900, 1700), (0, 1))
o.append(f'<text class="ttl" x="80" y="22">a · The 256 slots of the inverse FFT across the IQ band (±400 MHz around the LO)</text>')
segs = [(-128, -33, 'fillfg', 'rows 129–224: not kept'), (-32, -1, 'fillwhy', 'rows 225–256 → ch 1–32'),
        (0, 95, 'fillwhy', 'rows 1–96 → channels 33–128'), (96, 127, 'fillfg', 'rows 97–128: not kept')]
for k0, k1, cls, lab in segs:
    x0, x1 = fLO + k0*dF, fLO + (k1 + 1)*dF
    o.append(f'<rect class="{cls}" x="{A.X(x0):.1f}" y="{A.y0}" width="{A.X(x1) - A.X(x0):.1f}" height="{A.h}"/>')
ticks = ' '.join(f'M{A.X(fLO + k*dF):.1f} {A.y0 + A.h - 9} v9' for k in range(-128, 129))
o.append(f'<path class="thin" d="{ticks}" stroke-width="0.6"/>')
o.append(f'<path class="axis" d="M{A.X(900):.1f} {A.y0 + A.h} H{A.X(1700):.1f}"/>')
o.append(f'<path class="dash" d="M{A.X(fLO):.1f} {A.y0 - 14} V{A.y0 + A.h + 6}"/>')
o.append(f'<text class="sub" x="{A.X(fLO):.1f}" y="{A.y0 - 18}" text-anchor="middle">LO 1.3 GHz</text>')
for k0, k1, cls, lab in segs:
    xm = fLO + (k0 + k1 + 1)/2*dF
    o.append(f'<text class="sub" x="{A.X(xm):.1f}" y="{A.y0 + A.h/2 + 4:.1f}" text-anchor="middle">{lab}</text>')
for t in range(900, 1701, 100):
    o.append(f'<path class="axis" d="M{A.X(t):.1f} {A.y0 + A.h} v5"/>')
    o.append(f'<text class="sub" x="{A.X(t):.1f}" y="{A.y0 + A.h + 18}" text-anchor="middle">{t}</text>')
o.append(f'<text class="sub" x="{A.X(1300):.1f}" y="{A.y0 + A.h + 34}" text-anchor="middle">RF frequency [MHz]; ticks inside the strip: one per slot (3.125 MHz)</text>')

# (b) six channels around the LO
B = Ax(80, 250, 750, 200, (1288.5, 1311.5), (0, 1.15))
f = np.arange(1288.5, 1311.5001, 0.04)
o.append(f'<text class="ttl" x="80" y="212">b · Six channels around the LO: the same block, shifted by 3.125 MHz each time</text>')
o.append(B.frame(xticks=[1290, 1295, 1300, 1305, 1310], yticks=[0, 0.5, 1],
                 xlab='RF frequency [MHz]', ylab='gain (linear)'))
for i, kj in enumerate(range(-3, 3)):
    fk = fLO + (kj + beta)*dF
    j = kj + 33
    row = kj % K + 1
    cls = 's1' if i % 2 == 0 else 's2'
    o.append(f'<rect class="{cls}band" x="{B.X(fk - dF/2):.1f}" y="{B.Y(1):.1f}" width="{B.X(fk + dF/2) - B.X(fk - dF/2):.1f}" height="{B.Y(0) - B.Y(1):.1f}"/>')
    fj = f[np.abs(f - fk) <= 3.5]             # the stopband beyond is zero on this scale
    o.append(f'<path class="{cls}" d="{B.path(fj, response(h, n, fj - fk))}"/>')
    o.append(f'<text class="lbl" x="{B.X(fk):.1f}" y="{B.Y(0.62):.1f}" text-anchor="middle">ch {j}</text>')
    o.append(f'<text class="sub" x="{B.X(fk):.1f}" y="{B.Y(0.62) + 15:.1f}" text-anchor="middle">row {row}</text>')
    o.append(f'<text class="sub" x="{B.X(fk):.1f}" y="{B.Y(0.62) + 29:.1f}" text-anchor="middle">{fk:.2f}</text>')
o.append(f'<path class="dash" d="M{B.X(fLO):.1f} {B.y0 - 6} V{B.Y(0):.1f}"/>')
o.append(f'<text class="sub" x="{B.X(fLO):.1f}" y="{B.y0 - 10}" text-anchor="middle">LO: on the edge between channels 32 and 33</text>')
o.append(f'<text class="sub" x="{B.X(1296.9):.1f}" y="{B.Y(1.06):.1f}" text-anchor="middle">slopes overlap</text>')

svg = '''<figure class="scroll">
<svg class="wide" viewBox="0 0 860 500" role="img" aria-label="The 256 FFT slots across the IQ band with the 128 kept channels, and the responses of six channels around the LO as shifted copies of one block">
  ''' + '\n  '.join(o) + '''
</svg>
<figcaption><b>Figure 5a.2.</b> One prototype for every channel. (a) The 256 DFT slots span the 800 MHz IQ band. The 128 slots
inside 1.2–1.6 GHz are kept, and those below the LO come from the top rows. (b) Six computed channel responses: the same block,
3.125 MHz apart, labelled with each channel's DFT row and centre frequency.</figcaption>
</figure>'''
open('pfbbank.svg', 'w').write(svg)
