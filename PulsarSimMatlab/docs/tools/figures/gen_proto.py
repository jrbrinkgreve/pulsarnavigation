import numpy as np
from plotsvg import Ax
from pfb import design, response, db, envelope
# The prototype low-pass of channelizeIQ in time (a tapered sinc) and in frequency (a block with sloped edges).
h, n, w, I = design(192)
G, fc, fsOut, dF = I['G'], I['fc'], I['fsOut'], 3.125
stop = fsOut - dF/2
h2 = h**2
lobe = h2[np.abs(n) <= 192].sum()/h2.sum()
print(I, 'main lobe', lobe)

o = []
# (a) time
A = Ax(80, 58, 750, 170, (-2000, 2000), (-0.3, 1.1))
o.append(A.frame(xticks=[-1920, -960, 0, 960, 1920], yticks=[0, 0.5, 1],
                 xlab='tap r − G [samples at 800 MHz; 192 samples = 0.24 µs]', ylab='weight h (peak = 1)'))
o.append(f'<rect class="fillwhy" x="{A.X(-192):.1f}" y="{A.y0}" width="{A.X(192) - A.X(-192):.1f}" height="{A.h}"/>')
o.append(f'<path class="thin" d="M{A.X(-2000):.1f} {A.Y(0):.1f} H{A.X(2000):.1f}"/>')
o.append(f'<path class="dash" d="{A.path(n[::8], w[::8])}"/>')
o.append(f'<path class="curve" d="{A.path(n[::3], h[::3]/h.max())}"/>')
for j in range(1, 11):
    for s in (-1, 1):
        o.append(f'<circle class="dot2" cx="{A.X(s*192*j):.1f}" cy="{A.Y(0):.1f}" r="2.6"/>')
o.append(f'<text class="ttl" x="80" y="22">a · In time: a sinc, tapered by a Kaiser window</text>')
o.append(f'<text class="sub" x="{A.X(230):.1f}" y="{A.Y(0.93):.1f}">central ±192 samples: {100*lobe:.0f} % of Σh²</text>')
o.append(f'<text class="sub" x="{A.X(560):.1f}" y="{A.Y(0.55):.1f}">dashed: the Kaiser taper</text>')
o.append(f'<text class="sub" x="{A.X(560):.1f}" y="{A.Y(-0.2):.1f}">dots: zero crossings, every 192 samples</text>')
o.append(f'<text class="sub" x="{A.X(-1980):.1f}" y="{A.Y(0.93):.1f}">3853 taps = 4.8 µs</text>')

# (b) frequency
B = Ax(80, 320, 750, 190, (-7, 7), (-100, 6))
f = np.arange(-7, 7.0001, 0.002)
H = response(h, n, f)
fe, He = envelope(f, db(H), 0.04)
o.append(B.frame(xticks=[-6, -4, -2, 0, 2, 4, 6], yticks=[0, -20, -40, -60, -80, -100],
                 xlab='frequency offset from the channel centre [MHz]', ylab='gain [dB]'))
o.append(f'<rect class="fillwhy" x="{B.X(-dF/2):.1f}" y="{B.y0}" width="{B.X(dF/2) - B.X(-dF/2):.1f}" height="{B.h}"/>')
for a, b in ((-stop, -dF/2), (dF/2, stop)):
    o.append(f'<rect class="fillwarn" x="{B.X(a):.1f}" y="{B.y0}" width="{B.X(b) - B.X(a):.1f}" height="{B.h}"/>')
o.append(f'<path class="thin" d="M{B.X(-7):.1f} {B.Y(-80):.1f} H{B.X(7):.1f}" stroke-dasharray="2 3"/>')
ideal = f'M{B.X(-7):.1f} {B.y0 + B.h:.1f} H{B.X(-fc):.1f} V{B.Y(0):.1f} H{B.X(fc):.1f} V{B.y0 + B.h:.1f} H{B.X(7):.1f}'
o.append(f'<path class="dash" d="{ideal}"/>')
o.append(f'<path class="curve" d="{B.path(fe, He)}"/>')
for s in (-1, 1):
    o.append(f'<circle class="dot2" cx="{B.X(s*fc):.1f}" cy="{B.Y(-6.02):.1f}" r="3.5"/>')
o.append(f'<text class="ttl" x="80" y="292">b · In frequency: a block with sloped edges</text>')
o.append(f'<text class="sub" x="{B.X(0):.1f}" y="{B.Y(-30):.1f}" text-anchor="middle">useful band</text>')
o.append(f'<text class="sub" x="{B.X(0):.1f}" y="{B.Y(-30) + 14:.1f}" text-anchor="middle">±1.5625 MHz,</text>')
o.append(f'<text class="sub" x="{B.X(0):.1f}" y="{B.Y(-30) + 28:.1f}" text-anchor="middle">flat to 10⁻⁴</text>')
o.append(f'<text class="sub" x="{B.X(fc) + 8:.1f}" y="{B.Y(-6) + 4:.1f}">−6 dB at the block edge, ±2.083 MHz</text>')
o.append(f'<text class="sub" x="{B.X(stop) + 6:.1f}" y="{B.Y(-50):.1f}">slope: 1.04 MHz wide</text>')
o.append(f'<text class="sub" x="{B.X(3.2):.1f}" y="{B.Y(-80) - 6:.1f}">−80 dB or lower from ±2.604 MHz</text>')
o.append(f'<text class="sub" x="{B.X(-6.9):.1f}" y="{B.Y(-6) + 4:.1f}">dashed: the ideal block</text>')

svg = '''<figure class="scroll">
<svg class="wide" viewBox="0 0 860 560" role="img" aria-label="The prototype filter: a Kaiser-tapered sinc in time with zero crossings every 192 samples, and a block with sloped edges in frequency">
  ''' + '\n  '.join(o) + '''
</svg>
<figcaption><b>Figure 5a.1.</b> The prototype of <code>channelizeIQ</code>, from its design formula. (a) In time: a sinc with zero
crossings every 192 samples, tapered by a Kaiser window (dashed). (b) In frequency: the ideal block (dashed) gets slopes 1.04 MHz wide
(orange), with −6 dB at the block edge. The useful band (green) is flat.</figcaption>
</figure>'''
open('proto.svg', 'w').write(svg)
