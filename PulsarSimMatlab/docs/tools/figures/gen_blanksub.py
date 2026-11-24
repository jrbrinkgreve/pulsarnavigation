import numpy as np
from plotsvg import Ax, pow10
# Why blank instead of subtract (excision.html section 1): the energy of one radar pulse in units of one pulsar
# pulse grows as 1/rho (at -5 dB they are about equal), and any subtraction leaves a fixed fraction of it.
# Marked: J0437-4715 at the design point of 23 Nov 2026 (-44.2 dB): 10^3.92 = 8,300 and 8.3.
dB = np.linspace(-60, 0, 241)
radar = 10**((-5 - dB)/10)
o = []
A = Ax(90, 40, 720, 250, (-60, 0), (1e-5, 1e6), logy=True)
o.append(A.frame(xticks=[-60, -50, -40, -30, -20, -10, -5, 0], yticks=[1e-4, 1e-2, 1, 1e2, 1e4, 1e6], xlab='snrDB',
                 ylab='energy, in units of one pulsar pulse', yfmt=lambda v: f'{v:g}' if 1e-2 < v < 1e4 else pow10(v)))
o.append(f'<path class="fillwhy" d="M{A.x0} {A.Y(1):.1f} H{A.x0 + A.w} V{A.y0 + A.h} H{A.x0} Z"/>')
o.append(f'<path class="dash" d="M{A.x0} {A.Y(1):.1f} H{A.x0 + A.w}"/>')
o.append(f'<text class="sub" x="{A.X(-59):.1f}" y="{A.Y(0.05):.1f}">below one pulsar pulse: harmless</text>')
for i, (frac, cls, lab) in enumerate(((1, 'curve4', 'one radar pulse'), (1e-3, 'curve2', 'left after 99.9 % subtraction'), (1e-4, 'curve', 'left after 99.99 % subtraction'))):
    o.append(f'<path class="{cls}" d="{A.path(dB, radar*frac)}"/>')
    ly = float(A.Y(3e5)) + 16*i
    o.append(f'<path class="{cls}" d="M{A.X(-28):.1f} {ly - 4:.1f} h18"/><text class="sub" x="{A.X(-28) + 24:.1f}" y="{ly:.1f}">{lab}</text>')
xd = float(A.X(-44.2))
o.append(f'<path class="thin" d="M{xd:.1f} {A.y0} V{A.y0 + A.h}"/>')
o.append(f'<text class="sub" x="{xd + 6:.1f}" y="{A.y0 + A.h - 28}">J0437−4715 at the design point (−44.2 dB)</text>')
for v in (10**3.92, 10**0.92):
    o.append(f'<path class="dot" d="M{xd - 3.5:.1f} {A.Y(v):.1f} a3.5 3.5 0 1 0 7 0 a3.5 3.5 0 1 0 -7 0"/>')
o.append(f'<text class="sub" x="{xd + 8:.1f}" y="{A.Y(10**3.92) - 6:.1f}">~8,300</text>')
o.append(f'<text class="sub" x="{xd + 8:.1f}" y="{A.Y(10**0.92) - 6:.1f}">~8</text>')
o.append(f'<path class="curve3" stroke-width="3" d="M{A.x0} {A.y0 + A.h - 2} H{A.x0 + A.w}"/>')
o.append(f'<text class="sub" x="{A.X(-59):.1f}" y="{A.y0 + A.h - 8}">blanked: exactly 0 left, at every SNR (drawn at the bottom)</text>')
svg = '''<figure class="scroll">
<svg class="wide" viewBox="0 0 860 330" role="img" aria-label="Energy of one radar pulse and what subtraction would leave, in units of one pulsar pulse, against SNR">
  ''' + '\n  '.join(o) + '''
</svg>
<figcaption><b>Figure E.2.</b> The energy of one radar pulse in units of one pulsar pulse (top line): equal at −5 dB, then growing as 1/ρ as the pulsar gets
weaker. The lines below show what a subtraction of 99.9 % or 99.99 % would leave. For J0437−4715 at the design point (−44.2 dB, dots) a
99.9 % subtraction still leaves about 8 pulsar pulses of radar energy, and a 99.99 % one about as much as the pulsar itself.</figcaption>
</figure>'''
open('blanksub.svg', 'w').write(svg)
