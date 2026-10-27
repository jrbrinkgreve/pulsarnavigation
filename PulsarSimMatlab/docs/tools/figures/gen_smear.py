import numpy as np
from plotsvg import Ax, pow10
# Why dedisperse the voltage (coherently): the smear left inside one channel if each channel were only
# shifted by its delay (incoherent dedispersion), at the bottom of the band (worst case), against the
# channel width. tau(f) = K*DM/f^2; smear = tau(fLow) - tau(fLow + width).
K, fL = 4.148808e3, 1200.0                             # s MHz^2, MHz
w = np.logspace(-2, np.log10(400), 300)                # channel width [MHz]
smear = lambda DM: (K*DM/fL**2 - K*DM/(fL + w)**2)*1e6 # us
smearC = lambda DM: (K*DM/fL**2 - K*DM/(fL + 3.125)**2)*1e6   # us, exactly at the pipeline's 3.125 MHz (27 Oct 2026)
o = []
A = Ax(90, 40, 520, 250, (0.01, 400), (0.1, 1e5), logy=True)
lx = lambda v: A.x0 + (np.log10(v) - np.log10(0.01))/(np.log10(400) - np.log10(0.01))*A.w
o.append(A.frame(yticks=[0.1, 1, 10, 100, 1e3, 1e4, 1e5], ylab='smear inside the channel [µs]',
                 yfmt=lambda v: f'{v:g}' if v < 1e4 else pow10(v)))
for v, lab in [(0.01, '0.01'), (0.1, '0.1'), (1, '1'), (10, '10'), (100, '100')]:
    o.append(f'<path class="axis" d="M{lx(v):.1f} {A.y0+A.h} v5"/><text class="sub" x="{lx(v):.1f}" y="{A.y0+A.h+18}" text-anchor="middle">{lab}</text>')
o.append(f'<text class="sub" x="{A.x0 + A.w/2:.1f}" y="{A.y0+A.h+34}" text-anchor="middle">channel width at 1.2 GHz [MHz]</text>')
for DM, cls, lab in ((25, 'curve', 'DM 25 (the simulation)'), (100, 'curve2', 'DM 100 (the design target)')):
    s = smear(DM)
    pts = ' L'.join(f'{lx(a):.1f} {float(A.Y(b)):.1f}' for a, b in zip(w, s) if 0.1 <= b <= 1e5)
    o.append(f'<path class="{cls}" d="M{pts}"/>')
    print(DM, 'smear at 3.125 MHz', smearC(DM), 'full band', smear(DM)[-1])
ly = float(A.Y(5e4))
o.append(f'<path class="curve2" d="M{A.x0+10} {ly-4:.1f} h18"/><text class="sub" x="{A.x0+34}" y="{ly:.1f}">DM 100 (the design target)</text>')
o.append(f'<path class="curve" d="M{A.x0+10} {ly+12:.1f} h18"/><text class="sub" x="{A.x0+34}" y="{ly+16:.1f}">DM 25 (the simulation)</text>')
for val, lab in ((500, 'simulated pulse, FWHM 0.5 ms'), (141, 'J0437, W50 0.141 ms')):
    y = float(A.Y(val))
    o.append(f'<path class="dash" d="M{A.x0} {y:.1f} H{A.x0 + A.w}"/>')
    o.append(f'<text class="sub" x="{A.x0 + 6}" y="{y - 5:.1f}">{lab}</text>')
xc = lx(3.125)
o.append(f'<path class="thin" d="M{xc:.1f} {A.y0} V{A.y0+A.h}"/>')
o.append(f'<text class="sub" x="{xc + 4:.1f}" y="{A.y0 + 12}">the pipeline\'s channels: 3.125 MHz</text>')
for DM, cls in ((25, 'dot'), (100, 'dot2')):
    v = smearC(DM)
    o.append(f'<path class="{cls}" d="M{xc-4:.1f} {A.Y(v):.1f} a4 4 0 1 0 8 0 a4 4 0 1 0 -8 0"/>')
nx = 632
notes = [('ttl', 'If channels were only shifted'), ('sub', 'Each channel would keep its own'),
         ('sub', f'sweep: {smearC(25):.0f} µs in a 3.125 MHz'), ('sub', f'channel at DM 25, {smearC(100)/1e3:.1f} ms at DM 100,'),
         ('sub', 'more than a millisecond pulsar\'s'), ('sub', 'whole pulse.'), ('sub', ''),
         ('ttl', 'Dedispersing the voltage'), ('sub', 'removes the sweep at every'),
         ('sub', 'frequency: no smear at all, in any'), ('sub', 'channel width (stages 5 and 5b).')]
for i, (c, line) in enumerate(notes):
    o.append(f'<text class="{c}" x="{nx}" y="{58 + 17*i}">{line}</text>')
svg = '''<figure class="scroll">
<svg class="wide" viewBox="0 0 860 330" role="img" aria-label="Smear left inside one channel by incoherent dedispersion, against channel width, for DM 25 and DM 100">
  ''' + '\n  '.join(o) + '''
</svg>
<figcaption>What incoherent dedispersion would leave behind: shifting a channel by one delay leaves the sweep <em>inside</em> it, here at the
bottom of the band, where it is largest (τ(1.2 GHz) − τ(1.2 GHz + width), the code's constant). The dots are the pipeline's 3.125 MHz
channels. Filtering the voltage with the inverse chirp, as <code>applyInverseDispersion</code> does, leaves no smear at any width.</figcaption>
</figure>'''
open('smear.svg', 'w').write(svg)
