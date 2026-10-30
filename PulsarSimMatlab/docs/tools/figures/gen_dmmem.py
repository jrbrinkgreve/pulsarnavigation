import numpy as np
from plotsvg import Ax, pow10
# Memory of the dispersion filters against DM, with the code's sizing rules (smallest allowed FFT):
#  stage 2 (simulator, real 4 GHz): Nfft = 2^nextpow2(2*Nk - 2), 44 bytes per FFT sample (applyDispersionStream)
#  stage 5 (full band, complex 800 MHz): Nfft = 2^nextpow2(2*(nPast + nFuture)), 56 bytes (applyInverseDispersion)
#  stage 5b (one channel, complex 4.17 MHz, bottom channel, reference 1.6 GHz): same rule, 56 bytes
K = 4.148808e3*1e12                                  # s Hz^2 per DM unit
fL, fH = 1.2e9, 1.6e9
sweep = lambda DM: K*DM*(1/fL**2 - 1/fH**2)          # s
p2 = lambda x: 2.0**np.ceil(np.log2(x))
def fwd(DM):
    g = np.ceil(20/8e6*4e9); Nk = g + np.ceil(sweep(DM)*4e9) + g + 1
    return p2(2*Nk - 2)*44
def full(DM):
    fs = 800e6; g = np.ceil(20/8e6*fs)
    return p2(2*(np.ceil(sweep(DM)*fs) + 2*g))*56
def chan(DM):
    fs = 800e6/192; g = np.ceil(20/62.5e3*fs)
    return p2(2*(np.ceil(K*DM*(1/fL**2 - 1/fH**2)*fs) + 2*g))*56
DMs = np.logspace(0, np.log10(200), 400)
for DM in (5, 25, 30, 60, 100):
    print(DM, 'fwd %.1f GB' % (fwd(DM)/1e9), 'full %.2f GB' % (full(DM)/1e9), 'chan %.1f MB' % (chan(DM)/1e6))
o = []
A = Ax(90, 40, 520, 250, (1, 200), (1e6, 3e11), logy=True)
lx = lambda v: A.x0 + np.log10(v)/np.log10(200)*A.w
o.append(A.frame(yticks=[1e6, 1e7, 1e8, 1e9, 1e10, 1e11], ylab='memory of the dispersion filter [bytes]', yfmt=lambda v: {1e6: '1 MB', 1e7: '10 MB', 1e8: '100 MB', 1e9: '1 GB', 1e10: '10 GB', 1e11: '100 GB'}[v]))
for v in (1, 2, 5, 10, 20, 50, 100, 200):
    o.append(f'<path class="axis" d="M{lx(v):.1f} {A.y0+A.h} v5"/><text class="sub" x="{lx(v):.1f}" y="{A.y0+A.h+18}" text-anchor="middle">{v}</text>')
o.append(f'<text class="sub" x="{A.x0 + A.w/2:.1f}" y="{A.y0+A.h+34}" text-anchor="middle">DM [pc cm⁻³]</text>')
for y, lab in ((24e9, '24 GB: the machine'), (4e9, '4 GB: default budget of stage 5')):
    o.append(f'<path class="dash" d="M{A.x0} {A.Y(y):.1f} H{A.x0 + A.w}"/><text class="sub" x="{A.x0 + 6}" y="{A.Y(y) - 5:.1f}">{lab}</text>')
for fn, cls, lab in ((fwd, 'curve4', 'stage 2: simulator, real 4 GHz'), (full, 'curve2', 'stage 5: full band, 800 MHz'), (chan, 'curve', 'stage 5b: one channel, 4.17 MHz')):
    v = np.array([fn(d) for d in DMs])
    pts = ' L'.join(f'{lx(a):.1f} {float(A.Y(b)):.1f}' for a, b in zip(DMs, v))
    o.append(f'<path class="{cls}" d="M{pts}"/>')
for x in (25, 100):
    o.append(f'<path class="thin" d="M{lx(x):.1f} {A.y0} V{A.y0 + A.h}"/>')
o.append(f'<text class="sub" x="{lx(25) + 4:.1f}" y="{A.y0 + 12}">DM 25: now</text>')
o.append(f'<text class="sub" x="{lx(100) - 4:.1f}" y="{A.y0 + 12}" text-anchor="end">DM 100: target</text>')
ly = float(A.Y(5e7))
for i, (cls, lab) in enumerate((('curve4', 'stage 2: simulator, real 4 GHz'), ('curve2', 'stage 5: full band, 800 MHz'), ('curve', 'stage 5b: one channel, 4.17 MHz'))):
    o.append(f'<path class="{cls}" d="M{A.x0 + 10} {ly + 16*i - 4:.1f} h18"/><text class="sub" x="{A.x0 + 34}" y="{ly + 16*i:.1f}">{lab}</text>')
nx = 632
notes = [('ttl', 'At DM 25 (now)'), ('sub', 'simulator 11.8 GB, full band 3.8 GB'), ('sub', ''), ('ttl', 'At DM 100'), ('sub', 'simulator ~47 GB (open)'), ('sub', 'full band ~15 GB'), ('sub', 'one channel ~120 MB'), ('sub', ''),
         ('sub', 'A channel holds 1/128 of the'), ('sub', 'band at 1/192 of the rate, so its'),
         ('sub', 'filter spans its sweep in ~200×'), ('sub', 'fewer samples. Channels run one'), ('sub', 'after another.')]
for i, (c, line) in enumerate(notes):
    o.append(f'<text class="{c}" x="{nx}" y="{58 + 17*i}">{line}</text>')
svg = '''<figure class="scroll">
<svg class="wide" viewBox="0 0 860 330" role="img" aria-label="Memory of the dispersion filter against DM: simulator, full-band dedispersion and one channel">
  ''' + '\n  '.join(o) + '''
</svg>
<figcaption><b>Figure 5b.1.</b> Memory of one dispersion filter against DM, with each function's own sizing rule (the smallest FFT it may use: twice the reach of
the filter, rounded up to a power of two; 44 or 56 bytes per FFT sample). The steps are the powers of two. One full-band filter at
800 MHz passes the 4 GB default near DM 33 and needs ~15 GB at DM 100; a channel's filter stays near 100 MB. The simulator's
forward filter (stage 2, real samples at 4 GHz) is the remaining problem (<a href="roadmap.html#dm">roadmap</a>).</figcaption>
</figure>'''
open('dmmem.svg', 'w').write(svg)
