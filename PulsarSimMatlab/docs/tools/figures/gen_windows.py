import numpy as np
from scipy.stats import gamma
from plotsvg import Ax
# detectRFI's thresholds against window length (PFA 1e-6): exact values for n = 1..16 (excision.html), and a Gamma
# with the exact variance for longer windows (within 0.5 % of exact where both exist). Channel correlation rho(d)
# from the channelizer prototype (Kaiser, 3853 taps, D = 192). Against: the pulsar's peak power at -5 and -20 dB.
M, D, A_db = 3853, 192, 80.0
fc = 1/(2*D); beta = 0.1102*(A_db - 8.7)
n = np.arange(M) - (M - 1)/2
h = np.where(n == 0, 2*fc, np.sin(2*np.pi*fc*n)/(np.pi*np.where(n == 0, 1, n)))
h *= np.i0(beta*np.sqrt(1 - (2*np.arange(M)/(M - 1) - 1)**2))/np.i0(beta)
rho = np.array([np.sum(h[:M - d*D]*h[d*D:])/np.sum(h**2) for d in range(20)])
def eta(nn):
    r = np.array([[rho[abs(a - b)] if abs(a - b) < 20 else 0 for b in range(nn)] for a in range(nn)])
    var = np.sum(np.linalg.eigvalsh(r)**2)/nn**2
    return gamma.isf(1e-6, 1/var, scale=var)
exact = {1: 13.82, 2: 8.39, 4: 5.39, 8: 3.68, 16: 2.69}
longer = {nn: eta(nn) for nn in (32, 64, 128, 256)}
print({k: round(v, 3) for k, v in longer.items()})
o = []
A = Ax(90, 40, 720, 250, (0, 8), (0.9, 16), logy=True)
lx = lambda nn: float(A.X(np.log2(nn)))
o.append(A.frame(yticks=[1, 2, 4, 8, 16], ylab='threshold ηₙ (window mean ÷ baseline)', yfmt=lambda v: f'{v:g}'))
for nn in (1, 2, 4, 8, 16, 32, 64, 128, 256):
    o.append(f'<path class="axis" d="M{lx(nn):.1f} {A.y0 + A.h} v5"/><text class="sub" x="{lx(nn):.1f}" y="{A.y0 + A.h + 18}" text-anchor="middle">{nn}</text>')
o.append(f'<text class="sub" x="{A.x0 + A.w/2:.1f}" y="{A.y0 + A.h + 34}" text-anchor="middle">window length [channel samples, 0.24 µs each]</text>')
o.append(f'<rect class="fillwhy" x="{A.x0}" y="{A.y0}" width="{lx(16) - A.x0 + 10:.1f}" height="{A.h}"/>')
o.append(f'<text class="sub" x="{A.x0 + 6}" y="{A.Y(1.6):.1f}">windows used: 1–16</text>')
for val, lab in ((1.316, 'pulsar peak at −5 dB: 1.32'), (1.01, 'pulsar peak at −20 dB: 1.01')):
    o.append(f'<path class="dash" d="M{A.x0} {A.Y(val):.1f} H{A.x0 + A.w}"/>')
o.append(f'<text class="sub" x="{lx(16) + 16:.1f}" y="{A.Y(1.316) - 6:.1f}">pulsar peak at −5 dB: 1.32 × baseline</text>')
o.append(f'<text class="sub" x="{A.x0 + 6}" y="{A.Y(1.01) - 6:.1f}">at −20 dB: 1.01</text>')
xs = [lx(k) for k in exact] + [lx(k) for k in longer]
ys = [float(A.Y(v)) for v in exact.values()] + [float(A.Y(v)) for v in longer.values()]
o.append('<path class="curve" d="M' + ' L'.join(f'{x:.1f} {y:.1f}' for x, y in zip(xs[:5], ys[:5])) + '"/>')
o.append('<path class="curve2" stroke-dasharray="5 3" d="M' + ' L'.join(f'{x:.1f} {y:.1f}' for x, y in zip(xs[4:], ys[4:])) + '"/>')
for i, (x, y) in enumerate(zip(xs, ys)):
    o.append(f'<path class="{"dot" if i < 5 else "dot2"}" d="M{x - 3.5:.1f} {y:.1f} a3.5 3.5 0 1 0 7 0 a3.5 3.5 0 1 0 -7 0"/>')
o.append(f'<text class="sub" x="{lx(1) + 8:.1f}" y="{A.Y(13.82) + 4:.1f}">13.82</text>')
o.append(f'<text class="sub" x="{lx(16) + 8:.1f}" y="{A.Y(2.69) - 6:.1f}">2.69</text>')
o.append(f'<text class="sub" x="{lx(256) - 4:.1f}" y="{A.Y(1.11):.1f}" text-anchor="end">{longer[256]:.2f} at 256: the −5 dB pulse plus noise often exceeds it</text>')
svg = '''<figure class="scroll">
<svg class="wide" viewBox="0 0 860 330" role="img" aria-label="detectRFI thresholds against window length, compared with the pulsar's peak power at -5 and -20 dB">
  ''' + '\n  '.join(o) + '''
</svg>
<figcaption>Why the windows stop at 16 samples. The detection threshold of <code>detectRFI</code> at a false-flag probability of 10⁻⁶ against the window length. Blue: the exact
values the code uses (1–16 samples). Orange: longer windows, from a Gamma distribution with the exact variance (it agrees with the exact
form to 0.5 % where both exist). The dashed lines are the pulsar's own peak power: at −5 dB it reaches the 256-sample threshold, at
−20 dB it is far below all of them.</figcaption>
</figure>'''
open('windows.svg', 'w').write(svg)
