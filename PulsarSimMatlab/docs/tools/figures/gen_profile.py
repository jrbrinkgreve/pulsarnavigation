import numpy as np
from plotsvg import Ax, pow10
T = 10.0  # ms
def G(t, duty, neighbours=True):
    sig = duty/100*T/(2*np.sqrt(2*np.log(2)))
    n = int(np.floor(np.sqrt(2*np.log(1e30))*sig/T + 0.5)) if neighbours else 0
    k0 = np.round(t/T - 0.5)
    g = sum(np.exp(-(t - (k0 + j + 0.5)*T)**2/(2*sig**2)) for j in range(-n, n+1))
    peak = sum(np.exp(-(j*T)**2/(2*sig**2)) for j in range(-n, n+1))
    return g/peak, sig, n
o = []
# panel A
A = Ax(70, 46, 330, 150, (4.0, 6.0), (0, 1.08))
t = np.linspace(4, 6, 801); g, sig, n = G(t, 5)
print('5%: sigma', sig, 'neighbours', n)
o.append(A.frame(xticks=[4, 4.5, 5, 5.5, 6], yticks=[0, 0.5, 1], xlab='time [ms]', ylab='G(t)'))
o.append(f'<path class="fillsoft" d="{A.area(t, g)}"/>')
o.append(f'<path class="curve" d="{A.path(t, g)}"/>')
xa, xb = A.X(5 - 0.25), A.X(5 + 0.25); y5 = A.Y(0.5)
o.append(f'<path class="line" marker-start="url(#ar-pp)" marker-end="url(#ar-pp)" d="M{xa:.1f} {y5:.1f} H{xb:.1f}"/>')
o.append(f'<text class="lbl" x="{A.X(5.42):.1f}" y="{y5+4:.1f}">FWHM 0.5 ms</text>')
o.append(f'<text class="sub" x="{A.X(5.42):.1f}" y="{y5+19:.1f}">= 5 % of T</text>')
ys = A.Y(np.exp(-0.5)); xs = A.X(5 + sig)
o.append(f'<path class="dash" d="M{A.X(5):.1f} {A.Y(0):.1f} V{A.Y(1):.1f}"/>')
o.append(f'<path class="thin" d="M{A.X(5):.1f} {ys:.1f} H{xs:.1f}"/>')
o.append(f'<text class="sub" x="{xs+6:.1f}" y="{ys-6:.1f}">σ = 212.3 µs</text>')
o.append(f'<text class="sub" x="{A.X(5)-6:.1f}" y="{A.Y(1)-6:.1f}" text-anchor="end">centre (k + ½)·T</text>')
o.append(f'<text class="ttl" x="{A.x0+A.w/2}" y="20" text-anchor="middle">5 % duty (default): one pulse</text>')
# panel B: the tails halfway between the pulses at 5 and 15 ms
B = Ax(520, 46, 320, 150, (9, 11), (1e-8, 2e-5), logy=True)
t = np.linspace(9, 11, 2001); gs, sig2, n2 = G(t, 20); gn, _, _ = G(t, 20, neighbours=False)
gm = G(np.array([10.0]), 20)[0][0]; gnm = G(np.array([10.0]), 20, False)[0][0]
print('20%: sigma', sig2, 'neighbours', n2, 'mid', gm, gnm)
o.append(B.frame(xticks=[9, 9.5, 10, 10.5, 11], yticks=[1e-8, 1e-7, 1e-6, 1e-5], xlab='time [ms] (pulse centres at 5 and 15 ms)', ylab='G(t), log scale', yfmt=pow10))
o.append(f'<path class="curve2" stroke-dasharray="6 4" d="{B.path(t, gn)}"/>')
o.append(f'<path class="curve" d="{B.path(t, gs)}"/>')
xm = B.X(10)
o.append(f'<path class="dot" d="M{xm-3.5:.1f} {B.Y(gm):.1f} a3.5 3.5 0 1 0 7 0 a3.5 3.5 0 1 0 -7 0"/>')
o.append(f'<path class="dot2" d="M{xm-3.5:.1f} {B.Y(gnm):.1f} a3.5 3.5 0 1 0 7 0 a3.5 3.5 0 1 0 -7 0"/>')
xl = B.X(10.32)
o.append(f'<path class="thin" d="M{xm+4:.1f} {B.Y(gm):.1f} H{xl-4:.1f}"/><text class="sub" x="{xl:.1f}" y="{B.Y(gm)+4:.1f}">5.96·10⁻⁸ (the sum)</text>')
o.append(f'<path class="thin" d="M{xm+4:.1f} {B.Y(gnm):.1f} H{xl-4:.1f}"/><text class="sub" x="{xl:.1f}" y="{B.Y(gnm)+4:.1f}">2.98·10⁻⁸ (nearest)</text>')
o.append(f'<text class="ttl" x="{B.x0+B.w/2}" y="20" text-anchor="middle">20 % duty (widest real pulses): the tails</text>')
lx = float(B.X(9.5))
o.append(f'<path class="curve" d="M{lx} 60 h22"/><text class="sub" x="{lx+28}" y="64">± 1 neighbour (the code)</text>')
o.append(f'<path class="curve2" stroke-dasharray="6 4" d="M{lx} 78 h22"/><text class="sub" x="{lx+28}" y="82">nearest pulse only</text>')
svg = '''<figure class="scroll">
<svg class="wide" viewBox="0 0 860 245" role="img" aria-label="The power profile G(t): a Gaussian per pulse, and the periodic sum with neighbouring pulses">
  <defs>
    <marker id="ar-pp" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
      <path class="ah" d="M0 0 L10 5 L0 10 z"/>
    </marker>
  </defs>
  ''' + '\n  '.join(o) + '''
</svg>
<figcaption>Left: the power profile <var>G</var> of one pulse at the default 5 % duty, computed as in the code. Right: at 20 % duty
the generator adds one neighbouring pulse on each side (<code>nNeighbors</code> = 1). It matters only far out in the tails: halfway
between two pulses both contribute equally, so the nearest pulse alone would give half the true value (and a kink). At 5 % the neighbours are below 10⁻¹²⁰ there,
so the code uses the nearest pulse alone. Test 4 of <code>testGeneratePulsarSignal</code> checks exactly these tail bins.</figcaption>
</figure>'''
open('profile.svg', 'w').write(svg)
