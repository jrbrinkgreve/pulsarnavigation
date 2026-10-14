import numpy as np
from plotsvg import Ax
# Where the timing information is (physics.html, section 7). a: the template p (Gaussian, FWHM 5 % of the turn) and the timing
# information per unit time, proportional to p'^2 in the weak-signal limit (1/sigma_TOA^2 = B*integral (rho*p')^2 dt), both
# scaled to peak 1. b: the weak-signal bound sigma_TOA = sqrt(2)*sigma_t/SNR for the simulated pulsar (sigma_t = 212.3 us) and
# for J0437-4715 (W50 = 0.141 ms, sigma_t = 59.9 us), with the design threshold SNR 7.
st = 0.05/(2*np.sqrt(2*np.log(2)))
phi = np.linspace(-0.12, 0.12, 241)
p = np.exp(-0.5*(phi/st)**2)
d2 = (phi/st**2*p)**2
d2 /= d2.max()
o = []
A = Ax(60, 54, 340, 230, (-0.12, 0.12), (0, 1.12))
o.append(A.frame(xticks=[-0.1, -0.05, 0, 0.05, 0.1], yticks=[0, 0.5, 1], xlab='rotational phase [turns]', ylab='relative'))
o.append(f'<path class="fillwhy" d="{A.area(phi, d2)}"/>')
o.append(f'<path class="curve3" d="{A.path(phi, d2)}"/>')
o.append(f'<path class="curve" d="{A.path(phi, p)}"/>')
o.append(f'<text class="sub" x="{float(A.X(0)):.1f}" y="{float(A.Y(1.04)):.1f}" text-anchor="middle">pulse (template)</text>')
for s in (-1, 1):
    x = float(A.X(s*st))
    o.append(f'<path class="dash" d="M{x:.1f} {A.y0 + A.h} V{float(A.Y(1)):.1f}"/>')
o.append(f'<text class="sub" x="{float(A.X(-0.105)):.1f}" y="{float(A.Y(0.8)):.1f}">timing information</text>')
o.append(f'<text class="sub" x="{float(A.X(-0.105)):.1f}" y="{float(A.Y(0.8)) + 14:.1f}">(green): on the flanks,</text>')
o.append(f'<text class="sub" x="{float(A.X(-0.105)):.1f}" y="{float(A.Y(0.8)) + 28:.1f}">none at the peak</text>')
o.append(f'<text class="ttl" x="60" y="20">a · where a TOA gets its precision</text>')
B = Ax(530, 54, 300, 230, (0, 3), (0.1, 1000), logy=True)
lx = lambda s: B.x0 + np.log10(s)/3*B.w
o.append(f'<path class="axis" d="M{B.x0} {B.y0 + B.h} H{B.x0 + B.w} M{B.x0} {B.y0} V{B.y0 + B.h}"/>')
for s in (1, 3, 7, 10, 30, 100, 300, 1000):
    o.append(f'<path class="axis" d="M{lx(s):.1f} {B.y0 + B.h} v5"/><text class="sub" x="{lx(s):.1f}" y="{B.y0 + B.h + 18}" text-anchor="middle">{s}</text>')
for v in (0.1, 1, 10, 100, 1000):
    o.append(f'<path class="axis" d="M{B.x0} {float(B.Y(v)):.1f} h-5"/><text class="sub" x="{B.x0 - 8}" y="{float(B.Y(v)) + 4:.1f}" text-anchor="end">{v:g}</text>')
o.append(f'<text class="sub" x="{B.x0 + B.w/2:.1f}" y="{B.y0 + B.h + 34}" text-anchor="middle">SNR of the folded profile</text>')
o.append(f'<text class="sub" x="{B.x0}" y="{B.y0 - 8}">best TOA error σ<tspan baseline-shift="sub" font-size="9">TOA</tspan> [µs]</text>')
x7 = lx(7)
o.append(f'<path class="dash" d="M{x7:.1f} {B.y0} V{B.y0 + B.h}"/>')
ss = np.logspace(0, 3, 50)
for sig_t, cls, lab in ((212.3, 'curve', 'simulated pulsar, σ<tspan baseline-shift="sub" font-size="9">t</tspan> 212 µs'), (59.9, 'curve2', 'J0437−4715, σ<tspan baseline-shift="sub" font-size="9">t</tspan> 60 µs')):
    o.append(f'<path class="{cls}" d="M{" L".join(f"{lx(s):.1f} {float(B.Y(np.sqrt(2)*sig_t/s)):.1f}" for s in ss)}"/>')
    v7 = np.sqrt(2)*sig_t/7
    y7 = float(B.Y(v7))
    o.append(f'<path class="{"dot" if cls == "curve" else "dot2"}" d="M{x7 - 4:.1f} {y7:.1f} a4 4 0 1 0 8 0 a4 4 0 1 0 -8 0"/>')
    o.append(f'<text class="sub" x="{x7 + 9:.1f}" y="{y7 + (-6 if cls == "curve" else 16):.1f}">{v7:.0f} µs</text>')
    xl, yl = lx(40), float(B.Y(np.sqrt(2)*sig_t/40)) - 8
    ang = np.degrees(np.arctan((B.h/4)/(B.w/3)))
    o.append(f'<text class="sub" x="{xl:.1f}" y="{yl:.1f}" transform="rotate({ang:.1f} {xl:.1f} {yl:.1f})">{lab}</text>')
o.append(f'<text class="sub" x="{x7 + 4:.1f}" y="{B.y0 + 12}">SNR 7</text>')
o.append(f'<text class="ttl" x="530" y="20">b · the bound σ<tspan baseline-shift="sub" font-size="10">TOA</tspan> ≈ √2·σ<tspan baseline-shift="sub" font-size="10">t</tspan> / SNR</text>')
svg = '''<figure class="scroll">
<svg class="wide" viewBox="0 0 860 344" role="img" aria-label="Timing information sits on the flanks of the pulse; the best TOA error falls as pulse width over SNR">
  ''' + '\n  '.join(o) + '''
</svg>
<figcaption>Where the precision of a TOA comes from (weak-signal limit). a: a small shift of the pulse changes the profile most where
it is steepest, so the timing information (∝ <var>p</var>′², green) sits on the flanks, one σ<sub>t</sub> either side of the peak
(dashed); the peak itself says nothing about time. b: the resulting bound for a Gaussian pulse: proportional to the pulse width and
inversely to the SNR. At the design threshold SNR 7 the simulated pulsar gives ~43 µs and the narrower J0437 ~12 µs.</figcaption>
</figure>'''
open('flanks.svg', 'w').write(svg)
