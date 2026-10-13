import numpy as np
from plotsvg import Ax
# Where the simulated power sits in frequency: generator output (white up to f_in/2) and what the
# dispersion stage passes (|H|^2 = W^2, the receiver band). Constants as in the code.
fin, fL, fH = 4000.0, 1200.0, 1600.0             # MHz
edge = 0.02*(fH - fL)
f = np.unique(np.r_[np.linspace(0, 2400, 25), np.linspace(fL, fL + 0.02*(fH - fL), 25),
                   np.linspace(fH - 0.02*(fH - fL), fH, 25), fL - 1e-3, fH + 1e-3, fin/2 + 1e-3])
gen = np.where(f <= fin/2, 1.0, 0.0)
W = np.where((f >= fL) & (f <= fH), 1.0, 0.0)
lo = (f >= fL) & (f < fL + edge); hi = (f > fH - edge) & (f <= fH)
W[lo] = np.sin(np.pi/2*(f[lo] - fL)/edge)**2; W[hi] = np.sin(np.pi/2*(fH - f[hi])/edge)**2
frac = 2*(fH - fL - 1.25*edge)/fin
print('kernel energy (fraction of white power passed)', frac)
o = []
P = Ax(70, 58, 760, 130, (0, 2400), (0, 1.35))
o.append(P.frame(xticks=[0, 400, 800, 1200, 1600, 2000, 2400], xlab='frequency [MHz]', ylab='power spectral density (relative)'))
o.append(f'<path class="fillfg" d="{P.area(f, gen)}"/>')
o.append(f'<path class="fillsoft" d="{P.area(f, W**2)}"/>')
o.append(f'<path class="curve" d="{P.path(f, W**2)}"/>')
o.append(f'<path class="line" d="{P.path(f, gen)}"/>')
xn = P.X(fin/2)
o.append(f'<path class="dash" d="M{xn:.1f} {P.y0-30:.1f} V{P.y0+P.h:.1f}"/>')
o.append(f'<text class="lbl" x="{xn+6:.1f}" y="{P.y0-20:.1f}">f_in/2 = 2 GHz</text>')
o.append(f'<text class="sub" x="{xn+6:.1f}" y="{P.y0-6:.1f}">Nyquist limit of a</text>')
o.append(f'<text class="sub" x="{xn+6:.1f}" y="{P.y0+8:.1f}">real 4 GHz stream</text>')
o.append(f'<text class="sub" x="{P.X(600):.1f}" y="{P.Y(1)+18:.1f}" text-anchor="middle">generator output: white from 0 to 2 GHz (grey)</text>')
o.append(f'<text class="lbl" x="{P.X(1400):.1f}" y="{P.Y(1)-10:.1f}" text-anchor="middle">receiver band 1.2–1.6 GHz</text>')
o.append(f'<text class="sub" x="{P.X(1400):.1f}" y="{P.Y(0.55):.1f}" text-anchor="middle">kept by the</text>')
o.append(f'<text class="sub" x="{P.X(1400):.1f}" y="{P.Y(0.55)+14:.1f}" text-anchor="middle">dispersion stage</text>')
o.append(f'<text class="sub" x="{P.X(1400):.1f}" y="{P.Y(0.55)+28:.1f}" text-anchor="middle">({100*frac:.1f} % of the power)</text>')
# margin arrow between 1.6 and 2 GHz
ym = P.Y(0.2)
o.append(f'<path class="line" marker-start="url(#ar-rf)" marker-end="url(#ar-rf)" d="M{P.X(fH)+2:.1f} {ym:.1f} H{xn-2:.1f}"/>')
o.append(f'<text class="sub" x="{P.X(1800):.1f}" y="{ym-6:.1f}" text-anchor="middle">margin 0.4 GHz</text>')
svg = '''<figure class="scroll">
<svg class="wide" viewBox="0 0 860 240" role="img" aria-label="Spectrum of the simulated RF stream: white up to 2 GHz from the generator, the 1.2 to 1.6 GHz band kept by the dispersion stage">
  <defs>
    <marker id="ar-rf" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
      <path class="ah" d="M0 0 L10 5 L0 10 z"/>
    </marker>
  </defs>
  ''' + '\n  '.join(o) + '''
</svg>
<figcaption>The simulated RF stream in frequency, with the code's constants. A real (not complex) stream sampled at <var>f</var><sub>in</sub> = 4 GHz
holds frequencies up to 2 GHz. The generator fills all of it with white noise; the dispersion stage keeps only the receiver band,
1.2–1.6 GHz, with 8 MHz tapered edges (too narrow to see here). So 19.5 % of the generated power passes, the kernel energy Σ<var>h</var>²
that <code>applyDispersionStream</code> reports.</figcaption>
</figure>'''
open('rfspec.svg', 'w').write(svg)
