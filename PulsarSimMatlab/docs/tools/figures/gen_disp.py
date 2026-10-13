import numpy as np
from plotsvg import Ax
K, DM, fL, fH = 4.148808e3, 5.0, 1200.0, 1600.0   # MHz units
edge = 0.02*(fH - fL)
f = np.linspace(1150, 1650, 5001)
W = np.where((f >= fL) & (f <= fH), 1.0, 0.0)
lo = (f >= fL) & (f < fL + edge); hi = (f > fH - edge) & (f <= fH)
W[lo] = np.sin(np.pi/2*(f[lo] - fL)/edge)**2; W[hi] = np.sin(np.pi/2*(fH - f[hi])/edge)**2
tau = (K*DM/f**2 - K*DM/fH**2)*1e3   # ms
print('delay at 1.2 GHz', (K*DM/fL**2 - K*DM/fH**2)*1e3, 'bulk', K*DM/fH**2*1e3)
o = []
A = Ax(90, 34, 720, 50, (1150, 1650), (0, 1.15))
o.append(A.frame(yticks=[0, 1], ylab='W(f)'))
o.append(f'<path class="fillsoft" d="{A.area(f, W)}"/>')
o.append(f'<path class="curve" d="{A.path(f, W)}"/>')
o.append(f'<text class="sub" x="{A.X(1400):.1f}" y="{A.Y(1)+16:.1f}" text-anchor="middle">receiver band: W = 1 inside, 0 outside; sin² edges 8 MHz wide, just inside 1.2 and 1.6 GHz</text>')
B = Ax(90, 120, 720, 150, (1150, 1650), (0, 7))
inb = (f >= fL) & (f <= fH)
o.append(B.frame(xticks=[1200, 1300, 1400, 1500, 1600], yticks=[0, 2, 4, 6], xlab='RF frequency [MHz]', ylab='delay relative to 1.6 GHz [ms]'))
for fx in (fL, fH):
    o.append(f'<path class="dash" d="M{B.X(fx):.1f} {A.y0:.1f} V{A.y0+A.h:.1f} M{B.X(fx):.1f} {B.y0:.1f} V{B.y0+B.h:.1f}"/>')
o.append(f'<path class="curve2" d="{B.path(f[inb], tau[inb])}"/>')
x12, y12 = B.X(fL), B.Y(6.302)
o.append(f'<path class="dot2" d="M{x12-4:.1f} {y12:.1f} a4 4 0 1 0 8 0 a4 4 0 1 0 -8 0"/>')
o.append(f'<path class="thin" d="M{x12+5:.1f} {y12+5:.1f} L{x12+16:.1f} {B.Y(2.3)-4:.1f}"/>')
o.append(f'<text class="lbl" x="{x12+18:.1f}" y="{B.Y(2.3)+8:.1f}">6.302 ms at 1.2 GHz</text>')
o.append(f'<text class="sub" x="{x12+18:.1f}" y="{B.Y(2.3)+24:.1f}">→ a kernel of 25.2 million taps at 4 GHz</text>')
xr = B.X(fH)
o.append(f'<path class="dot2" d="M{xr-4:.1f} {B.Y(0):.1f} a4 4 0 1 0 8 0 a4 4 0 1 0 -8 0"/>')
o.append(f'<text class="lbl" x="{xr-14:.1f}" y="{B.Y(4.5):.1f}" text-anchor="end">RefFreq = 1.6 GHz keeps its time</text>')
o.append(f'<text class="sub" x="{xr-14:.1f}" y="{B.Y(4.5)+16:.1f}" text-anchor="end">(the bulk delay τ(1.6 GHz) = 8.103 ms is not applied)</text>')
o.append(f'<path class="thin" d="M{xr-8:.1f} {B.Y(4.5)-4:.1f} H{xr:.1f} V{B.Y(0)-6:.1f}"/>')
o.append(f'<text class="sub" x="{B.X(1310):.1f}" y="{B.Y(6.4):.1f}">τ(f) − τ(f_ref) = K·DM·(1/f² − 1/f_ref²), K = 4.148808·10³ s MHz², DM 5</text>')
svg = '''<figure class="scroll">
<svg class="wide" viewBox="0 0 860 310" role="img" aria-label="Dispersion delay across the band and the band taper W(f)">
  ''' + '\n  '.join(o) + '''
</svg>
<figcaption>What <code>applyDispersionStream</code> applies, computed with the code's constants. Bottom: the delay relative to 1.6 GHz;
low frequencies arrive up to 6.302 ms later. Top: the taper <var>W</var>(<var>f</var>) that limits the filter to the band (it models the
receiver bandpass); outside it nothing is passed, so no finite kernel has to represent the delays that diverge towards 0 Hz.</figcaption>
</figure>'''
open('dispersion.svg', 'w').write(svg)
