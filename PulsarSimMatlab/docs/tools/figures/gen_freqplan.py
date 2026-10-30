import numpy as np
from plotsvg import Ax
# The frequency plan (parameters.html): a: the real RF stream (f_in = 4 GHz, so up to 2 GHz), the band fLow-fHigh and the LO.
# b: complex baseband at fs = 800 MHz (-400..+400 MHz): the band moved by -fLO, the IQ low-pass (cutoff 0.9*fs/2 = 360 MHz,
# drawn as flat to 355 and zero at 365 MHz: transition ~9.8 MHz at 80 dB) and the 128 channel edges every 3.125 MHz
# (channel centres at fLow + (j - 1/2)*chanWidth, so the edges fall on the band edges).
f_in, fLow, fHigh, fLO, fs, dF = 4e3, 1200, 1600, 1300, 800, 3.125
o = []
o.append('<defs><marker id="ar-fp" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse">'
         '<path class="ah" d="M0 0 L10 5 L0 10 z"/></marker></defs>')
A = Ax(60, 46, 760, 90, (0, 2200), (0, 1.3))
o.append(A.frame(xticks=[0, 400, 800, 1200, 1600, 2000], yticks=[], xlab='', ylab='', xfmt='{:g}', left=False))
o.append(f'<text class="sub" x="{A.x0 + A.w}" y="{A.y0 + A.h + 18}" text-anchor="end">MHz</text>')
o.append(f'<rect class="fillfg" x="{float(A.X(0)):.1f}" y="{float(A.Y(0.55)):.1f}" width="{float(A.X(2000) - A.X(0)):.1f}" height="{float(A.Y(0) - A.Y(0.55)):.1f}"/>')
o.append(f'<rect class="fillsoft" x="{float(A.X(fLow)):.1f}" y="{float(A.Y(1)):.1f}" width="{float(A.X(fHigh) - A.X(fLow)):.1f}" height="{float(A.Y(0) - A.Y(1)):.1f}"/>')
o.append(f'<path class="curve" d="M{float(A.X(fLow)):.1f} {float(A.Y(0)):.1f} V{float(A.Y(1)):.1f} H{float(A.X(fHigh)):.1f} V{float(A.Y(0)):.1f}"/>')
o.append(f'<text class="lbl" x="{float(A.X(1460)):.1f}" y="{float(A.Y(0.45)):.1f}" text-anchor="middle">band 1.2–1.6 GHz</text>')
o.append(f'<text class="sub" x="{float(A.X(600)):.1f}" y="{float(A.Y(0.2)):.1f}" text-anchor="middle">simulated stream: white up to f_in/2</text>')
for f, lab, anc, dx in ((fLO, 'f_LO 1.3 GHz', 'end', -6), (2000, 'f_in/2 = 2 GHz', 'start', 6)):
    x = float(A.X(f))
    o.append(f'<path class="dash" d="M{x:.1f} {A.y0 - 6} V{A.y0 + A.h}"/>')
    o.append(f'<text class="sub" x="{x + dx:.1f}" y="{A.y0 + 2}" text-anchor="{anc}">{lab}</text>')
o.append(f'<text class="ttl" x="60" y="20">a · RF: a real stream sampled at f_in = 4 GHz</text>')
B = Ax(60, 236, 760, 110, (-450, 450), (0, 1.3))
o.append(f'<path class="line" marker-end="url(#ar-fp)" d="M{float(A.X(1400)):.1f} {A.y0 + A.h + 26} L{float(B.X(100)):.1f} {B.y0 - 14}"/>')
o.append(f'<text class="sub" x="{float(A.X(1400)) + 8:.1f}" y="{A.y0 + A.h + 50}">mix with the LO (shift by −1.3 GHz),</text>')
o.append(f'<text class="sub" x="{float(A.X(1400)) + 8:.1f}" y="{A.y0 + A.h + 64}">low-pass, keep every 5th sample</text>')
o.append(B.frame(xticks=[-400, -300, -200, -100, 0, 100, 200, 300, 400], yticks=[], xlab='', ylab='', left=False))
o.append(f'<text class="sub" x="{B.x0 + B.w}" y="{B.y0 + B.h + 18}" text-anchor="end">MHz</text>')
ff = np.array([-450, -365, -355, 355, 365, 450])
hh = np.array([0, 0, 1.15, 1.15, 0, 0])
o.append(f'<path class="curve2" d="{B.path(ff, hh)}"/>')
o.append(f'<text class="sub" x="{float(B.X(-350)):.1f}" y="{float(B.Y(1.2)):.1f}">IQ low-pass: flat to ±355 MHz, cutoff 360</text>')
o.append(f'<rect class="fillsoft" x="{float(B.X(-100)):.1f}" y="{float(B.Y(1)):.1f}" width="{float(B.X(300) - B.X(-100)):.1f}" height="{float(B.Y(0) - B.Y(1)):.1f}"/>')
edges = np.arange(-100, 300 + 1e-9, dF)
o.append('<path class="thin" d="' + ' '.join(f'M{float(B.X(e)):.1f} {float(B.Y(0)):.1f} V{float(B.Y(1)):.1f}' for e in edges) + '" opacity="0.6"/>')
o.append(f'<path class="curve" d="M{float(B.X(-100)):.1f} {float(B.Y(0)):.1f} V{float(B.Y(1)):.1f} H{float(B.X(300)):.1f} V{float(B.Y(0)):.1f}"/>')
o.append(f'<rect x="{float(B.X(30)):.1f}" y="{float(B.Y(0.66)):.1f}" width="{float(B.X(170) - B.X(30)):.1f}" height="34" rx="4" fill="var(--bg)" opacity="0.9"/>')
o.append(f'<text class="lbl" x="{float(B.X(100)):.1f}" y="{float(B.Y(0.66)) + 14:.1f}" text-anchor="middle">band −100…+300 MHz</text>')
o.append(f'<text class="sub" x="{float(B.X(100)):.1f}" y="{float(B.Y(0.66)) + 28:.1f}" text-anchor="middle">128 channels of 3.125 MHz</text>')
for f in (-400, 400):
    x = float(B.X(f))
    o.append(f'<path class="dash" d="M{x:.1f} {B.y0 - 6} V{B.y0 + B.h}"/>')
o.append(f'<text class="sub" x="{float(B.X(-400)) - 6:.1f}" y="{B.y0 + 2}" text-anchor="end">−fs/2</text>')
o.append(f'<text class="sub" x="{float(B.X(400)) + 6:.1f}" y="{B.y0 + 2}">fs/2</text>')
o.append(f'<path class="dash" d="M{float(B.X(0)):.1f} {B.y0 + 10} V{B.y0 + B.h}"/>')
o.append(f'<text class="sub" x="{float(B.X(0)):.1f}" y="{B.y0 + B.h + 32}" text-anchor="middle">(f_LO)</text>')
o.append(f'<text class="ttl" x="60" y="{B.y0 - 26}">b · complex baseband at fs = 800 MHz</text>')
svg = '''<figure class="scroll">
<svg class="wide" viewBox="0 0 860 380" role="img" aria-label="Frequency plan: the RF band 1.2 to 1.6 GHz and the LO; at baseband the band sits at -100 to +300 MHz inside the IQ low-pass and the 800 MHz window, split into 128 channels">
  ''' + '\n  '.join(o) + '''
</svg>
<figcaption><b>Figure 1.</b> The frequency plan, from the parameters. (a) The simulated RF stream is real and sampled at 4 GHz, so it holds frequencies up to
2 GHz; the receiver band is 1.2–1.6 GHz and the LO sits at 1.3 GHz, off centre on purpose. (b) After the IQ stage the band lies at
−100…+300 MHz, lopsided around zero, so a mirrored frequency mapping would put it at −300…+100 MHz and fail visibly. The +300 MHz edge is
why the sample rate is 800 MHz: it must sit inside the flat part of the IQ low-pass (cutoff 0.9·fs/2 = 360 MHz). The channelizer then cuts
the band into 128 channels whose edges fall on the band edges (fine lines).</figcaption>
</figure>'''
open('freqplan.svg', 'w').write(svg)
