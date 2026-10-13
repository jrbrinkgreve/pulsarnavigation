import numpy as np
from plotsvg import Ax
# The IQ low-pass of applyIQmodulation (kaiserLowpass: 2049 taps, cutoff 0.9*fs/2 = 360 MHz, 80 dB, at
# f_in = 4 GHz) against baseband frequency, with what it has to keep and what it has to remove.
fin, fs, fLO = 4000.0, 800.0, 1300.0                # MHz
M, A_db = 2049, 80.0
fc = 0.9*fs/2
beta = 0.1102*(A_db - 8.7)
G = (M - 1)//2; n = np.arange(M) - G
h = np.where(n == 0, 2*fc/fin, np.sin(2*np.pi*fc/fin*n)/(np.pi*np.where(n == 0, 1, n)))
w = np.i0(beta*np.sqrt(1 - (2*np.arange(M)/(M - 1) - 1)**2))/np.i0(beta)
h = h*w; h /= h.sum()
Nf = 2**20
H = np.fft.fftshift(np.abs(np.fft.fft(h, Nf)))
f = np.fft.fftshift(np.fft.fftfreq(Nf, 1/fin))
HdB = 20*np.log10(np.maximum(H, 1e-12))
transW = (A_db - 8)/(2.285*(M - 1))*fin/(2*np.pi)
print('transition', transW, 'flat to', fc - transW/2, 'stopband from', fc + transW/2,
      'max beyond stop', HdB[np.abs(f) > fc + transW/2].max())
sel = np.abs(f) <= 1650
fp, Hp = f[sel], HdB[sel]
# sidelobe envelope per 2 MHz bin keeps the path small
bins = np.arange(-1650, 1652, 2.0)
idx = np.digitize(fp, bins)
fx = np.array([fp[idx == i].mean() for i in np.unique(idx)])
hx = np.array([Hp[idx == i].max() for i in np.unique(idx)])
o = []
A = Ax(80, 60, 750, 210, (-1650, 1650), (-120, 8))
o.append(A.frame(xticks=[-1600, -1200, -800, -400, 0, 400, 800, 1200, 1600], yticks=[-120, -80, -40, 0],
                 xlab='baseband frequency [MHz]', ylab='low-pass gain [dB]'))
def zone(a, b, cls, y0=None, y1=None):
    return f'<rect class="{cls}" x="{A.X(a):.1f}" y="{A.y0 if y0 is None else y0:.1f}" width="{A.X(b) - A.X(a):.1f}" height="{(A.h if y1 is None else y1 - y0):.1f}"/>'
o.append(zone(-100, 300, 'fillwhy'))
o.append(zone(700, 1100, 'fillwarn')); o.append(zone(-900, -500, 'fillwarn'))
o.append(zone(1100, 1500, 'fillfg'))
o.append(f'<path class="curve" d="{A.path(fx, hx)}"/>')
for x in (-400, 400):
    o.append(f'<path class="dash" d="M{A.X(x):.1f} {A.y0 - 14} V{A.y0 + A.h}"/>')
o.append(f'<text class="sub" x="{A.X(400) + 4:.1f}" y="{A.y0 - 4}">output Nyquist ±400 MHz</text>')
o.append(f'<text class="lbl" x="{A.X(100):.1f}" y="{A.Y(-30):.1f}" text-anchor="middle">the band</text>')
o.append(f'<text class="sub" x="{A.X(100):.1f}" y="{A.Y(-30) + 14:.1f}" text-anchor="middle">−100…+300 MHz</text>')
o.append(f'<text class="sub" x="{A.X(100):.1f}" y="{A.Y(-30) + 28:.1f}" text-anchor="middle">flat (to ±355)</text>')
o.append(f'<text class="sub" x="{A.X(900):.1f}" y="{A.Y(-30):.1f}" text-anchor="middle">would fold onto</text>')
o.append(f'<text class="sub" x="{A.X(900):.1f}" y="{A.Y(-30) + 14:.1f}" text-anchor="middle">the band after</text>')
o.append(f'<text class="sub" x="{A.X(900):.1f}" y="{A.Y(-30) + 28:.1f}" text-anchor="middle">keeping every 5th</text>')
o.append(f'<text class="sub" x="{A.X(-700):.1f}" y="{A.Y(-30):.1f}" text-anchor="middle">folds onto</text>')
o.append(f'<text class="sub" x="{A.X(-700):.1f}" y="{A.Y(-30) + 14:.1f}" text-anchor="middle">the band too</text>')
o.append(f'<text class="sub" x="{A.X(1300):.1f}" y="{A.Y(-30):.1f}" text-anchor="middle">mirror image</text>')
o.append(f'<text class="sub" x="{A.X(1300):.1f}" y="{A.Y(-30) + 14:.1f}" text-anchor="middle">of the band</text>')
o.append(f'<text class="sub" x="{A.X(1300):.1f}" y="{A.Y(-30) + 28:.1f}" text-anchor="middle">(1.1–1.5 GHz)</text>')
o.append(f'<text class="sub" x="{A.X(-1600):.1f}" y="{A.Y(-80) - 6:.1f}">stopband: −80 dB or lower beyond ±365 MHz</text>')
o.append(f'<text class="ttl" x="80" y="22">The I/Q low-pass: 2049 taps, cut-off 360 MHz, 80 dB stopband (at 4 GHz, after mixing with the 1.3 GHz LO)</text>')
svg = '''<figure class="scroll">
<svg class="wide" viewBox="0 0 860 322" role="img" aria-label="Gain of the IQ low-pass filter against baseband frequency, with the band, the zones that would alias onto the band and the mirror image">
  ''' + '\n  '.join(o) + '''
</svg>
<figcaption>The low-pass of <code>applyIQmodulation</code>, computed with its design formula (Kaiser-windowed sinc). After mixing, the
band sits at −100…+300 MHz (green): inside the flat passband. Keeping every 5th sample (800 MHz) folds frequencies 800 MHz apart onto
each other, so the orange zones would land on the band; the filter suppresses them by about 80 dB or more. The mixer's mirror image of the
band (grey) and the receiver noise at all other frequencies are removed the same way.</figcaption>
</figure>'''
open('iqfilter.svg', 'w').write(svg)
