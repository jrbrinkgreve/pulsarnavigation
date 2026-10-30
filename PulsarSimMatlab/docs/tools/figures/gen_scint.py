import numpy as np
from plotsvg import Ax
# Interstellar scintillation as simulateFold makes it (C3b): F = |g|^2, g a complex Gaussian field with Gaussian
# correlations - intensity ACF 1/e at the scintle time, 1/2 at the scintle bandwidth (same construction as
# simulateFold's scintField: white complex noise convolved with Gaussian kernels of std sigma/sqrt(2)).
# Left: J0437 at 1.4 GHz (IAR 2020, rescaled from Keith et al. 2013: bandwidth ~740 MHz, time ~2290 s) over 8 h and
# the 400 MHz band; middle: the same with 50 MHz scintles (a pulsar further away / a lower frequency) for contrast;
# bottom: J0437's band-averaged flux factor against time.
rng = np.random.default_rng(3)
def field(T, nT, band, nF, bw, ts):
    dt = T / nT; dnu = band / nF
    sT, sF = ts, bw / np.sqrt(np.log(2))
    pT, pF = int(np.ceil(3*sT/dt)), int(np.ceil(3*sF/dnu))
    kT = np.exp(-(np.arange(-pT, pT+1)*dt)**2 / (2*(sT/np.sqrt(2))**2)); kT /= np.sqrt((kT**2).sum())
    kF = np.exp(-(np.arange(-pF, pF+1)*dnu)**2 / (2*(sF/np.sqrt(2))**2)); kF /= np.sqrt((kF**2).sum())
    e = (rng.standard_normal((nT + 2*pT, nF + 2*pF)) + 1j*rng.standard_normal((nT + 2*pT, nF + 2*pF))) / np.sqrt(2)
    g = np.apply_along_axis(lambda c: np.convolve(c, kT, 'valid'), 0, e)
    g = np.apply_along_axis(lambda r: np.convolve(r, kF, 'valid'), 1, g)
    return np.abs(g)**2
T = 8*3600; nT = 64; nF = 24
Fj = field(T, nT, 400e6, nF, 740e6, 2290.0)
Fn = field(T, nT, 400e6, nF, 50e6, 2290.0)
o = []
def heat(x0, y0, w, h, F, title):
    nT, nF = F.shape; cw, ch = w / nT, h / nF
    q = [f'<text class="ttl" x="{x0}" y="{y0 - 8}">{title}</text>']
    for i in range(nT):
        for j in range(nF):
            op = min(F[i, j] / 3.0, 1.0)
            q.append(f'<rect class="heat" x="{x0 + i*cw:.1f}" y="{y0 + h - (j+1)*ch:.1f}" width="{cw + 0.3:.1f}" height="{ch + 0.3:.1f}" fill-opacity="{op:.2f}"/>')
    q.append(f'<rect class="frame" x="{x0}" y="{y0}" width="{w}" height="{h}" fill="none"/>')
    for hh in (0, 2, 4, 6, 8):
        x = x0 + w * hh / 8
        q.append(f'<text class="sub" x="{x:.1f}" y="{y0 + h + 15}" text-anchor="middle">{hh}</text>')
    q.append(f'<text class="sub" x="{x0 + w/2:.1f}" y="{y0 + h + 30}" text-anchor="middle">time [h]</text>')
    q.append(f'<text class="sub" x="{x0 - 6}" y="{y0 + h + 4}" text-anchor="end">1.2</text>')
    q.append(f'<text class="sub" x="{x0 - 6}" y="{y0 + 8}" text-anchor="end">1.6</text>')
    q.append(f'<text class="sub" x="{x0 - 6}" y="{y0 + h/2 + 4}" text-anchor="end">GHz</text>')
    return '\n  '.join(q)
o.append(heat(60, 34, 340, 120, Fj, 'a · J0437 at 1.4 GHz: scintles ~740 MHz × ~38 min'))
o.append(heat(480, 34, 340, 120, Fn, 'b · the same with 50 MHz scintles'))
yTop = max(4, int(np.ceil(Fj.mean(axis=1).max())))
A = Ax(60, 230, 760, 70, (0, 8), (0, yTop))
o.append(A.frame(xticks=[0, 1, 2, 3, 4, 5, 6, 7, 8], yticks=list(range(0, yTop + 1)), xlab='time [h]',
                 ylab='c · J0437, flux factor averaged over the band (mean 1)'))
t = np.linspace(0, 8, nT)
o.append(f'<path class="dash" d="M{A.x0} {float(A.Y(1)):.1f} H{A.x0 + A.w}"/>')
o.append(f'<path class="curve" d="{A.path(t, Fj.mean(axis=1))}"/>')
svg = '''<figure class="scroll">
<svg class="wide" viewBox="0 0 860 340" role="img" aria-label="Simulated interstellar scintillation: dynamic spectra for J0437 and for narrower scintles, and J0437's band-averaged flux against time">
  ''' + '\n  '.join(o) + '''
</svg>
<figcaption><b>Figure P.7.</b> Interstellar scintillation as the fast simulator draws it (one random realization; stronger colour = more flux, full colour at
3× the mean). (a) J0437 at 1.4 GHz: its scintles (~740 MHz) are wider than our 400 MHz band, so the whole band brightens and fades together,
about every 40 minutes. (b) The same with 50 MHz scintles: the band holds several scintles, and at any moment some channels are bright and
others faint. (c) J0437's flux over the band against time: long stretches well below its catalogue average, and some well above.</figcaption>
</figure>'''
open('scint.svg', 'w').write(svg)
