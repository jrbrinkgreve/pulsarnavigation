import numpy as np
from plotsvg import Ax
# The pulse seen through 16 sub-bands of 25 MHz, before and after applyDispersionStream.
# Expected power (no noise). After dispersion each frequency f is delayed by tau(f) - tau(1.6 GHz);
# a sub-band's power is the average over its frequencies, so the pulse also smears inside a sub-band.
K, DM = 4.148808e3, 5.0                            # s MHz^2 pc^-1 cm^3
fL, fH, T, duty = 1200.0, 1600.0, 10.0, 5.0        # MHz, MHz, ms, %
sig = duty/100*T/(2*np.sqrt(2*np.log(2)))
dtau = lambda f: (K*DM/f**2 - K*DM/fH**2)*1e3      # ms
def G(t):                                          # pulses at 5, 15, 25 ms (the file starts at 0)
    return sum(np.exp(-(t - c)**2/(2*sig**2)) for c in (5.0, 15.0, 25.0))
nSub, bw = 16, 25.0
t = np.linspace(0, 23, 2301)
o = []
rowH, top = 15.0, 52
def panel(x0, w, after, title):
    A = Ax(x0, top, w, rowH*nSub, (0, 23), (0, 1))
    o.append(A.frame(xticks=[0, 5, 10, 15, 20], xlab='time [ms]', left=False))
    o.append(f'<text class="ttl" x="{x0}" y="22">{title}</text>')
    for i in range(nSub):
        fc = fH - (i + 0.5)*bw                     # top row = highest frequency
        yb = top + (i + 1)*rowH                    # baseline of this row
        if after:
            fs = fc + np.linspace(-bw/2, bw/2, 41)
            p = np.mean([G(t - dtau(f)) for f in fs], axis=0)
        else:
            p = G(t)
        X = np.round(A.X(t)[::2], 1); Y = np.round(yb - 0.95*rowH*p[::2], 1)
        keep = np.r_[True, (Y[1:-1] != Y[:-2]) | (Y[1:-1] != Y[2:]), True]   # drop points inside flat runs
        pts = ' L'.join(f'{a:.1f} {b:.1f}' for a, b in zip(X[keep], Y[keep]))
        o.append(f'<path class="curve" stroke-width="1.4" d="M{pts}"/>')
        if after and i == nSub - 1:
            print('bottom sub-band peak', p.max())
    if after:
        ff = np.linspace(fL + bw/2, fH - bw/2, 200)
        for c in (5.0, 15.0):
            xs = A.X(c + dtau(ff)); ys = top + (fH - ff)/bw*rowH
            pts = ' L'.join(f'{a:.1f} {b:.1f}' for a, b in zip(xs, ys))
            o.append(f'<path class="curve2" stroke-dasharray="5 3" stroke-width="1.4" d="M{pts}"/>')
    return A
A = panel(110, 330, False, 'a · as generated: every frequency at the same time')
B = panel(520, 330, True, 'b · after dispersion (DM 5): up to 6.3 ms late')
for f in (1600, 1500, 1400, 1300, 1200):
    y = top + (fH - f)/bw*rowH
    o.append(f'<text class="sub" x="{A.x0-8}" y="{y+4:.1f}" text-anchor="end">{f} MHz</text>')
    o.append(f'<path class="axis" d="M{A.x0} {y:.1f} h-4"/>')
svg = '''<figure class="scroll">
<svg class="wide" viewBox="0 0 860 336" role="img" aria-label="The pulse in 16 sub-bands before and after dispersion: after dispersion lower frequencies arrive later, along a 1/f-squared curve">
  ''' + '\n  '.join(o) + '''
</svg>
<figcaption>The expected pulse power (no noise) in 16 sub-bands of 25 MHz, with the code's constants. a: the generator's pulses arrive at all
frequencies at once. b: after <code>applyDispersionStream</code> each frequency <var>f</var> is late by τ(<var>f</var>) − τ(1.6 GHz) (dashed), up to
6.3 ms at the bottom: more than half a turn and 12 pulse widths. Inside each sub-band the delay changes too, so the pulse is also
smeared (by 0.6 ms in the lowest 25 MHz, which is why its peak is lower). Shifting whole sub-bands back cannot undo that smear;
filtering the voltage with the inverse chirp can (stages 5 and 5b).</figcaption>
</figure>'''
open('waterfall.svg', 'w').write(svg)
