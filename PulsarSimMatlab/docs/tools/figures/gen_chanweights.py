import numpy as np
from plotsvg import Ax
# The built-in notch of optimal weighting (toa.html section 6): estimateTOA weights each channel by 1/variance with the
# radiometer model, i.e. by 1/level^2. Expected level of each 3.125 MHz channel with the illustrative scenario's continuous
# sources (pipelineParams.m: GNSS L1 bpsk, GNSS L2 bpsk, carrier), from their INR: level = 1 + 128 * INR * (share of the
# source's power in that channel); BPSK power spectrum sinc^2 with the chip rate. Expected values, not a measurement.
fL, cw, nch = 1200.0, 3.125, 128
edges = fL + cw*np.arange(nch + 1)
srcs = [('bpsk', 1575.42, 1.023, -10, 'GNSS L1 (−10 dB)'), ('bpsk', 1227.60, 10.23, -15, 'GNSS L2 (−15 dB)'),
        ('cw', 1351.3, 0, -5, 'carrier (−5 dB)')]
level = np.ones(nch)
for typ, f0, rc, inr, _ in srcs:
    if typ == 'cw':
        share = np.zeros(nch); share[int((f0 - fL)//cw)] = 1.0
    else:
        f = np.linspace(f0 - 40*rc, f0 + 40*rc, 400001)
        psd = np.sinc((f - f0)/rc)**2/rc
        share = np.array([np.trapezoid(psd[(f >= a) & (f < b)], f[(f >= a) & (f < b)]) if np.any((f >= a) & (f < b)) else 0 for a, b in zip(edges[:-1], edges[1:])])
    level += nch*10**(inr/10)*share
w = 1/level**2
print('min weights', np.sort(w)[:8])
print('SNR with optimal weights / clean:', np.sqrt(w.sum()/nch))
o = []
A = Ax(90, 40, 720, 220, (1200, 1600), (1e-4, 60), logy=True)
o.append(A.frame(xticks=[1200, 1300, 1400, 1500, 1600], yticks=[1e-4, 1e-3, 1e-2, 0.1, 1], xlab='channel frequency [MHz]',
                 ylab='weight of the channel in the TOA fit (clean channel = 1)', yfmt=lambda v: f'{v:g}'))
for c in range(nch):
    x0, x1 = float(A.X(edges[c])) + 0.4, float(A.X(edges[c + 1])) - 0.4
    y = float(A.Y(w[c]))
    o.append(f'<rect class="box sky" x="{x0:.1f}" y="{y:.1f}" width="{x1 - x0:.1f}" height="{A.y0 + A.h - y:.1f}"/>')
for typ, f0, rc, inr, lab in srcs:
    c = int((f0 - fL)//cw)
    y = float(A.Y(w[c]))
    yl = float(A.Y(15))
    o.append(f'<path class="thin" d="M{A.X(f0):.1f} {yl + 4:.1f} V{y - 2:.1f}"/>')
    anc = 'end' if f0 > 1500 else 'start'
    xt = float(A.X(f0)) + (-6 if anc == 'start' else (6 if anc == 'end' else 0))
    o.append(f'<text class="sub" x="{xt:.1f}" y="{yl:.1f}" text-anchor="{anc}">{lab}: weight {w[c]:.1g}</text>')
svg = '''<figure class="scroll">
<svg class="wide" viewBox="0 0 860 310" role="img" aria-label="Weight of each channel in the optimal TOA fit with the illustrative RFI: the carrier and GNSS channels get almost no weight">
  ''' + '\n  '.join(o) + '''
</svg>
<figcaption><b>Figure 8.3.</b> The built-in notch of optimal weighting, for the illustrative scenario's three continuous sources (expected levels from
their INRs, not a measurement). The fit weights each channel by 1/level², so the carrier's channel gets about 1/1700 of a clean channel's
weight and is in effect left out, which costs only ~0.4 % of the SNR (one channel of 128); with all three sources the loss is ~3 %. This is
why constant-envelope RFI needs no excision on the channel path (<a href="excision.html#results">excision results</a>: SNR 91.9 with the
carrier, against 74 on the full band).</figcaption>
</figure>'''
open('chanweights.svg', 'w').write(svg)
