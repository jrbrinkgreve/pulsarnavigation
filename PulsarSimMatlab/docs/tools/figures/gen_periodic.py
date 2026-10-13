import numpy as np
from plotsvg import Ax
# periodicRFI in one picture (schematic, with the scenario radar's PRF of 373 Hz over 50 ms): weak sidelobe pulses,
# the ~42 % that detectRFI finds one by one (testPeriodicRFI's rotating-radar case), the period fitted from the gaps,
# and the predicted pulses that are blanked unconditionally.
P = 1e3/373.0                                        # ms
t0 = 0.9
tp = t0 + P*np.arange(19)
tp = tp[tp < 50]
rng = np.random.default_rng(16)
found = rng.random(tp.size) < 0.42
print(found.sum(), 'of', tp.size, 'found')
o = []
A = Ax(150, 40, 680, 170, (0, 50), (0, 1))
o.append(A.frame(xticks=[0, 10, 20, 30, 40, 50], xlab='time [ms]', left=False))
rows = [('true radar pulses', 'sidelobes, −35 dB'), ('found one by one', 'detectRFI: ~42 %'), ('predicted, blanked', 'periodicRFI: all')]
for i, (lab, sub) in enumerate(rows):
    y = A.y0 + 22 + i*52
    o.append(f'<text class="lbl" x="{A.x0 - 10}" y="{y}" text-anchor="end">{lab}</text>')
    o.append(f'<text class="sub" x="{A.x0 - 10}" y="{y + 14}" text-anchor="end">{sub}</text>')
    o.append(f'<path class="grid" d="M{A.x0} {y + 4} H{A.x0 + A.w}"/>')
    for k, t in enumerate(tp):
        x = float(A.X(t))
        if i == 0:
            o.append(f'<rect class="fillfg" x="{x - 1.5:.1f}" y="{y - 10}" width="3" height="20"/>')
        elif i == 1 and found[k]:
            o.append(f'<rect class="bad" x="{x - 1.5:.1f}" y="{y - 10}" width="3" height="20"/>')
        elif i == 2:
            o.append(f'<rect class="fillwarn" x="{x - 4:.1f}" y="{y - 12}" width="8" height="24"/>')
            o.append(f'<rect class="bad" x="{x - 1.5:.1f}" y="{y - 10}" width="3" height="20"/>')
# gaps between found pulses, as multiples of P
idx = np.flatnonzero(found)
yb = A.y0 + 22 + 52 + 22
for a, b in zip(idx[:-1], idx[1:]):
    xa, xb = float(A.X(tp[a])), float(A.X(tp[b]))
    o.append(f'<path class="thin" d="M{xa:.1f} {yb} H{xb:.1f}"/>')
    o.append(f'<text class="sub" x="{(xa + xb)/2:.1f}" y="{yb + 12}" text-anchor="middle">{b - a}P</text>')
svg = '''<figure class="scroll">
<svg class="wide" viewBox="0 0 860 260" role="img" aria-label="periodicRFI: a few detected radar pulses give the period, and every predicted pulse is blanked">
  ''' + '\n  '.join(o) + '''
</svg>
<figcaption>The idea of <code>periodicRFI</code> (schematic: the scenario radar's PRF, 373 Hz, over 50 ms; which pulses are found is random). Weak
sidelobe pulses are found only now and then. The gaps between the found ones are whole multiples of the period <var>P</var>, so dividing
the gaps by 1–8 gives candidate periods; a significance test picks the one the found pulses fit far better than chance. From
<var>P</var> and one pulse time every pulse is known, and each predicted pulse is blanked with a margin (orange), found or not.</figcaption>
</figure>'''
open('periodic.svg', 'w').write(svg)
