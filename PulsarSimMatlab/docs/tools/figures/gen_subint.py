import numpy as np
from plotsvg import Ax, pow10
# Choosing the sub-int length (fold.html): turns per sub-int for a target SNR, N = (SNR_target / SNR_turn)^2,
# SNR_turn = 373*rho (flat band, incl. the baseline-estimation loss), for the simulated 10 ms pulsar; at least 1 turn.
# Marked: the levels of the design point of 23 Nov 2026 (J0437-4715 -44.2 dB, B0329+54 -38.0 dB), -20 dB and -5 dB.
dB = np.linspace(-60, 0, 301)
rho = 10**(dB/10)
o = []
A = Ax(90, 40, 560, 250, (-60, 0), (0.5, 3e8), logy=True)
o.append(A.frame(xticks=[-60, -50, -40, -30, -20, -10, -5, 0], yticks=[1, 1e2, 1e4, 1e6, 1e8], xlab='snrDB',
                 ylab='turns per sub-int', yfmt=lambda v: f'{v:g}' if v < 1e4 else pow10(v)))
for snr, cls in ((7, 'curve'), (6, 'curve2')):
    N = np.maximum((snr/(373*rho))**2, 1)
    o.append(f'<path class="{cls}" d="{A.path(dB, N)}"/>')
for turns, lab in ((100, '1 s'), (6000, '1 min'), (360000, '1 h'), (8640000, '1 day')):
    y = float(A.Y(turns))
    o.append(f'<path class="grid" d="M{A.x0} {y:.1f} H{A.x0 + A.w}"/>')
    o.append(f'<text class="sub" x="{A.x0 + A.w + 6}" y="{y + 4:.1f}">{lab}</text>')
pts = [(-44.2, 'J0437−4715 (−44.2 dB): 2.4·10⁵ turns ≈ 41 min'), (-38, 'B0329+54 (−38.0 dB): 1.4·10⁴ turns ≈ 2.3 min'),
       (-20, '−20 dB: 3.5 turns'), (-5, '−5 dB: 1 turn')]
for d, lab in pts:
    N = max((7/(373*10**(d/10)))**2, 1)
    x, y = float(A.X(d)), float(A.Y(N))
    o.append(f'<path class="dot" d="M{x - 3.5:.1f} {y:.1f} a3.5 3.5 0 1 0 7 0 a3.5 3.5 0 1 0 -7 0"/>')
    o.append(f'<text class="sub" x="{x + 8:.1f}" y="{y - 6:.1f}">{lab}</text>')
    print(d, N, N*0.01/3600, 'h')
lx, ly = float(A.X(-14)), float(A.Y(3e7))
o.append(f'<path class="curve" d="M{lx} {ly} h18"/><text class="sub" x="{lx + 24}" y="{ly + 4}">SNR 7 per sub-int</text>')
o.append(f'<path class="curve2" d="M{lx} {ly + 16} h18"/><text class="sub" x="{lx + 24}" y="{ly + 20}">SNR 6</text>')
svg = '''<figure class="scroll">
<svg class="wide" viewBox="0 0 860 330" role="img" aria-label="Turns per sub-int needed for SNR 6 or 7 against snrDB, with the time for the 10 ms pulsar">
  ''' + '\n  '.join(o) + '''
</svg>
<figcaption><b>Figure 7.4.</b> The number of turns for an SNR of 6–7 per sub-int, <var>N</var> = (SNR/373·ρ)², for the simulated 10 ms pulsar (grid lines:
the time this takes). At −5 dB every turn is enough, which is why <code>main.m</code> folds one turn per sub-int. At the levels of the array's design point
(dots) one sub-int of this test pulsar, with one polarization in 400 MHz, takes about 2.3 minutes for B0329+54's −38.0 dB and 41 minutes
for J0437−4715's −44.2 dB.</figcaption>
</figure>'''
open('subint.svg', 'w').write(svg)
