import numpy as np
from plotsvg import Ax
# Which parts of the 0.1 s main.m file are fully supported after dedispersion, and which turns are complete,
# for both front ends (notes section 4). Turn k = [10k, 10k + 10) ms, pulse centres at 5, 15, ..., 95 ms.
# Full band: the last 6.302 ms (+ the 2.5 us guard) lack part of their sweep. Channel path: also the first
# 0.32 ms (the per-channel guard, 20 / 62.5 kHz).
L, T = 100.0, 10.0
sweep = 6.302442708
rows = [('full band', 0.0025, sweep + 0.0025, 9, 'turns 0–8 complete: 9 TOAs'),
        ('channel path', 0.32, sweep + 0.32, 8, 'turns 1–8 complete: 8 TOAs')]
o = []
A = Ax(130, 40, 690, 150, (0, 100), (0, 1))
o.append(A.frame(xticks=[0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100], xlab='time in the file [ms]', left=False))
for k in range(10):
    x0, x1 = A.X(k*T), A.X((k + 1)*T)
    o.append(f'<text class="sub" x="{(x0 + x1)/2:.1f}" y="{A.y0 - 6}" text-anchor="middle">turn {k}</text>')
    o.append(f'<path class="grid" d="M{x0:.1f} {A.y0} V{A.y0 + A.h}"/>')
    o.append(f'<path class="dot" d="M{A.X(k*T + 5) - 3:.1f} {A.y0 + 10:.1f} a3 3 0 1 0 6 0 a3 3 0 1 0 -6 0"/>')
o.append(f'<text class="sub" x="{A.x0 - 8}" y="{A.y0 + 14}" text-anchor="end">pulse centres</text>')
for i, (name, start, end, nT, lab) in enumerate(rows):
    y = A.y0 + 32 + i*58
    o.append(f'<text class="lbl" x="{A.x0 - 8}" y="{y + 15}" text-anchor="end">{name}</text>')
    o.append(f'<rect class="box proc" x="{A.X(start):.1f}" y="{y}" width="{A.X(L - end) - A.X(start):.1f}" height="22"/>')
    o.append(f'<rect class="fillwarn" x="{A.x0}" y="{y}" width="{max(A.X(start) - A.x0, 1.5):.1f}" height="22"/>')
    o.append(f'<rect class="fillwarn" x="{A.X(L - end):.1f}" y="{y}" width="{A.X(L) - A.X(L - end):.1f}" height="22"/>')
    first = 0 if name == 'full band' else 1
    for k in range(first, first + nT):
        o.append(f'<path class="line" d="M{A.X(k*T) + 2:.1f} {y + 30} H{A.X((k + 1)*T) - 2:.1f}"/>')
    o.append(f'<text class="sub" x="{A.X(L - end) - 4:.1f}" y="{y + 44}" text-anchor="end">{lab}</text>')
ly = A.y0 + A.h + 52
o.append(f'<rect class="box proc" x="{A.x0}" y="{ly-9}" width="18" height="10"/><text class="sub" x="{A.x0+24}" y="{ly}">fully supported</text>')
o.append(f'<rect class="fillwarn" x="{A.x0+130}" y="{ly-9}" width="18" height="10"/><text class="sub" x="{A.x0+154}" y="{ly}">not: the last 6.3 ms (sweep incomplete); on the channel path also the first 0.32 ms (guard)</text>')
o.append(f'<path class="line" d="M{A.x0} {ly+14} h18"/><text class="sub" x="{A.x0+24}" y="{ly+18}">complete turns: one TOA each</text>')
svg = '''<figure class="scroll">
<svg class="wide" viewBox="0 0 860 268" role="img" aria-label="Fully supported parts of the 0.1 s file after dedispersion, for the full band and the channel path, and the complete turns they contain">
  ''' + '\n  '.join(o) + '''
</svg>
<figcaption>The 0.1 s file of <code>main.m</code> after dedispersion (notes §4). Green: output samples that saw their complete dispersion
sweep (<code>fullySupported</code>). Orange: samples that did not. At the end the low frequencies of the last 6.3 ms were delayed past the
end of the file in stage 2; at the start the filter's guard reaches before the first sample (2.5 µs on the full band, too narrow to see; 0.32 ms per
channel). The fold uses only green data and the TOA fit only complete turns (lines), so the full band gives 9 TOAs
and the channel path 8: the file starts exactly on a turn boundary, so the 0.32 ms channel guard leaves turn 0 with only 96.8 % of its
phase bins covered, below <code>MinCoverage</code> = 1.</figcaption>
</figure>'''
open('support.svg', 'w').write(svg)
