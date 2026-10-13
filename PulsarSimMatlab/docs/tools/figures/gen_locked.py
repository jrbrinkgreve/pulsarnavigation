import numpy as np
from plotsvg import Ax, pow10
# runLockedRadar (B5b / B6, 8 Oct; excision.html section 6 table): TOA bias from a radar locked to the pulsar
# (400 Hz = 4 f0, worst pulse phase), against snrDB, for each excision setting. |bias| in us; the -20 dB value of
# the periodic-mask case is 0 (not drawn on a log axis).
snr = np.array([-20, -30, -40, -54])
cases = [('excision off: locks on the radar', None, 'curve4'),
         ('instant edges, excision on', [5.5, 64, 662, 1453], 'curve4'),
         ('excision on, no guard', [0.41, 3.9, 31, 147], 'curve2'),
         ('frequency guard 3', [0.02, 0.18, 1.8, 41], 'curve2'),
         ('frequency guard 5', [0.003, 0.025, 0.25, 6.4], 'curve'),
         ('frequency guard 7 (default)', [0.001, 0.008, 0.08, 2.0], 'curve'),
         ('guard 7 + periodic mask', [None, 0.002, 0.016, 0.41], 'curve3')]
o = []
A = Ax(90, 40, 470, 250, (-57, -17), (5e-4, 5e3), logy=True)
o.append(A.frame(xticks=[-54, -50, -40, -30, -20], yticks=[1e-3, 1e-2, 0.1, 1, 10, 100, 1e3],
                 xlab='snrDB', ylab='|TOA bias| from the locked radar [µs]', yfmt=lambda v: f'{v:g}'))
o.append(f'<path class="dash" d="M{A.x0} {A.Y(1):.1f} H{A.x0 + A.w}"/>')
o.append(f'<text class="sub" x="{A.x0 + A.w - 4}" y="{A.Y(1) - 5:.1f}" text-anchor="end">1 µs</text>')
o.append(f'<path class="thin" d="M{A.X(-54):.1f} {A.y0} V{A.y0 + A.h}"/>')
o.append(f'<text class="sub" x="{A.X(-54) + 4:.1f}" y="{A.y0 + 12}">design target</text>')
labels = []
for name, vals, cls in cases:
    if vals is None:
        continue
    pts = [(s, v) for s, v in zip(snr, vals) if v]
    d = ' L'.join(f'{A.X(s):.1f} {A.Y(v):.1f}' for s, v in pts)
    dash = ' stroke-dasharray="5 3"' if name.startswith('instant') else ''
    o.append(f'<path class="{cls}"{dash} d="M{d}"/>')
    for s, v in pts:
        o.append(f'<path class="dot" d="M{A.X(s) - 3:.1f} {A.Y(v):.1f} a3 3 0 1 0 6 0 a3 3 0 1 0 -6 0"/>')
    labels.append((float(A.Y(pts[0][1])), name, float(A.X(pts[0][0]))))
# labels at the right end of each line (the -20 dB side, or -30 dB), stacked to avoid overlap
labels.sort()
y_prev = -1e9
for y0, name, x0 in labels:
    y = max(y0, y_prev + 14); y_prev = y
    o.append(f'<path class="thin" d="M{x0 + 4:.1f} {y0:.1f} L{A.x0 + A.w + 6} {y:.1f}"/>')
    o.append(f'<text class="sub" x="{A.x0 + A.w + 10}" y="{y + 4:.1f}">{name}</text>')
o.append(f'<text class="sub" x="{A.x0 + A.w + 8}" y="{A.y0 - 10}">excision off: 172–2288 µs</text>')
o.append(f'<text class="sub" x="{A.x0 + A.w + 8}" y="{A.y0 + 4}">(locks on the radar)</text>')
svg = '''<figure class="scroll">
<svg class="wide" viewBox="0 0 860 330" role="img" aria-label="TOA bias from a radar locked to the pulsar against SNR, for each frequency-guard width">
  ''' + '\n  '.join(o) + '''
</svg>
<figcaption>A radar locked to the pulsar's period (PRF 400 Hz = 4 <var>f</var><sub>0</sub>, worst pulse phase; <code>runLockedRadar.m</code>,
8 Oct): the TOA bias left after excision, from the table in <a href="#locked">section 6</a>. A fixed leftover in the fold grows as 1/ρ when the pulsar
gets weaker, so the lines rise roughly ten-fold per 10 dB. Each two extra guard channels cut it 3–4×; guard 7 (the default) keeps it below
1 µs down to −51 dB, and the periodic mask (B6) to −58 dB. With instant radar edges (dashed) the guard does not help.</figcaption>
</figure>'''
open('locked.svg', 'w').write(svg)
