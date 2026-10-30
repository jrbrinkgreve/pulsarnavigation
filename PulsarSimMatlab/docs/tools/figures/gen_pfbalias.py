import numpy as np
from plotsvg import Ax
from pfb import design, response, db, envelope, fs, dF
# Sliding the frame by D samples samples each channel at fs/D, so its response repeats every fs/D.
# With D = 192 the copies stay out of the useful band; with D = 256 (same prototype) they reach 1 MHz into it.
h, n, w, I = design(192)
stop = I['fsOut'] - dF/2
f = np.arange(-8, 8.0001, 0.002)

o = []
panels = [(192, 58, 'a · Slide by D = 192: copies every 4.17 MHz'),
          (256, 330, 'b · Slide by D = 256, the same prototype: copies every 3.125 MHz')]
for D, y0, title in panels:
    rate = fs/D
    A = Ax(80, y0 + 30, 750, 170, (-8, 8), (-100, 6))
    o.append(f'<text class="ttl" x="80" y="{y0}">{title}</text>')
    o.append(A.frame(xticks=[-8, -6, -4, -2, 0, 2, 4, 6, 8], yticks=[0, -20, -40, -60, -80, -100],
                     xlab='frequency offset from the channel centre [MHz]', ylab='gain [dB]'))
    o.append(f'<rect class="fillwhy" x="{A.X(-dF/2):.1f}" y="{A.y0}" width="{A.X(dF/2) - A.X(-dF/2):.1f}" height="{A.h}"/>')
    worst = -200.0
    for m in (-2, -1, 1, 2):
        Hm = db(response(h, n, f - m*rate))
        fe, He = envelope(f, Hm, 0.07)
        o.append(f'<path class="curve2" opacity="0.75" d="{A.path(fe, He)}"/>')
        worst = max(worst, Hm[np.abs(f) <= dF/2].max())
    fe, He = envelope(f, db(response(h, n, f)), 0.07)
    o.append(f'<path class="curve" d="{A.path(fe, He)}"/>')
    o.append(f'<path class="thin" d="M{A.X(-8):.1f} {A.Y(-80):.1f} H{A.X(8):.1f}" stroke-dasharray="2 3"/>')
    print(D, 'largest copy inside the useful band [dB]', worst)
    if D == 192:
        o.append(f'<text class="sub" x="{A.X(0):.1f}" y="{A.Y(-45):.1f}" text-anchor="middle">copies below</text>')
        o.append(f'<text class="sub" x="{A.X(0):.1f}" y="{A.Y(-45) + 14:.1f}" text-anchor="middle">−80 dB here</text>')
        o.append(f'<text class="sub" x="830" y="{y0}" text-anchor="end">the nearest copy (+4.17 MHz) reaches −80 dB exactly at +1.5625 MHz</text>')
    else:
        bad = 3.125 - stop                    # where the copy rises above -80 dB inside the band
        for a, b in ((-dF/2, -bad), (bad, dF/2)):
            o.append(f'<rect class="fillwarn" x="{A.X(a):.1f}" y="{A.y0}" width="{A.X(b) - A.X(a):.1f}" height="{A.h}" opacity="0.45"/>')
        o.append(f'<text class="sub" x="830" y="{y0}" text-anchor="end">the copies reach 1 MHz into the useful band (shaded)</text>')

svg = '''<figure class="scroll">
<svg class="wide" viewBox="0 0 860 590" role="img" aria-label="A channel response and its copies after decimation: with D = 192 the copies stay below -80 dB in the useful band; with D = 256 they reach into it">
  ''' + '\n  '.join(o) + '''
</svg>
<figcaption><b>Figure 5a.4.</b> Sliding by <var>D</var> samples gives the channel response copies every
<var>f</var><sub>s</sub>/<var>D</var> (orange). (a) With <var>D</var> = 192, the copies stay below −80 dB inside the useful band
(green). (b) With <var>D</var> = 256 and the same prototype, they reach 1 MHz into it (shaded).</figcaption>
</figure>'''
open('pfbalias.svg', 'w').write(svg)
