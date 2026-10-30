import numpy as np
from plotsvg import Ax
# The f_out sweep of 9 Oct (logs/2026-10-09_fout-study.md, -5 dB table, matched nBin): channel path,
# 5 noise seeds, the same dedispersed channels for every f_out. SNR per sub-int and honesty of the
# error bars (rms error / stated error) against the width of a detected time bin.
binus = np.array([500, 200, 100, 33, 10, 3.4, 0.96])
snr = np.array([81.9, 91.3, 92.0, 91.5, 91.9, 92.0, 92.1])
ratio = np.array([1.68, 2.04, 1.08, 1.03, 1.04, 1.04, 1.04])
lo, hi = np.log10(0.5), np.log10(1000)
def lxf(A):
    return lambda v: A.x0 + (np.log10(v) - lo)/(hi - lo)*A.w
o = []
for k, (y, ylim, yt, data, ylab, ttl) in enumerate([
        (snr, (78, 95), [80, 85, 90, 95], snr, 'SNR per sub-int', 'a · sensitivity: flat once the bins resolve the pulse'),
        (ratio, (0.8, 2.2), [1, 1.5, 2], ratio, 'rms error ÷ stated error', 'b · error bars: honest only for fine bins')]):
    A = Ax(80 + 420*k, 46, 330, 190, (0, 1), ylim)
    lx = lxf(A)
    o.append(A.frame(yticks=yt, ylab=ylab))
    for v, lab in [(1, '1'), (10, '10'), (100, '100'), (1000, '1000')]:
        o.append(f'<path class="axis" d="M{lx(v):.1f} {A.y0+A.h} v5"/><text class="sub" x="{lx(v):.1f}" y="{A.y0+A.h+18}" text-anchor="middle">{lab}</text>')
    o.append(f'<text class="sub" x="{A.x0 + A.w/2:.1f}" y="{A.y0+A.h+34}" text-anchor="middle">width of a detected time bin [µs]</text>')
    o.append(f'<text class="ttl" x="{A.x0}" y="20">{ttl}</text>')
    # zones: bins coarser than W/5 (pulse FWHM 0.5 ms), and main's choice
    o.append(f'<rect class="fillwarn" x="{lx(100):.1f}" y="{A.y0}" width="{lx(1000) - lx(100):.1f}" height="{A.h}"/>')
    o.append(f'<path class="dash" d="M{lx(100):.1f} {A.y0} V{A.y0 + A.h}"/>')
    o.append(f'<text class="sub" x="{lx(105):.1f}" y="{A.y0 + A.h - 6}">coarser than W/5</text>')
    o.append(f'<path class="dash" d="M{lx(4.88):.1f} {A.y0} V{A.y0 + A.h}"/>')
    pts = ' L'.join(f'{lx(a):.1f} {float(A.Y(b)):.1f}' for a, b in zip(binus, data))
    o.append(f'<path class="curve" d="M{pts}"/>')
    for a, b in zip(binus, data):
        cls = 'dot2' if a == 0.96 else 'dot'
        o.append(f'<path class="{cls}" d="M{lx(a) - 3.5:.1f} {A.Y(b):.1f} a3.5 3.5 0 1 0 7 0 a3.5 3.5 0 1 0 -7 0"/>')
    if k == 0:
        o.append(f'<text class="sub" x="{lx(0.96) + 6:.1f}" y="{A.Y(92.1) + 16:.1f}">main.m: 0.96 µs</text>')
        o.append(f'<text class="sub" x="{lx(500) - 6:.1f}" y="{A.Y(81.9) + 4:.1f}" text-anchor="end">−11 %</text>')
    else:
        o.append(f'<text class="sub" x="{lx(4.88) + 4:.1f}" y="{A.Y(2.0):.1f}">≤ 4.9 µs: every one of</text>')
        o.append(f'<text class="sub" x="{lx(4.88) + 4:.1f}" y="{A.Y(2.0) + 13:.1f}">2048 phase bins filled</text>')
        o.append(f'<text class="sub" x="{lx(4.88) + 4:.1f}" y="{A.Y(2.0) + 26:.1f}">in one turn</text>')
        o.append(f'<text class="sub" x="{lx(200) + 8:.1f}" y="{A.Y(2.04) + 4:.1f}">2× too small</text>')
svg = '''<figure class="scroll">
<svg class="wide" viewBox="0 0 860 300" role="img" aria-label="Measured SNR and error-bar honesty against detected bin width: flat for fine bins, worse for bins coarser than a fifth of the pulse">
  ''' + '\n  '.join(o) + '''
</svg>
<figcaption><b>Figure 6.2.</b> The <var>f</var><sub>out</sub> sweep (channel path, −5 dB, 5 noise seeds, the same dedispersed data for every
<var>f</var><sub>out</sub>; <code>logs/2026-10-09_fout-study.md</code>). (a) Once a bin is finer than about W/5 = 100 µs, finer bins add no
sensitivity. (b) Coarse bins also smear the pulse, which the template does not know about, so the error bars become too small. The second
dashed line is the coverage limit for <code>NBin</code> = 2048 with one turn per sub-int. <code>main.m</code> uses 0.96 µs (orange dot).</figcaption>
</figure>'''
open('fout.svg', 'w').write(svg)
