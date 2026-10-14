import numpy as np
from plotsvg import Ax
# Why times are double precision, from absolute sample indices (conventions.html, section 5).
# a: rounding error of a time value t (half the spacing of representable numbers) in float32 and double, against t, with the
# TOA precisions that matter. b: a sample clock accumulated as t += 1/fs (fs = 800 MHz) in float32 and in double, against the
# exact k/fs; error after k samples (simulated with sequential cumulative sums).
o = []
A = Ax(80, 46, 320, 240, (-3, 5), (-18, 0))
lx = lambda t: A.x0 + (np.log10(t) + 3)/8*A.w
ly = lambda v: A.y0 + A.h - (np.log10(v) + 18)/18*A.h
o.append(f'<path class="axis" d="M{A.x0} {A.y0 + A.h} H{A.x0 + A.w} M{A.x0} {A.y0} V{A.y0 + A.h}"/>')
for t, lab in ((1e-3, '1 ms'), (1e-1, '0.1 s'), (1, '1 s'), (60, '1 min'), (3600, '1 h'), (86400, '1 d')):
    o.append(f'<path class="axis" d="M{lx(t):.1f} {A.y0 + A.h} v5"/><text class="sub" x="{lx(t):.1f}" y="{A.y0 + A.h + 18}" text-anchor="middle">{lab}</text>')
for v, lab in ((1e-18, '1 as'), (1e-15, '1 fs'), (1e-12, '1 ps'), (1e-9, '1 ns'), (1e-6, '1 µs'), (1e-3, '1 ms'), (1, '1 s')):
    o.append(f'<path class="axis" d="M{A.x0} {ly(v):.1f} h-5"/><text class="sub" x="{A.x0 - 8}" y="{ly(v) + 4:.1f}" text-anchor="end">{lab}</text>')
o.append(f'<text class="sub" x="{A.x0 + A.w/2:.1f}" y="{A.y0 + A.h + 34}" text-anchor="middle">time value t</text>')
o.append(f'<text class="sub" x="{A.x0}" y="{A.y0 - 8}">rounding error of t</text>')
tt = np.logspace(-3, 5, 400)
for dt, cls, lab in ((np.float32, 'curve4', 'float32 (single)'), (np.float64, 'curve3', 'double')):
    err = np.maximum(np.array([np.spacing(dt(t)) for t in tt], float)/2, 1e-18)
    o.append(f'<path class="{cls}" d="M{" L".join(f"{lx(t):.1f} {ly(e):.1f}" for t, e in zip(tt, err))}"/>')
    i = 300
    o.append(f'<text class="sub" x="{lx(tt[i]) - 4:.1f}" y="{ly(err[i]) - 8:.1f}" text-anchor="end">{lab}</text>')
for v, lab in ((12e-6, 'J0437 TOA (12 µs)'), (1e-9, 'ns goal (30 cm)')):
    o.append(f'<path class="dash" d="M{A.x0} {ly(v):.1f} H{A.x0 + A.w}"/>')
    if v > 1e-6:
        o.append(f'<text class="sub" x="{A.x0 + 6}" y="{ly(v) - 5:.1f}">{lab}</text>')
    else:
        o.append(f'<text class="sub" x="{A.x0 + A.w - 4}" y="{ly(v) + 15:.1f}" text-anchor="end">{lab}</text>')
e1 = np.spacing(np.float32(1.0))/2
o.append(f'<path class="dot2" d="M{lx(1) - 4:.1f} {ly(e1):.1f} a4 4 0 1 0 8 0 a4 4 0 1 0 -8 0"/>')
o.append(f'<text class="sub" x="{lx(1) + 8:.1f}" y="{ly(e1) + 16:.1f}">60 ns at 1 s</text>')
o.append(f'<text class="ttl" x="80" y="20">a · how finely a number can hold a time</text>')

# ---------- b ----------
fs = 800e6
N = 20_000_000
k = np.unique(np.round(np.logspace(1, np.log10(N), 200)).astype(int))
c32 = np.cumsum(np.full(N, 1/fs, np.float32), dtype=np.float32)
c64 = np.cumsum(np.full(N, 1/fs, np.float64))
exact = np.arange(1, N + 1)/fs
B = Ax(520, 46, 300, 240, (1, np.log10(N)), (-18, 0))
bx = lambda n: B.x0 + (np.log10(n) - 1)/(np.log10(N) - 1)*B.w
o.append(f'<path class="axis" d="M{B.x0} {B.y0 + B.h} H{B.x0 + B.w} M{B.x0} {B.y0} V{B.y0 + B.h}"/>')
for n, lab in ((10, '10'), (1e3, '10³'), (1e5, '10⁵'), (1e7, '10⁷')):
    o.append(f'<path class="axis" d="M{bx(n):.1f} {B.y0 + B.h} v5"/><text class="sub" x="{bx(n):.1f}" y="{B.y0 + B.h + 18}" text-anchor="middle">{lab}</text>')
for v, lab in ((1e-18, '1 as'), (1e-15, '1 fs'), (1e-12, '1 ps'), (1e-9, '1 ns'), (1e-6, '1 µs'), (1e-3, '1 ms'), (1, '1 s')):
    o.append(f'<path class="axis" d="M{B.x0} {ly(v):.1f} h-5"/><text class="sub" x="{B.x0 - 8}" y="{ly(v) + 4:.1f}" text-anchor="end">{lab}</text>')
o.append(f'<text class="sub" x="{B.x0 + B.w/2:.1f}" y="{B.y0 + B.h + 34}" text-anchor="middle">samples k at 800 MHz (10⁷ = 12.5 ms)</text>')
o.append(f'<text class="sub" x="{B.x0}" y="{B.y0 - 8}">error of t after k steps</text>')
for arr, cls, lab in ((c32, 'curve4', 't += 1/fs in float32'), (c64, 'curve3', 't += 1/fs in double')):
    err = np.abs(arr[k - 1].astype(float) - exact[k - 1])
    i0 = np.argmax(err > 1e-18)
    o.append(f'<path class="{cls}" d="M{" L".join(f"{bx(n):.1f} {ly(max(e, 1e-18)):.1f}" for n, e in zip(k[i0:], err[i0:]))}"/>')
    print(lab, 'error at 1e7:', err[np.searchsorted(k, 10_000_000)])
o.append(f'<text class="sub" x="{bx(3e4):.1f}" y="{ly(1e-6) + 2:.1f}" text-anchor="middle">float32 sum</text>')
o.append(f'<text class="sub" x="{bx(3e4):.1f}" y="{ly(1e-14) - 10:.1f}" text-anchor="middle">double sum</text>')
o.append(f'<text class="sub" x="{B.x0 + 8}" y="{ly(1e-17):.1f}">k/fs: exact</text>')
o.append(f'<text class="ttl" x="520" y="20">b · a clock built by adding 1/fs drifts</text>')
svg = '''<figure class="scroll">
<svg class="wide" viewBox="0 0 860 340" role="img" aria-label="Why times are double precision from sample indices: float32 rounds a time to 60 ns at 1 s, and a clock built by adding 1/fs drifts">
  ''' + '\n  '.join(o) + '''
</svg>
<figcaption>Why every time is computed in double precision from the absolute sample index. a: a float32 number can only hold a time to
about 60 ns at 1 s, and to milliseconds after a day; double stays below 10 ps even after a day. b: building a sample clock by adding 1/<var>f</var><sub>s</sub>
step by step lets rounding errors pile up; in float32 the clock is off by about 1 µs after 1 ms and by 1 ms after 12 ms. Computing
<var>k</var>/<var>f</var><sub>s</sub> from the index has no such drift, and it gives the same time whatever block size a stage streams with.</figcaption>
</figure>'''
open('timeprec.svg', 'w').write(svg)
