import numpy as np
# How big each effect that the simulation leaves out is (limits.html), as a timing error, on a log scale, against the J0437 TOA
# precision at SNR 7 (12 us) and the ns goal (1 ns = 30 cm). Sizes from the limits page: Earth's orbital motion v/c ~ 1e-4
# (100 us after 1 s, ~7 s after 20 h uncorrected); J0437's binary orbit ~3 s within one 20 h fold; clock rate off by 1e-9 over
# 20 h: 72 us; station position error 4-100 km: delta_r/c; DM change of 1e-3 pc/cm3 at 1.4 GHz: 2.1 us; residual of an exactly
# locked radar in the main beam at -54 dB with guard and periodic mask: 0.41 us; ionosphere 3-70 ns; scattering 0.4 ns at DM 5,
# about a microsecond at DM ~100.
c = 299792458.0
rows = [('Earth\'s orbital motion (uncorrected)', 100e-6, 7.2, 'after 1 s … after 20 h'),
        ('J0437\'s binary orbit, within one fold', 3.0, None, '~3 s in 20 h'),
        ('clock rate off by 10⁻⁹, over 20 h', 72e-6, None, '72 µs'),
        ('station position error 4–100 km', 4e3/c, 100e3/c, '13–330 µs'),
        ('DM change of 0.001 pc cm⁻³', 2.1e-6, None, '2.1 µs at 1.4 GHz'),
        ('locked radar left after excision (main beam)', 0.41e-6, None, '0.41 µs at −54 dB'),
        ('ionosphere', 3e-9, 70e-9, '3–70 ns'),
        ('scattering', 0.4e-9, 1e-6, '0.4 ns (DM 5) … µs (DM ~100)')]
x0, x1, lo, hi = 300, 820, -10, 1
X = lambda t: x0 + (np.log10(t) - lo)/(hi - lo)*(x1 - x0)
y0, dy = 50, 30
o = []
yb = y0 + len(rows)*dy
toa, ns = 12e-6, 1e-9
o.append(f'<rect class="fillwarn" x="{X(ns):.1f}" y="{y0 - 14}" width="{X(toa) - X(ns):.1f}" height="{yb - y0 + 4}"/>')
o.append(f'<rect class="fillfg" x="{X(toa):.1f}" y="{y0 - 14}" width="{x1 - X(toa):.1f}" height="{yb - y0 + 4}" opacity="0.18"/>')
for t, lab in ((toa, 'J0437 TOA at SNR 7: 12 µs'), (ns, '1 ns = 30 cm')):
    o.append(f'<path class="dash" d="M{X(t):.1f} {y0 - 22} V{y0 - 14}"/>')
    o.append(f'<text class="sub" x="{X(t):.1f}" y="{y0 - 26}" text-anchor="middle">{lab}</text>')
for i, (name, a, b, val) in enumerate(rows):
    y = y0 + i*dy
    o.append(f'<text class="sub" x="10" y="{y + 4}">{name}</text>')
    big = (b or a) > toa
    cls = 'bad' if big else ('dot2' if (b or a) > ns else 'dot')
    lcls = 'curve4' if big else ('curve2' if (b or a) > ns else 'curve')
    if b:
        o.append(f'<path class="{lcls}" style="stroke-width:5" stroke-linecap="round" d="M{X(a):.1f} {y} H{X(b):.1f}"/>')
        xe = X(b)
    else:
        o.append(f'<circle cx="{X(a):.1f}" cy="{y}" r="5" class="{cls}"/>')
        xe = X(a)
    if xe > x1 - 150:
        o.append(f'<text class="sub" x="{(X(a) if b else xe) - 10:.1f}" y="{y + 4}" text-anchor="end">{val}</text>')
    else:
        o.append(f'<text class="sub" x="{xe + 10:.1f}" y="{y + 4}">{val}</text>')
o.append(f'<path class="axis" d="M{x0} {yb} H{x1}"/>')
for t, lab in ((1e-9, '1 ns'), (1e-6, '1 µs'), (1e-3, '1 ms'), (1, '1 s')):
    o.append(f'<path class="axis" d="M{X(t):.1f} {yb} v5"/><text class="sub" x="{X(t):.1f}" y="{yb + 18}" text-anchor="middle">{lab}</text>')
for e in range(lo, hi + 1):
    o.append(f'<path class="axis" d="M{X(10.0**e):.1f} {yb} v3"/>')
o.append(f'<text class="sub" x="{(x0 + x1)/2:.1f}" y="{yb + 34}" text-anchor="middle">timing error if left out</text>')
o.append(f'<text class="sub" x="10" y="{yb + 34}">red: larger than a J0437 TOA · orange: between 1 ns and a TOA</text>')
svg = f'''<figure class="scroll">
<svg class="wide" viewBox="0 0 860 {yb + 44}" role="img" aria-label="Size of each effect the simulation leaves out, as a timing error on a log scale, against the J0437 TOA precision and the nanosecond goal">
  ''' + '\n  '.join(o) + '''
</svg>
<figcaption>How big each effect left out of the simulation is, as a timing error (log scale; numbers from the sections below). Right of the
12 µs line (grey) an effect is larger than a J0437 TOA and must be modelled before folding for hours: Earth's motion, the pulsar's binary
orbit, the clock and the station position, all part of block D. Between 1 ns and 12 µs (orange) an effect does not change today's TOAs but
matters for nanosecond navigation: DM changes, the ionosphere, scattering at high DM and the residual of a radar locked to the pulsar.</figcaption>
</figure>'''
open('effects.svg', 'w').write(svg)
