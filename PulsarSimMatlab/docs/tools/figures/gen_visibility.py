import numpy as np
from plotsvg import Ax
# Visibility (scenario.html): hours per sidereal day a pulsar is above 20 deg elevation, against station latitude, from
# sin(h) = sin(lat) sin(dec) + cos(lat) cos(dec) cos(H). The 20 deg mask is the scenario's visibility rule (max elevation >= 20 deg,
# i.e. latitude within dec +- 70 deg; notes section 11, deep analysis section 17). Declinations from the J2000 names.
psr = [('J0437−4715', -47.25, 'curve2', 2.6), ('Vela', -45.17, 'curve3', 2.6), ('B0329+54', 54.57, 'curve', 2.6),
       ('B0950+08', 7.92, 'thin', 1.3), ('B1749−28', -28.10, 'thin', 1.3), ('B1642−03', -3.28, 'thin', 1.3),
       ('B1929+10', 10.98, 'thin', 1.3)]
mask = np.radians(20)
lat = np.linspace(-89.5, 89.5, 359)


def hours(dec):
    la, de = np.radians(lat), np.radians(dec)
    c = (np.sin(mask) - np.sin(la)*np.sin(de))/(np.cos(la)*np.cos(de))
    return 24*np.degrees(np.arccos(np.clip(c, -1, 1)))/180


o = []
A = Ax(70, 40, 740, 250, (-90, 90), (0, 24))
o.append(A.frame(xticks=[-90, -60, -30, 0, 30, 60, 90], yticks=[0, 6, 12, 18, 24], xlab='station latitude',
                 ylab='hours per day above 20° elevation', xfmt='{:g}°'))
o.append(f'<text class="sub" x="{float(A.X(-75)):.1f}" y="{A.y0 + A.h + 34}" text-anchor="middle">south</text>')
o.append(f'<text class="sub" x="{float(A.X(75)):.1f}" y="{A.y0 + A.h + 34}" text-anchor="middle">north</text>')
for la, lab in ((-30, '30° S'), (52, '52° N (Europe)')):
    x = float(A.X(la))
    o.append(f'<path class="dash" d="M{x:.1f} {A.y0} V{A.y0 + A.h}"/>')
    o.append(f'<text class="sub" x="{x + 4:.1f}" y="{float(A.Y(1.2)) if la > 0 else A.y0 + 10:.1f}">{lab}</text>')
TC = {'curve': 'var(--accent)', 'curve2': 'var(--warn)', 'curve3': 'var(--why)', 'thin': 'var(--fg-muted)'}
for name, dec, cls, w in psr:
    o.append(f'<path class="{cls}" stroke-width="{w}" d="{A.path(lat, hours(dec))}"/>')
o.append(f'<text class="sub" x="{float(A.X(-26)):.1f}" y="{float(A.Y(19.5)):.1f}"><tspan style="fill:var(--warn)">J0437−4715</tspan> and <tspan style="fill:var(--why)">Vela</tspan></text>')
o.append(f'<text class="sub" x="{float(A.X(-26)):.1f}" y="{float(A.Y(19.5)) + 14:.1f}">(almost the same declination)</text>')
o.append(f'<text class="sub" x="{float(A.X(36)):.1f}" y="{float(A.Y(16.5)):.1f}" text-anchor="end" style="fill:var(--accent)">B0329+54</text>')
o.append(f'<text class="sub" x="{float(A.X(89)):.1f}" y="{float(A.Y(16)):.1f}" text-anchor="end">grey: B0950+08, B1749−28,</text>')
o.append(f'<text class="sub" x="{float(A.X(89)):.1f}" y="{float(A.Y(16)) + 14:.1f}" text-anchor="end">B1642−03, B1929+10</text>')
svg = '''<figure class="scroll">
<svg class="wide" viewBox="0 0 860 340" role="img" aria-label="Hours per day each candidate pulsar is above 20 degrees elevation, against the latitude of the station">
  ''' + '\n  '.join(o) + '''
</svg>
<figcaption>Where you are decides which pulsars you can use: hours per day each candidate is more than 20° above the horizon,
against the station's latitude (the scenario's visibility rule). J0437 and Vela, the two best sources, are southern: from 30° S
each is up about 12 h a day, from 52° N never. B0329+54 is the only fast northern source. A TOA that needs 21 hours on source
therefore takes about two days of calendar time for J0437 from 30° S (<a href="limits.html">limitations</a>).</figcaption>
</figure>'''
open('visibility.svg', 'w').write(svg)
