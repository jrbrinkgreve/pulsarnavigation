import numpy as np
from plotsvg import Ax
# The realistic RFI list (pipelineParams.m, rfiRealistic) sketched as a dynamic spectrum: where each
# source sits in frequency and when it is on. Widths: main lobe of a BPSK signal = +-chip rate, the
# occupied band of a 'noise' source; radar pulses (2 us) and the impulse (200 ns) drawn wider so they
# can be seen. Colour = how the pipeline deals with it.
K, DM, fH = 4.148808e3, 5.0, 1600.0
cont = [  # (centre MHz, half width MHz, label)
    (1575.42, 1.023, 'GNSS L1 band'), (1561.098, 2.046, 'BeiDou B1I'),
    (1278.75, 5.115, 'Galileo E6'), (1268.52, 10.23, 'BeiDou B3I'), (1246.0, 5.11, 'GLONASS G2'),
    (1227.60, 10.23, 'GPS L2'),
    (1459.5, 6.75, 'LTE'), (1472.0, 4.5, 'LTE'), (1484.5, 6.75, 'LTE'),
    (1528.0, 2.0, 'Inmarsat'), (1541.4, 3.2, 'Inmarsat'), (1550.5, 1.4, 'Inmarsat')]
cw = 1388.9
radars = [(1332.9, 973.0, 0.0), (1257.3, 351.0, 0.7)]   # MHz, PRF, StartTime [ms]
tImp = 13.2
o = []
A = Ax(80, 70, 540, 300, (0, 20), (1200, 1600))
o.append(A.frame(xticks=[0, 5, 10, 15, 20], yticks=[1200, 1300, 1400, 1500, 1600], xlab='time [ms]', ylab='frequency [MHz]'))
for fc, hw, _ in cont:
    o.append(f'<rect class="fillwarn" x="{A.x0}" y="{A.Y(fc + hw):.1f}" width="{A.w}" height="{max(A.Y(fc - hw) - A.Y(fc + hw), 1.5):.1f}"/>')
o.append(f'<path class="curve2" stroke-width="1.5" d="M{A.x0} {A.Y(cw):.1f} H{A.x0 + A.w}"/>')
for f, prf, t0 in radars:
    for tp in np.arange(t0, 20, 1e3/prf):
        x = float(A.X(tp))
        o.append(f'<rect class="bad" x="{x - 1:.1f}" y="{A.Y(f + 4):.1f}" width="2.5" height="{A.Y(f - 4) - A.Y(f + 4):.1f}"/>')
o.append(f'<path class="curve4" stroke-width="2" d="M{A.X(tImp):.1f} {A.y0} V{A.y0 + A.h}"/>')
# the pulsar, for orientation
ff = np.linspace(1200, 1600, 200)
for c in (5.0, 15.0):
    tt = c + (K*DM/ff**2 - K*DM/fH**2)*1e3
    ok = tt <= 20
    o.append(f'<path class="dash" d="{A.path(tt[ok], ff[ok])}"/>')
o.append(f'<text class="sub" x="{A.X(6.6):.1f}" y="{A.Y(1590):.1f}">pulsar sweep, far below the noise</text>')
# labels on the right
lx = A.x0 + A.w + 12
labels = [(1575.42, 1582, 'GNSS L1 band (noise)'), (1561.1, 1566, 'BeiDou B1I (bpsk)'),
          (1540, 1543, 'Inmarsat, 3 blocks (noise)'), (1472, 1472, 'LTE band 32, 3 operators (noise)'),
          (cw, 1389, 'local carrier (cw)'), (1332.9, 1333, 'radar 1, beam passage (pulsed)'),
          (1268, 1285, 'Galileo E6, BeiDou B3I (bpsk)'), (1257.3, 1258, 'radar 2, sidelobes (pulsed)'),
          (1236, 1231, 'GLONASS G2, GPS L2 (bpsk)')]
for fsrc, flab, txt in labels:
    o.append(f'<path class="thin" d="M{A.x0 + A.w + 2} {A.Y(fsrc):.1f} L{lx - 3} {A.Y(flab):.1f}"/>')
    o.append(f'<text class="sub" x="{lx}" y="{A.Y(flab) + 4:.1f}">{txt}</text>')
o.append(f'<text class="sub" x="{A.X(tImp) + 5:.1f}" y="{A.Y(1415):.1f}">impulse</text>')
o.append(f'<text class="sub" x="{A.X(tImp) + 5:.1f}" y="{A.Y(1415) + 13:.1f}">(all frequencies)</text>')
# legend
o.append(f'<rect class="fillwarn" x="80" y="16" width="22" height="10"/>')
o.append('<text class="sub" x="108" y="25">continuous: always on. Not blanked; the TOA fit gives these channels less weight.</text>')
o.append(f'<rect class="bad" x="80" y="31" width="22" height="10"/>')
o.append('<text class="sub" x="108" y="40">intermittent: short and strong. Found and blanked by excision.</text>')
svg = '''<figure class="scroll">
<svg class="wide" viewBox="0 0 860 412" role="img" aria-label="The realistic RFI scenario as a time-frequency sketch: continuous sources as horizontal bands, radar pulses and an impulse as short marks">
  ''' + '\n  '.join(o) + '''
</svg>
<figcaption>The realistic RFI list (<code>rfiRealistic</code> in <code>pipelineParams.m</code>) sketched in time and frequency over two pulsar
turns. Frequencies and repetition rates are the scenario's; widths are each source's occupied band (for BPSK the main lobe,
± the chip rate). Radar pulses (2 µs) and the impulse (200 ns) are drawn much wider than they are. Orange sources are on all the time
and look like extra noise in their channels; red ones are short and strong, which is what the excision stage looks for (see
<a href="excision.html">RFI excision</a>).</figcaption>
</figure>'''
open('rfiscene.svg', 'w').write(svg)
