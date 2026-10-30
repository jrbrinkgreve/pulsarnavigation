import numpy as np
from plotsvg import Ax, pow10
# Two radar details of rfiSource / addNoiseAndRFI:
# a: the spectrum of one 2 us pulse with instant edges (RiseTime 0) and with 0.1 us raised-cosine edges,
#    envelope as in rfiBlock ('pulsed'), no chirp; the sidelobe envelope (maximum per log bin) is drawn;
# b: the rotating antenna's power gain over one scan, g(t) = max(exp(-4 ln2 (dt/BeamWidth)^2), 10^(SidelobeDB/10)),
#    for the realistic scenario's radar 1 (10 s scan, 39 ms beam, -30 dB sidelobes).
fs, W = 4e9, 2e-6
N = 2**19                                          # 131 us window, 7.6 kHz resolution
t = (np.arange(N) - N//4)/fs                       # pulse starts at 0
def envelope(tr):
    u = t
    e = ((u >= -tr/2) & (u < W + tr/2)).astype(float)
    if tr > 0:
        up = (u >= -tr/2) & (u < tr/2); dn = (u > W - tr/2) & (u < W + tr/2)
        e[up] = 0.5*(1 - np.cos(np.pi*(u[up] + tr/2)/tr))
        e[dn] = 0.5*(1 + np.cos(np.pi*(u[dn] - W + tr/2)/tr))
    return e
f = np.fft.rfftfreq(N, 1/fs)
edges = np.logspace(np.log10(0.1e6), np.log10(300e6), 120)
def sidelobes(tr):
    P = np.abs(np.fft.rfft(envelope(tr)))**2
    P /= P[0]
    fc, pk = [], []
    for a, b in zip(edges[:-1], edges[1:]):
        m = (f >= a) & (f < b)
        if m.any():
            fc.append(np.sqrt(a*b)); pk.append(P[m].max())
    return np.array(fc), np.array(pk)
o = []
A = Ax(80, 46, 330, 220, (0.1, 300), (1e-10, 2), logy=True)
lx = lambda v: A.x0 + (np.log10(v) - np.log10(0.1))/(np.log10(300) - np.log10(0.1))*A.w
o.append(A.frame(yticks=[1e-10, 1e-8, 1e-6, 1e-4, 1e-2, 1], ylab='power, relative to the carrier', yfmt=pow10))
for v, lab in [(0.1, '0.1'), (1, '1'), (10, '10'), (100, '100')]:
    o.append(f'<path class="axis" d="M{lx(v):.1f} {A.y0+A.h} v5"/><text class="sub" x="{lx(v):.1f}" y="{A.y0+A.h+18}" text-anchor="middle">{lab}</text>')
o.append(f'<text class="sub" x="{A.x0 + A.w/2:.1f}" y="{A.y0+A.h+34}" text-anchor="middle">distance from the carrier [MHz]</text>')
for tr, cls in ((0.0, 'curve2'), (0.1e-6, 'curve')):
    fc, pk = sidelobes(tr)
    ok = pk >= 1e-10
    pts = ' L'.join(f'{lx(a/1e6):.1f} {float(A.Y(b)):.1f}' for a, b in zip(fc[ok], pk[ok]))
    o.append(f'<path class="{cls}" d="M{pts}"/>')
    print('RiseTime', tr, 'at 30 MHz', pk[np.argmin(abs(fc - 30e6))], 'at 100 MHz', pk[np.argmin(abs(fc - 100e6))])
xc = lx(3.125)
o.append(f'<path class="dash" d="M{xc:.1f} {A.y0} V{A.Y(1e-6):.1f}"/>')
o.append(f'<text class="sub" x="{xc+4:.1f}" y="{A.y0+12}">one channel (3.125 MHz)</text>')
x5 = lx(5)
ly1, ly2 = float(A.Y(1e-8)), float(A.Y(1e-8)) + 16
o.append(f'<path class="curve2" d="M{A.x0+10} {ly1-4:.1f} h18"/><text class="sub" x="{A.x0+34}" y="{ly1:.1f}">instant edges: falls as 1/f²</text>')
o.append(f'<path class="curve" d="M{A.x0+10} {ly2-4:.1f} h18"/><text class="sub" x="{A.x0+34}" y="{ly2:.1f}">0.1 µs edges (default): 1/f⁶</text>')
o.append(f'<text class="ttl" x="{A.x0}" y="20">a · one 2 µs radar pulse: how far its power spreads</text>')
# panel b: rotating antenna gain
B = Ax(520, 46, 310, 220, (0, 10), (-33, 3))
ts = np.linspace(0, 10, 20001)
dt = np.mod(ts - 5 + 5, 10) - 5
g = np.maximum(np.exp(-4*np.log(2)*(dt/0.039)**2), 10**(-30/10))
o.append(B.frame(xticks=[0, 2, 4, 6, 8, 10], yticks=[-30, -20, -10, 0], xlab='time in one scan [s]', ylab='antenna gain towards us [dB]'))
gd = 10*np.log10(g)
keep = np.r_[True, np.abs(np.diff(gd)) > 1e-9] | np.r_[np.abs(np.diff(gd)) > 1e-9, True]
o.append(f'<path class="curve4" d="{B.path(ts[keep], gd[keep])}"/>')
frac = np.mean(g > 10**(-30/10) * 1.0000001)
print('fraction of the scan above the sidelobe floor', frac)
o.append(f'<text class="sub" x="{B.X(5.3):.1f}" y="{B.Y(-2):.1f}">main beam: 39 ms</text>')
o.append(f'<text class="sub" x="{B.X(5.3):.1f}" y="{B.Y(-2)+14:.1f}">at −3 dB, once per scan</text>')
o.append(f'<text class="sub" x="{B.X(0.3):.1f}" y="{B.Y(-24):.1f}">sidelobes, −30 dB:</text>')
o.append(f'<text class="sub" x="{B.X(0.3):.1f}" y="{B.Y(-24)+14:.1f}">{100*(1-frac):.1f} % of the time</text>')
o.append(f'<text class="ttl" x="{B.x0}" y="20">b · a rotating radar (10 s scan)</text>')
svg = '''<figure class="scroll">
<svg class="wide" viewBox="0 0 860 310" role="img" aria-label="Radar details: the spectrum of one pulse with instant and with smooth edges, and the gain of a rotating radar antenna over one scan">
  ''' + '\n  '.join(o) + '''
</svg>
<figcaption><b>Figure 3.5.</b> Two radar details, computed with the code's formulas. (a) The power spectrum of one 2 µs pulse (envelope only, sidelobe
peaks drawn). With instant edges (orange) it falls as 1/<var>f</var>²: still 2.5·10⁻⁴ of the carrier 10 MHz away. With the default
0.1 µs raised-cosine edges (blue) it falls as 1/<var>f</var>⁶ beyond about 5 MHz. (b) The power gain of the realistic scenario's radar 1
towards us over one 10 s scan: most of the time only its −30 dB sidelobes reach us.</figcaption>
</figure>'''
open('radar.svg', 'w').write(svg)
