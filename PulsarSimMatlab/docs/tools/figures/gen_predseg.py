import numpy as np
from plotsvg import Ax
from delaymodel import delay, track, fixed, EPOCH
# makePredictor's segments around a ship's turn (fold.html#prediction #predictor): testMakePredictor test 5 replayed with a port of
# makePredictor's loop (Chebyshev degree 12 on one-hour segments, checked at the extrema between the nodes, split in two until
# the check is within 1 ns, shortest half-segment 1 s) and of delayModel (delaymodel.py). A 100 Hz pulsar in B0329+54's
# direction; the ship sails east at about 9 m/s, turns north after 3.3 h, and the track is linear between its waypoints.
# The port gives the test's numbers: 15 segments, the shortest 7 s.
psr = dict(raDeg=53.24738, decDeg=54.57878); f0 = 100.0; t1 = EPOCH
dlon = 10 * 3 * 3600 / (111.32e3 * np.cos(np.radians(52))); dlat = 10 * 3 * 3600 / 111.2e3
ship = track(t1 + np.array([0, 3.3, 6]) * 3600, [52, 52, 52 + dlat], [5, 5 + dlon, 5 + dlon])
d0 = delay(t1, ship, psr)['total'][0]
phi = lambda tau: f0 * (tau + (delay(t1 + tau, ship, psr)['total'] - d0))
K = 12; L = 6 * 3600.0; tol = 1e-9; minseg = 1.0
xn = np.cos(np.pi * (np.arange(K + 1) + 0.5) / (K + 1)); xc = np.cos(np.pi * np.arange(K + 1) / K)
Tn = np.cos(np.outer(np.arange(K + 1), np.arccos(xn)))
def cheb(c, x):
    b1 = np.zeros_like(x); b2 = b1.copy()
    for k in range(len(c) - 1, 0, -1):
        b0 = 2 * x * b1 - b2 + c[k]; b2 = b1; b1 = b0
    return x * b1 - b2 + c[0]
edges = np.unique(np.r_[np.arange(0, L, 3600.0), L]); todo = list(zip(edges[:-1], edges[1:])); seg = []
while todo:
    a, b = todo.pop(0); m = (a + b) / 2; h = (b - a) / 2
    c = (2 / (K + 1)) * (Tn @ phi(m + h * xn)); c[0] /= 2
    err = np.max(np.abs(cheb(c, xc) - phi(m + h * xc))) / f0
    if err > tol and h >= minseg:
        todo += [(a, m), (m, b)]; continue
    seg.append((a, b, err))
seg.sort()
tTurn = 3.3 * 3600
o = []
# a: the track
A = Ax(60, 60, 130, 130, (-10, 120), (-10, 120))
o.append(A.frame(xticks=[0, 50, 100], yticks=[0, 50, 100], xlab='east [km]', ylab='north [km]'))
o.append(f'<path class="s1" d="{A.path([0, 108, 108], [0, 0, 108])}"/>')
for x, y, lab, dx, anc in ((0, 0, '0 h', 0, 'middle'), (108, 0, 'turn, 3.3 h', -4, 'end'), (108, 108, '6 h', -8, 'end')):
    o.append(f'<circle class="s1f" cx="{float(A.X(x)):.1f}" cy="{float(A.Y(y)):.1f}" r="4"/>')
    o.append(f'<text class="sub" x="{float(A.X(x)) + dx:.1f}" y="{float(A.Y(y)) - 9:.1f}" text-anchor="{anc}">{lab}</text>')
o.append('<text class="ttl" x="30" y="22">a · the ship</text>')
# b: the delay the ship adds, with the segment boundaries
t = np.arange(0, L + 1, 60.0)
x = (delay(t1 + t, ship, psr)['total'] - delay(t1 + t, fixed(52, 5), psr)['total']) * 1e6
B = Ax(290, 60, 270, 200, (0, 6), (-150, 450))
for v in (-100, 0, 100, 200, 300, 400):
    o.append(f'<path class="grid" d="M{B.x0} {float(B.Y(v)):.1f} H{B.x0 + B.w}"/>')
o.append(B.frame(xticks=[0, 1, 2, 3, 4, 5, 6], yticks=[0, 200, 400], xlab='hours', ylab='delay the ship adds [µs]'))
o.append(f'<path class="dash" d="M{float(B.X(3.3)):.1f} {B.y0} V{B.y0 + B.h}"/>')
o.append(f'<path class="s1" d="{B.path(t / 3600, x)}"/>')
for a, b, e in seg:
    xa = float(B.X(a / 3600)); o.append(f'<path class="s2" style="stroke-width:1.2" d="M{xa:.1f} {B.y0 + B.h - 14} v14"/>')
o.append(f'<path class="s2" style="stroke-width:1.2" d="M{float(B.X(6)):.1f} {B.y0 + B.h - 14} v14"/>')
o.append(f'<text class="sub" x="{float(B.X(3.3)) + 5:.1f}" y="{B.y0 + 12}">the turn: a kink</text>')
o.append(f'<text class="sub" x="{float(B.X(3.55)):.1f}" y="{B.y0 + B.h - 20}" style="fill:var(--series-2)">segment boundaries</text>')
o.append('<text class="ttl" x="290" y="22">b · a kink in the delay</text>')
# c: zoom on the turn, segment lengths
w0, w1 = tTurn - 300, tTurn + 300
Cz = Ax(640, 60, 200, 200, (-300, 300), (0, np.log10(3600)))
lt = [(1, '1 s'), (10, '10 s'), (60, '1 min'), (600, '10 min'), (3600, '1 h')]
for v, lab in lt:
    y = float(Cz.Y(np.log10(v))); o.append(f'<path class="grid" d="M{Cz.x0} {y:.1f} H{Cz.x0 + Cz.w}"/>')
    o.append(f'<path class="axis" d="M{Cz.x0} {y:.1f} h-5"/><text class="sub" x="{Cz.x0 - 8}" y="{y + 4:.1f}" text-anchor="end">{lab}</text>')
o.append(Cz.frame(xticks=[-300, -150, 0, 150, 300], xlab='seconds from the turn', ylab='segment length'))
o.append(f'<path class="dash" d="M{float(Cz.X(0)):.1f} {Cz.y0} V{Cz.y0 + Cz.h}"/>')
for a, b, e in seg:
    if b < w0 or a > w1:
        continue
    xa, xb = float(Cz.X(max(a, w0) - tTurn)), float(Cz.X(min(b, w1) - tTurn)); y = float(Cz.Y(np.log10(b - a)))
    o.append(f'<path class="s2" style="stroke-width:5" stroke-linecap="butt" d="M{xa + 0.6:.1f} {y:.1f} H{xb - 0.6:.1f}">'
             f'<title>{b - a:.1f} s, checked to {e * 1e9:.2f} ns</title></path>')
short = min(b - a for a, b, e in seg)
o.append(f'<text class="sub" x="{float(Cz.X(12)):.1f}" y="{float(Cz.Y(np.log10(short))) + 4:.1f}" style="fill:var(--series-2)">{short:.0f} s</text>')
o.append('<text class="ttl" x="640" y="22">c · the segments near the turn</text>')

nS = len(seg); short = min(b - a for a, b, e in seg); eTurn = max(e for a, b, e in seg)
svg = '''<figure class="scroll">
<svg class="wide" viewBox="0 0 860 310" role="img" aria-label="makePredictor's segments around a ship's turn: the ship's track, the delay it adds with a kink at the turn, and the segment lengths, which halve from one hour down to 7 seconds at the turn">
  ''' + '\n  '.join(o) + '''
</svg>
<figcaption>How <code>makePredictor</code> handles a ship's sharp turn, replayed from its test 5. a: the ship sails east, and after
3.3 h turns north. b: the delay its position adds to a 100 Hz pulsar in B0329+54's direction. At the turn the delay has a kink, which no
smooth polynomial can follow. The orange ticks mark where segments begin. c: the segment lengths near the turn. Every segment that
misses the 1 ns check is split in two, so the segments halve towards the turn, from one hour down to ''' + f'{short:.0f}' + ''' s. Away
from the turn one hour passes. The result is ''' + f'{nS}' + ''' segments, and every one is within 1 ns. Generated by
<code>docs/tools/figures/gen_predseg.py</code>, with ports of <code>makePredictor</code>'s loop and of <code>delayModel</code>.</figcaption>
</figure>'''
open('predseg.svg', 'w').write(svg)
print(nS, 'segments, shortest', short, 'largest check', eTurn)
