import numpy as np
from plotsvg import Ax
# Why 'nearest' phase assignment was biased (fold.html): with 1 us time bins and NBin = 1024 phase bins of a 10 ms turn
# (9.765625 us), every turn has exactly 10,000 time bins, so the rounding error of each time bin's phase repeats
# identically turn after turn. Bin centroids at binTime0 + k*1 us (binTime0 = 0.499 us at 500 MHz, the run of 30 Sep).
T, dt, nb, t0 = 10000.0, 1.0, 1024, 0.499          # us
dph = T/nb
o = []
A = Ax(90, 46, 720, 200, (-60, 60), (-6, 6))
o.append(A.frame(xticks=[-60, -40, -20, 0, 20, 40, 60], yticks=[-5, 0, 5], xlab='time bin relative to the pulse centre [µs]',
                 ylab='rounding error of its phase position [µs]'))
o.append(f'<path class="dash" d="M{A.x0} {A.Y(0):.1f} H{A.x0 + A.w}"/>')
for turn, cls, r in ((0, 'dot', 3.5), (7, 'dot2', 2.0)):
    k = np.arange(-70, 71)
    t = turn*T + 5000.0 + k*dt + (t0 - 0.5)          # around the pulse centre of this turn (TRef = T/2)
    x = ((t - 5000.0)/T % 1.0)
    x = np.where(x >= 0.5, x - 1, x)*nb              # position in phase bins, in [-nb/2, nb/2)
    err = (np.round(x) - x)*dph                      # nearest bin minus true position, in us
    rel = (t - turn*T - 5000.0)
    sel = np.abs(rel) <= 60
    for a, b in zip(rel[sel], err[sel]):
        o.append(f'<path class="{cls}" d="M{A.X(a) - r:.1f} {A.Y(b):.1f} a{r} {r} 0 1 0 {2*r} 0 a{r} {r} 0 1 0 {-2*r} 0"/>')
o.append(f'<text class="ttl" x="90" y="20">Nearest-bin assignment: the same rounding errors in every turn</text>')
o.append(f'<path class="dot" d="M{A.x0 + 470} 31 a3.5 3.5 0 1 0 7 0 a3.5 3.5 0 1 0 -7 0"/><text class="sub" x="{A.x0 + 482}" y="35">turn 0</text>')
o.append(f'<path class="dot2" d="M{A.x0 + 540} 31 a2 2 0 1 0 4 0 a2 2 0 1 0 -4 0"/><text class="sub" x="{A.x0 + 552}" y="35">turn 7: on top of turn 0</text>')
svg = '''<figure class="scroll">
<svg class="wide" viewBox="0 0 860 300" role="img" aria-label="Nearest-bin assignment: the rounding errors of the time bins near the pulse are identical in every turn">
  ''' + '\n  '.join(o) + '''
</svg>
<figcaption>Why <code>'nearest'</code> was biased. Each 1 µs time bin is put in the nearest of 1024 phase bins (9.77 µs), which moves it by up to
±4.9 µs. With 10,000 time bins per turn the pattern repeats exactly every turn (turn 7, small dots, lies on top of turn 0), so the
errors never average out; weighted by the pulse they left a TOA bias of 0.5 ns. Linear assignment splits each time bin between the two
neighbouring phase bins so that its average position is exact: no rounding error at all.</figcaption>
</figure>'''
open('nearest.svg', 'w').write(svg)
