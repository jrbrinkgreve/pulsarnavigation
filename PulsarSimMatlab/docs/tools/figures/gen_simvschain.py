import numpy as np
from plotsvg import Ax, pow10
# The fast simulator (simulateFold, block C2) against the voltage chain, full band, 1 turn per sub-int: median TOA
# error and the unknown-phase detection probability against rho (snrDB). Numbers from tests/expSimVsChain.m
# (log 2026-10-20_blockC.md, C2): chain = the saved sweeps snrSweep_phaseDE.mat (90 TOAs per point) and
# snrSweep_threshold_6oct.mat (180), simulated = 2000 sub-ints per point.
dE = np.array([-25, -20, -15, -10, -5, 0, 10, 20])
sChainE = np.array([102.8, 68.03, 24.53, 8.216, 2.907, 1.225, 0.5377, 0.4721])
sSimE   = np.array([99.31, 71.07, 24.91, 8.257, 2.916, 1.225, 0.5378, 0.4722])
dT = np.array([-20, -19, -18, -17, -16, -15])
sChainT = np.array([73.37, 60.32, 50.07, 40.33, 32.21, 25.66])
sSimT   = np.array([70.62, 59.23, 47.88, 38.85, 31.35, 25.12])
pChainT = np.array([0.400, 0.717, 0.922, 1.0, 1.0, 1.0]); nT = 180
pSimT   = np.array([0.377, 0.725, 0.961, 1.0, 1.0, 1.0])
o = []
A = Ax(70, 40, 330, 220, (-27, 22), (0.3, 200), logy=True)
o.append(A.frame(xticks=[-25, -20, -15, -10, -5, 0, 5, 10, 15, 20], yticks=[1, 10, 100], yfmt=lambda v: f'{v:g}',
                 xlab='ρ [dB] (pulse peak / receiver noise)', ylab='median σ of one TOA [µs]'))
def dot(Ax_, x, y, cls, r=4.5):
    return f'<path class="{cls}" d="M{float(Ax_.X(x)) - r:.1f} {float(Ax_.Y(y)):.1f} a{r} {r} 0 1 0 {2*r} 0 a{r} {r} 0 1 0 {-2*r} 0"/>'
for x, y in zip(dE, sChainE): o.append(dot(A, x, y, 'dot'))
for x, y in zip(dT, sChainT): o.append(dot(A, x, y, 'dot', 3.5))
for x, y in zip(np.r_[dE, dT], np.r_[sSimE, sSimT]): o.append(dot(A, x, y, 'dot2', 2.5))
lx, ly = A.x0 + 150, A.y0 + 14
o.append(f'<path class="dot" d="M{lx - 4.5:.1f} {ly:.1f} a4.5 4.5 0 1 0 9 0 a4.5 4.5 0 1 0 -9 0"/>')
o.append(f'<text class="sub" x="{lx + 10}" y="{ly + 4}">voltage chain (6 Oct sweeps)</text>')
o.append(f'<path class="dot2" d="M{lx - 2.5:.1f} {ly + 18:.1f} a2.5 2.5 0 1 0 5 0 a2.5 2.5 0 1 0 -5 0"/>')
o.append(f'<text class="sub" x="{lx + 10}" y="{ly + 22}">simulateFold</text>')
B = Ax(500, 40, 300, 220, (-20.5, -14.5), (0, 1.05))
o.append(B.frame(xticks=[-20, -19, -18, -17, -16, -15], yticks=[0, 0.25, 0.5, 0.75, 1], xlab='ρ [dB]',
                 ylab='detected, unknown phase (good TOAs)'))
for x, p in zip(dT, pChainT):
    e = np.sqrt(max(p*(1-p), 1/nT)/nT)
    o.append(f'<path class="thin" d="M{float(B.X(x)):.1f} {float(B.Y(min(p+e, 1.05))):.1f} V{float(B.Y(max(p-e, 0))):.1f}"/>')
    o.append(dot(B, x, p, 'dot', 4.5))
for x, p in zip(dT, pSimT): o.append(dot(B, x, p, 'dot2', 2.5))
o.append(f'<text class="sub" x="{float(B.X(-17.3)):.1f}" y="{float(B.Y(0.35)):.1f}">bars: chain ±1σ (180 sub-ints)</text>')
svg = '''<figure class="scroll">
<svg class="wide" viewBox="0 0 860 310" role="img" aria-label="Median TOA error and detection probability against rho: voltage chain and fast simulator agree">
  ''' + '\n  '.join(o) + '''
</svg>
<figcaption>The fast simulator against the voltage chain (full band, one turn per sub-int; <code>expSimVsChain</code>). Left: the median
error bar of one TOA from −25 to +20 dB; right: the fraction of sub-ints detected at the unknown phase near the threshold. Large dots:
the saved voltage-chain sweeps (90 and 180 TOAs per point); small: <code>simulateFold</code> (2000 sub-ints per point). The differences
are within the chain's statistics; the two chain sweeps differ from each other by as much near the threshold.</figcaption>
</figure>'''
open('simvschain.svg', 'w').write(svg)
