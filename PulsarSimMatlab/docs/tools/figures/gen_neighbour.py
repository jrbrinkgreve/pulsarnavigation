import numpy as np
from plotsvg import Ax
# Bug #8 (lessons.html): linear assignment splits each detected time bin between the two nearest phase bins, so neighbouring
# phase bins share noise. With the asymptotic fold weights W : W2 : WX = 1 : 2/3 : 1/6 (fold.html) the neighbour correlation is
# (1/6)/(2/3) = 0.25. a: schematic of the split. b: variance of a sum of k adjacent phase bins relative to the value that
# ignores the covariance: 1 + 2*rho*(k-1)/k, which tends to 1.5.
rho = (1/6)/(2/3)
o = []
o.append('<defs><marker id="ar-nb" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse">'
         '<path class="ah" d="M0 0 L10 5 L0 10 z"/></marker></defs>')
# ---------- a: schematic ----------
o.append('<text class="ttl" x="10" y="20">a · one time bin feeds two phase bins</text>')
bx = [40, 140, 240, 340]
for i, x in enumerate(bx):
    o.append(f'<rect class="box" x="{x}" y="170" width="96" height="60" rx="5"/>')
    o.append(f'<text class="sub" x="{x + 48}" y="248" text-anchor="middle">phase bin {i + 1}</text>')
o.append('<rect class="box rx" x="160" y="60" width="96" height="44" rx="5"/>')
o.append('<text class="sub" x="208" y="78" text-anchor="middle">one detected</text>')
o.append('<text class="sub" x="208" y="94" text-anchor="middle">time bin</text>')
o.append('<path class="line" marker-end="url(#ar-nb)" d="M196 104 L188 166"/>')
o.append('<path class="line" marker-end="url(#ar-nb)" d="M222 104 L280 166"/>')
o.append('<text class="sub" x="150" y="140" text-anchor="end">70 %</text>')
o.append('<text class="sub" x="262" y="140">30 %</text>')
o.append('<path class="curve2" d="M188 200 H288" stroke-dasharray="4 3"/>')
o.append('<text class="sub" x="238" y="194" text-anchor="middle">shared noise</text>')
o.append('<text class="sub" x="40" y="280">Its phase lies between two bin centres, so its power is split</text>')
o.append('<text class="sub" x="40" y="296">by distance (unbiased). Both bins now carry part of the same</text>')
o.append('<text class="sub" x="40" y="312">noise: neighbours are correlated, 0.25 for a fold of many turns.</text>')
# ---------- b ----------
B = Ax(520, 50, 300, 230, (1, 40), (0.9, 1.6))
o.append(B.frame(xticks=[1, 10, 20, 30, 40], yticks=[1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6], xlab='adjacent phase bins summed, k',
                 ylab='true variance ÷ assumed variance', yfmt=lambda v: f'{v:.1f}'))
k = np.arange(1, 41)
r = 1 + 2*rho*(k - 1)/k
o.append(f'<path class="dash" d="M{B.x0} {float(B.Y(1)):.1f} H{B.x0 + B.w}"/>')
o.append(f'<path class="dash" d="M{B.x0} {float(B.Y(1.5)):.1f} H{B.x0 + B.w}"/>')
o.append(f'<path class="curve4" d="{B.path(k, r)}"/>')
for kk, rr in zip(k[::3], r[::3]):
    o.append(f'<circle cx="{float(B.X(kk)):.1f}" cy="{float(B.Y(rr)):.1f}" r="3" class="bad"/>')
o.append(f'<text class="sub" x="{B.x0 + B.w - 4}" y="{float(B.Y(1)) + 16:.1f}" text-anchor="end">old error bars: neighbours independent</text>')
o.append(f'<text class="sub" x="{B.x0 + B.w - 4}" y="{float(B.Y(1.5)) - 8:.1f}" text-anchor="end">→ 1.5: variance 1.5× too small</text>')
o.append(f'<text class="ttl" x="520" y="20">b · summing neighbours: the missing covariance</text>')
svg = '''<figure class="scroll">
<svg class="wide" viewBox="0 0 860 330" role="img" aria-label="Bug 8: linear assignment makes neighbouring phase bins share noise; ignoring it makes the variance of a sum of bins up to 1.5 times too small">
  ''' + '\n  '.join(o) + '''
</svg>
<figcaption><b>Figure 1.</b> Lesson #8 in one picture. a: linear assignment splits a time bin's power between the two nearest phase bins, which removes
a timing bias but makes neighbours share noise (schematic, split 70/30 as an example). b: any sum over bins, such as the template fit,
then has more variance than a model with independent bins predicts: 1 + 2·0.25·(<var>k</var>−1)/<var>k</var>, up to 1.5× for a wide
pulse. The fold now tracks this covariance (<code>weightX</code>), so the error bars include it.</figcaption>
</figure>'''
open('neighbour.svg', 'w').write(svg)
