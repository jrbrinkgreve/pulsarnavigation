import numpy as np
# Error-bar honesty across the validation record (results.html): measured rms / predicted rms of the TOA errors, with the 1-sigma
# uncertainty of each ratio (stated in the record, or ~1/sqrt(2n) for n TOAs). Numbers from results.html sections 1, 2, 3, 4, 9
# and the H0 run (std of the detection statistic T0, expected 1).
rows = [('header', 'Noise-free chain and −5 dB reference'),
        ('noise-free, 1 s (99 TOAs)', 1.053, 1/np.sqrt(2*99)),
        ('−5 dB reference run (9 TOAs)', 0.868, 0.24),
        ('Monte Carlo, 10 noise seeds (90 TOAs)', 0.950, 0.075),
        ('header', 'SNR sweep, phase D (90 TOAs per point)'),
        ('−15 dB (SNR 12)', 1.030, 0.075), ('−10 dB (SNR 35)', 1.040, 0.075), ('−5 dB (SNR 97)', 1.053, 0.075),
        ('0 dB (SNR 216)', 1.061, 0.075), ('+10 dB (SNR 496)', 1.044, 0.075), ('+20 dB (SNR 657)', 1.036, 0.075),
        ('header', 'Threshold sweep, good TOAs only (per SNR per sub-int)'),
        ('SNR 3.8 (72 good TOAs)', 1.20, 0.08), ('SNR 4.8', 0.97, 1/np.sqrt(2*130)), ('SNR 6.0', 1.01, 1/np.sqrt(2*166)),
        ('SNR 7.5', 1.01, 1/np.sqrt(2*180)), ('SNR 9.4', 1.02, 1/np.sqrt(2*180)), ('SNR 11.8', 1.02, 1/np.sqrt(2*180)),
        ('header', 'Noise only (runH0, 900 profiles)'),
        ('std of the detection statistic T0', 1.0076, 0.024)]
x0, x1, lo, hi = 330, 820, 0.6, 1.4
X = lambda v: x0 + (v - lo)/(hi - lo)*(x1 - x0)
o = []
y = 40
ys = []
for r in rows:
    if r[0] == 'header':
        y += 10
        o.append(f'<text class="ttl" x="10" y="{y + 4}">{r[1]}</text>')
        y += 18
        continue
    name, v, s = r
    o.append(f'<text class="sub" x="22" y="{y + 4}">{name}</text>')
    o.append(f'<path class="line" d="M{X(max(v - s, lo)):.1f} {y} H{X(min(v + s, hi)):.1f}"/>')
    o.append(f'<path class="dot2" d="M{X(v) - 4.5:.1f} {y} a4.5 4.5 0 1 0 9 0 a4.5 4.5 0 1 0 -9 0"/>')
    o.append(f'<text class="sub" x="{x1 + 8}" y="{y + 4}">{v:.3f}</text>')
    ys.append(y)
    y += 18
ytop, ybot = 30, y
o.insert(0, f'<rect class="fillwhy" x="{X(0.95):.1f}" y="{ytop}" width="{X(1.05) - X(0.95):.1f}" height="{ybot - ytop}"/>')
o.insert(1, f'<path class="dash" d="M{X(1):.1f} {ytop} V{ybot}"/>')
o.append(f'<path class="axis" d="M{x0} {ybot} H{x1}"/>')
for v in (0.6, 0.8, 1.0, 1.2, 1.4):
    o.append(f'<path class="axis" d="M{X(v):.1f} {ybot} v5"/><text class="sub" x="{X(v):.1f}" y="{ybot + 18}" text-anchor="middle">{v:g}</text>')
o.append(f'<text class="sub" x="{(x0 + x1)/2:.1f}" y="{ybot + 34}" text-anchor="middle">measured scatter ÷ stated error bar (1 = honest; green band ±5 %)</text>')
o.append(f'<text class="ttl" x="10" y="18">Are the error bars honest? Every check in the record</text>')
svg = f'''<figure class="scroll">
<svg class="wide" viewBox="0 0 860 {ybot + 46}" role="img" aria-label="Error-bar honesty across the validation record: measured over predicted scatter for every check, all near 1">
  ''' + '\n  '.join(o) + '''
</svg>
<figcaption>Every honesty check of the record on one axis: the measured scatter of (TOA − truth) divided by the stated error bar, with its
own 1σ uncertainty (from the record, or ≈ 1/√(2<var>n</var>) for <var>n</var> TOAs). All agree with 1 within a few percent, except
the good TOAs at SNR 3.8, below the design threshold of 6–7. The phase D points at −25 and −20 dB are off the scale (11 and 12) because
of outliers; with the good-TOA rule those outliers are removed (threshold sweep). The last row is the noise-only run: the detection
statistic has the predicted spread.</figcaption>
</figure>'''
open('honesty.svg', 'w').write(svg)
