import numpy as np
# Error-bar honesty across the validation record (results.html): measured rms / predicted rms of the TOA errors, with the 1-sigma
# uncertainty of each ratio (stated in the record, or ~1/sqrt(2n) for n TOAs or pulls). Numbers from results.html sections 1, 2, 3, 4, 9,
# the H0 run (std of the detection statistic T0, expected 1), and the block C and D checks (results.html 12-15).
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
        ('std of the detection statistic T0', 1.0076, 0.024),
        # added 22 Oct 2026: the channel path, the fast simulator and navigation (logs/2026-10-20_blockC.md; results.html 12-15)
        ('header', 'Channel path with blanking (A3b, 50 seeds)'),
        ('Monte Carlo on processIQ (400 TOAs)', 0.956, 0.034),
        ('header', 'Fast simulator against its truth (block C)'),
        ('−5 dB, 2000 sub-ints (expSimVsChain)', 1.002, 1/np.sqrt(2*2000)),
        ('RFI: noise only (expSimRFI, 400 TOAs)', 1.004, 0.035),
        ('RFI: realistic scenario', 0.957, 0.035),
        ('RFI: narrow carriers', 1.069, 0.035),
        ('B0329+54, 100 days, all 978 TOAs', 1.156, 1/np.sqrt(2*978)),
        ('J0437, 100 days, 76 TOAs', 1.129, 1/np.sqrt(2*76)),
        ('header', 'Navigation fit and filter (block D)'),
        ('navFit, east, 400 realizations', 0.991, 1/np.sqrt(800)),
        ('navFit, radial', 0.973, 1/np.sqrt(800)),
        ('navFit, 2 % outliers, robust', 1.045, 1/np.sqrt(800)),
        ('navFit, clock estimate (crystal, 1000)', 0.981, 1/np.sqrt(2000)),
        ('navFilter, 40 ships, east (160 pulls)', 1.029, 1/np.sqrt(320)),
        ('navFilter, north', 0.911, 1/np.sqrt(320)),
        ('runShip, calibrated, east (32 pulls)', 0.86, 1/np.sqrt(64)),
        ('runShip, calibrated, north', 0.84, 1/np.sqrt(64)),
        ('its ideal twin, east', 0.81, 1/np.sqrt(64)),
        ('its ideal twin, north', 0.81, 1/np.sqrt(64))]
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
own 1σ uncertainty (from the record, or ≈ 1/√(2<var>n</var>) for <var>n</var> TOAs or pulls). The voltage-chain checks agree with 1
within a few percent, except the good TOAs at SNR 3.8, below the design threshold of 6–7. The phase D points at −25 and −20 dB are off the scale (11 and 12) because
of outliers; with the good-TOA rule those outliers are removed (threshold sweep). The noise-only row shows that the detection
statistic has the predicted spread. The fast simulator reproduces honest error bars; its 100-day runs show B0329+54's heavy tail (1.9 %
of the TOAs beyond 3σ), while the core of the distribution has a width of 1.01. The navigation rows are pulls of the fitted position
(and clock) against the truth. The fit agrees with 1; the filter's north pulls are 0.91. The ship's pulls (32 per axis) are 0.84–0.86,
and its ideal twin's on the same TOAs 0.81: the TOA chain adds nothing measurable.</figcaption>
</figure>'''
open('honesty.svg', 'w').write(svg)
