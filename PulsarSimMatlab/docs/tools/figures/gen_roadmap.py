# The roadmap at a glance (roadmap.html): the blocks on the calendar, block C's steps with their status,
# and what D and E hold. Statuses from roadmap.html / logs (update the lists below and rerun).
# Status codes: V validated (Jasper's run), C done (Claude's run), R running / under investigation, N next / open.
import datetime as dt

DAY0 = dt.date(2026, 9, 29)
TODAY = dt.date(2026, 10, 22)
blocks = [  # (label, first day, last day, lane, status)
    ('Block 1 · full-band chain', dt.date(2026, 9, 29), dt.date(2026, 10, 6), 0, 'V'),
    ('A · channels', dt.date(2026, 10, 6), dt.date(2026, 10, 8), 1, 'V'),
    ('B · excision', dt.date(2026, 10, 8), dt.date(2026, 10, 9), 0, 'V'),
    ('review week', dt.date(2026, 10, 12), dt.date(2026, 10, 16), 1, 'V'),
    ('processIQ', dt.date(2026, 10, 19), dt.date(2026, 10, 19), 0, 'V'),
    ('C · fast simulator', dt.date(2026, 10, 20), dt.date(2026, 10, 22), 1, 'V'),     # closed 22 Oct
    ('D · motion', dt.date(2026, 10, 21), dt.date(2026, 10, 22), 0, 'V', 'right'),   # closed 22 Oct; label beside the bar
]
cCols = [  # block C in detail: (column title, [(item, status)])
    ('simulator (simulateFold)', [('C1 noise model', 'V'), ('C2 against the voltage chain', 'C'),
                                  ('C4a channels, band-edge taper', 'V'), ('C3a passes and gaps', 'V'),
                                  ('C3b scintillation (+ grid fix)', 'V'), ('C3c adaptive sub-ints', 'V'),
                                  ('dual pol, groups over gaps', 'V'), ('multi-component profiles', 'V')]),
    ('RFI model', [('C4b RFI profiles', 'V'), ('testRFIModel (synthetic)', 'V'),
                   ('expSimRFI: cost of RFI', 'C'), ('testRFIChain (voltage chain)', 'V'),
                   ('narrow carriers: fixed, passes', 'V'), ('noise calibration: simulator', 'V'),
                   ('noise calibration: chain check', 'C'), ('channel 61: blanked bins', 'C')]),
    ('scenario (runScenario)', [('default DM 25, L 0.2 s', 'C'), ('clean, 30 days', 'C'),
                                ('B0329+54 five components', 'C'), ('RFI variants, calibrated', 'C'),
                                ('J0437 measured profile', 'C'), ('adopt J0437 profile', 'C'),
                                ('100-day final runs', 'C')]),
    ('open questions', [('error bars at SNR ~7: explained', 'C'), ('J0437 narrow pulls: noise', 'C'),
                        ('detectRFI thresholds per channel', 'V'), ('rfiProfile: lags beyond 7', 'V'),
                        ('low-SNR outliers: robust fit (D)', 'V')]),
]
after = [  # blocks D and E: (title, now/next, [(item, status)])
    ('D · motion and phase', 'closed', [('D1 trajectory, delays to the barycentre', 'V'),
                                     ('D2 phase predictor, folding with it', 'V'),
                                     ('D3 simulator truth, receiver clock', 'V'),
                                     ('D4 navigation fit, ship filter', 'V'),
                                     ('D5 DE440, TDB − TT, precession', 'V')]),
    ('E · full system, simulated', 'next', [('several pulsars at full scale', 'N'),
                                ('navigation solution: position + clock', 'N'),
                                ('3×3 array, second polarization', 'N')]),
    ('F · real implementation', 'later', [('IQ from the antenna, mixer and SDR', 'N'),
                                ('processIQ on the recordings', 'N'),
                                ('GPU / FPGA port, stage by stage', 'N'),
                                ('real time', 'N')]),
]
cls = {'V': 'box val', 'C': 'box proc', 'R': 'box rx', 'N': 'box'}
o = []
# ---- 1. calendar -----------------------------------------------------------------------------
xa, xb, nd = 20, 590, 25                     # axis from DAY0 over nd days
px = (xb - xa) / nd
X = lambda d: xa + (d - DAY0).days * px
o.append('<text class="ttl" x="10" y="16">1 · the blocks on the calendar</text>')
yAx = 112
o.append(f'<path class="axis" d="M{xa} {yAx} H{xb}"/>')
for d in [dt.date(2026, 9, 29), dt.date(2026, 10, 6), dt.date(2026, 10, 13), dt.date(2026, 10, 20)]:
    o.append(f'<path class="axis" d="M{X(d):.1f} {yAx} v4"/>')
    o.append(f'<text class="sub" x="{X(d):.1f}" y="{yAx + 16}" text-anchor="middle">{d.day} {d.strftime("%b")}</text>')
for lab, a, b, lane, st, *where in blocks:
    y = 36 + lane * 36
    x0, x1 = X(a), X(b) + px
    dash = ' stroke-dasharray="5 4"' if st == 'N' else ''
    o.append(f'<rect class="{cls[st]}" x="{x0:.1f}" y="{y}" width="{x1 - x0:.1f}" height="18" rx="4"{dash}/>')
    if where == ['right']:
        o.append(f'<text class="lbl" x="{x1 + 4:.1f}" y="{y + 13}">{lab}</text>')
    else:
        o.append(f'<text class="lbl" x="{x0:.1f}" y="{y - 4}">{lab}</text>')
xt = X(TODAY) + px
o.append(f'<path class="dash" d="M{xt:.1f} 76 V{yAx}"/><text class="sub" x="{xt + 4:.1f}" y="{yAx - 4}">today</text>')
o.append(f'<path class="line" d="M{X(TODAY) + px:.1f} 81 h18" marker-end="url(#ar-rd)"/>')
for k, (lab, when, _) in enumerate(after):
    y = 36 + k * 36
    st = {'now': 'R', 'closed': 'V'}.get(when, 'N')
    dash = ' stroke-dasharray="5 4"' if st == 'N' else ''
    o.append(f'<rect class="{cls[st]}" x="640" y="{y}" width="210" height="18" rx="4"{dash}/>')
    o.append(f'<text class="lbl" x="745" y="{y + 13}" text-anchor="middle">{lab} ({when})</text>')
# ---- 2. block C in detail --------------------------------------------------------------------
y2 = 156
o.append(f'<text class="ttl" x="10" y="{y2}">2 · block C, step by step</text>')
cw, ch, gap = 205, 22, 6
for j, (title, items) in enumerate(cCols):
    x = 10 + j * (cw + 10)
    o.append(f'<text class="lbl" x="{x}" y="{y2 + 22}" font-weight="600">{title}</text>')
    for i, (it, st) in enumerate(items):
        y = y2 + 30 + i * (ch + gap)
        dash = ' stroke-dasharray="5 4"' if st == 'N' else ''
        o.append(f'<rect class="{cls[st]}" x="{x}" y="{y}" width="{cw}" height="{ch}" rx="5"{dash}/>')
        o.append(f'<text class="sub" x="{x + 8}" y="{y + 15}">{it}</text>')
nMax = max(len(items) for _, items in cCols)
y3 = y2 + 30 + nMax * (ch + gap) + 26
# ---- 3. after C ------------------------------------------------------------------------------
o.append(f'<text class="ttl" x="10" y="{y3}">3 · blocks D, E and F</text>')
for j, (title, when, items) in enumerate(after):
    x = 10 + j * 283
    o.append(f'<text class="lbl" x="{x}" y="{y3 + 22}" font-weight="600">{title} ({when})</text>')
    for i, (it, st) in enumerate(items):
        y = y3 + 30 + i * (ch + gap)
        dash = ' stroke-dasharray="5 4"' if st == 'N' else ''
        o.append(f'<rect class="{cls[st]}" x="{x}" y="{y}" width="275" height="{ch}" rx="5"{dash}/>')
        o.append(f'<text class="sub" x="{x + 8}" y="{y + 15}">{it}</text>')
yL = y3 + 30 + max(len(i) for _, _, i in after) * (ch + gap) + 16
leg = [('V', 'validated (Jasper\'s run)'), ('C', 'done (Claude\'s run)'), ('R', 'running / under investigation'),
       ('N', 'next / open')]
xl = 10
for st, txt in leg:
    dash = ' stroke-dasharray="4 3"' if st == 'N' else ''
    o.append(f'<rect class="{cls[st]}" x="{xl}" y="{yL - 11}" width="22" height="14" rx="3"{dash}/>')
    o.append(f'<text class="sub" x="{xl + 30}" y="{yL}">{txt}</text>')
    xl += 210
svg = f'''<figure class="scroll">
<svg class="wide" viewBox="0 0 860 {yL + 10}" role="img" aria-label="Roadmap: the blocks on the calendar from 29 September, block C's steps with their status, and what blocks D and E hold">
  <defs>
    <marker id="ar-rd" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse">
      <path class="ah" d="M0 0 L10 5 L0 10 z"/>
    </marker>
  </defs>
  ''' + '\n  '.join(o) + f'''
</svg>
<figcaption>The roadmap on {TODAY.day} {TODAY.strftime("%B")}. (1) The blocks on the calendar. C, the fast simulator, ran on 20–22 Oct; it came before D because a
TOA takes hours to tens of hours of data, far beyond the voltage-level simulation. D, motion, phase and the navigation fit, ran on 21–22 Oct and is closed. (2) Block C's
steps, all done. The simulator, the voltage-chain test of the RFI model, the per-channel detectRFI thresholds and the extra lags in the
RFI profiles are validated (Jasper's run of 22 Oct). The experiments, the RFI variants, J0437's measured profile and the 100-day runs are
done in Claude's runs; the robust fit against the low-SNR tail came with D4. (3) Blocks D and E: D1–D5 are validated, a few items are
parked (navigation page); E brings the full system in simulation, F the real implementation on antenna data and a GPU or FPGA.
Generated by <code>docs/tools/figures/gen_roadmap.py</code> (update its lists).</figcaption>
</figure>'''
open('roadmap.svg', 'w').write(svg)
