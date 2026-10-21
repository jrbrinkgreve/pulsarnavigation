import numpy as np
# Project timeline (history.html): the blocks as bars over the calendar, 29 Sep - 21 Oct 2026, with the session days
# (the dated logs) as dots. Dates from history.html and roadmap.html.
day0 = 0  # 29 Sep
dates = ['29', '30'] + [str(k) for k in range(1, 22)]
sessions = [0, 1, 2, 3, 6, 7, 8, 9, 10, 13, 14, 15, 16, 17, 20, 21, 22]
x0, x1 = 40, 830
nd = len(dates)
X = lambda d: x0 + d/nd*(x1 - x0)
blocks = [('Block 1 · full-band chain', 0, 8, 'box sky', 'noise-free chain, −5 dB reference, Monte Carlo, SNR sweep, detection; closed 6 Oct'),
          ('A · channelized front end', 7, 10, 'box proc', 'filterbank, exact noise, data weights, optimal fit; complete 8 Oct'),
          ('B · RFI excision', 9, 11, 'box rx', 'blank before dedispersion, realistic scenario; closed for now 9 Oct'),
          ('review week', 13, 18, 'box val', 'each stage explained; code review, frozen test data, χ² report, docs'),
          ('processIQ', 20, 21, 'box proc', 'the processing as one function, IQ in, fold out'),
          ('C · fast simulator', 21, 23, 'box sky', 'simulateFold, RFI profiles, scintillation, scenario runs')]
o = []
o.append('<defs><marker id="ar-tl" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse">'
         '<path class="ah" d="M0 0 L10 5 L0 10 z"/></marker></defs>')
y = 34
for name, a, b, cls, txt in blocks:
    o.append(f'<rect class="{cls}" x="{X(a):.1f}" y="{y}" width="{X(b) - X(a):.1f}" height="20" rx="4"/>')
    if b - a >= 4:
        o.append(f'<text class="lbl" x="{X(a) + 8:.1f}" y="{y + 14}">{name}</text>')
        o.append(f'<text class="sub" x="{X(a) + 8:.1f}" y="{y + 36}">{txt}</text>')
    elif a < 12:
        o.append(f'<text class="lbl" x="{X(b) + 8:.1f}" y="{y + 14}">{name}</text>')
        o.append(f'<text class="sub" x="{x1:.1f}" y="{y + 36}" text-anchor="end">{txt}</text>')
    else:
        o.append(f'<text class="lbl" x="{X(a) - 8:.1f}" y="{y + 14}" text-anchor="end">{name}</text>')
        o.append(f'<text class="sub" x="{X(b):.1f}" y="{y + 36}" text-anchor="end">{txt}</text>')
    y += 52
ya = y + 6
o.append(f'<path class="axis" d="M{x0} {ya} H{x1}"/>')
for d in range(nd + 1):
    o.append(f'<path class="axis" d="M{X(d):.1f} {ya} v4"/>')
for d in range(nd):
    xc = (X(d) + X(d + 1))/2
    if d in sessions:
        o.append(f'<circle cx="{xc:.1f}" cy="{ya - 8}" r="4" class="dot"/>')
    o.append(f'<text class="sub" x="{xc:.1f}" y="{ya + 18}" text-anchor="middle">{dates[d]}</text>')
o.append(f'<text class="sub" x="{(X(0) + X(2))/2:.1f}" y="{ya + 34}" text-anchor="middle">September</text>')
o.append(f'<text class="sub" x="{(X(2) + X(nd))/2:.1f}" y="{ya + 34}" text-anchor="middle">October 2026</text>')
o.append(f'<path class="thin" d="M{X(2):.1f} {ya + 22} V{ya + 38}"/>')
o.append(f'<text class="sub" x="{x0}" y="{ya + 54}">dots: days with a session (a dated log in logs/)</text>')
o.append('<text class="ttl" x="40" y="20">From the first chain to the fast simulator</text>')
svg = f'''<figure class="scroll">
<svg class="wide" viewBox="0 0 860 {ya + 64}" role="img" aria-label="Project timeline: block 1, A, B, the review week, processIQ and block C as bars over the calendar from 29 September to 21 October 2026, with the session days">
  ''' + '\n  '.join(o) + '''
</svg>
<figcaption>The project so far on one calendar. Block 1 took a week, most of it validation against the truth. Blocks A and B followed,
from 6 to 9 October. The review week (12–16 Oct) added no new stages: each stage was explained, reviewed and refined. Then the
processing became one function (19 Oct), and block C, the fast simulator, started on 20 Oct.</figcaption>
</figure>'''
open('timeline.svg', 'w').write(svg)
