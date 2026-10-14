# What each run script keeps fixed and what it varies (scripts.html), from the script headers and scripts.html:
# runMonteCarlo: main's sky, noise seeds seed+(1:10); runSNRSweep: 10 pulsar+noise realizations (seed+100(j-1), +1),
# 8 snrDB values with common random numbers, plus a zero-sky H0 pass; runH0: zero sky, 100 noise seeds seed+10000+i;
# runRFITest: main's sky and noise seed, 15 RFI cases, channel path without and with excision; runLockedRadar: main's sky and
# noise seed, locked-radar cases, excision variants, weak SNR from the leftover in baseline units.
cols = ['pulsar signal', 'receiver noise', 'signal strength', 'interference', 'front end']
rows = [('main.m', [('seed 42', 0), ('seed 43', 0), ('−5 dB', 0), ('off (rfiOn)', 0), ('frontEnd', 0)]),
        ('runMonteCarlo', [("main's", 0), ('10 seeds', 1), ('−5 dB', 0), ('none', 0), ('full band', 0)]),
        ('runSNRSweep', [('10 realizations', 1), ('new with each', 1), ('8 values, −25…+20 dB', 1), ('none', 0), ('full band', 0)]),
        ('runH0', [('none (zero sky)', 1), ('100 seeds', 1), ('no pulsar', 0), ('none', 0), ('full band', 0)]),
        ('runRFITest', [("main's", 0), ("main's seed", 0), ('−5 dB', 0), ('15 cases', 1), ('channels, ± excision', 1)]),
        ('runLockedRadar', [("main's", 0), ("main's seed", 0), ('−5 dB, then −5…−60 dB*', 1), ('locked radar cases', 1), ('channels, excision options', 1)])]
x0, w0, cw, y0, rh = 10, 140, 140, 40, 34
o = []
for j, c in enumerate(cols):
    o.append(f'<text class="ttl" x="{x0 + w0 + j*cw + cw/2:.0f}" y="{y0 - 10}" text-anchor="middle">{c}</text>')
for i, (name, cells) in enumerate(rows):
    y = y0 + i*rh
    o.append(f'<text class="mono" x="{x0}" y="{y + 21}">{name}</text>')
    for j, (txt, varied) in enumerate(cells):
        x = x0 + w0 + j*cw
        cls = 'box proc' if varied else 'box'
        o.append(f'<rect class="{cls}" x="{x + 3}" y="{y + 3}" width="{cw - 6}" height="{rh - 6}" rx="5"/>')
        o.append(f'<text class="sub" x="{x + cw/2:.0f}" y="{y + 21}" text-anchor="middle">{txt}</text>')
yb = y0 + len(rows)*rh + 18
o.append(f'<rect class="box proc" x="{x0}" y="{yb - 11}" width="22" height="14" rx="3"/><text class="sub" x="{x0 + 30}" y="{yb}">varied: what the experiment is about</text>')
o.append(f'<rect class="box" x="{x0 + 290}" y="{yb - 11}" width="22" height="14" rx="3"/><text class="sub" x="{x0 + 320}" y="{yb}">kept fixed</text>')
o.append(f'<text class="sub" x="{x0 + 430}" y="{yb}">* weak signals from noise-free folds plus the measured leftover</text>')
svg = f'''<figure class="scroll">
<svg class="wide" viewBox="0 0 860 {yb + 14}" role="img" aria-label="What each run script keeps fixed and what it varies: pulsar signal, receiver noise, signal strength, interference and front end">
  ''' + '\n  '.join(o) + '''
</svg>
<figcaption>Each experiment changes one or two things and keeps the rest as in <code>main.m</code>, so its result can be traced to that
change. Varying only the receiver noise tests the error bars; varying the pulsar too adds its self-noise; a zero sky is the
true "no pulsar" case; the RFI experiments keep the noise of <code>main.m</code>, so a difference between two cases is the
interference alone.</figcaption>
</figure>'''
open('experiments.svg', 'w').write(svg)
