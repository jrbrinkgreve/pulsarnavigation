import numpy as np
# Time scales of the pipeline on one log ruler (parameters.html), from the parameters: RF sample 1/f_in, IQ sample 1/fs,
# channel sample 1/(4/3 * 3.125 MHz), detected bin 4 channel samples, phase bin T/nBin, pulse sigma_t and FWHM, dispersion
# sweep at DM 25 (since 20 Oct; 6.3 ms at DM 5), period T, file length L (0.2 s since 20 Oct); RFI: impulse 200 ns, radar pulse 2 us; time per TOA: J0437 worst case ~21 h
# (scenario.html), the simulated pulsar at -54 dB ~62 h (SNR 7, ~373*rho per turn).
T = 10e-3
items = [  # (seconds, name, value text, category, row)
    (1/4e9, 'RF sample', '0.25 ns', 'grid', 'a1'),
    (1/800e6, 'IQ sample', '1.25 ns', 'grid', 'a2'),
    (200e-9, 'impulse burst', '200 ns', 'rfi', 'a3'),
    (1/(3.125e6*4/3), 'channel sample', '0.24 µs', 'grid', 'b1'),
    (4/(3.125e6*4/3), 'detected bin', '0.96 µs', 'grid', 'a1'),
    (2e-6, 'radar pulse', '2 µs', 'rfi', 'b2'),
    (T/2048, 'phase bin', '4.9 µs', 'grid', 'a2'),
    (0.05*T/(2*np.sqrt(2*np.log(2))), 'pulse σ_t', '212 µs', 'psr', 'a1'),
    (0.05*T, 'pulse FWHM', '0.5 ms', 'psr', 'a2'),
    (31.51e-3, 'dispersion sweep (DM 25)', '31.5 ms', 'psr', 'b1'),
    (T, 'period', '10 ms', 'psr', 'a1'),
    (0.2, 'file length L', '0.2 s', 'run', 'a2'),
    (1.0, 'reference runs', '1 s', 'run', 'a1'),
    (21*3600, 'J0437 TOA, worst case', '21 h', 'run', 'a1'),
    (62*3600, 'TOA at −54 dB (sim. pulsar)', '62 h', 'run', 'b1')]
col = {'grid': 'var(--accent)', 'psr': 'var(--why)', 'rfi': 'var(--warn)', 'run': 'var(--fg)'}
x0, x1, lo, hi = 40, 820, -10, 6
X = lambda t: x0 + (np.log10(t) - lo)/(hi - lo)*(x1 - x0)
ya = 150
rowy = {'a1': ya - 26, 'a2': ya - 62, 'a3': ya - 98, 'b1': ya + 52, 'b2': ya + 88}
o = []
o.append(f'<path class="axis" d="M{x0} {ya} H{x1}"/>')
for e in range(lo, hi + 1):
    o.append(f'<path class="axis" d="M{X(10**e):.1f} {ya} v4"/>')
for t, lab in ((1e-9, '1 ns'), (1e-6, '1 µs'), (1e-3, '1 ms'), (1, '1 s'), (3600, '1 h'), (86400, '1 d')):
    o.append(f'<path class="axis" d="M{X(t):.1f} {ya - 4} v12"/><text class="sub" x="{X(t):.1f}" y="{ya + 20}" text-anchor="middle">{lab}</text>')
for t, name, val, cat, row in items:
    x = X(t)
    y = rowy[row]
    above = row[0] == 'a'
    yl = y + 6 if above else y - 28
    o.append(f'<path class="thin" d="M{x:.1f} {ya} V{yl:.1f}" opacity="0.7"/>')
    o.append(f'<circle cx="{x:.1f}" cy="{ya}" r="3.5" style="fill:{col[cat]}"/>')
    name = name.replace('σ_t', 'σ<tspan baseline-shift="sub" font-size="9">t</tspan>')
    xt = x - (16 if name == 'RF sample' else 0)
    o.append(f'<text class="sub" x="{xt:.1f}" y="{y - 8:.1f}" text-anchor="middle" style="fill:{col[cat]}">{name}</text>')
    o.append(f'<text class="lbl" x="{xt:.1f}" y="{y + 5:.1f}" text-anchor="middle">{val}</text>')
lx = 40
for cat, lab in (('grid', 'sampling and processing grid'), ('psr', 'the pulsar'), ('rfi', 'interference'), ('run', 'runs and time per TOA')):
    o.append(f'<circle cx="{lx + 4}" cy="{ya + 126}" r="4" style="fill:{col[cat]}"/><text class="sub" x="{lx + 13}" y="{ya + 130}">{lab}</text>')
    lx += 190
o.append(f'<text class="ttl" x="40" y="20">time scales, from a quarter nanosecond to days (log scale)</text>')
svg = '''<figure class="scroll">
<svg class="wide" viewBox="0 0 860 296" role="img" aria-label="Time scales of the pipeline on one logarithmic ruler, from the 0.25 ns RF sample to the 62 hours a TOA takes at -54 dB">
  ''' + '\n  '.join(o) + '''
</svg>
<figcaption>Every time scale of the pipeline on one ruler, computed from the parameters. Each step of the processing moves to a coarser grid
that is still much finer than the next scale it must resolve: the 0.96 µs detected bins resolve a 2 µs radar pulse; the 4.9 µs phase
bins resolve the 212 µs pulse. The 0.2 s file holds twenty turns (16 complete after the 31.5 ms sweep); a TOA at the design point needs about a million times longer, which
is why block C plans a fast simulator.</figcaption>
</figure>'''
open('timescales.svg', 'w').write(svg)
