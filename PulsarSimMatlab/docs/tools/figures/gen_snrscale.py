import numpy as np
from plotsvg import Ax, pow10
# What snrDB means: the matched-filter SNR of one pulse against snrDB, from the radiometer equation that
# addNoiseAndRFI prints (predictPerformance: flat band B = 400 MHz, Gaussian power profile, self-noise).
B, T, duty = 400e6, 10e-3, 5.0
sp = duty/100*T/(2*np.sqrt(2*np.log(2)))            # 212.3 us
tt = np.linspace(-8*sp, 8*sp, 40001); p = np.exp(-0.5*(tt/sp)**2)
def snrPulse(rho):
    return np.sqrt(B*np.trapezoid((rho*p)**2/(1 + rho*p)**2, tt))
dB = np.linspace(-60, 25, 341)
s = np.array([snrPulse(10**(d/10)) for d in dB])
weak = np.sqrt(B*sp*np.sqrt(np.pi))                  # SNR = weak*rho for weak pulses
print('weak-signal factor', weak)
for d in (-54, -20, -5, 20):
    print(d, snrPulse(10**(d/10)))
o = []
A = Ax(90, 40, 560, 250, (-60, 25), (1e-4, 3e3), logy=True)
o.append(A.frame(xticks=[-60, -50, -40, -30, -20, -10, 0, 10, 20], yticks=[1e-3, 1e-2, 1e-1, 1, 10, 100, 1000],
                 xlab='snrDB = 10·log₁₀(S_peak / N)', ylab='SNR of one pulse (matched filter)', yfmt=lambda v: f'{v:g}'))
dw = dB[dB <= 5]
o.append(f'<path class="curve2" stroke-dasharray="6 4" d="{A.path(dw, weak*10**(dw/10))}"/>')
o.append(f'<path class="curve" d="{A.path(dB, s)}"/>')
y7 = A.Y(7)
o.append(f'<path class="dash" d="M{A.x0} {y7:.1f} H{A.x0+A.w}"/>')
o.append(f'<text class="sub" x="{A.X(-59):.1f}" y="{y7-6:.1f}">SNR 7: one good TOA</text>')
pts = [(-54, 'design target', 'right'), (-20, 'original plan (30 Sep)', 'right'), (-5, 'test runs', 'left')]
for d, lab, side in pts:
    v = snrPulse(10**(d/10)); x, y = float(A.X(d)), float(A.Y(v))
    o.append(f'<path class="dot2" d="M{x-4:.1f} {y:.1f} a4 4 0 1 0 8 0 a4 4 0 1 0 -8 0"/>')
    txt = f'{d} dB: {v:.2g} per pulse' if v < 1 else f'{d} dB: {v:.3g} per pulse'
    if side == 'right':
        o.append(f'<text class="lbl" x="{x+9:.1f}" y="{y+4:.1f}">{txt}</text>')
        o.append(f'<text class="sub" x="{x+9:.1f}" y="{y+18:.1f}">{lab}</text>')
    else:
        o.append(f'<text class="lbl" x="{x-9:.1f}" y="{y-14:.1f}" text-anchor="end">{txt}</text>')
        o.append(f'<text class="sub" x="{x-9:.1f}" y="{y}" text-anchor="end">{lab}</text>')
o.append(f'<text class="sub" x="{A.X(-36):.1f}" y="{A.Y(weak*10**(-3.0)):.1f}" transform="rotate(-30 {A.X(-36):.1f} {A.Y(weak*10**(-3.0)):.1f})">weak pulses: SNR = 388·ρ</text>')
o.append(f'<text class="sub" x="{A.X(24):.1f}" y="{A.Y(s[-1])+18:.1f}" text-anchor="end">self-noise: the pulsar</text>')
o.append(f'<text class="sub" x="{A.X(24):.1f}" y="{A.Y(s[-1])+32:.1f}" text-anchor="end">limits itself</text>')
# side notes
nx = 672
notes = [('ttl', 'At the design target'), ('sub', 'ρ = 10⁻⁵·⁴ = 4.0·10⁻⁶'), ('sub', 'One pulse: SNR 0.0015.'),
         ('sub', 'SNR grows as √(turns folded), so'), ('sub', 'one TOA at SNR 7 needs ~2·10⁷ turns:'),
         ('sub', '~62 h for this 10 ms pulsar.'), ('sub', ''), ('ttl', 'The test runs'),
         ('sub', '−5 dB is ~80,000× brighter:'), ('sub', 'every turn gives a usable TOA,'),
         ('sub', 'so a short run (0.2 s, 15–16 TOAs)'), ('sub', 'tests the whole chain.')]
for i, (c, line) in enumerate(notes):
    o.append(f'<text class="{c}" x="{nx}" y="{58 + 17*i}">{line}</text>')
svg = '''<figure class="scroll">
<svg class="wide" viewBox="0 0 860 330" role="img" aria-label="SNR of one pulse against snrDB: linear for weak pulses, limited by self-noise for strong ones; the design target and the test runs marked">
  ''' + '\n  '.join(o) + '''
</svg>
<figcaption>What a value of <code>snrDB</code> means, from the radiometer equation that <code>addNoiseAndRFI</code> prints (flat 400 MHz band,
the 0.5 ms Gaussian pulse; ~1.5 % optimistic). For weak pulses the SNR of one pulse is 388·ρ, a straight line here (dashed). At high
<code>snrDB</code> the pulsar's own noise (self-noise) takes over and more brightness hardly helps. The time per TOA (side note) uses the
fold's 373·ρ per turn, which includes the loss from estimating the baseline; the 62 h is the <a href="scenario.html">scenario page</a>'s
number for this pulsar.</figcaption>
</figure>'''
open('snrscale.svg', 'w').write(svg)
