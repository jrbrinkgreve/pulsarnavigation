import numpy as np
from plotsvg import Ax
# Why square first (physics.html, section 2): a: two turns of a toy voltage x = sqrt(1 + rho*G)*n (receiver noise 1,
# pulsar rho = 2 so it shows in one turn; 300 samples per turn, power profile G with 5 % FWHM). b: the average voltage over
# 500 turns at each phase: zero, no pulse. c: the average power x^2 over the same turns: 1 + rho*G, the pulse.
rng = np.random.default_rng(7)
ns, nt, rho = 300, 500, 2.0
ph = (np.arange(ns) + 0.5)/ns
sig = 0.05/(2*np.sqrt(2*np.log(2)))
G = sum(np.exp(-0.5*((ph - 0.5 - k)/sig)**2) for k in (-1, 0, 1))
x = np.sqrt(1 + rho*G)[None, :]*rng.standard_normal((nt, ns))
o = []
A = Ax(60, 40, 770, 120, (0, 2), (-4.6, 4.6))
o.append(A.frame(xticks=[0, 0.5, 1, 1.5, 2], yticks=[-4, 0, 4], xlab='', ylab='voltage x', xfmt='{:g}'))
tt = np.concatenate([ph, ph + 1])
xx = np.concatenate([x[0], x[1]])
env = np.sqrt(1 + rho*np.concatenate([G, G]))
o.append(f'<path class="curve2" d="{A.path(tt, env)}" opacity="0.8"/>')
o.append(f'<path class="curve2" d="{A.path(tt, -env)}" opacity="0.8"/>')
o.append(f'<path d="{A.path(tt, xx)}" stroke="var(--accent)" stroke-width="0.9" fill="none"/>')
o.append(f'<text class="sub" x="{A.x0 + A.w/2:.1f}" y="{A.y0 + A.h + 34}" text-anchor="middle">time [turns]</text>')
o.append(f'<text class="ttl" x="60" y="20">a · two turns of the voltage: noise that gets louder at the pulse (orange: ±√(1 + ρG), its standard deviation)</text>')
for i, (lab, y, cls, ttl) in enumerate([
        ('average voltage', x.mean(0), 'curve', 'b · average voltage over 500 turns'),
        ('average power x²', (x**2).mean(0), 'curve3', 'c · average power over the same 500 turns')]):
    B = Ax(60 + i*420, 250, 350, 140, (0, 1), (-0.5, 3.5))
    o.append(B.frame(xticks=[0, 0.25, 0.5, 0.75, 1], yticks=[0, 1, 2, 3], xlab='rotational phase [turns]', ylab=lab))
    o.append(f'<path class="dash" d="M{B.x0} {float(B.Y(0)):.1f} H{B.x0 + B.w}"/>')
    o.append(f'<path class="{cls}" d="{B.path(ph, y)}"/>')
    o.append(f'<text class="ttl" x="{B.x0}" y="218">{ttl}</text>')
    if i == 0:
        o.append(f'<text class="sub" x="{B.x0 + B.w/2:.1f}" y="{float(B.Y(0.75)):.1f}" text-anchor="middle">zero everywhere: the pulse is gone</text>')
    else:
        o.append(f'<text class="sub" x="{float(B.X(0.56)):.1f}" y="{float(B.Y(2.9)):.1f}">1 + ρ·G: the pulse</text>')
        o.append(f'<text class="sub" x="{float(B.X(0.03)):.1f}" y="{float(B.Y(1.25)):.1f}">receiver noise level</text>')
svg = '''<figure class="scroll">
<svg class="wide" viewBox="0 0 860 440" role="img" aria-label="Why the processing squares the voltage: the average voltage over many turns is zero, the average power shows the pulse">
  ''' + '\n  '.join(o) + '''
</svg>
<figcaption>Why the processing squares the voltage (toy numbers: 300 samples per turn, a pulsar twice as strong as the receiver noise
at its peak so it shows in one turn). Every turn is new noise; only its loudness repeats. Averaging the voltage over many turns gives
zero, pulse or not (b). Averaging the power gives the profile (c). The real signal is ~50 dB weaker, which is why it takes millions of
turns instead of 500.</figcaption>
</figure>'''
open('whysquare.svg', 'w').write(svg)
