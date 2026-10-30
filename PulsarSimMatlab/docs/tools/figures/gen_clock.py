import numpy as np
from plotsvg import Ax, pow10
# The receiver clock (receiver.html#clock): makeClock's "crystal" preset, h = [1.8e-21 6.492e-22 1.37e-24] (functions/makeClock.m), and its
# Allan deviation from the code's own formula (clk.sigmaY; Van Dierendonck, McGraw & Brown 1984, eqs. 6-8):
#   sigma_y^2(tau) = h0 / (2 tau) + 2 ln2 h_1 + (2 pi)^2 / 6 h_2 tau.
# b: tau * sigma_y(tau), the typical time error built up over tau, against the oscillator classes on the clock page (a TCXO ~3.6 us
# per hour, an OCXO 4-36 ns per hour: logs/2026-10-20_blockC.md, 21 Oct) and the two pulse widths (B0329+54 W50 6.6 ms, J0437 0.141 ms).
h0, h1, h2 = 1.8e-21, 6.492e-22, 1.37e-24
white = lambda t: np.sqrt(h0 / (2 * t))
flick = lambda t: np.sqrt(2 * np.log(2) * h1) * np.ones_like(t)
rwalk = lambda t: np.sqrt((2 * np.pi)**2 / 6 * h2 * t)
sig = lambda t: np.sqrt(white(t)**2 + flick(t)**2 + rwalk(t)**2)
H, DAY = 3600.0, 86400.0
o = []
# a: sigma_y(tau)
A = Ax(90, 50, 320, 220, (0, np.log10(10 * DAY)), (1e-12, 1e-8), logy=True)
tt = np.logspace(0, np.log10(10 * DAY), 300); lt = np.log10(tt)
xt = [(0, '1 s'), (np.log10(60), '1 min'), (np.log10(H), '1 h'), (np.log10(DAY), '1 d'), (np.log10(10 * DAY), '10 d')]
for e in (-11, -10, -9):
    o.append(f'<path class="grid" d="M{A.x0} {float(A.Y(10.0**e)):.1f} H{A.x0 + A.w}"/>')
o.append(A.frame(xticks=[], yticks=[1e-12, 1e-11, 1e-10, 1e-9, 1e-8], yfmt=pow10, ylab='Allan deviation σ<tspan baseline-shift="sub" font-size="9">y</tspan>(τ)'))
for v, lab in xt:
    x = float(A.X(v)); o.append(f'<path class="axis" d="M{x:.1f} {A.y0 + A.h} v5"/>')
    o.append(f'<text class="sub" x="{x:.1f}" y="{A.y0 + A.h + 18}" text-anchor="middle">{lab}</text>')
o.append(f'<text class="sub" x="{A.x0 + A.w / 2:.1f}" y="{A.y0 + A.h + 34}" text-anchor="middle">averaging time τ</text>')
for f, lab, tl, dy in ((white, 'white', 2.0, -7), (flick, 'flicker', 3.4, 14), (rwalk, 'random walk', 4.3, 18)):
    o.append(f'<path class="dash" d="{A.path(lt, f(tt))}"/>')
    y = float(A.Y(f(np.array([10**tl]))[0]))
    o.append(f'<text class="sub" x="{float(A.X(tl)):.1f}" y="{y + dy:.1f}" text-anchor="middle">{lab}</text>')
o.append(f'<path class="s1" d="{A.path(lt, sig(tt))}"><title>the crystal preset, all three parts</title></path>')
o.append(f'<text class="sub" x="{float(A.X(0.9)):.1f}" y="{float(A.Y(sig(np.array([10**0.9]))[0])) - 10:.1f}" style="fill:var(--series-1)">crystal preset</text>')
o.append('<text class="ttl" x="90" y="24">a · how stable the frequency is</text>')
# b: tau * sigma_y(tau)
B = Ax(530, 50, 290, 220, (np.log10(60), np.log10(10 * DAY)), (1e-9, 1e-1), logy=True)
tb = np.logspace(np.log10(60), np.log10(10 * DAY), 300)
for e in (-8, -7, -6, -5, -4, -3, -2):
    o.append(f'<path class="grid" d="M{B.x0} {float(B.Y(10.0**e)):.1f} H{B.x0 + B.w}"/>')
o.append(B.frame(xticks=[], yticks=[1e-9, 1e-6, 1e-3], yfmt=lambda v: {1e-9: '1 ns', 1e-6: '1 µs', 1e-3: '1 ms'}[v],
                 ylab='time error τ·σ<tspan baseline-shift="sub" font-size="9">y</tspan>(τ)'))
for v, lab in xt[1:]:
    x = float(B.X(v)); o.append(f'<path class="axis" d="M{x:.1f} {B.y0 + B.h} v5"/>')
    o.append(f'<text class="sub" x="{x:.1f}" y="{B.y0 + B.h + 18}" text-anchor="middle">{lab}</text>')
o.append(f'<text class="sub" x="{B.x0 + B.w / 2:.1f}" y="{B.y0 + B.h + 34}" text-anchor="middle">time since the clock was last known, τ</text>')
for w, lab in ((6.6e-3, 'B0329+54 pulse, 6.6 ms'), (0.141e-3, 'J0437 pulse, 0.14 ms')):
    y = float(B.Y(w))
    o.append(f'<path class="dash" d="M{B.x0} {y:.1f} H{B.x0 + B.w}"/>')
    o.append(f'<text class="sub" x="{B.x0 + 6}" y="{y - 5:.1f}">{lab}</text>')
o.append(f'<path class="s1" d="{B.path(np.log10(tb), tb * sig(tb))}"><title>the crystal preset</title></path>')
for t, lab in ((H, '0.66 µs'), (DAY, '76 µs')):
    v = t * sig(np.array([t]))[0]
    o.append(f'<circle class="s1f" cx="{float(B.X(np.log10(t))):.1f}" cy="{float(B.Y(v)):.1f}" r="4.5"><title>crystal preset: {lab} after {"1 h" if t == H else "1 day"}</title></circle>')
    o.append(f'<text class="sub" x="{float(B.X(np.log10(t))) + 8:.1f}" y="{float(B.Y(v)) + 14:.1f}" style="fill:var(--series-1)">{lab}</text>')
x = float(B.X(np.log10(H)))
o.append(f'<circle class="s2f" cx="{x:.1f}" cy="{float(B.Y(3.6e-6)):.1f}" r="4.5"><title>TCXO estimate: 3.6 µs per hour</title></circle>')
o.append(f'<text class="sub" x="{x - 8:.1f}" y="{float(B.Y(3.6e-6)) - 6:.1f}" text-anchor="end" style="fill:var(--series-2)">TCXO, 3.6 µs</text>')
o.append(f'<path class="s3" style="stroke-width:5" stroke-linecap="round" d="M{x:.1f} {float(B.Y(4e-9)):.1f} V{float(B.Y(36e-9)):.1f}"><title>OCXO estimate: 4-36 ns per hour</title></path>')
o.append(f'<text class="sub" x="{x + 8:.1f}" y="{float(B.Y(12e-9)) + 4:.1f}" style="fill:var(--series-3)">OCXO, 4–36 ns</text>')
o.append('<text class="ttl" x="530" y="24">b · how far the clock wanders</text>')
svg = '''<figure class="scroll">
<svg class="wide" viewBox="0 0 860 330" role="img" aria-label="The crystal clock preset: its Allan deviation against averaging time with its white, flicker and random-walk parts, and the time error it builds up, against TCXO and OCXO estimates and the pulse widths">
  ''' + '\n  '.join(o) + '''
</svg>
<figcaption><b>Figure 4a.1.</b> The receiver clock's <code>"crystal"</code> preset, from <code>makeClock</code>'s own formula. (a) Its Allan deviation
σ<sub><var>y</var></sub>(τ) and its three noises. Flicker is the largest from about 1 s to 100 s; beyond that the random walk takes over.
(b) The time error this builds up, τ·σ<sub><var>y</var></sub>(τ): 0.66 µs after an hour and 76 µs after a day. It reaches J0437's 0.14 ms
pulse after about a day and a half, and B0329+54's 6.6 ms only after about 20 days. For comparison, the TCXO (3.6 µs per hour) and OCXO
(4–36 ns per hour) estimates. Generated by <code>docs/tools/figures/gen_clock.py</code>.</figcaption>
</figure>'''
open('clock.svg', 'w').write(svg)
print('1 h', H * sig(np.array([H]))[0], '1 d', DAY * sig(np.array([DAY]))[0], 'crossovers', h0 / (4 * np.log(2) * h1), 2 * np.log(2) * h1 / ((2 * np.pi)**2 / 6 * h2))
