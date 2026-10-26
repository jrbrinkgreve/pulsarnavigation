import numpy as np
from plotsvg import Ax
# What makeClock draws (receiver.html#clock-how): a port of its three noise generators on the 60 s grid (functions/makeClock.m), with the
# "crystal" preset h = [1.8e-21 6.492e-22 1.37e-24]. White frequency noise: independent steps of the time error, variance h0 D/2.
# Flicker: the step-averaged frequency by FFT synthesis with the spectrum h_1/f over twice the span, integrated. Random walk: the
# exact two-state steps, q2 = 2 pi^2 h_2. numpy's random numbers, so these are other draws than makeClock's for the same seed.
# a: one draw of each noise over a day, each on its own scale. b: eight draws of the whole preset over 30 days, against
# +-tau sigma_y(tau) from the preset's formula (gen_clock.py).
h0, h1, h2 = 1.8e-21, 6.492e-22, 1.37e-24
Dt = 60.0; DAY = 86400.0

def white(K, rng):
    return np.r_[0, np.cumsum(np.sqrt(h0 * Dt / 2) * rng.standard_normal(K))]

def flicker(K, rng):
    N = int(2**np.ceil(np.log2(2 * K)))
    f = np.arange(1, N // 2 + 1) / (N * Dt); A = np.sqrt(N * (h1 / (2 * f)) / Dt)
    Y = np.zeros(N, complex)
    Y[1:N // 2] = A[:-1] * (rng.standard_normal(N // 2 - 1) + 1j * rng.standard_normal(N // 2 - 1)) / np.sqrt(2)
    Y[N // 2] = A[-1] * rng.standard_normal()
    Y[N // 2 + 1:] = np.conj(Y[N // 2 - 1:0:-1])
    yF = np.real(np.fft.ifft(Y))
    return np.r_[0, np.cumsum(Dt * yF[:K])]

def rwalk(K, rng):
    q2 = 2 * np.pi**2 * h2
    L = np.linalg.cholesky(q2 * np.array([[Dt**3 / 3, Dt**2 / 2], [Dt**2 / 2, Dt]]))
    w = L @ rng.standard_normal((2, K)); x = np.zeros(K + 1); y = 0.0
    for k in range(K):
        x[k + 1] = x[k] + Dt * y + w[0, k]; y += w[1, k]
    return x

sig = lambda t: np.sqrt(h0 / (2 * t) + 2 * np.log(2) * h1 + (2 * np.pi)**2 / 6 * h2 * t)
o = []
# a: one day of each noise
K1 = int(DAY / Dt); t1 = np.arange(K1 + 1) * Dt / 3600
rng = np.random.default_rng(7)
strips = [('white', white(K1, rng), 'ns', 1e9, 'independent steps'), ('flicker', flicker(K1, rng), 'µs', 1e6, 'in between'),
          ('random walk', rwalk(K1, rng), 'µs', 1e6, 'the frequency wanders')]
for i, (name, x, unit, sc, note) in enumerate(strips):
    v = x * sc; m = np.max(np.abs(v)) * 1.15
    A = Ax(90, 52 + i * 92, 300, 64, (0, 24), (-m, m))
    o.append(f'<path class="grid" d="M{A.x0} {float(A.Y(0)):.1f} H{A.x0 + A.w}"/>')
    tk = float(np.round(m / 1.15, 1 if m < 10 else 0)) or round(m / 1.15, 2)
    o.append(A.frame(xticks=[0, 6, 12, 18, 24] if i == 2 else [], yticks=[-tk, 0, tk], yfmt=lambda u: f'{u:g}', xlab='hours' if i == 2 else ''))
    o.append(f'<text class="sub" x="{A.x0 - 8}" y="{A.y0 - 6}" text-anchor="end">{unit}</text>')
    o.append(f'<path class="s1" style="stroke-width:1.4" d="{A.path(t1, v)}"/>')
    o.append(f'<text class="lbl" x="{A.x0 + A.w + 12}" y="{A.y0 + 28}">{name}</text>')
    o.append(f'<text class="sub" x="{A.x0 + A.w + 12}" y="{A.y0 + 44}">{note}</text>')
o.append('<text class="ttl" x="90" y="24">a · one draw of each noise over a day, each on its own scale</text>')
# b: the preset over 30 days, eight draws
K30 = int(30 * DAY / Dt); t30 = np.arange(K30 + 1) * Dt
B = Ax(560, 52, 270, 248, (0, 30), (-30, 30))
for v in (-20, -10, 10, 20):
    o.append(f'<path class="grid" d="M{B.x0} {float(B.Y(v)):.1f} H{B.x0 + B.w}"/>')
o.append(B.frame(xticks=[0, 10, 20, 30], yticks=[-20, -10, 0, 10, 20], xlab='days', ylab='time error [ms]'))
env = t30[1:] * sig(t30[1:]) * 1e3
o.append(f'<path class="dash" d="{B.path(np.r_[0, t30[1:] / DAY], np.r_[0, env])}"/>')
o.append(f'<path class="dash" d="{B.path(np.r_[0, t30[1:] / DAY], np.r_[0, -env])}"/>')
ends = []
for s in range(8):
    r = np.random.default_rng(100 + s)
    x = (white(K30, r) + flicker(K30, r) + rwalk(K30, r)) * 1e3
    ends.append(x[-1])
    o.append(f'<path class="s1" style="stroke-width:1.2;opacity:0.75" d="{B.path(t30[::30] / DAY, x[::30])}"/>')
e30 = 30 * DAY * sig(np.array([30 * DAY]))[0] * 1e3
o.append(f'<text class="sub" x="{B.x0 + 6}" y="{float(B.Y(26)):.1f}">dashed: ±τ·σ<tspan baseline-shift="sub" font-size="9">y</tspan>(τ), {e30:.0f} ms at 30 days</text>')
o.append('<text class="ttl" x="530" y="24">b · eight crystal clocks over 30 days</text>')

svg = '''<figure class="scroll">
<svg class="wide" viewBox="0 0 860 340" role="img" aria-label="What makeClock draws: one day of white, flicker and random-walk frequency noise as time error, each on its own scale; and eight draws of the crystal preset over 30 days, spreading to about 13 ms, inside the envelope from the Allan deviation">
  ''' + '\n  '.join(o) + '''
</svg>
<figcaption>What <code>makeClock</code> draws for the <code>"crystal"</code> preset, with a port of its three generators. a: the time error that
each noise builds up over a day. White frequency noise gives independent steps, so the time error is rough, but it stays at nanoseconds.
Flicker is in between. The random walk is smooth, because the frequency itself wanders, and it is the largest: tens of microseconds after a
day. b: eight clocks with all three noises over 30 days. The dashed lines are their typical size, τ·σ<sub><var>y</var></sub>(τ) from the
figure above: ''' + f'{e30:.0f}' + ''' ms after 30 days. A typical size is one standard deviation, so some draws end outside it. That is why the clock must be part of the navigation fit. An SDR's rate of up to
2·10⁻⁶ is not shown: it alone adds 5.2 s in 30 days, so it must be estimated first. The draws use Python's random numbers, so they are
other draws than <code>makeClock</code>'s. Generated by <code>docs/tools/figures/gen_clockdraw.py</code>.</figcaption>
</figure>'''
open('clockdraw.svg', 'w').write(svg)
print('1-day max |x|:', [f'{n} {np.max(np.abs(x)):.3g}' for n, x, *_ in strips], '30 d envelope ms', e30, 'ends ms', np.round(ends, 1))
