import numpy as np
from plotsvg import Ax
# makeDispersionKernel at the reference settings (DM 5, 1.2-1.6 GHz, f_ref = 1.6 GHz, 4 GHz):
# a: which frequency the kernel carries at each time (a descending chirp);
# b: the kernel inside its design grid, and where the wrapped ringing goes.
K, DM, fs = 4.148808e3, 5.0, 4e9
KDM = K*1e12*DM
fL, fH = 1.2e9, 1.6e9
edge = 0.02*(fH - fL); g = int(np.ceil(20/edge*fs))
s0 = g; post = KDM/fL**2 - KDM/fH**2
Nk = s0 + int(np.ceil(post*fs)) + g + 1
Nf = 2**int(np.ceil(np.log2(2*Nk)))
print('g', g, 'end of sweep', s0 + int(np.ceil(post*fs)), 'Nk', Nk, 'Nf', Nf, 'empty', Nf - Nk, 'Nk/Nf', Nk/Nf)
tk = np.linspace(0, post, 400)                       # s after the zero lag s0
fk = 1/np.sqrt(tk/KDM + 1/fH**2)
print('cycles in the sweep', np.trapezoid(fk, tk))
o = []
# panel a: frequency vs kernel time
A = Ax(90, 120, 330, 130, (0, 6.6), (1150, 1650))
o.append(A.frame(xticks=[0, 2, 4, 6], yticks=[1200, 1400, 1600], xlab='time inside the kernel [ms]', ylab='frequency [MHz]'))
o.append(f'<path class="curve2" d="{A.path(tk*1e3, fk/1e6)}"/>')
o.append(f'<path class="dot2" d="M{A.X(0)-4:.1f} {A.Y(1600):.1f} a4 4 0 1 0 8 0 a4 4 0 1 0 -8 0"/>')
o.append(f'<path class="dot2" d="M{A.X(post*1e3)-4:.1f} {A.Y(1200):.1f} a4 4 0 1 0 8 0 a4 4 0 1 0 -8 0"/>')
o.append(f'<text class="sub" x="{A.X(0.3):.1f}" y="{A.Y(1625):.1f}">← 1.6 GHz first (sample s0 = 10,000)</text>')
o.append(f'<text class="sub" x="{A.X(0.3):.1f}" y="{A.Y(1225):.1f}">1.2 GHz last, 6.302 ms later →</text>')
# a toy chirp above it (schematic): fast oscillation first, slower later, like the real kernel
tt = np.linspace(0, post, 3000); ft = 1/np.sqrt(tt/KDM + 1/fH**2)
ph = np.cumsum(ft)*(tt[1] - tt[0]); ph = ph/ph[-1]*34      # 34 schematic cycles
S = Ax(90, 46, 330, 40, (0, 6.6), (-1.2, 1.2))
o.append(f'<path class="thin" d="{S.path(tt*1e3, np.cos(2*np.pi*ph))}"/>')
o.append(f'<text class="ttl" x="90" y="22">a · the kernel is a descending chirp</text>')
o.append(f'<text class="sub" x="90" y="38">schematic: the real kernel runs at 1.2–1.6 GHz (~9 M cycles)</text>')
# panel b: the design grid
x0, w, yb, hb = 470, 370, 92, 34
sc = w/Nf
xk = x0 + Nk*sc
o.append(f'<text class="ttl" x="{x0}" y="22">b · why nothing wraps: the design grid</text>')
o.append(f'<text class="sub" x="{x0}" y="38">one period of the design grid: N_f = 2²⁶ = 67.1 M (to scale)</text>')
o.append(f'<rect class="box sky" x="{x0}" y="{yb}" width="{xk-x0:.1f}" height="{hb}"/>')
o.append(f'<rect class="box" x="{xk:.1f}" y="{yb}" width="{x0+w-xk:.1f}" height="{hb}"/>')
o.append(f'<text class="lbl" x="{(x0+xk)/2:.1f}" y="{yb+15}" text-anchor="middle">kernel, kept</text>')
o.append(f'<text class="sub" x="{(x0+xk)/2:.1f}" y="{yb+29}" text-anchor="middle">N_k = 25.2 M (37.6 %)</text>')
o.append(f'<text class="lbl" x="{(xk+x0+w)/2:.1f}" y="{yb+15}" text-anchor="middle">empty zone, cut</text>')
o.append(f'<text class="sub" x="{(xk+x0+w)/2:.1f}" y="{yb+29}" text-anchor="middle">41.9 M samples</text>')
# ringing: late tail just after N_k, wrapped pre-ringing at the end of the grid (schematic)
r = np.linspace(0, 1, 60)
tail = xk + 3 + r*34; ytail = yb + hb + 14 - 8*np.exp(-4*r)*np.cos(2*np.pi*5*r)
o.append('<path class="curve4" stroke-width="1.3" d="M' + ' L'.join(f'{a:.1f} {b:.1f}' for a, b in zip(tail, ytail)) + '"/>')
pre = x0 + w - 3 - r*34; ypre = yb + hb + 14 - 8*np.exp(-4*r)*np.cos(2*np.pi*5*r)
o.append('<path class="curve4" stroke-width="1.3" d="M' + ' L'.join(f'{a:.1f} {b:.1f}' for a, b in zip(pre, ypre)) + '"/>')
o.append(f'<text class="sub" x="{xk+2:.1f}" y="{yb+hb+40}">late ringing</text>')
o.append(f'<text class="sub" x="{x0+w:.1f}" y="{yb+hb+40}" text-anchor="end">ringing from before sample 0</text>')
# wrap arrow: from the left edge round to the right end
ya = yb - 10
o.append(f'<path class="line" marker-end="url(#ar-kn)" d="M{x0+4} {yb-2} C{x0+4} {yb-30} {x0+w-4} {yb-30} {x0+w-6} {yb-2}"/>')
o.append(f'<text class="sub" x="{x0+w/2:.1f}" y="{yb-28:.1f}" text-anchor="middle">periodic grid: what falls before 0 reappears at the end</text>')
# the three facts
fy = yb + hb + 70
for i, line in enumerate([
        'The true response lives in [0, N_k). With N_f ≥ 2·N_k,',
        'both ringing tails land in the empty zone, never on',
        'the kernel, and the cut removes them.',
        'leak = energy cut = 3.5·10⁻¹³ (warning above 10⁻⁸).',
        'Frequency view: H changes by at most N_k/N_f = 0.376',
        'turn of phase between neighbouring bins (< ½).']):
    o.append(f'<text class="sub" x="{x0}" y="{fy + 16*i + (8 if i >= 4 else 0)}">{line}</text>')
svg = '''<figure class="scroll">
<svg class="wide" viewBox="0 0 860 302" role="img" aria-label="The dispersion kernel: a descending chirp from 1.6 to 1.2 GHz, and its place in the design grid with an empty zone where wrapped ringing is cut">
  <defs>
    <marker id="ar-kn" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
      <path class="ah" d="M0 0 L10 5 L0 10 z"/>
    </marker>
  </defs>
  ''' + '\n  '.join(o) + '''
</svg>
<figcaption><code>makeDispersionKernel</code> at the reference settings. a: the kernel's impulse response is a chirp. High frequencies come out
first, low frequencies 6.302 ms later; that delay difference <em>is</em> dispersion. b: the kernel is designed by sampling <var>H</var>(<var>f</var>) on
<var>N</var><sub>f</sub> frequency bins. That makes its time response periodic with period <var>N</var><sub>f</sub>, so anything outside the kernel
"wraps around". Making the grid at least twice the kernel leaves an empty zone for the wrapped ringing (schematic, magnified), which the
truncation to <var>N</var><sub>k</sub> taps then removes.</figcaption>
</figure>'''
open('kernel.svg', 'w').write(svg)
