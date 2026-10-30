import numpy as np
from plotsvg import Ax
from pfb import design, response, fs, K, dF
# Folding keeps every 16th bin: the spectrum of one weighted frame on the 4096-point grid with three tones,
# the 256-point FFT's bins marked; and how alike two outputs of one frame are against their spacing.
fLO = 1300.0
h, n, w, I = design(192)
L = I['P']*K                                  # 4096
fine = fs/L                                   # 0.195 MHz
tones = [fLO + 1.5*dF, fLO + 4.5*dF + 1.0, fLO + 7*dF]     # ch 34 centre; 1 MHz above ch 37's centre; edge ch 39/40
nf = np.arange(0, 160 + 1)                    # fine bins 0..160 = 0..31.25 MHz above the LO
ff = fLO + nf*fine
S = np.abs(sum(response(h, n, ff - t) for t in tones))
fc_ = np.arange(fLO, fLO + 31.2501, 0.04)
Sc = np.abs(sum(response(h, n, fc_ - t) for t in tones))
kept = nf % 16 == 8                           # channel centres (k + 1/2)*3.125 MHz
for t in tones:
    print(t, [(33 + k, round(float(np.abs(response(h, n, fLO + (k + .5)*dF - t)[0])), 3)) for k in range(10)
              if abs(response(h, n, fLO + (k + .5)*dF - t)[0]) > 1e-3])

o = []
A = Ax(80, 64, 750, 210, (fLO, fLO + 31.25), (0, 1.5))
o.append(f'<text class="ttl" x="80" y="22">a · One weighted frame on the 4096-point grid (0.195 MHz) and the 256-point FFT of the folded frame</text>')
o.append(A.frame(xticks=[1300, 1305, 1310, 1315, 1320, 1325, 1330], yticks=[0, 0.5, 1],
                 xlab='RF frequency [MHz]', ylab='|spectrum| of the weighted frame'))
o.append(f'<path class="grid" d="' + ' '.join(f'M{A.X(fLO + k*dF):.1f} {A.y0} V{A.y0 + A.h}' for k in range(1, 10)) + '"/>')
for k in range(10):
    o.append(f'<text class="sub" x="{A.X(fLO + (k + .5)*dF):.1f}" y="{A.Y(1.12):.1f}" text-anchor="middle">ch {33 + k}</text>')
o.append(f'<path class="thin" d="{A.path(fc_, Sc)}"/>')
for x, y, k in zip(ff, S, kept):
    if not k:
        o.append(f'<circle class="dot" cx="{A.X(x):.1f}" cy="{A.Y(y):.1f}" r="1.7" opacity="0.55"/>')
for x, y, k in zip(ff, S, kept):
    if k:
        o.append(f'<circle class="dot2" cx="{A.X(x):.1f}" cy="{A.Y(y):.1f}" r="4.6"/>')
lab = [('tone at the centre', 'of ch 34'), ('tone 1 MHz above', 'the centre of ch 37'), ('tone on the edge', 'of ch 39 and 40')]
for t, (l1, l2) in zip(tones, lab):
    o.append(f'<text class="sub" x="{A.X(t):.1f}" y="{A.Y(1.4):.1f}" text-anchor="middle">{l1}</text>')
    o.append(f'<text class="sub" x="{A.X(t):.1f}" y="{A.Y(1.4) + 13:.1f}" text-anchor="middle">{l2}</text>')
o.append(f'<circle class="dot" cx="96" cy="{A.y0 + A.h + 48}" r="2" opacity="0.55"/>')
o.append(f'<text class="sub" x="104" y="{A.y0 + A.h + 52}">fine bins of a 4096-point FFT (never computed)</text>')
o.append(f'<circle class="dot2" cx="396" cy="{A.y0 + A.h + 48}" r="4.6"/>')
o.append(f'<text class="sub" x="406" y="{A.y0 + A.h + 52}">every 16th: the 256-point FFT of the folded frame = the channel outputs</text>')

# (b) likeness of two outputs of one frame
B = Ax(80, 430, 750, 150, (0, 6), (0, 1.05))
d = np.linspace(0, 6, 241)
h2 = (h**2)/np.sum(h**2)
rho = np.abs(np.cos(2*np.pi*np.outer(d, n)/fs) @ h2)
r1 = np.abs(np.cos(2*np.pi*fine*n/fs) @ h2)
r16 = np.abs(np.cos(2*np.pi*dF*n/fs) @ h2)
print('rho', r1, r16)
o.append(f'<text class="ttl" x="80" y="392">b · How alike are two outputs of the same frame? Correlation against their spacing (white noise)</text>')
o.append(B.frame(xticks=[0, 1, 2, 3, 4, 5, 6], yticks=[0, 0.5, 1], xlab='spacing of the two bins [MHz]', ylab='correlation'))
o.append(f'<path class="curve" d="{B.path(d, rho)}"/>')
o.append(f'<circle class="dot2" cx="{B.X(fine):.1f}" cy="{B.Y(r1):.1f}" r="4"/>')
o.append(f'<text class="sub" x="{B.X(fine) + 9:.1f}" y="{B.Y(r1) - 6:.1f}">neighbouring fine bins (0.195 MHz): {r1:.3f}</text>')
o.append(f'<circle class="dot2" cx="{B.X(dF):.1f}" cy="{B.Y(r16):.1f}" r="4"/>')
o.append(f'<text class="sub" x="{B.X(dF) - 10:.1f}" y="{B.Y(r16) + 20:.1f}" text-anchor="end">neighbouring channels (3.125 MHz): {r16:.2f}</text>')
o.append(f'<text class="sub" x="{B.X(4.4):.1f}" y="{B.Y(0.62):.1f}">nearly independent beyond</text>')
o.append(f'<text class="sub" x="{B.X(4.4):.1f}" y="{B.Y(0.62) + 14:.1f}">about 4 MHz: one frame</text>')
o.append(f'<text class="sub" x="{B.X(4.4):.1f}" y="{B.Y(0.62) + 28:.1f}">resolves about 4 MHz</text>')

svg = '''<figure class="scroll">
<svg class="wide" viewBox="0 0 860 630" role="img" aria-label="The spectrum of one weighted frame with three tones on the fine grid, with every 16th bin marked as a channel output; and the correlation of two outputs against their spacing">
  ''' + '\n  '.join(o) + '''
</svg>
<figcaption><b>Figure 5a.3.</b> Folding keeps every 16th bin. (a) A weighted frame with three tones. Each tone becomes a 4 MHz
block. A 4096-point DFT would sample it every 0.195 MHz (small dots), and the folded 256-point DFT gives every 16th sample (large dots).
(b) The correlation of two outputs of one frame against their spacing.</figcaption>
</figure>'''
open('pfbfold.svg', 'w').write(svg)
