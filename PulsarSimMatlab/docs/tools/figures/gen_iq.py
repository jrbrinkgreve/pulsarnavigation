from plotsvg import Ax
o = []
def band(ax, f1, f2, cls='fillsoft', h=0.8, extra=''):
    x1, x2 = float(ax.X(f1)), float(ax.X(f2)); yb = float(ax.Y(0)); yt = float(ax.Y(h))
    return f'<rect class="{cls}" x="{x1:.1f}" y="{yt:.1f}" width="{x2-x1:.1f}" height="{yb-yt:.1f}" {extra}/>' + \
           f'<path class="curve" d="M{x1:.1f} {yb:.1f} V{yt:.1f} H{x2:.1f} V{yb:.1f}" {extra}/>'
def label(ax, f, y, txt, cls='sub', anchor='middle'):
    return f'<text class="{cls}" x="{float(ax.X(f)):.1f}" y="{y:.1f}" text-anchor="{anchor}">{txt}</text>'
# panel 1: RF
A = Ax(60, 40, 760, 60, (-2.0, 2.0), (0, 1))
o.append(f'<text class="ttl" x="60" y="22">1 · RF: real samples at f_in = 4 GHz</text>')
o.append(A.frame(xticks=[-2, -1.6, -1.2, 0, 1.2, 1.6, 2], xfmt='{:g}', left=False))
o.append(f'<text class="sub" x="826" y="{A.Y(0)+4:.1f}">GHz</text>')
o.append(band(A, 1.2, 1.6)); o.append(band(A, -1.6, -1.2))
o.append(label(A, 1.17, float(A.Y(0.45)), 'the band', anchor='end'))
o.append(label(A, -1.4, float(A.Y(0.8)) - 6, 'its mirror (a real signal)'))
xl = float(A.X(1.3)); o.append(f'<path class="curve4" d="M{xl:.1f} {A.Y(0):.1f} V{A.Y(1.12):.1f}"/>')
o.append(f'<text class="sub" x="{xl+6:.1f}" y="{A.Y(1.12)+8:.1f}">f_LO = 1.3</text>')
# panel 2: after mixing, before decimation
B = Ax(60, 172, 760, 60, (-2.0, 2.0), (0, 1))
o.append(f'<text class="ttl" x="60" y="150">2 · × 2·e^(−i2π f_LO t): every frequency moves down by 1.3 GHz, then the low-pass</text>')
o.append(B.frame(xticks=[-2, -0.36, 0.36, 1.1, 1.5, 2], xfmt='{:g}', left=False))
o.append(f'<text class="sub" x="826" y="{B.Y(0)+4:.1f}">GHz</text>')
x1, x2 = float(B.X(-0.36)), float(B.X(0.36))
o.append(f'<path class="curve3" stroke-dasharray="5 3" d="M{float(B.X(-0.37)):.1f} {B.Y(0):.1f} L{x1:.1f} {B.Y(1):.1f} H{x2:.1f} L{float(B.X(0.37)):.1f} {B.Y(0):.1f}"/>')
o.append(label(B, -0.4, float(B.Y(0.75)), 'low-pass', anchor='end'))
o.append(label(B, -0.4, float(B.Y(0.75)) + 15, '±360 MHz', anchor='end'))
o.append(band(B, -0.1, 0.3))
o.append(band(B, 1.1, 1.5, cls='fillfg', extra='opacity="0.5"'))
o.append(label(B, 1.3, float(B.Y(0.8)) - 6, 'image: −2.9…−2.5 GHz ≡ +1.1…+1.5'))
o.append(label(B, 1.3, float(B.Y(0)) + 34, '(4 GHz sampling) → removed'))
o.append(label(B, 0.1, float(B.Y(0.5)) + 4, '−0.1…+0.3', cls='lbl'))
# panel 3: complex baseband at 800 MHz
C = Ax(60, 322, 760, 60, (-0.45, 0.45), (0, 1))
o.append(f'<text class="ttl" x="60" y="290">3 · keep every 5th sample: complex baseband at f_s = 800 MHz (zoomed)</text>')
o.append(C.frame(xticks=[-0.4, -0.3, -0.1, 0, 0.1, 0.3, 0.4], xfmt='{:g}', left=False))
o.append(f'<text class="sub" x="826" y="{C.Y(0)+4:.1f}">GHz</text>')
for fx in (-0.4, 0.4):
    o.append(f'<path class="dash" d="M{float(C.X(fx)):.1f} {C.Y(0):.1f} V{C.Y(1.05):.1f}"/>')
o.append(f'<text class="sub" x="{float(C.X(-0.4))-6:.1f}" y="{float(C.Y(0.5)):.1f}" text-anchor="end">−f_s/2</text>')
o.append(f'<text class="sub" x="{float(C.X(0.4))+6:.1f}" y="{float(C.Y(0.5)):.1f}">+f_s/2</text>')
o.append(f'<path class="curve4" stroke-dasharray="6 4" d="M{float(C.X(-0.3)):.1f} {C.Y(0):.1f} V{C.Y(0.62):.1f} H{float(C.X(0.1)):.1f} V{C.Y(0):.1f}"/>')
o.append(band(C, -0.1, 0.3))
o.append(label(C, 0.1, float(C.Y(0.8)) - 6, 'the band: f_bb = f_RF − f_LO, −100…+300 MHz'))
o.append(label(C, -0.2, float(C.Y(0.62)) + 16, 'sign error:', cls='sub'))
o.append(label(C, -0.2, float(C.Y(0.62)) + 30, '−300…+100', cls='sub'))
svg = '''<figure class="scroll">
<svg class="wide" viewBox="0 0 860 410" role="img" aria-label="Spectra through the IQ stage: RF band, mixing with the LO, low-pass and decimation to complex baseband">
  ''' + '\n  '.join(o) + '''
</svg>
<figcaption>Where the band goes in <code>applyIQmodulation</code> (to scale; amplitudes schematic). 1: a real RF signal has its band at
±(1.2–1.6) GHz. 2: the complex LO shifts everything down by 1.3 GHz; the band lands at −0.1…+0.3 GHz and the mirror far away (at 4 GHz
sampling, −2.9…−2.5 GHz is the same as +1.1…+1.5 GHz), where the low-pass (cut-off 360 MHz, 2049 taps) removes it. 3: after keeping every
5th sample the band sits inside ±400 MHz, off-centre on purpose: a mirrored mapping would move it to −300…+100 MHz and fail visibly.</figcaption>
</figure>'''
open('iqspec.svg', 'w').write(svg)
