import numpy as np
from plotsvg import Ax
# Folding (physics.html, section 6). a: simulated folded profiles (128 phase bins, Gaussian pulse of 5 % FWHM) after 1, 16, 256
# and 4096 turns, at a matched-filter SNR of 0.25 per turn, shown in units of the per-bin noise: the noise stays, the pulse
# grows as sqrt(N). b: SNR against observing time for the simulated pulsar (10 ms period), SNR per turn ~373*rho (results.html,
# threshold sweep: flat band x baseline loss 0.962), with the design threshold SNR 7.
rng = np.random.default_rng(5)
nb = 128
ph = (np.arange(nb) + 0.5)/nb - 0.5
sb = 0.05*nb/(2*np.sqrt(2*np.log(2)))
g = np.exp(-0.5*(ph*nb/sb)**2)
a1 = 0.25/np.sqrt(np.sum(g**2))
o = []
o.append('<text class="ttl" x="60" y="20">a · the folded profile after N turns (SNR 0.25 per turn)</text>')
H = 64
for i, N in enumerate((1, 16, 256, 4096)):
    S = Ax(60, 46 + i*(H + 12), 330, H, (-0.5, 0.5), (-3.2, 8.2))
    prof = a1*np.sqrt(N)*g + rng.standard_normal(nb)
    o.append(f'<path class="axis" d="M{S.x0} {S.y0 + S.h} H{S.x0 + S.w}"/>')
    o.append(f'<path class="fillfg" d="M{S.x0} {float(S.Y(1)):.1f} H{S.x0 + S.w} V{float(S.Y(-1)):.1f} H{S.x0} Z"/>')
    o.append(f'<path class="curve" d="{S.path(ph, prof)}"/>')
    o.append(f'<path class="curve2" d="{S.path(ph, a1*np.sqrt(N)*g)}" opacity="0.8"/>')
    o.append(f'<text class="lbl" x="{S.x0 + S.w + 10}" y="{S.y0 + 26}">N = {N}</text>')
    o.append(f'<text class="sub" x="{S.x0 + S.w + 10}" y="{S.y0 + 42}">SNR {0.25*np.sqrt(N):g}</text>')
o.append(f'<text class="sub" x="{60 + 165}" y="{46 + 4*(H + 12) + 4}" text-anchor="middle">rotational phase (grey band: ±1σ noise per bin; orange: the true pulse)</text>')
# ---------- b ----------
B = Ax(530, 46, 300, 250, (np.log10(1e-2), np.log10(1e6)), (1e-3, 1e3), logy=True)
lx = lambda t: B.x0 + (np.log10(t) + 2)/8*B.w
o.append(f'<path class="axis" d="M{B.x0} {B.y0 + B.h} H{B.x0 + B.w} M{B.x0} {B.y0} V{B.y0 + B.h}"/>')
for t, lab in ((1e-2, '10 ms'), (1, '1 s'), (60, '1 min'), (3600, '1 h'), (86400, '1 d')):
    o.append(f'<path class="axis" d="M{lx(t):.1f} {B.y0 + B.h} v5"/><text class="sub" x="{lx(t):.1f}" y="{B.y0 + B.h + 18}" text-anchor="middle">{lab}</text>')
for v, lab in ((1e-3, '0.001'), (1e-1, '0.1'), (10, '10'), (1e3, '1000')):
    o.append(f'<path class="axis" d="M{B.x0} {float(B.Y(v)):.1f} h-5"/><text class="sub" x="{B.x0 - 8}" y="{float(B.Y(v)) + 4:.1f}" text-anchor="end">{lab}</text>')
o.append(f'<text class="sub" x="{B.x0 + B.w/2:.1f}" y="{B.y0 + B.h + 34}" text-anchor="middle">observing time</text>')
o.append(f'<text class="sub" x="{B.x0}" y="{B.y0 - 8}">SNR of the folded profile</text>')
yt = float(B.Y(7))
o.append(f'<path class="fillwhy" d="M{B.x0} {B.y0} H{B.x0 + B.w} V{yt:.1f} H{B.x0} Z"/>')
o.append(f'<path class="dash" d="M{B.x0} {yt:.1f} H{B.x0 + B.w}"/>')
o.append(f'<rect class="fillwhy" x="{B.x0 + 196}" y="{B.y0 + 161}" width="20" height="10"/><text class="sub" x="{B.x0 + 222}" y="{B.y0 + 170}">SNR ≥ 7: usable TOA</text>')
for rdb, cls in ((-5, 'curve'), (-20, 'curve3'), (-40, 'curve2'), (-54, 'curve4')):
    s1 = 373*10**(rdb/10)
    tend = min(1e6, 0.01*(1e3/s1)**2)
    tt = np.logspace(-2, np.log10(tend), 60)
    snr = s1*np.sqrt(tt/0.01)
    o.append(f'<path class="{cls}" d="M{" L".join(f"{lx(t):.1f} {float(B.Y(s)):.1f}" for t, s in zip(tt, snr))}"/>')
    tc = 0.01*(7/s1)**2
    ly = B.y0 + 186 + 15*[-5, -20, -40, -54].index(rdb)
    o.append(f'<path class="{cls}" d="M{B.x0 + 196} {ly - 4} h20"/><text class="sub" x="{B.x0 + 222}" y="{ly}">ρ = {rdb} dB'.replace('-', '−') + '</text>')
    if tc > 0.01:
        o.append(f'<path class="dot2" d="M{lx(tc) - 3.5:.1f} {yt:.1f} a3.5 3.5 0 1 0 7 0 a3.5 3.5 0 1 0 -7 0"/>')
        lab_t = f'{tc*1e3:.0f} ms' if tc < 1 else (f'{tc/60:.1f} min' if tc < 3600 else f'{tc/3600:.0f} h')
        o.append(f'<text class="sub" x="{lx(tc) + 6:.1f}" y="{yt + 16:.1f}">{lab_t}</text>')
o.append(f'<text class="ttl" x="530" y="20">b · SNR grows as √time (10 ms pulsar)</text>')
svg = '''<figure class="scroll">
<svg class="wide" viewBox="0 0 860 360" role="img" aria-label="Folding: the pulse emerges from the noise as more turns are added; SNR against observing time for several signal strengths">
  ''' + '\n  '.join(o) + '''
</svg>
<figcaption><b>Figure P.5.</b> Folding. (a) Simulated folded profiles at an SNR of 0.25 per turn, scaled so the noise per bin stays the same: the pulse grows
as √<var>N</var> and is clear after a few thousand turns. (b) SNR against observing time for the simulated 10 ms pulsar (SNR per turn ≈
373·ρ): each 10 dB weaker signal needs 100× longer. At −54 dB a TOA with SNR 7 takes about 62 hours; J0437, with its narrower pulse and
both polarizations, needs about 20 in the catalogue estimate; the fast simulator with its measured profile gives 82–95 h (16 h at nominal
flux; <a href="scenario.html#hundred">scenario</a>).</figcaption>
</figure>'''
open('foldgrow.svg', 'w').write(svg)
