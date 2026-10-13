import numpy as np
from plotsvg import Ax
# Why noise and RFI are added after dispersion: what dedispersion does to the pulsar, to a broadband
# impulse and to a carrier. 16 sub-bands of 25 MHz, expected power, the code's constants (DM 5,
# f_ref = 1.6 GHz). Dedispersion moves frequency f earlier by tau(f) - tau(f_ref).
K, DM = 4.148808e3, 5.0
fL, fH, T, duty = 1200.0, 1600.0, 10.0, 5.0
sig = duty/100*T/(2*np.sqrt(2*np.log(2)))
dtau = lambda f: (K*DM/f**2 - K*DM/fH**2)*1e3       # ms
G = lambda t: sum(np.exp(-(t - c)**2/(2*sig**2)) for c in (5.0, 15.0, 25.0))
nSub, bw, rowH, top = 16, 25.0, 15.0, 60
t = np.linspace(0, 23, 1151)
tImp, fCw = 11.0, 1351.3
iCw = int((fH - fCw)//bw)
o = []
def simplify(X, Y):
    X = np.round(X, 1); Y = np.round(Y, 1)
    keep = np.r_[True, (Y[1:-1] != Y[:-2]) | (Y[1:-1] != Y[2:]), True]
    return 'M' + ' L'.join(f'{a:.1f} {b:.1f}' for a, b in zip(X[keep], Y[keep]))
def panel(x0, w, dedispersed, title):
    A = Ax(x0, top, w, rowH*nSub, (0, 23), (0, 1))
    o.append(A.frame(xticks=[0, 5, 10, 15, 20], xlab='time [ms]', left=False))
    o.append(f'<text class="ttl" x="{x0}" y="22">{title}</text>')
    # the carrier: on all the time, in one sub-band
    o.append(f'<rect class="fillwarn" x="{A.x0}" y="{top + iCw*rowH + 2:.1f}" width="{A.w}" height="{rowH - 4:.1f}"/>')
    for i in range(nSub):
        fc = fH - (i + 0.5)*bw
        yb = top + (i + 1)*rowH
        if dedispersed:
            p = G(t)
        else:
            fs = fc + np.linspace(-bw/2, bw/2, 41)
            p = np.mean([G(t - dtau(f)) for f in fs], axis=0)
        o.append(f'<path class="curve" stroke-width="1.4" d="{simplify(A.X(t), yb - 0.95*rowH*p)}"/>')
        # the impulse: one instant at every frequency at the input; after dedispersion each
        # frequency moves earlier, so the sub-band holds it from t - dtau(top) to t - dtau(bottom)
    if dedispersed:
        ff = np.linspace(fL, fH, 200)
        o.append(f'<path class="curve4" stroke-width="2" d="{simplify(A.X(tImp - dtau(ff)), top + (fH - ff)/bw*rowH)}"/>')
    else:
        o.append(f'<path class="curve4" stroke-width="2" d="M{A.X(tImp):.1f} {top:.1f} V{top + nSub*rowH:.1f}"/>')
    return A
A = panel(110, 330, False, 'a · at the receiver input (output of stage 3)')
B = panel(520, 330, True, 'b · after dedispersion (stage 5)')
for f in (1600, 1500, 1400, 1300, 1200):
    y = top + (fH - f)/bw*rowH
    o.append(f'<text class="sub" x="{A.x0-8}" y="{y+4:.1f}" text-anchor="end">{f} MHz</text>')
    o.append(f'<path class="axis" d="M{A.x0} {y:.1f} h-4"/>')
# legend
ly = 40
o.append(f'<path class="curve" d="M110 {ly} h22"/><text class="sub" x="138" y="{ly+4}">pulsar (added at the sky, dispersed)</text>')
o.append(f'<path class="curve4" stroke-width="2" d="M360 {ly} h22"/><text class="sub" x="388" y="{ly+4}">broadband impulse (added at the receiver)</text>')
o.append(f'<rect class="fillwarn" x="652" y="{ly-5}" width="22" height="10"/><text class="sub" x="680" y="{ly+4}">carrier at 1351.3 MHz</text>')
svg = '''<figure class="scroll">
<svg class="wide" viewBox="0 0 860 352" role="img" aria-label="What dedispersion does: the pulsar sweep becomes aligned, a broadband impulse becomes an inverse sweep, a carrier stays a carrier">
  ''' + '\n  '.join(o) + '''
</svg>
<figcaption>Why noise and RFI are added <em>after</em> dispersion (expected power in 16 sub-bands of 25 MHz, the code's constants). a: the
receiver sees the pulsar already dispersed, but an impulse made at the receiver arrives at all frequencies at once, and a carrier is
simply on. b: dedispersion moves each frequency earlier by τ(<var>f</var>) − τ(1.6 GHz). The pulsar lines up; the impulse becomes an
<em>inverse</em> sweep, 6.3 ms long; the carrier stays where it was. Real data behaves like this, so the simulator must add RFI at this
point. It also shows why excision works <em>before</em> dedispersion: there the impulse is short and in one place.</figcaption>
</figure>'''
open('rxorder.svg', 'w').write(svg)
