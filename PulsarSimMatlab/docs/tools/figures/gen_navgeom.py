import numpy as np
from plotsvg import Ax
# Navigation geometry (physics.html, section 1): a: a pulse front reaches a position 600 m closer to the pulsar 2 µs earlier;
# a sideways offset changes nothing. b: three pulsars, each TOA puts the station on one line of position; they cross at
# the true position. Positions in metres, c = 299,792,458 m/s.
c = 299792458.0
o = []
o.append('<defs><marker id="ar-ng" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse">'
         '<path class="ah" d="M0 0 L10 5 L0 10 z"/></marker></defs>')


def seg(A, p, q, cls, extra=''):
    return f'<path class="{cls}" d="M{float(A.X(p[0])):.1f} {float(A.Y(p[1])):.1f} L{float(A.X(q[0])):.1f} {float(A.Y(q[1])):.1f}"{extra}/>'


def dot(A, p, cls, r=4.5):
    x, y = float(A.X(p[0])), float(A.Y(p[1]))
    return f'<path class="{cls}" d="M{x - r:.1f} {y:.1f} a{r} {r} 0 1 0 {2*r} 0 a{r} {r} 0 1 0 {-2*r} 0"/>'


def ring(A, p, r=5):
    x, y = float(A.X(p[0])), float(A.Y(p[1]))
    return f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{r}" fill="var(--bg)" stroke="var(--fg)" stroke-width="1.5"/>'


# ---------- a: one pulsar ----------
A = Ax(20, 40, 380, 280, (-1000, 1000), (-1000*280/380, 1000*280/380))
ang = np.radians(30)
n = np.array([np.cos(ang), np.sin(ang)])
t = np.array([-n[1], n[0]])
P = np.array([-200.0, -150.0])
Q = P + 600*n
R = P + 430*t
L = 760
for base, cls in ((P, 'dash'), (Q, 'curve')):
    o.append(seg(A, base - L*t, base + L*t, cls))
# arrow towards the pulsar
o.append(seg(A, Q + 260*n - 1*t, Q + 520*n - 1*t, 'line', ' marker-end="url(#ar-ng)"'))
xq, yq = float(A.X((Q + 540*n)[0])), float(A.Y((Q + 540*n)[1]))
o.append(f'<text class="sub" x="{xq - 6:.1f}" y="{yq - 6:.1f}" text-anchor="end">to the pulsar</text>')
o.append(seg(A, P, Q, 'curve2', ' stroke-width="3"'))
o.append(seg(A, P, R, 'thin', ' stroke-dasharray="3 3"'))
o.append(ring(A, P))
o.append(dot(A, Q, 'dot2'))
o.append(dot(A, R, 'dot'))
xp, yp = float(A.X(P[0])), float(A.Y(P[1]))
o.append(f'<text class="lbl" x="{xp - 12:.1f}" y="{yp + 5:.1f}" text-anchor="end">predicted position</text>')
xm, ym = float(A.X((P + 300*n)[0])), float(A.Y((P + 300*n)[1]))
o.append(f'<text class="lbl" x="{xm - 30:.1f}" y="{ym + 30:.1f}">600 m closer:</text>')
o.append(f'<text class="lbl" x="{xm - 30:.1f}" y="{ym + 45:.1f}">pulse {600/c*1e6:.0f} µs early</text>')
xr, yr = float(A.X(R[0])), float(A.Y(R[1]))
o.append(f'<text class="sub" x="{xr + 10:.1f}" y="{yr - 2:.1f}">sideways: on the same</text>')
o.append(f'<text class="sub" x="{xr + 10:.1f}" y="{yr + 12:.1f}">front, no change</text>')
o.append(f'<text class="ttl" x="20" y="20">a · one pulsar measures one direction</text>')
fa = Q - 0.62*L*t
o.append(f'<text class="sub" x="{float(A.X(fa[0])) + 8:.1f}" y="{float(A.Y(fa[1])) + 4:.1f}">pulse front at t</text>')
fb = P - 0.52*L*t
o.append(f'<text class="sub" x="{float(A.X(fb[0])) - 8:.1f}" y="{float(A.Y(fb[1])) + 4:.1f}" text-anchor="end">the same front 2 µs later</text>')

# ---------- b: three pulsars, lines of position ----------
B = Ax(460, 40, 380, 280, (-1000, 1000), (-1000*280/380, 1000*280/380))
P = np.array([0.0, 0.0])
d = np.array([430.0, -170.0])
cls3 = ['curve', 'curve2', 'curve3']
angs = [30, 140, 255]
TC = {'curve': 'var(--accent)', 'curve2': 'var(--warn)', 'curve3': 'var(--why)'}
for k, (a_deg, cls) in enumerate(zip(angs, cls3)):
    a = np.radians(a_deg)
    nk = np.array([np.cos(a), np.sin(a)])
    tk = np.array([-nk[1], nk[0]])
    proj = nk @ d
    foot = P + proj*nk
    o.append(seg(B, foot - 800*tk, foot + (650 if k == 0 else 800)*tk, cls))
    o.append(seg(B, P + 330*nk, P + 600*nk, 'line', ' marker-end="url(#ar-ng)"'))
    lp = P + 690*nk
    o.append(f'<text class="sub" x="{float(B.X(lp[0])):.1f}" y="{float(B.Y(lp[1])) + 4:.1f}" text-anchor="middle" style="fill:{TC[cls]}">to pulsar {k + 1}</text>')
    dt = proj/c*1e6
    word = 'early' if dt > 0 else 'late'
    ll = foot + {0: -720, 1: -520, 2: -700}[k]*tk
    anc = {0: 'start', 1: 'start', 2: 'middle'}[k]
    o.append(f'<text class="sub" x="{float(B.X(ll[0])) + (8 if anc == "start" else 0):.1f}" y="{float(B.Y(ll[1])) - 6:.1f}" text-anchor="{anc}" style="fill:{TC[cls]}">TOA {abs(dt):.1f} µs {word}</text>')
o.append(ring(B, P))
o.append(dot(B, d, 'dot2', 5.5))
o.append(f'<text class="lbl" x="{float(B.X(0)) - 8:.1f}" y="{float(B.Y(0)) + 20:.1f}" text-anchor="end">predicted</text>')
o.append(f'<text class="lbl" x="{float(B.X(d[0])) + 10:.1f}" y="{float(B.Y(d[1])) - 8:.1f}">true position</text>')
o.append(f'<text class="ttl" x="460" y="20">b · several pulsars: the lines of position cross</text>')

svg = '''<figure class="scroll">
<svg class="wide" viewBox="0 0 860 340" role="img" aria-label="Navigation geometry: a pulse arrives early when the station is closer to the pulsar; three pulsars give three lines of position that cross at the true position">
  ''' + '\n  '.join(o) + '''
</svg>
<figcaption>How arrival times give a position (to scale, in a plane). a: a pulse front moves at the speed of light, so a station
600 m closer to the pulsar meets it 2 µs earlier; a sideways offset changes nothing. b: each pulsar's TOA therefore puts the station
on one line, perpendicular to that pulsar's direction. Where the lines cross is the position. In 3-D the lines become planes and
three pulsars fix the position; a receiver clock that is off moves all arrival times by the same amount, a fourth unknown, so four
pulsars are needed.</figcaption>
</figure>'''
open('navgeom.svg', 'w').write(svg)
