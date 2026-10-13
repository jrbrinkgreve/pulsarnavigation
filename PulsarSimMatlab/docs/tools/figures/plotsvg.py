import numpy as np
class Ax:
    def __init__(s, x0, y0, w, h, xlim, ylim, logy=False):
        s.x0, s.y0, s.w, s.h, s.xlim, s.ylim, s.logy = x0, y0, w, h, xlim, ylim, logy
    def X(s, x):
        return s.x0 + (np.asarray(x, float) - s.xlim[0]) / (s.xlim[1] - s.xlim[0]) * s.w
    def Y(s, y):
        y = np.asarray(y, float)
        if s.logy:
            y = np.log10(np.maximum(y, 1e-300)); lo, hi = np.log10(s.ylim[0]), np.log10(s.ylim[1])
        else:
            lo, hi = s.ylim
        return s.y0 + s.h - (y - lo) / (hi - lo) * s.h
    def path(s, x, y, clip=True):
        X, Y = s.X(x), s.Y(y)
        if clip:
            Y = np.clip(Y, s.y0 - 2, s.y0 + s.h + 2)
        pts = [f'{a:.1f} {b:.1f}' for a, b in zip(X, Y)]
        # drop consecutive duplicates
        out = [pts[0]] + [p for i, p in enumerate(pts[1:], 1) if p != pts[i-1]]
        return 'M' + ' L'.join(out)
    def area(s, x, y, base=None):
        X, Y = s.X(x), s.Y(y)
        yb = s.y0 + s.h if base is None else s.Y(base)
        pts = ' L'.join(f'{a:.1f} {b:.1f}' for a, b in zip(X, Y))
        return f'M{X[0]:.1f} {yb:.1f} L{pts} L{X[-1]:.1f} {yb:.1f} Z'
    def frame(s, xticks=(), yticks=(), xlab='', ylab='', xfmt='{:g}', yfmt='{:g}', left=True):
        o = [f'<path class="axis" d="M{s.x0} {s.y0+s.h} H{s.x0+s.w}"/>']
        if left:
            o.append(f'<path class="axis" d="M{s.x0} {s.y0} V{s.y0+s.h}"/>')
        for t in xticks:
            x = float(s.X(t)); o.append(f'<path class="axis" d="M{x:.1f} {s.y0+s.h} v5"/>')
            o.append(f'<text class="sub" x="{x:.1f}" y="{s.y0+s.h+18}" text-anchor="middle">{xfmt.format(t) if not isinstance(t, str) else t}</text>')
        for t in yticks:
            y = float(s.Y(t)); o.append(f'<path class="axis" d="M{s.x0} {y:.1f} h-5"/>')
            lab = yfmt(t) if callable(yfmt) else yfmt.format(t)
            o.append(f'<text class="sub" x="{s.x0-8}" y="{y+4:.1f}" text-anchor="end">{lab}</text>')
        if xlab:
            o.append(f'<text class="sub" x="{s.x0+s.w/2:.1f}" y="{s.y0+s.h+34}" text-anchor="middle">{xlab}</text>')
        if ylab:
            o.append(f'<text class="sub" x="{s.x0}" y="{s.y0-8}" text-anchor="start">{ylab}</text>')
        return '\n  '.join(o)

SUP = str.maketrans('-0123456789', '⁻⁰¹²³⁴⁵⁶⁷⁸⁹')
def pow10(v):
    e = int(round(np.log10(v)))
    return '1' if e == 0 else '10' + str(e).translate(SUP)
