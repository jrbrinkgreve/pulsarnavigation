import numpy as np
# The prototype low-pass of channelizeIQ, with the same grid and design formulas (kaiserLowpass), in MHz and samples.
fs, dF, K = 800.0, 3.125, 256           # IQ rate, channel width [MHz], FFT size fs/dF


def design(D=192, A=80.0):
    """Prototype for decimation D: flat to dF/2, stopband (A dB) from fs/D - dF/2. Returns h, tap index n, window w, info."""
    fsOut = fs/D
    transW = fsOut - dF
    fc = fsOut/2
    M = int(np.ceil((A - 8)/(2.285*2*np.pi*transW/fs))) + 1
    if M % 2 == 0:
        M += 1
    beta = 0.1102*(A - 8.7)
    G = (M - 1)//2
    n = np.arange(M) - G
    h = np.where(n == 0, 2*fc/fs, np.sin(2*np.pi*fc/fs*n)/(np.pi*np.where(n == 0, 1, n)))
    w = np.i0(beta*np.sqrt(1 - (2*np.arange(M)/(M - 1) - 1)**2))/np.i0(beta)
    h = h*w
    h /= h.sum()
    return h, n, w, dict(D=D, fsOut=fsOut, transW=transW, fc=fc, M=M, G=G, P=int(np.ceil(M/K)))


def response(h, n, f):
    """Zero-phase gain H(f) at offsets f [MHz] from the channel centre (h is symmetric, so H is real)."""
    f = np.atleast_1d(np.asarray(f, float))
    out = np.empty(f.size)
    for i in range(0, f.size, 512):
        out[i:i + 512] = np.cos(2*np.pi*np.outer(f[i:i + 512], n)/fs) @ h
    return out


def db(x):
    return 20*np.log10(np.maximum(np.abs(x), 1e-12))


def envelope(f, y, step):
    """Largest y per bin of width step: keeps the stopband ripple's envelope with a short path."""
    edges = np.arange(f[0], f[-1] + step, step)
    idx = np.digitize(f, edges)
    keys = np.unique(idx)
    return (np.array([f[idx == k].mean() for k in keys]), np.array([y[idx == k].max() for k in keys]))
