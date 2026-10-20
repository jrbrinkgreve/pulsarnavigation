import numpy as np
# Data volume per stage (pipeline.html, sizes): bytes written per second of simulated signal, from the parameters:
# f_in 4 GHz real float32; fs 800 MHz complex float32; 128 channels of 3.125 MHz, 4/3 oversampled (4.1667 MHz complex);
# detected power float32 at f_out (full band 1 MHz; channels 4.1667 MHz / 4 = 1.0417 MHz per channel); data weights 9 float32
# per bin and channel (W, V, X(1..7)). Log scale, with the 24 GB memory of the machine.
f_in, fs, nch, fsch = 4e9, 800e6, 128, 3.125e6*4/3
rows = [('sky', 'test.dat', 'generator, real 4 GHz', f_in*4),
        ('sky', 'test_dispersed.dat', 'dispersed', f_in*4),
        ('rx', 'test_rx.dat', 'noise + RFI added', f_in*4),
        ('rx', 'test_rx_IQ.dat', 'complex, 800 MHz', fs*8),
        ('fb', 'test_IQ_dedispersed.dat', 'full band: dedispersed', fs*8),
        ('fb', 'test_envelope.dat', 'full band: power, 1 µs', 1e6*4),
        ('ch', '…_IQ_chan_ch###', 'channels: 128 × 4.17 MHz', nch*fsch*8),
        ('ch', '…_blanked_ch###', 'excision: blanked copies', nch*fsch*8),
        ('ch', 'dedispB_ch###', 'channels dedispersed', nch*fsch*8),
        ('ch', '…_power.dat', 'channels: power, 0.96 µs', nch*fsch/4*4),
        ('ch', '…_weights.dat', 'excision: W, V, X(1..7)', nch*fsch/4*4*9)]
cls = {'sky': 'box sky', 'rx': 'box rx', 'fb': 'box proc', 'ch': 'box proc'}
x0, x1 = 330, 800
lo, hi = 6, 11
X = lambda b: x0 + (np.log10(b) - lo)/(hi - lo)*(x1 - x0)
o = []
y0, dy = 50, 24
o.append(f'<text class="ttl" x="10" y="20">bytes written per second of simulated signal (the default L = 0.2 s writes a fifth)</text>')
for i, (g, f, what, b) in enumerate(rows):
    y = y0 + i*dy + (12 if i >= 4 else 0) + (12 if i >= 6 else 0)
    o.append(f'<rect class="{cls[g]}" x="{x0}" y="{y}" width="{X(b) - x0:.1f}" height="{dy - 7}" rx="3"/>')
    o.append(f'<text class="mono" x="10" y="{y + 13}">{f}</text>')
    o.append(f'<text class="sub" x="{x0 - 8}" y="{y + 13}" text-anchor="end">{what}</text>')
    lab = f'{b/1e9:.1f} GB' if b >= 1e9 else (f'{b/1e6:.0f} MB')
    if b > 1e9:
        o.append(f'<text class="sub" x="{X(b) - 6:.1f}" y="{y + 13}" text-anchor="end">{lab}</text>')
    else:
        o.append(f'<text class="sub" x="{X(b) + 6:.1f}" y="{y + 13}">{lab}</text>')
ybot = y0 + len(rows)*dy + 24
o.append(f'<path class="axis" d="M{x0} {ybot} H{x1}"/>')
for e, lab in ((6, '1 MB'), (7, '10 MB'), (8, '100 MB'), (9, '1 GB'), (10, '10 GB'), (11, '100 GB')):
    o.append(f'<path class="axis" d="M{X(10**e):.1f} {ybot} v5"/><text class="sub" x="{X(10**e):.1f}" y="{ybot + 18}" text-anchor="middle">{lab}</text>')
xm = X(24e9)
o.append(f'<path class="dash" d="M{xm:.1f} {y0 - 8} V{ybot}"/>')
o.append(f'<text class="sub" x="{xm + 4:.1f}" y="{y0 - 12}">24 GB memory</text>')
tot = {'full band': sum(b for g, *_, b in rows if g in ('sky', 'rx', 'fb')),
       'channels': sum(b for g, f, w, b in rows if g in ('sky', 'rx') or w in ('channels: 128 × 4.17 MHz', 'channels dedispersed', 'channels: power, 0.96 µs')),
       'channels + excision': sum(b for g, *_, b in rows if g in ('sky', 'rx', 'ch'))}
txt = ' · '.join(f'{k} {v/1e9:.1f} GB' for k, v in tot.items())
o.append(f'<text class="sub" x="10" y="{ybot + 40}">Total per second of signal: {txt}</text>')
H = ybot + 52
svg = f'''<figure class="scroll">
<svg class="wide" viewBox="0 0 860 {H}" role="img" aria-label="Bytes written per second of simulated signal by each stage, on a log scale, against the 24 GB memory">
  ''' + '\n  '.join(o) + '''
</svg>
<figcaption>How much each stage writes per second of simulated signal, computed from the parameters (log scale; blue: sky, orange:
receiver, green: processing). At <var>L</var> = 1 s each RF file is 16 GB and a run writes 61–73 GB, far more than the 24 GB of
memory, which is why every stage streams in blocks. Detection is where the data shrink: on the full band 800 complex samples become one power value (1600× fewer bytes); on the
channel path each 0.96 µs bin keeps all 128 channels, so it shrinks 8×. With excision the weight file adds 9 numbers per bin and channel.</figcaption>
</figure>'''
open('datarate.svg', 'w').write(svg)
print({k: round(v/1e9, 2) for k, v in tot.items()})
