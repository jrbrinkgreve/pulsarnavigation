# 2026-10-09 – investigation: fold and correlate in the dispersed domain (no dedispersion)

Jasper's idea: skip coherent dedispersion of the channel voltages; fold the detected power of
each channel as it is (dispersed) and do the correlation (template fit) there. Investigation
only: no code changed, no MATLAB run. Estimates by Claude with two Python scripts (session
scratchpad, not in the repo; method below so they can be redone).

## Physics

- Whole band without channels: impossible. The sweep over 1.2–1.6 GHz is 6.3 ms (simulated
  pulsar, DM 5) = 63 % of the period; J0437 3.3 ms = 58 %. The channelizer is the minimum.
- Per channel, two parts of the dispersion:
  1. the delay of the channel relative to `refFreq`: a pure time shift → exactly absorbed in
     the fold (per-channel phase offset) or in the template (phase ramp on its harmonics). Free.
  2. the sweep inside the channel (smear): 8.3 µs·DM·Δf[MHz]/f[GHz]³, e.g. 75 µs at 1.2 GHz
     for DM 5 and 3.125 MHz. This is the only real difference.
- Weak-signal limit (receiver noise dominates; the noise is the same in both paths): the
  expected detected power of channel c is the true power profile convolved with
  k_c = |h_c|² (power impulse response of channel bandpass × leftover chirp, Σk_c = 1).
  With the known k_c in the template the fit is exact (unbiased); the cost is information:
  detection SNR² ∝ Σ_c ∫ w_c², TOA Fisher information ∝ Σ_c ∫ (w_c')², w_c = profile ⊛ k_c.
  (Locally optimal quadratic detector: x^H H D H^H x for coherent vs x^H D_w x for dispersed;
  deflection ratio Σw²/Σa⁴.) k_c follows from ephemeris DM + channel bandpass → ground-truth
  rule satisfied.
- Computed exactly on the pulse harmonics: K_c(ν) = ∫ H_c(f) H_c*(f−ν) df / ∫|H_c|², with H_c
  the channel bandpass (flat, sin² edges 2 %) times the in-channel chirp; Gaussian profiles
  (W50 from the catalogue table, notes §10), 10 µs detection bins in both paths, equal flux per
  channel, channels combined optimally (as 'optimal' weighting).

## Results: loss of the dispersed domain (ratio dispersed / coherent, whole band)

| Pulsar (P, W50, DM) | channel | smear at 1.2 GHz | detection SNR | TOA error | worst channel TOA | bias if smear ignored |
|---|---|---|---|---|---|---|
| simulated (10 ms, 0.5 ms, 5) | 3.125 MHz | 75 µs | 0.9988 | ×1.0034 | ×1.0074 | 23 ns |
| | 0.781 MHz | 19 µs | 0.9999 | ×1.0002 | ×1.0005 | 1 ns |
| J0437 (5.76 ms, 0.141 ms, 2.64) | 3.125 MHz | 40 µs | 0.9959 | ×1.012 | ×1.026 | 12 ns |
| | 0.781 MHz | 10 µs | 0.9997 | ×1.0007 | ×1.0016 | 1 ns |
| Vela (89.3 ms, 1.7 ms, 67.8) | 3.125 MHz | 1.0 ms | 0.982 | ×1.054 | ×1.12 | 313 ns |
| | 0.781 MHz | 0.25 ms | 0.9988 | ×1.0034 | ×1.0074 | 20 ns |
| B0329+54 (714 ms, 6.6 ms, 26.8) | 3.125 MHz | 0.40 ms | 0.9998 | ×1.0006 | ×1.0012 | 124 ns |
| stress: J0437 shape at DM 100 | 3.125 / 0.781 / 0.195 MHz | 1.5 / 0.37 / 0.09 ms | 0.43 / 0.80 / 0.98 | ×7.3 / ×1.9 / ×1.06 | | |

"Bias if smear ignored" = centroid of k_c in the lowest channel: the TOA offset when the
channel is shifted by its centre delay and the plain template is used. Removed by putting
k_c in the template (or, to first order, shifting by the channel's power-weighted mean delay).
Rule of thumb: channel width such that the in-channel smear is ≲ ¼ of W50.

## Results: computing cost (400 MHz band, 128 × 3.125 MHz, per second of data)

| Stage | GFLOP/s |
|---|---|
| filterbank (`channelizeIQ`: 3853 taps, 256-pt FFT every 192 samples) | 110 |
| `dedisperseChannels` as now (reference = top of band) | 100 (DM 2.6) … 125 (DM 100) |
| – of which the band-edge trim alone (same filter with DM 0) | 91 |
| replacement in the dispersed domain: 16-point second FFT per channel (260 kHz channels) | 11 |
| detection / fold | 2 / 0.3 |

- Total arithmetic about halves (≈ 215–240 → ≈ 125 GFLOP/s); the filterbank becomes the main
  cost. Not an order of magnitude: FFT dedispersion per channel is already cheap (cost grows
  only with log of the filter length).
- Surprise: ~90 % of `dedisperseChannels`' cost at our DMs is the band-edge trim, not the
  dispersion: the sharp 62.5 kHz edge taper needs a 320 µs guard on each side (1334 samples),
  longer than the in-channel sweep. The dispersion itself adds only 10–35 GFLOP/s.
- Bigger gains than arithmetic: (1) memory — coherent dedispersion buffers the whole sweep to
  `refFreq` (bottom channel 0.53 M samples at DM 100), dispersed needs none; (2) files in the
  simulator — 4.3 GB of dedispersed channel files per second of data written and read again
  (plus blanked copies); (3) `blankingWeights` (≈ 20 s per 0.1 s of data in the B3 run, the
  slowest processing step) exists because dedispersion smears each blank over the filter's
  reach — without it a blank only touches its own detected bin; (4) DM becomes a fit parameter
  (DM variations, ionosphere has the same 1/f² law) without reprocessing.

## What would have to change (not done; needs Jasper's OK)

1. Channels must stay independent: today the dedispersion filter also trims each 4/3
   oversampled channel to exactly ±1.5625 MHz. Without it neighbouring channels share noise.
   → a cheap trim (second FFT, drop the edge bins) or model the inter-channel correlation.
2. Per-channel delay in the fold (`foldProfile`) or in the template; per-channel smeared
   template in `estimateTOA` / `detectPulsar` ('optimal' uses one template for all channels).
3. Higher DM needs narrower channels (Vela ≤ 0.8 MHz); cheap with the second FFT.
4. RFI excision is unchanged: `detectRFI` / `periodicRFI` already work on the dispersed channel
   voltages.

The coherent path stays the reference (exact); the dispersed path would be an option next to
it and validated against it on the same seed. A lossless fold-first method exists (cyclic
spectroscopy: fold cross-products of neighbouring frequencies, remove dispersion afterwards)
but it costs more than it saves here.

## Decision (Jasper, 9 Oct)

Just an idea to investigate: noted (notes §10, deferred list), **not to be implemented yet**.

## First step, if it is ever picked up (measure first)

`foldProfile` option for a per-channel time offset (from `ephem.DM` and the channel
frequencies), run with `dedisperseChannels(..., DM = 0)` so the validated band trim and noise
model stay unchanged and only the dispersion handling differs; compare TOAs with the coherent
run (same seed). Expected for the simulated pulsar: TOA error ×1.003, offset ≲ 23 ns (smear not
yet in the template). Saves no time yet (the trim costs the same); it checks the physics.
Then: smear in the template; then the cheap trim if the speed is wanted.
