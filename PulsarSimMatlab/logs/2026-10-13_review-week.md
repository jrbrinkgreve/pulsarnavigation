# 2026-10-13 – review week, continued: stages 2 and 3

Continues `2026-10-12_review-week.md` (stage 1). Same method: physics, code walk, tests,
refinements noted (not implemented) in notes §10 item 11 and the docs roadmap. The code-review
work of the same day is in `2026-10-13_code-review.md` (separate sessions).

## Stage 2 – dispersion (`applyDispersionStream`, notes §5.2)

Walk-through given (code unchanged since 12 Oct):
- Physics: τ(f) = K·DM/f², DM 5 → 8.103 ms at 1.6 GHz, 14.406 ms at 1.2 GHz, smear 6.3 ms
  (12× the pulse). A phase-only filter on the voltage, exact 1/f phase (no Taylor expansion):
  group delay of φ = 2π·K·DM·(f − fRef)²/(f·fRef²) is exactly τ(f) − τ(fRef); bulk τ(fRef)
  recorded, not applied. Band-limiting = the receiver bandpass (dispersion and bandpass commute);
  power 0.195 = 2·(400 − 1.25·8)/4000 (sin⁴ edges average 3/8).
- `makeDispersionKernel` in detail (Jasper's question "why is nothing wrapping?"): the impulse
  response is a descending chirp 1.6 → 1.2 GHz from sample 10,000 to 25,219,771 between two
  2.5 µs guards (Nk = 25,229,772). Sampling H on Nf bins makes the time response periodic
  (period Nf); h_true lives in [0, Nk), Nf = 2^26 ≥ 2·Nk leaves an empty zone of 41.9 M samples
  where the late tail and the wrapped pre-ringing land and are cut; `leak` (3.5e-13) measures
  exactly that energy, every run. Frequency view: largest phase step Nk/Nf = 0.376 turn per bin
  (< ½). Second guard: the streaming FFT, n + M − 1 ≤ Nfft (circular = linear convolution);
  the carry is overlap-add, not wrapping. Guard rule 20/edgeW is empirical.
- Flowcharts drawn in the chat (overview; one overlap-add pass).
- Tests: NumPy ports (not in the repo), `plotDispersionCheck`, end-to-end runs; no MATLAB unit test.

Noted (notes §10 item 11, roadmap row "Dispersion refinements"): **D1** unit test
`testApplyDispersionStream`; **D2** one `refFreq` for the simulation's bulk delay and the
observer's dedispersion (truth holds only because they are equal; linked to R3 and code-review
item 24). DM 100 (item 5): this stage is the blocker, nothing new.

## Stage 3 – receiver noise and RFI (`addNoiseAndRFI`, `rfiSource`, notes §5.3)

Walk-through given (current code, i.e. after today's `0499da7`: `:weak` warning removed,
`'noise'` spectrum cached):
- Physics: noise and RFI enter at the receiver, after dispersion, so they are not dispersed
  (dedispersion later smears impulses, shifts carriers). `snrDB` = peak PSD ratio, SNR = A²/σ_n².
  Numbers: −5 dB → σ_n 1.78, P_band 0.632 vs pulsar in-band peak 0.195 (ratio 0.309 = ρ minus the
  2.5 % edge tapers); −54 dB → σ_n 501, SNR per pulse 388·ρ = 0.0015, ~60 h per TOA. Self-noise
  (1 + ρp)². INR = total in-band power (a carrier is huge per Hz but in one channel only).
- Code: noise level and P_band; per-source setup (amplitude from INR, own `mrg32k3a` stream,
  noise filter, impulse event list); `predictPerformance` (flat band, ~1.5 % optimistic); stream
  loop; `rfiBlock` per type, phases from the absolute sample index and substreams by absolute
  position (block-size independent); radar edges; rotating antenna; `noiseBandpass`.
- Tests: `testRadarEdges`, `testRfiGating`, `testRfiNoise` (new features vs formulas; older types
  only bit-identical to the pre-B5a reference); noise level end to end (−5 dB Monte Carlo).

Noted (notes §10 item 11, roadmap row "Receiver-stage test"): **N1** unit test
`testAddNoiseAndRFI` (SNR definition, INR per type, block size for all types, 388·ρ). Links only:
receiver noise not band-passed, no ADC, `'noise'` speed (item 10b), `predictPerformance`
duplicate (R2, item 20e).

## Stage 4 – IQ conversion (`applyIQmodulation`, notes §5.4)

Walk-through given (code unchanged since 12 Oct), with today's parameters: complex envelope
z = LPF{2·x·e^(−i2π·fLO·t)}; band 1.2–1.6 GHz → −100…+300 MHz (no inversion), mirror at
+1.1…+1.5 GHz (mod 4 GHz) removed; factor 2 → mean|z|² = 2·mean(x²) (0.195 → 0.39); D = 5,
800 MHz complex covers ±400 MHz; Kaiser low-pass 2049 taps, cutoff 360 MHz, β 7.86, transition
9.8 MHz (flat to ±355, stopband from 365 < 400); G = 1024 (256 ns) removed exactly, output k ↔
input 5k on the global grid; LO phase from the absolute index (r = 13/40); Nfft 2^22, 955 blocks
per s, ~350 MB; 800,000,001 complex samples per s. It is a digital downconverter: ideal I/Q.
Tests: NumPy ports (not in the repo), `plotIQCheck`, end-to-end runs; no MATLAB unit test.

Noted (notes §10 item 11 "Stage 4", roadmap row "IQ-stage refinements"): **I1** unit test,
**I2** stale §5.4 numbers (text only; Jasper's choice fix-now vs later still open), **I3**
ideal-I/Q physics gap (imbalance, DC/LO leakage, phase noise; LO inside our band).

## Back to stage 2: visual guide to `makeDispersionKernel` (Jasper's request)

An interactive toy kernel in the chat (same recipe, band 0.3–0.4·fs = 1.2–1.6 GHz at 4 GHz,
sweep 400 samples, edge 0.01·fs): spectrum on the design grid (W and Re H per bin) and one
period of the inverse FFT in dB with kept window / empty zone; sliders Nf/Nk and guard g;
readouts leak (what the code measures) and folded-in error vs an 8× grid (what it cannot
see: tails wrapped into the kept window, the real reason for Nf ≥ 2·Nk). Things to try:
factor 2 clean, 1.1 tail folds back, < 1 the chirp wraps (Re H under-sampled), g = 0
pre-ringing at the period end.

## Next session (Jasper, 13 Oct: "let's continue tomorrow with the dispersion thing")

1. Continue with stage 2 / `makeDispersionKernel` (the visual guide; Jasper's questions).
2. Stage 4: Jasper's choice on I2 (fix §5.4 now, text only, or later).
3. Stage 5: channelizer (`channelizeIQ`, notes §5.15).

Not committed (Jasper decides): this log, notes §10 item 11 (D1, D2, N1, I1–I3), roadmap rows
(Dispersion refinements, Receiver-stage test, IQ-stage refinements). Other files in the
working tree (`docs/physics.html`, `docs/tools/figures/gen_*.py`, `2026-10-13_deep-analysis.md`
additions) are from the other sessions of the day.
