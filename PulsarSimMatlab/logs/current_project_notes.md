# Pulsar navigation – synthetic data & processing pipeline: project notes

Reference document for the MATLAB pipeline that simulates a pulsar signal through a
realistic receiver chain and processes it back into times of arrival (TOAs). It records
what exists, the conventions every stage follows, the physics and formulas behind each
step, what has been validated (with numbers), what went wrong along the way, and what is
still to do.

This file tracks the **current state** of the code (started 30 September 2026 as
`2026-09-30_project_notes.md`). Last updated: 5 October 2026 (§9: physical-fidelity and
scintillation assessment; bandpass review and exact radiometer optimum in
`expectedPowerModel`, §5.3/§5.10/§6/§10; phase D SNR sweep `runSNRSweep.m`, §5.12/§7;
RFI test, bpsk fix, `rfiSelect`, §5.3/§7/§8; see `2026-10-05_fidelity.md`). Then 6 October
2026: phase E (`detectPulsar`, `flagChi2`, good TOAs), §5.8/§5.12/§5.13/§7, and the long
H0 run (`runH0.m`, §5.14/§7), and the TOA threshold sweep (design rule SNR ≥ 6–7, baseline
loss, turns per sub-int, §6/§7/§10), and the centroid-offset test (Block 1 closed), reference scenario draft (§11); see
`2026-10-06_phase-e.md`. Status markers: **[validated]** = run in MATLAB and
checked against ground truth; **[written]** = code exists, not yet run;
**[todo]** = not implemented.

---

## 1. Goal

Pulsar signals arrive at a known location at extremely predictable times. A pulse
arriving earlier than predicted means the observer is closer to the pulsar along that
direction; TOAs from several pulsars (≥ 4 for 3-D position + clock offset) give a
position. This project builds, in MATLAB:

1. a **physics-realistic synthetic data generator** (pulsar → interstellar medium →
   receiver), with full ground truth, and
2. a **processing pipeline** that turns receiver data into TOAs with honest
   uncertainties, and later into a navigation solution,

and validates every processing stage against the ground truth.

---

## 2. Conventions (apply to every stage)

**File I/O**
- Every stage is **file in → file out, streamed in blocks**; memory is bounded by the
  block/FFT size, never by the file size.
- Real data: **float32, little-endian** (`fopen(..., 'ieee-le')`), no header.
- Complex data: **interleaved float32 I,Q** (I0, Q0, I1, Q1, …; the `.cf32`
  convention), little-endian.
- Detected power: float32 `[nChan × nBins]`, column-major (all channels of bin 1, then
  bin 2, …).
- Every `fwrite` count is checked (a full disk raises an error, not a truncated file).

**Info structs**
- Every stage returns an `info` struct and saves it as `<outputname>_info.mat` next to
  its data file. `main.m` chains stages through these (`info.file`,
  `info.actualFsOut`, …) and reloads them with `loadInfo` when a stage is skipped.
- Common fields: `file`, `inFile`, `precision`, `byteOrder`, `isComplex`, `N`, `fs`,
  `actualFsOut`, stage-specific parameters and conventions, `elapsed`.

**Arguments**
- Required inputs positional; everything else name-value (`arguments` blocks,
  R2019b+; `mustBeTextScalar` needs R2020b+). No Signal Processing Toolbox needed
  (filters are hand-built).

**Time**
- Times are computed in **double precision** from absolute sample indices (single
  precision would round to ~60 ns at t = 1 s). Phases of oscillators/chirps are
  computed from the absolute index, never accumulated.
- Sample k (0-based) of a stream at rate fs is at time k/fs (+ T0, currently 0).

**Ground-truth rule**
- **Processing stages only use what a real observer would know**: the ephemeris
  (`ephem` struct in main: f0, F1, TRef, DM, profile width) and receiver settings
  (band, fLO, fs, from `info_IQ`). Ground truth (`info_gen`, `info_disp`, `info_rx`) is
  used only by the synthetic stages and by the check/plot functions.

**Reproducibility**
- Private `RandStream`s: pulsar signal (`mt19937ar`, `seed`), receiver noise
  (`mt19937ar`, `noiseSeed`), RFI (`mrg32k3a` substreams). The global RNG is never
  touched. Same seed → identical data, independent of block sizes.

---

## 3. Pipeline map (`main.m`)

**Scripts** (all run from the `PulsarSimMatlab` folder; `dataDir` is relative):
- `pipelineParams.m`: all parameters (pulsar, simulation, receiver, noise/RFI,
  processing, the observer's `ephem` struct, file names). A script, so its variables
  land in the caller's workspace. Edit parameters here.
- `main.m`: run control (`runStage.*`, `plots.*`), `addpath(functions)`,
  `pipelineParams;`, then the stages below for one realization.
- `runMonteCarlo.m`: many receiver-noise realizations through the same stages (§5.11).
- `runSNRSweep.m`: phase D, TOA precision and error bars versus SNR, pulsar and noise
  realizations both varied (§5.12).
- `runH0.m`: long noise-only run, detection statistics and bin-bin noise correlation under
  H0 (§5.14).

```
SKY (synthetic)                      runStage.sky
  1. generatePulsarSignal   ─► test.dat              real, f_in = 4 GHz   [validated]
  2. applyDispersionStream  ─► test_dispersed.dat    real, 4 GHz          [validated]
RECEIVER (synthetic)                 runStage.receiver
  3. addNoiseAndRFI         ─► test_rx.dat           real, 4 GHz          [written]
  4. applyIQmodulation      ─► test_rx_IQ.dat        complex, 500 MHz     [validated*]
PROCESSING                           runStage.process / .fold / .toa
  5. applyInverseDispersion ─► test_IQ_dedispersed.dat  complex, 500 MHz  [validated]
  6. detectPower            ─► test_envelope.dat     power, 1 MHz bins    [validated]
  7. foldProfile            ─► test_fold.mat         profiles             [validated]
  8. estimateTOA            (in memory: toa, info_toa)                    [validated]
PROCESSING, channelized (in development, §10 item 1; not yet in main)
  4b. channelizeIQ          ─► data/chan/*_ch###.dat 128 × complex, 4.17 MHz [validated]
      (excision per channel goes here, step 2)
  5b. dedisperseChannels    ─► data/chan/*_ch###.dat 128 × complex, aligned at 1.6 GHz [validated]
  6b. detectChannels        ─► *_power.dat [128 × nBins], 0.96 µs bins + noise stats [validated]
      (next: foldProfile with per-channel weights and lag covariances, 3c)
VALIDATION
  9. validateTOA            (val struct + figure)                         [validated]
CHECKS (ground truth, optional via plots.*)
  plotDispersionCheck, plotIQCheck, plotDetectedPower, plotFoldCheck,
  expectedPowerModel (shared model)
```
\* validated on noise-free data. The receiver stage (3) and the noise-aware check
functions are validated with noise: reference run −5 dB, L = 0.1 s, fLO 1.3 GHz,
fs 800 MHz, noise seed 43 (SNR 93.4, σ 2.880 µs), Monte Carlo over noise seeds, SNR and
TOA threshold sweeps down to −25 dB, long H0 run (all in §7).

**Run control** (`main.m`): `runStage.sky` (generator + dispersion, the slow part),
`runStage.receiver` (noise/RFI + IQ; rerun alone to change SNR or RFI),
`runStage.process` (dedispersion + detection), `runStage.fold`, `runStage.toa`.
Skipped stages reload their `info` from disk; `checkConsistency` warns when files on disk
were made with different parameters (T, f_in, L, DM, band, SNR, receiver input, fLO, fs,
RFI sources; the RFI check since 5 Oct, `isequaln` on `info_rx.rfi.params` vs `rfi`).
`plots.*` switch the check figures; `closeFigures` runs `close all`.

**Helpers** (in `functions/`, shared by main and the MC script): `loadInfo(dataFile)`
(reload `<name>_info.mat`), `gaussianTemplate(nBin, fwhmTurns)` (periodic Gaussian,
peak 1 at bin 1), `noiseBandwidth(fLow, fHigh, edgeWidths)` (= (∫W²)²/∫W⁴ for a
product of raised-cosine edge tapers). **Local function in main**: `checkConsistency`.
**Detection** (processing, after `estimateTOA` in main and `runSNRSweep`):
`detectPulsar` (§5.13); main defines `toa.good = valid & detectedUnknown & ~flagChi2`
and runs `validateTOA` also on the good TOAs only.

**Disk use** ≈ L·(3·f_in·4 + 2·fs·8 + f_out·4) bytes: for L = 1 s at fs = 500 MHz about
56 GB (three 16 GB RF files, two 4 GB complex files, 4 MB power); at fs = 800 MHz about
61 GB. `test.dat` can be deleted after dispersion (only `plotDispersionCheck` needs it).
The MC adds rx + IQ + dedispersed files in `data/mc/` (≈ 2.9 GB at L = 0.1 s, 800 MHz).

---

## 4. Reference configuration

| Group | Parameter | Value | Notes |
|---|---|---|---|
| Pulsar | T | 10 ms | f0 = 100 Hz |
| | A | 1 | pulsar noise std at pulse peak (signal scale) |
| | dutycycle | 5 % | FWHM of the **power** profile = 0.5 ms, σ = 212.3 µs |
| | genEnvMode | `'power'` | power profile Gaussian with that FWHM |
| | DM | 5 pc cm⁻³ | |
| Simulation | f_in | 4 GHz | real RF sampling |
| | L | 1 s (0.1 s for quick tests) | 100 pulse centres at 5, 15, …, 995 ms |
| | seed | 42 | pulsar signal |
| Receiver | band | 1.2–1.6 GHz | B = 400 MHz |
| | fLO | 1.4 GHz | |
| | fs | 500 MHz | IQ rate, D = 8 |
| | filterOrder | 2048 | 2049-tap IQ low-pass |
| | snrDB | −20 dB | S_peak/SEFD, see §5.3 |
| | noiseSeed | seed + 1 | |
| | rfiOn | false (then true) | scenario in main |
| Processing | refFreq | fHigh = 1.6 GHz | dedispersion reference |
| | f_out | 1 MHz | 1 µs detection bins |
| | nBin | 1024 | 9.77 µs phase bins |
| | subintPeriods | 10 (1 for noise-free tests) | turns per sub-integration; with L = 0.1 s the fold spans turns 0–9 (data ends at 93.7 ms, turn 9 ~37 % covered) → only **one** sub-int / one TOA |

**Current `pipelineParams.m`** (1 Oct 2026) differs from the reference: L = 0.1 s,
fLO = 1.3 GHz, fs = 800 MHz (D = 5; band at −100…+300 MHz baseband, which also tests
the RF↔baseband mapping), snrDB = −5 (ρ = 0.316: SNR per pulse 97.9, best TOA error per
pulse 2.85 µs), subintPeriods = 1, nBin = 2048 (4.88 µs bins).

---

## 5. Stage reference

### 5.1 `generatePulsarSignal(outFile, T, f_in, A, L, dutycycle, ...)` [validated]

**Physics.** Pulsar radio emission is an incoherent sum of many emitters → the voltage
is Gaussian noise whose **variance** is modulated once per rotation:
x(t) = A·√G(t)·n(t), n white unit-variance, G(t) periodic train of unit-peak Gaussians.
E[x²] = A²·G(t). The noise is white over 0 – f_in/2 (a property of the simulation grid,
not of real pulsars, whose spectra fall as a power law).

**Code.**
- FWHM = dutycycle/100·T; σ = FWHM/(2√(2 ln 2)).
- `'power'` mode: amplitude envelope √G, so the power profile G has the stated FWHM
  (`'amplitude'` mode: amplitude G, power FWHM smaller by √2).
- Pulse centres at (k + 0.5)·T (first pulse at T/2, signal starts near zero) →
  `info.pulseCenters` = ground truth.
- Periodic sum over `nNeighbors = max(1, ceil(6σ/T + 0.5))` neighbours, nearest pulse
  **per sample** (`k0 = round(t/T − 0.5)`), normalized by the exact peak value
  (`peakNorm`) → no clipping for overlapping pulses.
- N = round(L·f_in) + 1 (inclusive endpoint → file lengths end in …001).
- Blocks of 4e6 samples (1001 blocks for 1 s).
- Options: `Seed`, `EnvelopeMode`, `BlockSize`, `SaveInfo`, `Verbose`.

**Info:** format, N, fs, T, A, L, dutycycle, envelopeMode, sigma, FWHM_power,
FWHM_amplitude, pulseCenters, nPulses, seed, rngType, bytesWritten, elapsed.

Progress (with `Verbose`) is printed about every 10 % of the blocks.

### 5.2 `applyDispersionStream(inFile, outFile, DM, fs, fLow, fHigh, ...)` [validated]

**Physics.** Cold plasma: **phase advance, group delay**. Delay relative to infinite
frequency τ(f) = K·DM/f², K = 4.148808×10³ s·MHz²·pc⁻¹·cm³. With DM = 5:
τ(1.6 GHz) = 8.103 ms, τ(1.2 GHz) = 14.406 ms, smear 6.302 ms. Only the **relative**
delay τ(f) − τ(fRef) is applied (fRef = fHigh by default: the top of the band keeps its
time, lower frequencies arrive later); the bulk τ(fRef) is recorded
(`info.bulkDelayRef`), not applied.

**Transfer function** (standard coherent-dispersion chirp):

H(f) = W(f) · exp(+i·2π·K·DM·(f − fRef)² / (f·fRef²)) · exp(−i·2π·f·s0/fs)

Group delay −(1/2π)·dφ/df = K·DM/f² − K·DM/fRef² + s0/fs. W(f): 1 in the band,
raised-cosine (sin²) edges of width `EdgeFrac`·B = 8 MHz **inside** the band, 0
outside.

**Why band-limited.** (1) Delays diverge as f → 0, no finite kernel can represent them.
(2) Dispersion is phase-only (|H| = 1) and frequency-preserving: energy never moves
between frequencies, so out-of-band content cannot affect in-band data. (3) Dispersion
and the receiver bandpass are both LTI filters and commute, so band-limiting before or
after dispersion gives identical in-band data. W therefore models the receiver
bandpass, not the sky. (Caveat: real receivers are not ideal; strong out-of-band RFI can
cause in-band products through non-linearity; not modelled.)

**Kernel design (`makeDispersionKernel`).**
- Guard g = ceil(GuardTime·fs), GuardTime = 20/edge width = 2.5 µs → 10,000 samples.
- s0 = g + ceil((τ(fRef) − τ(fHigh))·fs) (= 10,000 for fRef = fHigh): zero-lag index
  of fRef.
- Kernel length Nk = s0 + ceil((τ(fLow) − τ(fRef))·fs) + g + 1 = **25,229,772** at 4 GHz.
- Designed on a grid Nf = 2^nextpow2(2·Nk) (2^26), in-band bins only, phase in
  **double** (chirp spans ~1.08 million turns); linear shift term via
  `mod(kb·s0, Nf)` (exact). Hermitian spectrum → real kernel; truncated to Nk;
  `leakage` = energy lost (3.5e-13). Kernel energy Σh² ≈ 0.195 = fraction of white
  noise power passed ≈ 2·(B − 1.25·edge)/fs.

**Streaming (overlap-add).** Carry = M − 1 samples. `chooseBlockSize` requires
blockLen ≥ M − 1 (a carry only spills into the next block) and picks the power-of-two
Nfft with lowest total cost nBlocks·Nfft·log₂Nfft within `MaxMemoryGB` (default 8;
estimate 44 bytes per FFT sample). Output index = full-convolution index − s0 (fRef on
its original time). Final carry flushed; energy delayed past the file end is dropped.
Output length = input length.

**Options:** `RefFreq`, `EdgeFrac` (0.02), `GuardTime`, `MaxLeakage` (1e-8),
`BlockLen` ([] = auto), `MaxMemoryGB` (16), `SaveInfo`, `Verbose`.

**Info:** DM, fLow, fHigh, refFreq, dispersionConst, delayConvention, smearTime,
bulkDelayRef, edgeWidth, kernelLen, kernelZeroLag, kernelLeakage, kernelEnergy, Nfft,
blockLen, nBlocks, memEstimateGB, …

**Run (L = 1 s):** Nfft = 2^28, blockLen = 243,205,685 (91 %), 17 blocks, ~11.8 GB,
81.5 s. The `MaxMemoryGB` default is now 16 GB in the code (machine: 24 GB unified
memory), which covers this 11.8 GB run.

### 5.3 `addNoiseAndRFI(inFile, outFile, fs, ...)` + `rfiSource(type, ...)` [validated: noise −5 dB MC; all RFI types 5 Oct, §7]

**Physics / placement.** Noise and man-made interference are added **at the receiver
input (RF), after dispersion** → they are not dispersed. Dedispersion later applies the
inverse chirp to them, exactly as with real data: a broadband impulse gets smeared
(inverse sweep, ~6.3 ms), a narrowband carrier is only shifted by its frequency's delay.

**SNR definition (what `snrDB` in main means).**
SNR = S_peak / N: pulsar power spectral density at the pulse peak over receiver-noise
PSD, in the band (= S_peak/SEFD in radio-astronomy terms). Independent of bandwidth,
time resolution and dispersion. With the generator's white signal (std A at peak) and
white receiver noise (std σ_n), both at f_in: SNR = A²/σ_n², i.e.
σ_n = A / √(10^(SNRdB/10)). `Inf` → no noise.

**Predictions** (printed, stored in `info.prediction`; radiometer equation for a
Gaussian power profile p(t), incl. self-noise, flat band B = fHigh − fLow, ρ = 10^(SNRdB/10)):
- SNR_pulse² = B·∫(ρp)²/(1 + ρp)² dt; weak pulses: SNR_pulse ≈ ρ·√(B·σ_t·√π) ≈ 388·ρ.
- 1/σ_TOA² = B·∫(ρp′)²/(1 + ρp)² dt; weak: σ_TOA ≈ √2·σ_t / SNR_pulse.
- Folded over N pulses: SNR × √N, σ_TOA / √N.

| snrDB = −20 | predicted |
|---|---|
| SNR per pulse | 3.88 |
| SNR per 10-turn sub-integration | 12.3 |
| SNR all ~100 pulses | 38.8 |
| best TOA error per pulse / sub-int / file | 77 µs / 24.5 µs / 7.7 µs |

These come from `predictPerformance` with a **flat** band B = 400 MHz: a quick
approximation (labelled as such in the print and in `info.prediction.note`), useful for
choosing `snrDB` before anything runs. It is **~1.5 % optimistic** (SNR too high, σ_TOA
too low) at every SNR, because the band tapers pass fewer independent samples
(noise-dominated limit: effective B = (∫W_f²W_i²)²/∫W_i⁴ = 388.4 MHz). The exact
optimum (same Fisher integral with the real tapers) is `expectedPowerModel`
`snrPulse` / `toaErrPulse` (§5.10), printed by `main.m` and `runMonteCarlo.m` since
5 Oct 2026. At −20 dB: 3.79 per pulse, 79.0 µs; at −5 dB: 96.5, 2.892 µs.

**RFI model (`rfiSource`).** Power `INRdB` = RFI power while on, relative to **all**
receiver-noise power in the analysis band (0 dB = as much as the whole in-band noise).
- `'cw'`: carrier, `Freq`, `Drift` [Hz/s], `Phase`.
- `'bpsk'`: GNSS-like spread spectrum, `Freq`, `ChipRate`; random ±1 rectangular chips.
- `'pulsed'`: radar, `Freq`, `PulseWidth`, `PRF`, `ChirpBW` (linear chirp in the pulse),
  `StartTime`.
- `'impulse'`: broadband white bursts, Poisson `Rate`, `Duration`.
Amplitudes: carrier types a = √(2·INR·P_nb) (power a²/2); impulse σ = √INR·σ_n;
P_nb = σ_n²·B/(fs/2). Phases from absolute sample index; chips from `mrg32k3a`
substreams per 65,536-chip chunk; impulse events generated once for the whole file,
each burst from its own substream → block-size independent.

**Scenario in `pipelineParams`** (illustrative L-band; `rfiOn` switch; `rfiSelect` picks
sources by index, 1 L1, 2 L2, 3 radar, 4 carrier, 5 impulses, [] = all):

| Source | Type | Parameters | INR |
|---|---|---|---|
| GNSS L1-like | bpsk | 1575.42 MHz, 1.023 Mchip/s | −10 dB |
| GNSS L2-like | bpsk | 1227.60 MHz, 10.23 Mchip/s | −15 dB |
| L-band radar | pulsed | 1300 MHz, 2 µs, PRF 373 Hz, 1 MHz chirp | +20 dB when on |
| Spurious carrier | cw | 1350 MHz | −5 dB |
| Broadband impulses | impulse | 200 ns, 50 /s | +20 dB |

Radar PRF is deliberately **not** a harmonic of the pulsar frequency (400 Hz would fold
coherently into the profile and bias TOAs — a real hazard to test once excision exists).

**Options:** `Band`, `SNRdB`, `SignalInfo` (info_gen) or `SignalAmp`, `NoiseStd`, `RFI`,
`Seed`, `BlockSize`, `SaveInfo`, `Verbose`. Warns below −60 dB (single-precision limit).

**Info:** noiseStd, noisePSD, noisePowerBand, snrDB, snrDefinition, noiseAddedAt, rfi
(with computed amplitudes and parameters), prediction, seed.

### 5.4 `applyIQmodulation(inFile, outFile, f_in, fs, fLO, ...)` [validated noise-free]

**Physics.** Analog I/Q receiver: mix with LO, low-pass, sample. Complex envelope
z(t) = LPF{2·x(t)·e^(−i2π·fLO·t)}, so x(t) = Re{z(t)·e^(i2π·fLO·t)}.
- Mapping f_baseband = f_RF − fLO, **no spectral inversion** (1.2–1.6 GHz → ±200 MHz).
- The mirror image (−1.6…−1.2 GHz → −3.0…−2.6 GHz ≡ +1.0…+1.4 GHz at 4 GHz sampling)
  lies far outside the pass band and is removed.
- Factor 2 (`'Gain','envelope'`, default): compensates the discarded image; tone
  A·cos(2πft + φ) → A·e^(iφ)·e^(i2π(f−fLO)t); mean|z|² = 2·mean(x²) in band.
  (`'unity'` = old behaviour, amplitudes halved.)

**Code.**
- D = round(f_in/fs) = 8.
- Kaiser-windowed sinc, **odd length** (2049 taps) → integer group delay G = 1024 input
  samples (256 ns), **removed exactly**. Cutoff default 0.9·fs/2 = 225 MHz, stopband
  80 dB (β = 7.86), transition ≈ 9.8 MHz. Checks: stopband start (≈230 MHz) < output
  Nyquist (250 MHz); `'Band'` inside flat pass band (±200 < ±220 MHz).
- LO phase in cycles `mod(n0·r, 1) + (0:n−1)·r`, r = fLO/f_in = 0.35, from the absolute
  index (no drift).
- Overlap-add with carry 2048; Nfft = 2^22 (auto; ~350 MB), 955 blocks for 1 s.
- `writeDecimated`: removes G and keeps every D-th sample on the **global** grid →
  output sample k ↔ input sample k·D ↔ time k/fs exactly. Output length ceil(N/D)
  (500,000,001 for 1 s).

**Options:** `FilterOrder`, `Cutoff`, `StopbandDB`, `Band`, `Gain`, `BlockLen`,
`MaxMemoryGB` (2), `SaveInfo`, `Verbose`.

**Info:** format cf32, N, fs/actualFsOut, fsIn, decimationFactor, fLO, freqMapping,
loPhaseRef, gainConvention, gainFactor, cutoff, transitionWidth, stopbandDB, filterLen,
groupDelayRemoved, t0 = 0, Nfft, blockLen.

### 5.5 `applyInverseDispersion(inFile, outFile, fs, fLO, DM, fLow, fHigh, ...)` [validated]

**Physics.** Exact inverse of 5.2, applied at baseband: per bin with RF frequency
f = fLO + f_bb inside the band,
H_inv = W(f)·exp(−i·2π·K·DM·(f − fRef)²/(f·fRef²)); 0 outside. Group delay
−(τ(f) − τ(fRef)): lower frequencies advanced (up to 6.3 ms) → non-causal. Coherent
(acts on voltages): removes dispersion exactly at every frequency, unlike incoherent
(channelized) dedispersion, which leaves intra-channel smearing.

**Code.**
- Checks: band strictly inside fLO ± fs/2 (1.15–1.65 GHz at 500 MHz); fRef in band.
- Overlap: nFuture = ceil((τ(fLow) − τ(fRef))·fs) + g = 3,152,472;
  nPast = ceil((τ(fRef) − τ(fHigh))·fs) + g = 1250 (g = 2.5 µs·fs).
- H written directly on the FFT grid (in-band signed bins, chunks of 4M, double phase,
  `idx = mod(ks, Nfft) + 1`); complex data → no Hermitian symmetry.
- **Leakage check**: one `ifft(H)`; energy outside lags [−nFuture, +nPast] must be tiny
  (run: 0.0e+00; warn > 1e-6).
- FFT size: candidates from 2^nextpow2(2·overlap) to whole-file, cost
  nBlocks·Nfft·log₂Nfft, 56 bytes/sample, `MaxMemoryGB` (4). For the reference run:
  Nfft = 2^26, step = 63,955,142 (95 %), 8 blocks, 3.76 GB (2^26 is also the unconstrained
  optimum; costs 2^25: 1.43e10, 2^26: 1.40e10, 2^27: 1.45e10).
- **Overlap-save** with one preallocated buffer [past nPast | outputs step | future
  nFuture], shifted by `step` each block; zero start edge; zero fill after EOF.
- `fullySupported = [nPast+1, N − nFuture] = [1251, 496,847,529]`: outputs that saw
  their complete sweep (the last ~6.3 ms did not; the forward stage dropped that energy).

**Options:** `RefFreq`, `AllowRefOutsideBand` (6 Oct 2026), `EdgeFrac`, `GuardTime`,
`Nfft`, `MaxMemoryGB`, `CheckLeakage`, `MaxLeakage`, `SaveInfo`, `Verbose`.

**`AllowRefOutsideBand` (6 Oct 2026, for `dedisperseChannels`).** Default false: fRef must
lie in [fLow, fHigh] as before. True: any fRef > 0, so every narrow channel can be
referenced to the top of the whole band. The filter formula is exact for any fRef
(phase K·DM·(f/fRef² − 2/fRef + 1/f): the 1/f dispersion term, a linear term = pure
delay, a constant); a reference above the band only adds a pure advance
τ(f) − τ(fRef) > 0, so the overlap is then all on the future side. Code: nFuture/nPast
= max(ceil(…), 0) + g (inside the band both terms are ≥ 0 → unchanged); info field
`refOutsideBand`. Regression: the full-band path reproduces `test_IQ_dedispersed.dat`
bit for bit (`tests/testDedisperseChannels.m`, test 1).

**Not checked internally** (see §9): correctness of DM, reference-frequency convention
vs TOA definition, sideband/frequency mapping of the receiver, true fs, limits of the
cold-plasma model.

**Run (1 s):** 7.4 s.

### 5.6 `detectPower(inFile, outFile, fs, fLO, f_out, ...)` [validated]

**Physics.** After dedispersion the pulse is still noise; only the variance carries the
signal → square-law detection |z|² = I² + Q², averaged per time bin. For
variance-modulated Gaussian noise, a profile-weighted sum of |z|² is the optimal (NP)
statistic, so square-law + folding + matched filtering is the optimal order. Linear,
phase-sensitive steps (dedispersion, beamforming, RFI excision) must happen **before**
squaring.

**Code.**
- binLen = round(fs/f_out) = 500 → 1 µs bins, 1,000,000 bins for 1 s.
- Fast path (1 channel, no band): I² + Q² from the raw read (no complex arrays, no
  sqrt), per-bin sums in double.
- Channelized path (`NChan` > 1 + `Band`): rectangular non-overlapping frames whose length
  divides binLen (≥ `FineBinsPerChan` = 8 FFT bins per channel, FFT-friendly sizes),
  |Z|² summed into channels by one matrix (incl. Parseval 1/(L·binLen)); channels of a
  full-baseband band sum exactly to the total. At 1 MHz bins: up to 32 channels clean.
- Timing: bin k (1-based) time = centroid of its samples:
  t_k = T0 + ((k−1)·binLen + (binLen−1)/2)/fs = `binTime0 + (k−1)·binDt`
  (binTime0 = 0.499 µs). Boxcar → no bias.
- `fullySupportedBins` from `FullySupported` samples: [4, 993,695] for the reference run.

**Info:** layout, nChan, N (bins), fs = actualFout, fsIn, binLen, binDt, binTime0,
chanFreqs, chanEdges, chanFftBins, frameLen, chanMeanPower, samplesDropped,
fullySupportedBins.

**Noise per bin:** relative std 1/√(B_noise·Δt) ≈ 1/√390 ≈ 5.1 % (1 µs, noise-free
case; B_noise ≈ 389.6 MHz for both tapers). Run: 1.2 s for 1 s of data (I/O bound).

### 5.7 `foldProfile(info_det, outFile, f0, ...)` [validated]

**Idea.** Average all rotations at equal rotational phase; noise drops as √(number of
turns). Phase model φ(t) = Phi0 + f0·(t − TRef) + ½·F1·(t − TRef)². With TRef = first
true pulse centre (T/2) and Phi0 = 0, every true pulse peak is at an integer phase.

**Code.**
1. Only `fullySupportedBins`.
2. Phase per time bin (double).
3. `turn = round(φ)`, `frac = φ − turn` ∈ [−0.5, 0.5]: each turn runs from half a turn
   before to half a turn after its pulse → cuts fall between pulses.
4. 1024 phase bins; bin j (1-based) centred at phase (j−1)/1024; negative phases wrap to
   the end of the array. Position x = frac·1024.
5. **Linear assignment** (default): j0 = floor(x), a = x − j0; weights 1 − a and a to
   bins j0 and j0+1. Preserves each sample's mean position exactly; `'nearest'` gave a
   0.5 ns bias (1 µs time grid vs 9.765625 µs phase grid repeat identically every turn).
6. Accumulated per (phase bin, sub-integration): **sum** Σw·P, **weight** Σw,
   **weight2** Σw² (variance of a weighted mean = σ²·Σw²/(Σw)²), **weightX**
   Σw_left·w_right (covariance of neighbouring bins that share time bins).
7. Sub-integrations of `SubintPeriods` turns; `profTotal` = Σsum/Σweight (exact).
8. Per sub-integration: `turnRef` (middle turn), `tRef` = time where φ = turnRef
   (stable quadratic root; F1 = 0: tRef = TRef + (turnRef − Phi0)/f0), `fRef` = spin
   frequency at tRef, `tMean`, `nTimeBins`.
9. Saved to `.mat` (fold + info). Vectorized `accumarray`, chunked reading; 0.05 s.

**Why weightX matters:** with linear assignment a phase bin looks like it has
σ²/14.7 variance per turn (instead of σ²/9.77) because it shares data with its
neighbours; ignoring the positive covariance underestimates the variance of any sum over
bins by ~1.5× (fold check predicted 0.226 µs instead of the correct 0.277 µs until fixed).

### 5.8 `estimateTOA(fold, info_fold, template, ...)` [validated noise-free]

**Model.** p(φ) = a + b·s(φ − τ) + noise; s = template (peak 1 at phase 0); τ = phase
offset from the ephemeris prediction. **TOA = tRef + τ/fRef** = arrival time of template
phase 0 at the dedispersion reference frequency (1.6 GHz).

**FFTFIT** (Taylor 1992):
- X_k = P_k·conj(S_k), harmonics k = 1…floor(N/2)−1 (DC carries only the baseline;
  Nyquist ambiguous). C(τ) = Re Σ X_k·e^(i2πkτ): continuous cross-correlation (a shift is
  a phase ramp in the Fourier domain).
- Coarse: zero-padded inverse FFT on 8×N points (1/8192 turn = 1.22 µs).
- Refine: Newton on C′(τ) = 0 with analytic C′, C″; steps clamped to one coarse step;
  C″ < 0 required; stop at 1e-14 turns.
- b = C(τ̂)/Σ|S_k|²; a = (P₀ − b·S₀)/N; τ wrapped to [−0.5, 0.5).

**Uncertainty.** Noise model `'radiometer'` (default): per time bin var = m²/(B_noise·Δt),
m = a + b·s the **fitted** model (observer-usable; with receiver noise m includes the
system-power baseline). Per phase bin: var = σ_tb²·W2/W², covariance with neighbour via
weightX. τ solves C′(τ) = Σ d_j·p_j (d from the template slope, largest on the pulse
flanks) → var(τ) = dᵀ·Cov·d / C″². Same for σ_b via C(τ). `'offpulse'` model: variance
and lag-1 covariance from off-pulse residuals (needs receiver noise).

**Main passes** `Bnoise` = noiseBandwidth of the **observer's own dedispersion taper
only** (`info_dedisp.edgeWidth`) ≈ 391.6 MHz, always, regardless of whether noise is
present (ground-truth rule: the observer does not know the dispersion-stage taper).
Receiver noise sees only W_inv, so this is exact for it; when the pulsar dominates
(high SNR, in-pulse) the true value is ~0.5 % lower (both tapers, ≈ 389.6 MHz).
(Before commit f4ac189 main switched to both tapers in the noise-free case.)

**Outputs** per sub-int: valid, coverage (default `MinCoverage` = 1: complete turns
only), phase, phaseErr, toa, toaErr, amp, ampErr, baseline, snr (= b/σ_b), redChi2
(on-pulse bins, model > 1 % of peak), flagChi2, tRef, fRef, turnRef; `toa.total` (whole
fold): phase, phaseErr, timeOffset (= phase/f0, relative to the ephemeris, **not** an
absolute TOA), timeOffsetErr, amp, ampErr, snr, redChi2, flagChi2.

**Quality flag** (6 Oct 2026): option `MaxRedChi2` (default 2) → `flagChi2` = redChi2 > 2.
Red. χ² = Σ (p_j − m_j)²/v_j over the on-pulse bins / (n − 3): normalized residuals
squared and averaged, ≈ 1 ± √(2/n) (≈ ± 0.09 for ~260 bins) when the profile is template +
noise; ≫ 1 when something else is in it. The error bar only measures the sharpness of the
correlation peak and cannot see a wrong peak; χ² checks the fit itself (radar lock: χ²
939–1625 with σ_TOA 1 µs). Threshold 2 ≫ clean scatter and the known high-SNR excess
(1.05–1.1); not sensitive to moderate RFI (χ² ~1.3), by choice.

**Optimality note.** FFTFIT weights by the template slope: optimal for white
(receiver-dominated) noise. With pure self-noise (variance ∝ profile²) a flatter
weighting (centroid) is better: 0.28 µs vs 0.47 µs per pulse in the noise-free run, and
the reduced χ² comes out ~1.06–1.09 (reproduced by simulation with a perfect noise
model; excess from low-level tail bins). Both effects should disappear once receiver
noise dominates. A weighted fit using the radiometer model would be optimal in both
regimes (todo).

**One TOA vs many:** one profile of all pulses has the same precision as the weighted
mean of per-pulse TOAs (constant profile, exact phase model). Sub-integrations are needed
to validate error bars, see drift (phase-model errors, motion, clock, **position error**)
and catch bad data. For a single TOA set `subintPeriods` ≥ number of turns (e.g. 1e6;
`Inf` is rejected).

**Template.** Main uses `gaussianTemplate(nBin, ephem.profileFWHM)`. `estimateTOA`
accepts any [NBin × 1] profile with phase 0 at bin 1, so real data with other pulse
shapes needs only a different template source (multi-component model or high-SNR
observed profile); its phase 0 then defines the TOA.

### 5.9 `validateTOA(toa, info_gen, ...)` [validated]

Matches each valid TOA to the nearest true pulse centre; errors, normalized errors
z = err/toaErr; mean ± expected, rms vs predicted rms (ratio), χ² with both tail
probabilities (`gammainc`), fractions |z| < 1 and < 2 (expect 68.3 % / 95.4 %), total-fold
offset in σ. Figure: errors with ±1σ, histogram vs N(0,1), Q-Q plot, amplitude and SNR per
sub-int. **Monte Carlo:** pass struct arrays `[toa1 toa2 …]` and `[info_gen1 …]`; all
TOAs are pooled.

### 5.10 Check / plot functions (ground truth)

- **`expectedPowerModel(info_gen, info_disp, info_IQ, info_dedisp, info_rx)`** [written]:
  complex baseband PSD S = s_s·p(t)·Wf²Wi² + s_n·Wi² (noise added after dispersion sees
  only Wi). Mean: Ps = sigScale·p(t), Pn; variance of a bin mean over dt:
  (css·Ps² + 2·csn·Ps·Pn + cnn·Pn²)/dt with css = ∫(WfWi)⁴/Js², csn = ∫Wf²Wi⁴/(Js·Jn),
  cnn = ∫Wi⁴/Jn², Js = ∫(WfWi)², Jn = ∫Wi². sigScale = g²A²Js/fsIn, Pn = g²σ_n²Jn/fsIn.
  Variance formula verified numerically (0 ≤ ρ ≤ 100, within 1 %). RFI is **not** in the
  model (shows up as excess in the checks).
  **Best achievable** per pulse (added 5 Oct 2026; Fisher information of the detected
  power, total-power detection, optimal time weighting): with
  v = css·Ps² + 2·csn·Ps·Pn + cnn·Pn², SNR² = ∫Ps²/v dt, 1/σ_TOA² = ∫Ps′²/v dt
  (dt cancels) → `M.snrPulse`, `M.toaErrPulse` (NaN without noise: the integrals diverge
  in the tails). Flat band (all c = 1/B) reduces to `predictPerformance`. Same Gaussian
  profile and ±8σ grid as `predictPerformance`. Checked: 96.5 / 2.892 µs at −5 dB in
  MATLAB = Python port.
- **`plotDetectedPower`** [validated noise-free; noise version written]: overview,
  optional channel waterfall, zoom vs expected, normalized residual (all bins with noise;
  on-pulse only without). Centroid checks on (P − measured off-pulse baseline). Prints
  per-pulse offsets (≤ 20 pulses listed), summaries, normalized-residual stats on- and
  off-pulse. `DisplayAverage` (auto) averages bins for display only. Output `check` is a
  struct (`pulses`, `baseline`, `displayAverage`, `zAll`).
- **`plotFoldCheck`** [validated noise-free; noise version written]: total profile vs
  expected (kernel-averaged: time-bin boxcar + assignment kernel), zoom, normalized
  residual, sub-int stack (baseline-subtracted), per-sub-int centroid offsets with
  covariance-correct errors (`combNoise`).
- **`plotIQCheck`** [validated noise-free; patched for noise]: dynamic spectra of dispersed
  and dedispersed IQ with expected sweep, power vs expected, mean spectra, raw I/Q zoom,
  per-pulse centroid/energy.
- **`plotDispersionCheck`** [validated]: raw vs dispersed dynamic spectra with expected
  sweep, band-integrated power.

### 5.11 `runMonteCarlo.m` [validated, −5 dB, 10 seeds]

**Idea.** The pulsar signal is fixed (sky files from a `main.m` run, same `seed`); only
the receiver noise changes. The spread of TOAs between realizations is exactly what the
radiometer error bars of `estimateTOA` claim to describe, so the pooled normalized errors
z = (TOA − truth)/toaErr must be N(0,1), to within ~1/√(number of TOAs).

**Code.**
- `nSeeds` at the top (default 10); noise seeds `seed + (1:nSeeds)`. Seed 1 equals
  main's `noiseSeed = seed + 1`, so its TOAs must be identical to main's (built-in check).
- `pipelineParams;` for all parameters; `info_gen`, `info_disp` via `loadInfo`. Errors
  if the sky files do not match T, L, f_in, seed, DM or band.
- Per seed: the same stage calls as `main.m` (`addNoiseAndRFI` → `applyIQmodulation` →
  `applyInverseDispersion` → `detectPower` → `foldProfile` → `estimateTOA`), with
  `'Verbose', false`, fold not saved, intermediate files in `data/mc/` (overwritten each
  seed, so main's files are never touched). One progress line per seed.
- `validateTOA(toaAll, repmat(info_gen, 1, nSeeds))` pools all TOAs; the radiometer
  prediction is printed at the end.
- Cost at L = 0.1 s: roughly 15–20 s per seed.

### 5.12 `runSNRSweep.m` [validated, 5 Oct 2026]

**Idea.** Phase D: how TOA precision and error-bar honesty change with SNR, against the
exact optimum (`expectedPowerModel` `toaErrPulse`). Unlike `runMonteCarlo.m` the **pulsar
realization changes too**: at high SNR the self-noise would otherwise be identical in
every realization and missing from the scatter.

**Code.**
- Top: `nReal` (10), `snrList` ([−25 −20 −15 −10 −5 0 10 20] dB), `replotOnly` (true:
  only redraw the figure from `data/mc/snrSweep.mat`).
- Outer loop j: new sky (generator + dispersion) with pulsar seed `seed + 100·(j−1)`,
  noise seed = pulsar seed + 1 (j = 1 = main's 42/43: built-in check). Inner loop over
  `snrList`, same sky and same noise seed (noise only rescaled: common random numbers,
  so SNR points are correlated, TOAs within one point are independent across j).
- Same stage calls as `main.m`; files `data/mc/sweep_*` (main's / MC files untouched).
  Optimum per sub-int from `expectedPowerModel` (once per SNR, j = 1).
- Per SNR point: pooled `validateTOA` (`'Plot', false`), plus robust σ = 1.4826·median|err|
  and outlier fraction |z| > 5. Summary table printed; results saved to
  `data/mc/snrSweep.mat` (res, toaC, genAll, …).
- Figure (local `plotSweep`): σ_TOA vs SNR (rms, robust, predicted, optimum); ratios
  rms/pred and pred/opt on a linear 0.8–3 axis (off-scale points on the edge, labelled);
  |z| < 1 and outlier fraction vs best SNR per sub-int.
- Cost: ≈ 7 s sky + ≈ 8 s per SNR value per realization; full run 11.6 min.
- Phase E additions (6 Oct): `detectPulsar` per (realization, SNR) (`detC`); a true **H0
  pass** per realization: zero sky (`generatePulsarSignal` with A = 0, made once,
  `sweep_zero.dat`) + receiver noise (`'NoiseStd'` 1, the realization's noise seed),
  processing, only `detectPulsar` (`detH0C`, `res.h0`). (A very low SNR such as −50 dB is
  not H0: it is a real small-antenna operating point.) Second table: T0 mean/std, P_D
  known/unknown vs theory Q(η − SNRopt), χ² flags, good TOAs, outliers and rms/pred among
  good TOAs, H0 line. Figure 2×2 with a detection panel (H0 false alarms in its title).
  Full run 13.0 min.

### 5.13 `detectPulsar(fold, info_fold, template, ...)` [validated, 6 Oct 2026]

**Question.** Is the pulsar there? Neyman–Pearson test of H0 (baseline + noise) against
H1 (baseline + template-shaped pulse + noise) at a fixed false-alarm probability P_FA
(default 1e-3 per profile), for every sub-int and the total fold. Observer-only (fold,
template, B_n, P_FA). Separate from the TOA error bar, which assumes the found peak is the
pulse; the detector asks whether the best peak stands out more than noise ever would.

**Code / maths.**
- H0 noise normalization: baseline a = mean(p) (exact under H0), per time bin var
  a²/(B_n·Δt), per phase bin v_j = that·W2_j/W_j², neighbour covariance from weightX (as
  `estimateTOA`); `normProfile` = (p − a)/√v.
- Known phase (`'Phase'`, default 0 = ephemeris): T0 = cᵀ(p − a)/√(cᵀ·Cov·c), c = template
  at that phase minus its mean. H0: N(0,1) → η₀ = Φ⁻¹(1 − P_FA) = 3.09.
- Unknown phase: T(τ) for all bin shifts (numerator and variance by FFT correlation),
  Tmax, phaseMax. Look-elsewhere: Rice, P_FA(η) ≈ Q(η) + √λ₂/(2π)·e^(−η²/2),
  λ₂ = Σ(2πk)²|S_k|²/Σ|S_k|² (k ≥ 1) → Rice factor 5.51, η = 4.15.
- Noise check: off-pulse variance (template at phaseMax < 1 % of peak) / radiometer
  variance with the off-pulse baseline (`noiseRatio`; carrier ~1.65, radar ≫ 1).
- Output `detection` (not `det`: would shadow `det()`), `info` (PFA, lambda2, riceFactor).

**Validation.** NumPy (20,000 H0 profiles with neighbour correlation): T0 mean 0.014, std
0.996; P(Tmax > η) 0.0098 at 1e-2, 0.0009 at 1e-3 (Rice slightly conservative far from
the tail: 0.26 at 0.3). MATLAB −5 dB: 9/9 detected, noise ratio 0.992, T0 114.5 (> SNR
93.4 by ×1.23: T0 uses the H0 noise and ignores self-noise; irrelevant for detection and
for the small-antenna regime). Sweep results in §7. Long H0 run (§5.14): T0 std 1.008 ±
0.024, false alarms at the nominal rates down to 1e-2…1e-3.

### 5.14 `runH0.m` [validated, 6 Oct 2026]

**Idea.** Many independent noise-only (H0) profiles to test the H0 side of `detectPulsar`
and the noise covariance model behind it (and behind the TOA error bars): T0 ~ N(0,1),
exceedance rates at several P_FA, and the correlation of `normProfile` between phase bins
(model: lag 1 from weightX, 0 beyond).

**Code.** `nPass` (100) passes; zero sky (generator with A = 0, made once,
`h0_zero.dat`); per pass receiver noise (`'NoiseStd'` 1, seed `seed + 10000 + i`, no
overlap with the sweeps) → same processing as main → `detectPulsar` (9 complete turns per
pass at L = 0.1 s). Collects T0, Tmax, noiseRatio, the circular autocorrelation of each
`normProfile` (FFT) and the model lag-1 correlation WX/√(W2_j·W2_{j+1}). Implied T0 std
from the measured correlation: √(Σ R_c·r / Σ R_c·r_model), R_c = template
autocorrelation. Exceedances at P_FA 1e-1, 1e-2, 1e-3 (Φ⁻¹ and Rice thresholds). Saves
`data/mc/h0Run.mat`; figure: T0 histogram vs N(0,1), P(statistic > η) vs theory (log),
correlation vs lag. Files `data/mc/h0_*`; ~8 s per pass (100 passes ≈ 14 min).

### 5.15 `channelizeIQ(inFile, outBase, fs, fLO, fLow, fHigh, ...)` [validated 6 Oct 2026 (Jasper's run); channelized front end, unit 1]

**Idea.** First stage of the channelized front end (§10 agreed order item 1): split the IQ
stream into nChan channels that tile [fLow, fHigh] exactly, each a complex voltage stream
at a reduced rate, so RFI excision and coherent dedispersion can work per channel.
Channel j is what a separate receiver would give: mix down by its centre
fc_j = fLow + (j − ½)·ChanWidth (LO phase 0 at input sample 0), zero-phase low-pass,
keep every D-th sample: y_j(m) = Σ_n h(mD − n)·x(n)·e^{−i2π(fc_j − fLO)n/fs}.
- Sample m at t0 + mD/fs (group delay removed → no TOA offset); RF = fc_j + f_bb (no
  inversion) → each channel file is a normal IQ file with fLO = fc_j, fs = fs/D.
- Unit DC gain: tone amplitude and passband PSD preserved; the channels' parts inside
  ±ChanWidth/2 add up to the IQ band power.
- **Oversampled** (default 4/3): K = fs/ChanWidth = 256, D = K·3/4 = 192, channel rate
  4.17 MHz. Prototype flat to ±ChanWidth/2 (1.5625 MHz), stopband from fs/D − ChanWidth/2
  (2.604 MHz) → decimation aliases only into the transition region, never into the useful
  band; the transition is removed later (per-channel dedispersion keeps ±ChanWidth/2).
  A critically sampled filterbank (D = K) would alias into the channel edges.

**Method (polyphase form).** With fc_j − fLO = (k_j + β)·fs/K (β = ½ here, band edges on
the K-grid) and prototype h(r) centred at G = (M−1)/2, zero-padded to L = P·K:
y_j(m) = e^{−i2π(k_j+β)(mD+G)/K} · Σ_p e^{+i2πk_j p/K}·u_m(p),
u_m(p) = Σ_q hb(p+qK)·x(mD+G−p−qK), hb(r) = h(r)·e^{+i2πβr/K}: weight a frame of L
samples, fold into K bins, one K-point IFFT for all channels, phase correction from the
absolute index (double). Same arithmetic as K mix-filter-decimate receivers at ~1/K cost.
Prototype: Kaiser-windowed sinc (copy of `applyIQmodulation`'s `kaiserLowpass`), cutoff
fs/(2D), transition fs/D − ChanWidth → 3853 taps (16 per branch) at 80 dB.

**Files.** `<outBase>_ch001.dat` … (cf32 per channel), `<outBase>_info.mat`
(`loadInfo(outBase)` works). Info: nChan, chanFreqs, chanWidth, fs (channel rate), K,
decimation, oversampling, prototype (taps), passband/stopband edges, groupDelayRemoved,
t0, fullySupported (1-based outputs that see only real data: inputs mD ± G inside the
file), conventions. Edge outputs see zero-padding.

**Tests (`tests/testChannelizeIQ.m`, run by Claude and Jasper 6 Oct 2026, identical results).**
- Tones (4 tones, channels at, beside and far from them): every sample equals the exact
  prediction A·e^{iφ}·Hc(Δf)·e^{i2πΔf·t_m} to ≤ 2.7e-7 (single precision) → frequency
  mapping, LO phase, timing and gain exact. |Hc| = 1 ± 1.2e-4 to ±1.5625 MHz; stopband
  −80.2 dB.
- White noise (2^23 samples): channel power / Σh² = 1.0007 (per channel 0.988–1.016,
  1σ 0.55 %); passband PSD level 1.0004 of expected; flat within ±3 % (Welch scatter).
- Real data (`test_rx_IQ.dat`, seed 43, −5 dB): Σ channel band powers / IQ band power
  = 1.00007; 0.1 s → 128 × 416,667 samples in 4.1 s (`data/chan/`, 427 MB).

### 5.16 `dedisperseChannels(info_chan, outBase, DM, ...)` [validated 6 Oct 2026 (Jasper's run); channelized front end, unit 2]

**Idea.** Coherent dedispersion per channel of `channelizeIQ`, reusing the validated
`applyInverseDispersion` on each channel file (an ordinary IQ file, fLO = channel
centre, fs = 4.17 MHz):
- band = channel centre ± ChanWidth/2: the chirp filter is zero outside, which removes
  the channelizer's transition region → channels tile the band exactly;
- RefFreq = one common frequency (default top of the band, 1.6 GHz) with
  `AllowRefOutsideBand` (§5.5): each channel's chirp includes its delay to the reference
  → all channels aligned at 1.6 GHz, no separate shift step;
- taper (`EdgeFrac` 0.02 of the channel) inside each channel: neighbours do not overlap →
  independent noise; ∫W² = 390.0 MHz, the same as the full-band path (0.02 of 400 MHz).
- Memory: kernel per channel = sweep to the reference at 4.17 MHz (DM 100 bottom channel
  ~5e5 samples, ~100 MB) instead of ~1e8 at 800 MHz → processing side of §10 item 5
  solved. DM 5: Nfft 2^14…2^18, max leakage 6e-8, 128 channels in 1.7 s.

**Info.** chanFiles, chanFreqs, chanWidth, fs, N, t0, DM, refFreq, edgeWidth, per-channel
nPast / nFuture / Nfft / leakage, fullySupported (valid in every channel: channelizer
range shrunk by each channel's filter reach), BnoiseChan (3.0596 MHz) and
BnoiseTotal = nChan·BnoiseChan (391.6 MHz, as the full band), bulkDelayRef.

**Tests (`tests/testDedisperseChannels.m`, run by Claude and Jasper 6 Oct 2026, identical results).**
1. Regression: full-band `applyInverseDispersion` with defaults bit-identical to
   `test_IQ_dedispersed.dat` (nPast, nFuture, Nfft equal).
2. Commutation: channelize(full-band dedispersed) vs dedisperseChannels(channelize(IQ))
   inside each channel's flat band: |ρ| ≥ 0.999993, phase ≤ 1.2e-5 rad, power ratio
   1 ± 3e-5 (channels 4–125; the outer 3 see the full-band 8 MHz taper) → chirp
   segments, inter-channel delays and phase exact.
3. End to end (seed 43, −5 dB; per-channel `detectPower` at fs/4 = 1.0417 MHz, channel
   powers summed, fold, TOAs): 8 TOAs, channelized − full band mean +0.19 µs, rms
   0.51 µs (0.18 σ). Expected: the paths keep slightly different parts of the same
   noise (channel-edge tapers vs the full-band taper), detected-power noise correlation
   ρ = ∫(W_i·W_ch)²/√(∫W_i⁴·∫W_ch⁴) = 0.979 → each difference scatters with
   σ·√(2(1−ρ)) = 0.205 σ. Test: χ² of the normalized differences 5.98 for 8 (99.8 %
   range 0.86–26.1), mean +0.32 (limit ±1.16) → consistent; a processing error would
   show as an offset many standard errors large (or in test 2). Median SNR 91.79 vs
   93.36 (ratio 0.983, expected 0.979: in the
   *simulation* the pulsar also passed the forward-dispersion 8 MHz taper, which overlaps
   the full-band dedispersion taper, so the channel tapers cost 2.1 % of the signal; on
   real data both paths collect the same 390.0 MHz); fold noise ratio 0.955 vs 0.990
   (scatter ~0.033 each).
4. Time-bin noise of the summed channel power: var / radiometer 0.884 (predicted 0.888
   ± 0.005 from the channel spectrum), lag-1 correlation 0.054 (0.052 ± 0.004) → the
   narrow-channel correlation-time effect (§9) is real and modelled exactly.

### 5.17 `powerCovariance(chanWidth, edgeWidth, fs, nPerBin, ...)` + `detectChannels(info_dc, outFile, f_out, ...)` [validated 6 Oct 2026 (Jasper's run); unit 3b]

**`powerCovariance`: exact noise statistics of detected time bins.** A bin is the mean of
|y|² over n samples; for Gaussian y with spectrum S (normalized autocorrelation R), relative
to the radiometer value m²/(B·dt): variance V = C(0)/rad, covariance with bin k+L
X_L = C(L)/rad, C(L) = (1/n²)Σ_{a,b}|R(a−b+Ln)|², rad = Σ_l|R(l)|²/n. Summed over all lags
V + 2ΣX_L = 1 exactly → the old "independent bins" model is right for long sums, not bin
by bin when B·dt is small. S = W² of the per-channel taper (chirp is all-pass, drops out
without blanking). Lmax = smallest L with V + 2Σ_{≤L}X ≥ 1 − 0.002. Values (3.125 MHz,
62.5 kHz edges, 4.1667 MHz): n = 4 (0.96 µs): V 0.8877, X 0.0458 / 0.0051 / 0.0020 …,
Lmax 7 (captured 0.9984; slow 1/L² tail from the sharp channel edges); n = 16: V 0.9625,
X₁ 0.0173, Lmax 2; full band n = 800: V 0.9986, X₁ 7e-4, Lmax 1. `MaxLag` default from the
frequency grid (no wrap-around).

**`detectChannels`: detection of all channels into one file.** Validated `detectPower` per
dedispersed channel (rate fs/round(fs/f_out), one warning if it differs: 1 MHz →
1.0417 MHz, 4 samples), interleaved block-wise into the existing layout
float32 [nChan × nBins] (temporaries deleted) → `foldProfile` reads it unchanged.
Info: detectPower's fields used by the fold (file, nChan, N, binTime0, binDt,
fullySupportedBins, byteOrder, chanFreqs, …) + chanWidth, chanMeanPower, BnoiseChan,
BnoiseTotal, `noise` (powerCovariance struct, the same for every channel and bin without
blanking; per-bin streams from the mask in step 2). 0.1 s, 128 channels: 0.5 s.

**Tests (`tests/testDetectChannels.m`).** (1) powerCovariance vs the statistics of the actual
dedispersion impulse response of channel 1: V, X identical (max |ΔX| 6e-10). (2) Rows of
the interleaved file bit-identical to detectPower per channel. (3) foldProfile on the
128-channel file: channel sum = fold of the summed power (1.2e-8). (4) Measured on the
seed-43 data, 68,803 off-pulse bins × 128 channels: V 0.8866 (pred 0.8877 ± 0.0006),
X₁₋₃ 0.0455 / 0.0053 / 0.0018 (pred 0.0458 / 0.0051 / 0.0020 ± 0.0004); variance of 5-bin
sums (~ one 4.9 µs phase bin) 0.9685 × 5·rad (pred 0.9692 ± 0.0014; the independent-bin
model says 1).

---

## 6. Key formulas (quick reference)

| Quantity | Formula |
|---|---|
| Dispersion delay | τ(f) = K·DM/f², K = 4.148808e3 s·MHz²·pc⁻¹·cm³ (4.148808e15 s·Hz² per DM unit) |
| Forward chirp | exp(+i2π·K·DM·(f − fRef)²/(f·fRef²)) |
| Group delay | τ_g = −(1/2π)·dφ/df |
| Complex envelope | z = LPF{2·x·e^(−i2π·fLO·t)}, f_bb = f_RF − fLO |
| Pulsar power after IQ | sigScale = g²·A²·Js/f_in = 4·388.4 MHz/4 GHz ≈ 0.388 (A = 1) |
| Power bookkeeping | generator 1 → after dispersion ≈ 0.195 (in band) → after IQ ≈ 0.388 |
| Bin noise (flat band) | σ/P = 1/√(B_noise·Δt) |
| Noise bandwidth | B_noise = (∫W²)²/∫W⁴ (≈ 389.6 MHz pulsar = both tapers, ≈ 391.6 MHz noise-only = dedispersion taper) |
| Weighted mean variance | σ²·Σw²/(Σw)² |
| SNR definition | ρ = S_peak/N (per Hz, in band) = A²/σ_n² |
| Per-pulse SNR (weak) | ≈ ρ·√(B·σ_t·√π) ≈ 388·ρ (0.5 ms pulse, 400 MHz) |
| Per-pulse TOA error (weak) | ≈ √2·σ_t / SNR_pulse |
| Best achievable, exact tapers | SNR² = ∫Ps²/v dt, 1/σ_TOA² = ∫Ps′²/v dt, v = css·Ps² + 2·csn·Ps·Pn + cnn·Pn² |
| Power vs statistical bandwidth | ∫W² (mean power: 390.0 / 388.4 MHz) vs (∫W²)²/∫W⁴ (variance: 391.6 / 389.6 MHz); B_n ≥ ∫W² for W ≤ 1 |
| TOA | tRef + τ/fRef |
| Detection thresholds | known phase η₀ = Φ⁻¹(1 − P_FA); unknown phase Q(η) + √λ₂/(2π)·e^(−η²/2) = P_FA (3.09 / 4.15 at 1e-3) |
| Detection probability | P_D = Q(η − SNR) (known phase; unknown phase approx.) |
| Reduced χ² | Σ (p_j − m_j)²/v_j / (n − 3), ≈ 1 ± √(2/n) for template + noise |
| Turns per sub-int | N = (SNR_target/SNR_turn)², SNR_turn ≈ ρ·√(B·σ_t·√π)·loss (≈ 373·ρ here); target 6–7 |
| Baseline-estimation loss | √(1 − (Σs)²/(N_bin·Σs²)) = 0.962 for the 5 % Gaussian template |
| FFTFIT uncertainty | σ_τ = √(dᵀ·Cov·d)/|C″| |

---

## 7. Validation record

**Noise-free, L = 1 s, fs = 500 MHz, f_out = 1 MHz, 1024 bins, 1 turn per sub-int:**

| Check | Result | Expected |
|---|---|---|
| Detected-power normalized residual (127,512 on-pulse bins) | mean −0.0016, std 1.0015 | 0 ± 0.0028, 1 ± 0.0020 |
| Per-pulse centroid offsets (99 pulses) | mean +0.005 µs, rms 0.301 µs | rms 0.277 µs (χ² ≈ 117/99, p ≈ 0.1) |
| Per-pulse energy ratio | mean 1.0000, scatter ~0.24 % | ~0.19 % |
| Fold: total centroid | +0.005 ± 0.028 µs, energy ratio 0.99999 | 0 |
| Fold: normalized residual (131 bins) | mean −0.065, std 0.973 | 0 ± 0.087, 1 ± 0.062 |
| Fold: sub-int centroid rms | 0.301 µs vs expected 0.277 µs (after covariance fix) | identical to per-pulse offsets |
| TOAs | 99/100 fitted, median SNR 448.7, median σ 0.466 µs, median red. χ² 1.086 | σ ≈ 0.469 µs (NumPy MC) |
| validateTOA | mean −0.020 ± 0.047 µs, rms 0.491 vs 0.466 µs (ratio 1.053), χ² 109.9/99 (p = 0.21), |z|<1 63.6 %, |z|<2 97.0 % | 68.3 % ± 4.7 %, 95.4 % |
| Total fold | τ = −2.04e-6 turns = −0.020 ± 0.047 µs (−0.44σ), SNR 4464 (≈ 449·√99), amp 0.3883 ± 0.0001 | amp 0.3884 |
| Reproducibility | first 9 pulse offsets identical between 0.1 s and 1 s runs (same seed) | |

**With receiver noise, −5 dB (ρ = 0.316), L = 0.1 s, fLO = 1.3 GHz, fs = 800 MHz,
2048 bins, 1 turn per sub-int, noiseSeed 43 (1 Oct 2026, single run):**

| Check | Result | Expected |
|---|---|---|
| Noise baseline | 1.232 measured | 1.233 (g²σ_n²·∫W²/f_in) |
| Detected-power normalized residual, off-pulse (73,604 bins) | mean −0.0054, std 0.9995 | 0 ± 0.0037, 1 ± 0.0026 |
| Detected-power normalized residual, on-pulse (11,592 bins) | mean −0.0065, std 0.9959 | 0 ± 0.0093, 1 ± 0.0066 |
| Fold normalized residual on / off (263 / 1785 bins) | std 0.987 / 0.990 | 1 ± 0.044 / 1 ± 0.017 |
| TOAs | 9/10 fitted (turn 9 partial → rejected), median SNR 93.4, median σ 2.880 µs, median red. χ² 0.963 | SNR 97.9 (flat 400 MHz band), best σ 2.85 µs |
| validateTOA | mean −0.21 ± 0.96 µs, rms 2.50 vs 2.88 µs (ratio 0.868; ±0.24 for 9 TOAs), χ² 6.8/9 (p = 0.34), \|z\|<1 66.7 % | ratio 1 |
| Total fold | −0.21 ± 0.96 µs (−0.22σ), SNR 280.8 | 0 |
| Check-function centroids | per pulse mean +3.7 µs (noise ~6.2 µs each), total +3.68 ± 2.06 µs (1.8σ) | 0; **resolved 6 Oct: fluctuation of seed 43** (5 seeds, see below) |

Conclusions: noise model and B_noise correct; FFTFIT within ~1 % of the radiometer
optimum at this SNR (little self-noise); off-centre baseband (fLO 1.3 GHz) mapping
correct. Reference values for MC seed 43: median SNR 93.4, median σ 2.880 µs.

**Monte Carlo, same setup, 10 noise seeds (43–52), 90 TOAs (`runMonteCarlo.m`, 1 Oct 2026):**

| Check | Result | Expected |
|---|---|---|
| Mean error | −0.000 µs | 0 ± 0.31 µs |
| rms / predicted rms | 2.762 / 2.907 µs = 0.950 | 1 ± 0.075 (−0.7σ) |
| Reduced χ² | 0.91 (90 dof) | 1 ± 0.15 |
| \|z\| < 1 | 72 % | 68.3 ± 4.9 % |
| Q-Q plot | on the line incl. tails | |

→ `estimateTOA` uncertainties validated with receiver noise at −5 dB; the noisy chain
is validated at this SNR. Mean predicted σ 2.907 µs ≈ 2 % above the flat-band optimum
(2.85 µs), but only **0.5 %** above the exact optimum with the real tapers (2.892 µs,
§5.10; computed 5 Oct): FFTFIT is nearly optimal at −5 dB.
Caveat: the seeds share the pulsar realization; its pure self-noise share is ~4 % of the
per-bin variance at the peak here, so the TOAs are nearly independent. At high SNR the
MC must also vary the pulsar seed. The check-function centroid offset (+3.7 µs) is not
tested by the MC (no plot checks there).

**Phase D SNR sweep (`runSNRSweep.m`, 5 Oct 2026): 10 pulsar + noise realizations ×
8 SNR values, L = 0.1 s, 1 turn per sub-int, 90 TOAs per point (times per sub-int, µs):**

| snrDB | SNR opt | rms | robust | pred | opt | rms/pred | pred/opt | χ²_red | \|z\|<1 | \|z\|>5 |
|---|---|---|---|---|---|---|---|---|---|---|
| −25 | 1.2 | 2430 | 1470 | 217 | 249 | 11.2 | 0.87 | 561 | 18 % | 57 % |
| −20 | 3.8 | 995 | 83.9 | 82.2 | 79.0 | 12.1 | 1.04 | 116 | 56 % | 5.6 % |
| −15 | 11.8 | 26.3 | 27.6 | 25.6 | 25.3 | 1.030 | 1.012 | 1.07 | 62 % | 0 |
| −10 | 35.4 | 8.60 | 8.53 | 8.28 | 8.28 | 1.040 | 1.000 | 1.08 | 66 % | 0 |
| −5 | 96.5 | 3.07 | 2.80 | 2.91 | 2.89 | 1.053 | 1.007 | 1.11 | 66 % | 0 |
| 0 | 216 | 1.30 | 1.25 | 1.22 | 1.17 | 1.061 | 1.050 | 1.12 | 68 % | 0 |
| +10 | 496 | 0.561 | 0.537 | 0.538 | 0.358 | 1.044 | 1.502 | 1.09 | 63 % | 0 |
| +20 | 657 | 0.489 | 0.482 | 0.472 | 0.198 | 1.036 | 2.380 | 1.07 | 61 % | 0 |

Built-in check: j = 1 at −5 dB gives median SNR 93.4, σ 2.880 µs (= main / MC seed 43).

Conclusions:
1. **Error bars honest for SNR per sub-int ≳ 10**, from receiver-noise- to
   self-noise-dominated (−15…+20 dB): rms/pred 1.03–1.06 (± 0.075), χ²_red 1.07–1.12,
   no outliers. The self-noise part of the radiometer model is now tested (pulsar varied).
   Watch item: all ratios 3–6 % above 1, but the points share the realizations (≈ one
   measurement of ~1.04 ± 0.075); consistent with the noise-free run (1.053). Settle
   with nReal ≈ 40 at two SNRs if needed.
2. **Threshold:** below SNR ≈ 10 per sub-int outliers appear (locking on a noise peak,
   errors up to ±T/2): 5.6 % at SNR 3.8 (Gaussian core still correct: robust 83.9 vs
   82.2 µs), 57 % at 1.2 (TOAs ~random; fitted SNR biased up by peak picking, median
   ~2.1). Rule: sub-int length for SNR ≥ 10; navigation needs outlier rejection.
   Transition between SNR 3.8 and 11.8 not resolved (extra points −19…−16 dB).
3. **FFTFIT efficiency:** optimal (≤ 1 %) up to −5 dB; 1.05 at 0 dB, 1.50 at +10 dB,
   2.38 at +20 dB. FFTFIT saturates at ≈ 0.47 µs (noise-free self-noise value; median
   SNR 444.7 vs 448.7 noise-free) while the optimum keeps falling (information in the far
   tails where ρp ≈ 1; needs an exactly known profile shape). Realistic small-antenna
   navigation (≈ −25…−50 dB, §9) is weak-signal → FFTFIT optimal → **phase F low
   priority**; phase E first.
4. The quick test's +10 dB ratio 0.56 (2 realizations, 18 TOAs) was a fluke.

**RFI test (`main.m`, −5 dB, L = 0.1 s, noise seed 43, one source at a time via
`rfiSelect`, then all; 5 Oct 2026).** Reference: noise only, median SNR 93.4, σ 2.880 µs.

| Run | Generator check | Effect | TOAs |
|---|---|---|---|
| Carrier 1350 MHz, −5 dB | baseline 1.632 (pred. 1.633) | off-pulse residual std 1.2885 vs 1.285 predicted from the cross term 2·P_cw·P_n (a carrier adds less variance than noise of equal power, which would give 1.32) | SNR 74.1 (pred. ~75), σ 3.69 µs (×1.28); ratio 0.805, χ² 5.9/9, total +0.33σ: valid; error bars ~3 % pessimistic expected |
| Radar 1300 MHz, +20 dB | pulses at k/373 s; after dedispersion 4.17 ms earlier (τ(1.3) − τ(1.6)), ~20 µs streaks (1 MHz chirp × 18.9 µs/MHz; peaked, time–bandwidth product 2); median baseline unchanged (1.236) | spikes ~21 in 1 µs bins; energy ≈ one pulsar pulse each, 3.7 per turn | **destroyed**: rms 2013 µs, ratio 1941, total −40 µs (−28σ), red. χ² 1625, **with confident error bars** (median σ 0.96 µs) |
| Impulses 50/s, 200 ns, +20 dB (4 bursts) | vertical lines before dedispersion, reversed 6.3 ms sweeps after; ~0.3 % of baseline per µs (off-pulse residual mean +0.010) | coherent dedispersion spreads broadband impulses: natural suppression | unaffected: SNR 93.2, σ 2.866 µs, ratio 0.919 |
| All five | GNSS L1 peak, L2 ~20 MHz bump in the mean spectrum; baseline 1.795 / 1.808 (pred. 1.80 = 1.233·(1 + 1.026·(0.316 + 0.1 + 0.032))) | as radar, plus SNR loss from carrier/GNSS | destroyed: ratio 1674, total −40 µs (−20σ), red. χ² 939 |

Conclusions: all five RFI types verified (after the bpsk fix, §8 #10). Constant-envelope
RFI (carrier, GNSS) lowers SNR but keeps TOAs valid; impulses are harmless after
coherent dedispersion; the radar is catastrophic and only the reduced χ² flags it →
red. χ² flag in `estimateTOA` (cheap, phase E) and RFI excision before dedispersion
(Block 2; the radar is trivial to find there). Check figures with strong RFI: the
dynamic-spectrum colour scale is set by the carrier and the residual panel is pushed
off-axis by the baseline offset (cosmetic; percentile colour limits would help).
Check-function centroid offset again positive (+6.3 µs carrier run, +3.1 µs impulse run;
same noise seed as the +3.7 µs item).

**Check-function centroid offset: 5 noise seeds (`main.m`, −5 dB, pulsar seed 42, noise
seeds 43–47, 9 pulses each; 6 Oct 2026):**

| noise seed | 43 | 44 | 45 | 46 | 47 | pooled |
|---|---|---|---|---|---|---|
| mean per-pulse centroid offset [µs] | +3.72 | +0.45 | −3.46 | −1.19 | −2.43 | **−0.58 ± 0.93** (45 pulses) |
| total-fold centroid [µs] | +3.68 | +0.37 | −3.54 | −1.07 | −2.42 | −0.60 ± 0.92 |
| FFTFIT total fold [µs] | −0.21 | +0.63 | +0.28 | −1.03 | +0.82 | +0.10 ± 0.43 |

No bias: the +3.7 µs was a +1.8σ fluctuation of seed 43, repeated by the carrier and impulse
runs (same noise). Per-seed means scatter as expected (std 2.8 vs 2.07, 5 values). Centroid
and FFTFIT disagree per seed (correlation −0.25): different weighting; the centroid is
dominated by the window wings (lever arm t − c) and is ~2× noisier. Centroid rms over 45
pulses 5.24 µs vs predicted 6.21 (ratio 0.84 ± 0.11, 1.5σ; noise-free run went the other
way, 0.301 vs 0.277) → not significant, not pursued.

**Phase E sweep (`runSNRSweep.m` with detection, 6 Oct 2026; same seeds as phase D, TOA
table identical; P_FA 1e-3, η 3.09 / 4.15; 90 sub-ints per point):**

| snrDB | SNR opt | T0 mean | T0 std | P_D known (theory) | P_D unknown (theory) | χ² flags | good | \|z\|>5 good | rms/pred good |
|---|---|---|---|---|---|---|---|---|---|
| H0 (no pulsar) | 0 | −0.011 | 1.102 | 0/90 false alarms (0.09) | 0/90 (0.09) | – | – | – | – |
| −25 | 1.21 | 1.15 | 1.11 | 0.022 (0.030) | 0.011 (0.002) | 0 | 1 | 0 | 1.16 (n = 1) |
| −20 | 3.79 | 3.66 | 1.11 | 0.711 (0.759) | 0.433 (0.358) | 0 | 39 | **0 %** (all: 5.6 %) | 1.12 ± 0.11 |
| −15 | 11.79 | 11.60 | 1.13 | 1 (1) | 1 (1) | 0 | 90 | 0 | 1.030 |
| −10…+20 | ≥ 35 | ≥ SNR | 1.16–4.8 | 1 | 1 | 0 | 90 | 0 | = all |

Conclusions:
1. False alarms 0/90 under a true H0; T0 mean = matched-filter SNR in the weak regime
   (3.66 vs 3.79, 11.60 vs 11.79); P_D on the theory curves (−20 dB: known −1.1σ, unknown
   +1.5σ; the unknown-phase theory is approximate and slightly low).
2. **Good TOAs are outlier-free**: at −20 dB 39 of 90 kept, 0 % outliers (rms 995 → 73 µs,
   |z| < 2 94.9 %); at ≥ −15 dB nothing rejected; no χ² false flags without RFI.
3. T0 > SNR at high SNR (self-noise ignored by the H0 normalization; T0 std grows to ~4.8)
   — expected, irrelevant for the small-antenna target.
4. **Watch item:** T0 std 1.102 ± 0.075 under H0 (1.4σ), 1.11–1.13 at −25…−15 dB, and the
   phase D error-bar ratios 1.03–1.06 — all from the same 10 noise realizations, so ≈ one
   hint that the noise model underestimates the σ of template-weighted sums over many bins
   by ~5–10 % (per-bin variance is right: noise ratio 0.992, residual std 0.9995 → it would
   have to be correlation beyond neighbouring bins; no mechanism identified). If real, the
   effective threshold is ~3.8σ → ~4–5× the nominal P_FA. Test: long H0 run (no pulsar,
   ~900 sub-ints → T0 std ± 0.024, ~0.9 false alarms expected vs ~4 if real).
   **→ Resolved by the long H0 run (below): a fluctuation of the shared noise realizations.**

**Long H0 run (`runH0.m`, 6 Oct 2026; 100 passes, 900 noise-only profiles):**

| Check | Measured | Expected |
|---|---|---|
| T0 std / mean | **1.0076** / +0.034 | 1 ± 0.024 / 0 ± 0.033 |
| noise ratio | median 0.996 | 1 |
| normProfile correlation lag 0 / lag 1 | 0.999 / **0.2484** | 1 / model 0.2501 (weightX) |
| lags 2–200 | mean −0.0008, max \|r\| 0.0026 | 0 (noise 0.0007; mean-subtraction gives −(1 + 2ρ₁)/N ≈ −0.0007) |
| T0 std implied by the measured correlation | 0.997 | measured 1.008 |
| exceedances P_FA 1e-1 / 1e-2 / 1e-3 (expected 90 / 9 / 0.9) | T0: 96 / 7 / 3; Tmax: 88 / 13 / 1 | (σ 10 % low would give ~25 at 1e-2) |

Conclusions: the noise covariance model is complete for noise-dominated data (per-bin
variance, neighbour covariance, nothing beyond) → the detector thresholds give the stated
P_FA (tested down to 1e-2…1e-3) and the TOA error bars are right in the weak-signal regime.
The phase E T0 std 1.10 and the phase D ratios 1.03–1.06 were most likely the same
shared-noise fluctuation (any remaining few % at high SNR would be in the self-noise part,
irrelevant for the small-antenna target). Rice follows the measured Tmax tail for η ≳ 2.5
(below, the approximation exceeds 1; thresholds are at 3.5–4.2).

**TOA threshold sweep (`runSNRSweep.m`, 6 Oct 2026; nReal 20 → 180 sub-ints per point,
−20…−15 dB in 1 dB steps; phase D/E results kept in `data/mc/snrSweep_phaseDE.mat`):**

| SNR per sub-int (opt) | 3.8 | 4.8 | 6.0 | 7.5 | 9.4 | 11.8 |
|---|---|---|---|---|---|---|
| outliers \|z\| > 5, all TOAs | 5.0 % | 1.7 % | 0 (< 1.7 %, 95 %) | 0 | 0 | 0 |
| good TOAs (yield = P_D unknown phase) | 72 (40 %) | 129 (72 %) | 166 (92 %) | 180 (100 %) | 100 % | 100 % |
| outliers among good TOAs | 0 | 0 | 0 | 0 | 0 | 0 |
| rms/pred, good TOAs | 1.20 ± 0.08 | 0.97 | 1.01 | 1.01 | 1.02 | 1.02 |
| P_D known (theory) | 0.711 (0.759) | 0.883 (0.953) | 0.994 (0.998) | 1 | 1 | 1 |
| T0 mean / SNR opt | 0.950 | 0.956 | 0.962 | 0.965 | 0.972 | 0.978 |

H0: 0/180 false alarms, T0 std 1.020. Outlier rate close to the rough prediction
(√λ₂/2π)/√2·e^(−SNR²/4) (~1 % at 4.8, ~0.05 % at 6.0; overestimates ~2× at 3.8).

Conclusions:
1. **Design rule: SNR per sub-int ≥ 6–7** (was ≥ 10): 92 % usable TOAs at 6, 100 % at 7.5;
   the detector removes all outliers at any SNR; error bars of good TOAs honest from SNR ≈ 5
   (slight excess 1.20 at 3.8, where only 40 % survive and selection favours upward noise).
   10 → ~6.5 allows (10/6.5)² ≈ 2.4× shorter sub-ints.
2. **Baseline-estimation loss:** T0 ≈ 0.96 × SNRopt. SNRopt assumes a known noise level;
   the observer estimates the baseline from the same profile (zero-mean template), which
   discards the template's DC part: loss = √(1 − (Σs)²/(N·Σs²)) = 0.962 for the 5 %
   Gaussian (grows with pulse width; depends on the real profile). With it the P_D theory
   agrees (−19 dB known 0.926 vs 0.883, −1.8σ; −18 dB unknown 0.94 vs 0.922, −0.9σ). The
   sweep's theory curves do not include it yet (optional fix). Real physics: any observer
   without a known baseline pays it.
3. pred/opt rises at low SNR (1.03 → 1.19): the predicted σ uses the fitted (noisy)
   amplitude, E[1/b̂²] > 1/b²; rms/pred stays 1.01–1.02, so the error bars remain honest.

**Turns per sub-int (design).** Weak signal: SNR_sub = SNR_turn·√N →
N = (SNR_target / SNR_turn)², SNR_turn ≈ ρ·√(B·σ_t·√π)·0.962 ≈ 373·ρ for this pulsar and
band. Time per TOA ∝ (target/ρ)² (−10 dB in ρ → 100× longer; target 10 → 7 saves 2×):

| ρ | SNR per turn | target 6 | target 7 | target 10 |
|---|---|---|---|---|
| −20 dB | 3.7 | 3 turns (30 ms) | 4 turns (40 ms) | 8 turns (80 ms) |
| −30 dB | 0.37 | 260 (2.6 s) | 350 (3.5 s) | 720 (7.2 s) |
| −40 dB | 0.037 | 2.6e4 (4.3 min) | 3.5e4 (5.9 min) | 7.2e4 (12 min) |
| −50 dB | 0.0037 | 2.6e6 (7.2 h) | 3.5e6 (9.8 h) | 7.2e6 (20 h) |

Total precision of an observation does not depend on the sub-int length (more, less
precise TOAs); shorter sub-ints give more frequent TOAs (motion tracking, outlier
robustness). Hour-long sub-ints need Doppler inside the phase model.

**NumPy validations (algorithm ports):** dispersion kernel group delays exact to 1e-4 µs,
leakage ~1e-12, overlap-add = direct convolution (~1e-7); forward → IQ → inverse round
trip ~2e-6; IQ tones ~1.5e-6, burst centroid shift 0.04 ns, power ratio 2.000;
detection channel mapping + Parseval exact; fold indexing and nearest/linear bias;
FFTFIT MC: bias < 0.1 ps (noise-free), σ predicted 0.469 vs empirical 0.46–0.47 µs
(0.383 without covariance); red. χ² excess reproduced (1.065 median, tail-driven);
mixed signal+noise variance formula within 1 %.

**Runtimes (L = 1 s, M-series Mac, 24 GB):** dispersion 81.5 s, IQ 46 s, dedispersion
7.4 s, detection 1.2 s, fold 0.05 s; generator dominated by random-number generation.

**Experiment 3a — detected power after blanking + per-channel dedispersion
(`tests/expBlankingVariance.m`, run by Claude and Jasper 6 Oct 2026, identical results).** Noise only,
bottom channel (1.2015625 GHz, 4.1667 MHz), 4-sample bins (0.96 µs), DM 5 (filter span
715 µs incl. guards) and DM 100 (2135 µs); masks: radar 4 µs at 373 Hz (0.2 %), random
50 µs blanks (10 %), 1 ms gaps every 5 ms (20 %). Valid fraction per output sample
w = (|h|² ∗ keep)/Σ|h|²; W = bin mean.
- Mean: P/W = 1 in every W group (E|y|² = w·m exactly).
- Variance: **exact** for Gaussian input from the sample covariance
  C_{a,a+τ} = Σ_l h(l)h*(l+τ)·keep(a−l) — a convolution of the mask with
  q_τ = h·h*(·+τ), one FFT per lag τ = 0…2n−1 → cheap for every bin (matches the direct
  2n×2n matrix to ≤ 7e-8 of var0). Measured / exact over ~1e6 bins: 0.9998, 1.0001,
  1.0002 (DM 5), 1.0018, 1.0027, 1.0019 (DM 100), ± 0.0019. Unblanked bin variance 0.302
  of m² = ν₀/(B·dt) as in §5.16.
- **Neither simple model holds**: at W 0.25–0.5 exact 0.078–0.088 vs w² model 0.041–0.044
  and w model 0.110–0.113 (DM 5 random/gaps, DM 100 gaps); DM 100 random 50 µs (blanks
  ≪ sweep): w² model nearly right (0.226 vs 0.228 at W 0.75–0.9). Physics: in a dispersed
  channel lag ↔ frequency (chirp), so a blank removes a frequency slice from the affected
  output samples (fewer independent frequencies → towards w), while with B·dt ≈ 3 the
  samples in a bin are partly correlated (towards w²); short blanks spread thinly only
  scale the amplitude (w²). Errors of the simple models up to ±45 %.
- Next-bin covariance at blanked edges stays ≈ 0.05·var0 (exact 0.050–0.057), not the
  w²-scaled 0.016–0.020.
→ Design consequence (unit 3): no g(w) model; carry exact per-time-bin variance and
lag covariances (from the mask convolutions) where a mask exists, constants from the
channel spectrum where not.

---

## 8. Bugs found and fixed (lessons)

1. **Old dispersion kernel had the wrong sign**: τ(f) was put into the phase
   (exp(−i2πf·τ(f))), giving group delay −τ(f) (low frequencies arrived *earlier*),
   plus wrap-around inside the kernel. Hidden because the old inverse conjugated the
   same kernel (round trip still worked). Fix: standard chirp with correct group delay.
2. **fftshift offset**: the old causal mode added a bulk delay of Nfft/2 samples
   (8.4 ms) depending on kernel length; the "centered" mode was off by one sample.
3. **Hann window on the chirp** attenuated the low end of the band (time ↔ frequency
   mapping in a chirp). Fix: band-limit in frequency with raised-cosine edges, no time
   window.
4. **Out-of-band content** received the dispersive phase and wrapped in time. Fix:
   H = 0 outside the band.
5. **IQ low-pass group delay** (256 ns) was not removed. Fix: odd-length filter, delay
   removed exactly.
6. **IQ amplitude** halved (no factor 2). Fix: complex-envelope convention.
7. **Oversized blocks** (1e8–1e9 samples, ~10–30 GB). Fix: automatic block sizing.
8. **Missing covariance** of neighbouring phase bins (linear assignment) in fold
   uncertainties. Fix: `weightX`.
9. Main: `generate_new_data = false` crashed (info structs missing) → `loadInfo`;
   `mode` shadowed MATLAB's `mode()` → `detMode`; `run` shadowed `run()` → `runStage`;
   `fs = 500e9` typo → `500e6`.
10. **GNSS (`'bpsk'`) RFI crashed on first use** (5 Oct): `vals(sub2ind(...))` returns a
    column when `vals` has one column (one 65,536-chip chunk per 1 ms block, the normal
    case: chunks last 64 ms for L1, 6.4 ms for L2), because indexing a vector keeps the
    vector's orientation; `chips .* cos(...)` then expanded to 4e6 × 4e6. Fix:
    `reshape(..., size(nAbs))`. Lesson: code marked [written] is untested; NumPy ports
    cannot catch MATLAB indexing-shape bugs; test each branch once.
11. **Test scripts leaked state through the base workspace** (6 Oct, Jasper's runs of
    the unit tests one after another): `testDedisperseChannels` / `testDetectChannels`
    did `oldDir = cd(root); restoreDir = onCleanup(@() cd(oldDir))` as *scripts*; the
    next test's assignment to `restoreDir` destroyed the previous object, whose cleanup
    then ran at once and switched back to `tests/` → "Unrecognized function or variable
    'pipelineParams'". Claude's runs used a fresh MATLAB per test and never saw it. Fix:
    all tests are functions (own workspace; cleanup when the test ends; `run(...)`
    unchanged). Lessons: run tests the way Jasper does (several in one session); never
    run two test sessions at once — they share the files in `data/chan/` (an
    overlapping run gave garbage: noise 72× the radiometer value).

---

## 9. Known limitations and modelling gaps

- **Processing assumptions not checked**: DM correctness (a wrong DM leaves a residual
  sweep), fRef convention vs TOA definition (constant offset), receiver sideband /
  spectral inversion (assumed f_RF = fLO + f_bb), true fs (trusted via info chaining).
- **Propagation**: no scattering (irreversible pulse broadening), no DM variations, no
  higher-order frequency terms, no scintillation.
  - DM variations: δDM = 1e-3 pc cm⁻³ shifts TOAs by ≈ 2.1 µs at 1.4 GHz (chromatic;
    needs sub-band TOAs, Block 2).
  - Scintillation/scattering (multipath; Δν_d ≈ 1/(2π·τ_d)). Empirical τ_d (Bhat et al.
    2004, ±1 dex): ≈ 0.4 ns at DM 5, 1.4 GHz → Δν_d ~ 0.1–1 GHz → only a few scintles in
    400 MHz → ~100 % band-averaged flux fading on minutes–hours timescales. Effect on the
    **SNR budget** (fading margin, missed TOAs), not on timing (ns). At DM ≳ 100: τ_d of
    a few µs → time-variable profile distortion = µs-level timing noise; Δν_d ~ 50 kHz,
    little fading. With a DM error or profile evolution, scintillation's random frequency
    weighting turns the chromatic error into TOA jitter.
  - Modelling choice: within a ≤ 1 s voltage simulation diffractive scintillation is
    static → a random frequency gain g(f) (correlation bandwidth Δν_d) multiplied into
    the forward filter, plus optionally an exponential multipath tail (not implemented).
    Long-term fading statistics, Doppler and DM(t) belong at the TOA/SNR level of a later
    mission simulation. A scintillated (or power-law) pulsar spectrum S(f) has signal
    bandwidth (∫S)²/∫S² < the receiver-noise B_n; optimal detection then weights channels
    by pulsar S/N (`detectPower` `NChan` > 1).
- **Realistic SNR regime**: −5 dB is a bright pulsar on a large dish. Small antenna
  (3 m, T_sys ≈ 30 K, SEFD ≈ 2e4 Jy): ρ ≈ −27 dB (Vela, peak ~40 Jy), ≈ −35 dB (J0437,
  peak ~6 Jy), ≲ −50 dB (typical MSP) → per-pulse SNR ≪ 1, TOAs need minutes–hours of
  folding, where motion, DM(t) and scintillation matter. Phase D should sweep below −30 dB.
- **Spectrum**: generator is white; real pulsars S ∝ f^α, mean α ≈ −1.6 (factor ≈ 1.6
  across 1.2–1.6 GHz) → pulsar B_n below the noise B_n; combined with profile evolution
  and a DM error this biases wideband TOAs.
- **Receiver**: ideal linear filters, infinite dynamic range; no amplifier compression,
  ADC quantization/clipping (2-bit costs ≈ 12 % SNR), intermodulation, bandpass ripple,
  gain drift; single polarisation (dual pol: √2 in SNR).
- **Timing**: no clock jitter/drift, no relative motion (Doppler), single pulsar only.
  Earth's orbital v/c ≈ 1e-4 → ≈ 100 µs phase drift per second at f0 = 100 Hz: the
  navigation signal itself; needs Doppler terms in the fold phase model (Block 3).
- **Pulse**: stable Gaussian profile, no pulse-to-pulse jitter (a noise floor for bright
  MSPs, not in the low-SNR regime), no profile evolution with frequency, no polarisation. Gaussian-specific places to revisit for real data or other
  shapes: `ephem.profileFWHM` (→ template per pulsar), the generator and
  `expectedPowerModel` profile, the 1 %-of-peak on-pulse threshold in `estimateTOA`.
- **Estimator**: FFTFIT not optimal under self-noise (1.5× at +10 dB, 2.4× at +20 dB;
  §7 phase D); weighted (radiometer-model) fit not implemented. Below SNR ≈ 10 per
  sub-int TOAs have outliers; the good-TOA definition (detected + χ² ok, §5.13) removes
  them (§7 phase E). Detection P_FA validated down to 1e-2…1e-3 (900 H0 profiles, §7);
  smaller P_FA (e.g. 1e-6 for navigation) only by extrapolating Gaussian / Rice tails.
- **Edges**: last ~6.3 ms of each file not fully supported; partial first/last turns
  skipped (`MinCoverage` = 1).
- **Scale**: single-FFT coherent dedispersion becomes memory-limited for large DM
  (e.g. DM ≈ 71 → ~90 ms sweep); channelized (filterbank) dedispersion needed.
  Target DMs go up to ~100 → see §10 item 5 (memory-efficient FFT / alternatives).
- **Noise model in narrow channels (found 6 Oct 2026, §5.16)**: the radiometer model
  var = m²/(B·dt) per time bin with independent time bins assumes B·dt ≫ 1. In a
  3.125 MHz channel the detected power has a correlation time ~1/B_ch ≈ 0.33 µs: with
  0.96 µs bins the time-bin variance is 0.888 of m²/(B·dt) and neighbouring bins
  correlate (+0.052), measured 0.884 / 0.054. The two nearly cancel within a phase bin
  (predicted phase-bin variance 0.987 of the model at 10 ms / 2048 bins; more for
  shorter phase bins, e.g. J0437 at 2.8 µs). Full band: B·dt ≈ 400, no effect. Effect:
  error bars and P_FA slightly conservative. Fix with the weight plumbing (§10 item 1):
  include the time-bin autocovariance (known from the channel spectrum) in the fold
  noise model.

---

## 10. Open items / immediate next steps

*(Tidied 6 Oct 2026: open items first, then deferred, then a short done list with
pointers; details of finished work are in §7 and the session logs.)*

**Agreed order (6 Oct 2026, revised the same day: excision before the fast simulator).**
Why (Jasper): real data contains RFI, excision is mandatory, simulations without RFI are
not realistic. Physics: (a) at −54 dB one radar pulse ≈ 80,000 pulsar pulses (5 Oct: one
radar pulse ≈ one pulsar pulse at −5 dB) → even 99.9 % suppression leaves more than the
pulsar → excision must **blank** (weight 0), not subtract; (b) excision acts on voltages
before dedispersion (radar = 2 µs blip there; smeared over the 6–126 ms sweep after) →
the fast simulator (profile level) can only reproduce its *effects*, which must come from
a real excision stage; (c) blanking needs a data weight per time bin and channel → the
fold data product changes (today `fold.weight` is assignment-only and the same for all
channels) → fix the interface before writing a simulator that emulates it.
Decisions (Jasper): excision first; **channelized front end** (current full-band path
stays as validated reference); generic L-band RFI scenario until site measurements exist.
1. **Front-end architecture + data product** (next; design first, notes, for Jasper's
   OK). Channelized coherent dedispersion as separate functions: channelizer (polyphase
   filterbank of the IQ stream, ~3 MHz channels) → RFI detection and blanking per channel
   (spectral kurtosis + robust power threshold; whole-channel flags for persistent
   narrowband) → per-channel coherent dedispersion (short kernels: intra-channel sweep ~ms
   at DM 100 → solves the processing side of item 5) → inter-channel delay alignment →
   `detectPower` / `foldProfile` carry a data-weight stream w(t, chan); `fold.weight`,
   `weight2`, `weightX` become [NBin × nSub × nChan]; `detectPulsar` / `estimateTOA`
   combine channels with these weights; nChan = 1 must reproduce today's results exactly.
   Also enables sub-band TOAs (DM check, Block 2).
   Progress (tests pass in Claude's runs, Jasper's runs pending): unit 1 `channelizeIQ`
   (§5.15, 128 × 3.125 MHz, oversampling 4/3); unit 2 `dedisperseChannels` + option
   `AllowRefOutsideBand` in `applyInverseDispersion` (§5.16, §5.5; decisions: common
   reference via the option, taper inside each channel).
   **Unit 3 design (agreed 6 Oct): power, data weights and an exact noise model through
   the fold.** Physics: (1) blanked samples carry no signal or noise; after dedispersion
   each output sample has a valid fraction w = (|h|² ∗ mask); detected power has mean
   w·m, so the fold keeps Σa·P and Σa·w separately (unbiased even for phase-correlated
   blanking); variance σ²·g(w) with g between w² and w → decided by experiment (3a).
   (2) Narrow-channel time-bin covariance ν_l (ν₀ 0.888, ν₁ 0.046, Σ = 1): the fold
   accumulates time-bin pairs at lags l weighted by ν_l into weight2/weightX (exact;
   ν = [1] → unchanged). (3) Channels combined as the sum of channel profiles; variance =
   sum of per-channel variances (own baseline a_c, B_c); a channel missing in some phase
   bins of a sub-int is excluded from that sub-int (no dips). Decisions (Jasper):
   weights per phase bin × sub-int × channel (exact; sub-band grouping option later for
   long runs); equal-weight channel sum (SNR weighting later, with scintillation);
   exclusion rule; exact lag terms; per-channel b_c by a linear fit with τ fixed.
   Sub-units: 3a experiment `tests/expBlankingVariance.m` (g(w), ν at blanked edges);
   3b `detectChannels` + helper `powerCovariance` (power file [nChan × nBins] + weight
   file, ν from the channel spectrum); 3c `foldProfile` (weight file, per-channel
   weights, `NoiseCoeffs`); 3d `detectPulsar` + `estimateTOA` (per-channel variance sum,
   a_c/b_c, exclusion); 3e `tests/testChannelWeights.m` (synthetic masks: unbiased TOAs,
   honest error bars). Full-band path bit-identical throughout. Then a switch in main.
   **Revised after 3a (6 Oct, agreed):** no g(w) model (experiment §7: errors up to ±45 %);
   each time bin carries W (valid fraction), V (variance) and X(L) (lag covariances),
   relative to m²/(B·dt) — constants from the channel spectrum without blanking (3b),
   per-bin streams from the mask convolutions with blanking (computed in step 2, where
   the kernel h and the mask are known; the fold only has to accept them in unit 3).
   Progress: 3a done (§7), 3b done (§5.17, Lmax 7 at 0.96 µs bins). Next: 3c.
2. **RFI excision on short voltage data**, one function at a time. Characterize per RFI
   type: flagged fraction, residual noise ratio, TOA ratio and χ², false-flag rate on
   clean noise (the pulsar at weak SNR must never trigger flags); compare with the 5 Oct
   no-excision table (§7). Include a harmonic-PRF radar (400 Hz, folds coherently) and an
   extended generic L-band scenario: ATC radar with antenna rotation (new gating option in
   `rfiSource`), GNSS L1/L2/E6, Inmarsat, LTE (frequencies from memory → verify).
3. **Fast simulator for long observations** (`simulateFold`). Hours cannot be simulated
   at voltage level (~61 GB per second of data). Draws the fold struct directly (same
   fields as `foldProfile`, so `detectPulsar` / `estimateTOA` run unchanged), plus a
   separate `truth` struct. Design draft (6 Oct):
   - Units: SEFD units, baseline 1, pulse S(φ)/SEFD; SEFD = 2k(T_rx + T_sky)/A_eff.
   - Mean of bin j: 1 + (S_mean/SEFD)·⟨f(φ − Δφ)⟩, f = true shape with mean 1 over a turn
     (S_mean = catalogue flux, no W_eq rule); ⟨⟩ over the linear-assignment triangle
     kernel and the drift of Δφ(t) = φ_true(t) − φ_ephem(t) within the sub-int (hook for
     step 4).
   - Noise: per time bin σ² = m²/(n_pol·B·binDt) (dual pol = factor n_pol; observer
     passes Bnoise = n_pol·noiseBandwidth); per phase bin var σ²·W2/W², neighbour
     covariance σ_jσ_{j+1}·WX/(W_j W_{j+1}); drawn with the Cholesky factor of the circular
     tridiagonal correlation matrix; optional exact taper coefficients [css csn cnn].
   - Weights: 'asymptotic' W = n, W2 = 2n/3, WX = n/6 (lag-1 correlation 0.25, as `runH0`
     measured), n = sub-int time/(NBin·binDt), binDt bookkeeping only (default 1 µs); or
     copied from a voltage-chain fold (validation).
   - Time: `obs.passes` [N × 2] start/end, `obs.subintTime` ~60 s (whole turns);
     turnFirst/Last/Ref, tRef, fRef from the observer's ephemeris as in `foldProfile`.
   - Inputs: psr (shape, S_mean, true phase model or φ_true(t), optional FluxScale(t) for
     scintillation), rx (T_rx, T_sky, A_eff, B, n_pol), obs, ephem.
   - Validation: reproduce §7 numbers with weights copied from a voltage fold (H0 T0 std
     1.008, lag-1 0.2484; −20 dB P_D 0.711 / 0.433, TOA rms/pred 1.12; −5 dB SNR 93.4,
     σ 2.880 µs).
   - Added by the reorder: RFI/excision-effects model calibrated from step 2 (per-channel /
     time weight loss, excess noise, corrupted-sub-int rate, topocentric periodic RFI that
     smears under the barycentric fold); writes the per-channel fold struct of step 1.
   The adaptive sub-int length (item 6) fits here.
4. **Relative motion / barycentric phase prediction** (Block 3) — a **hard requirement**
   since §11: a TOA is folded across passes, so the phase must stay coherent over days.
   Doppler in the generator and in the fold phase model (~100 µs/s drift at Earth's
   orbital speed); Earth rotation; station position as input.
In parallel / later: generator-side DM memory (item 5); rest of Block 2 (DM check from
sub-band TOAs, clock jitter); Block 3 further (multiple pulsars from the §11 table,
barycentric corrections, navigation solution).

**Open items**
5. **[todo] DMs up to ~100 (Jasper, 5 Oct: 60; raised to 100 on 6 Oct for Vela, DM 67.8,
   and margin): FFT memory.** The single-FFT dispersion
   kernels grow with the sweep (∝ DM). Minimum FFT sizes (sweep 1.2–1.6 GHz + guards):

   | DM | sweep | forward dispersion (real, 4 GHz, 44 B/sample) | dedispersion (complex, 800 MHz, 56 B/sample) |
   |---|---|---|---|
   | 5 (now) | 6.3 ms | 2^26 → 3 GB (run uses 2^28, 11.8 GB) | 2^24 → 0.9 GB |
   | 30 | 37.8 ms | 2^29 → 24 GB | 2^26 → 3.8 GB |
   | 60 | 75.6 ms | 2^30 → **47 GB** (kernel design grid alone ~17 GB) | 2^27 → 7.5 GB (above the 4 GB default `MaxMemoryGB`) |
   | 100 | 126 ms | ~2^31 → **~94 GB** (extrapolated ×2) | ~2^28 → ~15 GB (extrapolated) |

   → the forward dispersion in `applyDispersionStream` is the blocker on the 24 GB
   machine; dedispersion is feasible with a raised `MaxMemoryGB`. Also: L must be well
   above the sweep (L = 0.1 s at DM 60 leaves only ~24 ms fully supported); Jasper: switch
   to 1–2 s test data when needed (disk ≈ 61 GB per second of data at fs = 800 MHz).
   Option (Jasper): a **memory-efficient FFT** (e.g. out-of-core / four-step FFT that
   works on disk-backed chunks, real-to-complex and single-precision transforms).
   Alternatives to weigh: forward dispersion at complex baseband (5× fewer samples than
   at 4 GHz RF); channelized (sub-band) dispersion and dedispersion with short per-channel
   kernels (standard in pulsar software). **Processing side decided (6 Oct): channelized
   front end** (agreed order item 1; per-channel kernels ~1e4 samples instead of ~1e8).
   **Generator side still open** (forward dispersion ~94 GB at DM 100): memory-efficient
   FFT, dispersion at complex baseband, or per-channel generation; only a simulation
   problem.
6. **[todo] Adaptive sub-int length (6 Oct).** `subintPeriods` is fixed by hand in
   `pipelineParams` (foldProfile also accepts `'SubintTime'`). Observer-only choice from the
   data: fold a first chunk or the whole observation, measure its SNR (`toa.total.snr` or
   `detection.total.T0`), SNR_turn ≈ SNR_total/√N_total, N = (SNR_target/SNR_turn)² with
   `SNR_target` (~7) as the parameter instead of `subintPeriods`. Adapts to scintillation
   (flux varies ~100 % over minutes–hours, §9). Helper e.g.
   `chooseSubintPeriods(...)`. Belongs with the fast simulator (item 3: long
   observations, varying SNR). Optional: include the baseline loss (0.962) in the
   `runSNRSweep` P_D theory curves.
7. **[todo] Factor the receiver/processing chain into one function.** The same stage
   calls appear in `main`, `runMonteCarlo`, `runSNRSweep` and `runH0` (noted 6 Oct).
8. **[todo] Second polarization in the voltage chain** (hardware is dual pol, §11): two
   independent noise streams per element, |X|² + |Y|² after detection. The fast simulator
   (item 3) covers it statistically first.
9. **[todo] Legacy `envelopeReconstruction` and `plotEnvelope`** (now in `old/`): retire
   or update to `binTime0` and `'ieee-le'`.

**Deferred / low priority**
- **Phase F**: weighted radiometer-model fit in `estimateTOA` — only gains at high SNR
  (§7), not for the weak-signal target.
- **Receiver noise does not pass W_f** in the simulation (real noise passes the same
  analog bandpass as the sky signal); < 1 % in SNR; fix by making the dispersion band
  slightly wider than the analysis band if needed (`2026-10-05_fidelity.md`).
- **Real data**: the receiver bandpass is not ideal (ripple, slopes) and not known
  exactly → measure W² from the off-pulse spectrum (or off-pulse var/mean²) and multiply
  it into W for the statistical bandwidth (∫W²)²/∫W⁴ (§6).
- Optional percentile colour limits in the RFI check figures.

**Done (details in §7 and the logs)**
- Phase C refactor (`pipelineParams`, helpers in `functions/`): `main.m` runs; bit-identity
  vs before not checkable (no earlier numbers kept), the first noisy run after it is the
  reference (`2026-10-01_session.md`).
- −5 dB / L = 1 s run (phase B): confirmed by Jasper, numbers never recorded; superseded
  by the L = 0.1 s reference (seed 43: SNR 93.4, σ 2.880 µs) and the Monte Carlo.
- `runMonteCarlo.m` (ratio 0.950, χ²_red 0.91, 72 % within 1σ); seed 43 = main
  (built-in check, §7).
- −20 dB reference test (`subintPeriods = 10`, L = 1 s): superseded by the phase D/E SNR
  sweep and the TOA threshold sweep at −20…−15 dB (§7).
- Phases D (SNR sweep) and E (noise normalization, NP detector, `flagChi2`, good TOAs,
  long H0 run); TOA threshold (SNR per sub-int ≥ 6–7); centroid offset (fluctuation of
  seed 43). Block 1 closed 6 Oct.
- Bandpass filtering / noise bandwidth review (5 Oct, `2026-10-05_fidelity.md`): all
  effects ≤ 1.5 %; flat-band prediction fixed in `expectedPowerModel`.
- RFI types observed and validated (5 Oct, §7); red. χ² flag done in phase E; RFI
  excision moved to Block 2.
- Reference scenario (6 Oct, §11): questions answered, psrcat table, DM target 100,
  T_sky in T_sys.
- `MaxMemoryGB` default 16 GB in `applyDispersionStream`; generator progress print per
  ~10 %; `*.asv` git-ignored.
per ~10 %; `*.asv` git-ignored.

---

## 11. Roadmap

**Block 1 – noise and detection (closed 6 Oct 2026; phase F and adaptive sub-int length deferred)**
- Receiver noise + RFI [noise validated at −5 dB (MC); all RFI types validated 5 Oct]
- TOA quality flag (red. χ² ≫ 1 → invalid; RFI-locked TOAs have small error bars)
  [validated: `flagChi2`, good TOAs, 6 Oct]
- Noise normalization (baseline, variance → SNR units) [validated: `detectPulsar`]
- Sub-integration length as a design parameter [validated 6 Oct: SNR per sub-int ≥ 6–7,
  N = (SNR_target/SNR_turn)²; adaptive choice from the data: §10 item 6]
- NP detector (matched filter on folded profile, threshold from false-alarm rate)
  [validated: `detectPulsar`; P_FA nominal down to 1e-2…1e-3, `runH0.m`]
- Re-evaluate TOA estimator under noise; weighted radiometer-model fit
- Monte Carlo harness: over noise seeds [validated: `runMonteCarlo.m`]; TOA error vs SNR
  vs exact optimum, pulsar seed varied [validated: `runSNRSweep.m`, phase D]

**Block 2 – processing robustness**
- DM check: sub-band TOAs (`NChan`), fit vs 1/f² → DM correction ± error
- Convention assertions: dedispersion takes `info_IQ` (fLO, fs, mapping, fRef)
- RFI excision before dedispersion (e.g. spectral kurtosis, time-frequency blanking of the
  radar); test harmonic-PRF radar. Baseline without excision: §7 RFI test.
  **Moved ahead (6 Oct 2026)**: now §10 agreed order items 1–2, with a channelized front
  end (filterbank → blanking per channel → per-channel coherent dedispersion → data
  weights per time bin and channel through the fold), before the fast simulator.
- Clock jitter in the generator and its TOA effect

**Block 3 – navigation physics**
- Relative motion: Doppler/period change in the generator; Doppler/orbital terms in the
  fold phase model; acceleration/period search
- Multiple pulsars (≥ 4) with realistic P, DM, flux, profiles
- Timing corrections: barycentric (Roemer, Shapiro, Einstein), clock corrections,
  infinite-frequency conversion (`bulkDelayRef`)
- Residuals vs timing model; navigation solution (least squares → Kalman)

**Block 4 – later**
- 3×3 array: element signals with geometric delays; beamforming at IQ level before
  dedispersion
- Scattering, DM variations; polarisation
- High DM (target up to ~100): processing side solved by the channelized front end
  (§10 agreed order item 1); generator side (forward dispersion) still open (§10 item 5;
  needed before realistic pulsars)
- Performance of generation for long simulations; optional fold-in-detection fusion for
  millisecond pulsars

**Open design decision:** a **reference scenario** (specific millisecond pulsar(s) with
P, DM, flux, profile; receiver: antenna/array gain, T_sys, bandwidth; ground vs
spacecraft; affordable observation time). It fixes realistic noise levels, sub-int
lengths, and turns Monte Carlo curves into "X ns per pulsar after Y minutes".

**Reference scenario – draft (6 Oct 2026, from Jasper's hardware; questions answered below).**

*Receiver.* 3×3 array, **0.61 m²** (effective or physical: open), **T_sys 60–80 K**, band
**1.2–1.6 GHz** (B = 400 MHz), **ground station**, observation time "as good as possible
for navigation". SEFD = 2kT_sys/A_eff = 2·1.38e-23·70/0.61 ≈ **3.2e5 Jy** (2.7–3.6e5 for
60–80 K; Parkes ~30 Jy, 25 m dish ~1000 Jy).

*Time per TOA* (radiometer equation SNR = (S_mean/SEFD)·√(n_pol·B·t)·√((P − W_eq)/W_eq),
n_pol 2, B 400 MHz, target SNR 7 from the threshold sweep). P, DM, W50, S₁₄₀₀ from the
**ATNF catalogue v2.6.5** (queried 6 Oct 2026; all pulsars with S₁₄₀₀ > 25 mJy plus MSPs
> 4 mJy). Rules: **W_eq = 2·W50** (Gaussian matched filter gives 1.5·W50; factor 2 as
margin for wings/components, e.g. J0437 W10 = 1.02 ms); ρ = S_peak/SEFD_nom with
S_peak = S·P/W_eq; σ_TOA = √2·(W50/2.355)/7; visible = max. elevation ≥ 20° → station
latitude within dec ± 70°. Catalogue S₁₄₀₀ are averages; low-DM pulsars scintillate
(factor ~2 or more either way). Uncertainty of t: ~×2.

| Pulsar | dec | P [ms] | DM | W50 [ms] | S₁₄₀₀ [mJy] | ρ [dB] | t nom / **worst** [h] | σ_TOA [µs] (c·σ) | latitudes |
|---|---|---|---|---|---|---|---|---|---|
| **J0437−4715** (MSP) | −47 | 5.76 | 2.6 | 0.141 | 150 | −50 | 3.9 / **21** | **12** (3.6 km) | < 23° N |
| **J0835−4510** (Vela) | −45 | 89.3 | 67.8 | 1.7 | 1050 | −41 | 0.1 / **0.3** | 146 (44 km) | < 25° N |
| **J0332+5434** (B0329+54) | +55 | 714.5 | 26.8 | 6.6 | 203 | −45 | 0.8 / **4.2** | 566 (170 km) | > 15° S |
| J0953+0755 (B0950+08) | +8 | 253.1 | 3.0 | 8.6 | 100 | −53 | 12 / 67 | 738 (221 km) | 62° S–78° N |
| J1752−2806 (B1749−28) | −28 | 562.6 | 50.3 | 6.6 | 48 | −52 | 18 / 96 | 566 (170 km) | < 42° N |
| J1645−0317 (B1642−03) | −3 | 387.7 | 35.8 | 3.8 | 26 | −54 | 52 / 278 | 326 (98 km) | 73° S–67° N |
| J1932+1059 (B1929+10) | +11 | 226.5 | 3.2 | 5.6 | 29 | −57 | 106 / 571 | 480 (144 km) | 59° S–81° N |
| J1939+2134 (B1937+21, MSP) | +22 | 1.56 | 71.0 | ~0.05 | 13.9 | −62 | 670 / 3600 | 5 (1.4 km) | > 48° S |
| J1713+0747 (MSP) | +8 | 4.57 | 16.0 | 0.30 | 8.3 | −67 | 3700 / 20000 | 26 | |
| J2145−0750 (MSP) | −8 | 16.05 | 9.0 | 0.40 | 5.5 | −65 | 3000 / 16000 | 34 | |

Bright but DM far too high for now: J1644−4559 (DM 479), J0738−4042 (161), J1935+1616
(159), J0837−4135 (147). Corrections to the earlier from-memory draft: B1929+10 is
29 mJy (not ~35–200), normal-pulsar times are 1–100 h nominal (not 0.5–10 h).

*Figure of merit for navigation* = σ_TOA·√t (µs·√h, smaller is better; σ after time T
is this/√T). Worst case: J0437 55, **Vela 80**, B0329+54 1160, B1749−28 5550, B0950+08
6040, B1642−03 5440. **Vela is as good as J0437** (bright: TOA in ~20 min; wider pulse) —
but DM 67.8 (now within the DM ≤ 100 target, §10 item 5) and glitches / timing noise
(ephemeris must be recent). Both are southern (visible < ~25° N); for northern sites
**B0329+54** is the only source with TOAs within a pass (~70 km per day of observation).
*Decided (Jasper, 6 Oct 2026):*
- **High-DM target 100** (was 60): includes Vela (67.8), B1937+21 (71.0), B0740−28 (73.8)
  and margin; needs the memory-aware FFT anyway (§10 item 5: sweep 126 ms at DM 100).
- **Sky temperature in T_sys**: T_sys = T_rx (60–80 K, assumed to exclude sky) +
  T_sky(direction). The array beam is wide (~λ/D ≈ 0.21/0.8 m ≈ 15°) → T_sky is the
  beam-averaged sky: ~3–5 K off the Galactic plane (CMB 2.7 K + Galactic), ~5–20 K for
  plane sources (Vela, B1749−28), more towards the Galactic centre. Per pulsar a T_sky
  value (from a 1.4 GHz sky map, e.g. Haslam 408 MHz scaled with spectral index ~−2.6).

t ∝ S⁻² (15× fainter → 225× longer) → **J0437 is the only realistic MSP**; it resembles the
simulated pulsar (P 5.8 ms vs 10 ms, DM 2.6 vs 5) at ρ ≈ −50 dB (deep weak-signal regime).
Normal pulsars are brighter but have wider pulses (worse TOA precision) and more timing
noise (glitches, red noise). J0437 TOA precision: main-peak FWHM ~0.14 ms → σ_t ≈ 60 µs →
σ_TOA ≈ √2·60/7 ≈ **12 µs per ~4 h TOA** (c·σ ≈ 3.6 km along the line of sight; ~1.5 km
per day).

*Consequences.*
1. A 3-D fix needs ≥ 4 pulsars (position + clock); with one MSP not directly possible →
   options: (a) J0437 + bright normal pulsars, (b) known ground position (demonstrate the
   timing chain), (c) more collecting area (10× area → 100× shorter integrations).
2. The fast simulator (§10 agreed order item 3) becomes essential (hour-long sub-ints).
3. Motion (§10 agreed order item 4) becomes urgent: within a 4 h sub-int, Earth rotation
   (~0.3 km/s at mid-latitudes) and orbit (30 km/s) shift the phase far beyond the pulse
   width → full barycentric phase prediction (TEMPO2-like) from the start.
4. Scintillation: J0437 (low DM) has scintles of hundreds of MHz and ~100 % flux
   variation over minutes–hours → SNR budget risk and opportunity (integrate longer when
   bright) → adaptive sub-int length (§10 item 6).
5. DM ≤ 100 (decided below) covers J0437, Vela and all feasible candidates; the bright
   pulsars above it (DM 147–479) stay out.
6. The simulation is single-polarization; dual pol gives √2 in SNR (2× in time). The
   hardware is dual pol (answer 2 below) → add the second polarization.

*Answers (Jasper, 6 Oct 2026) and design point.*
1. **Area unknown → design for the noisy case.** Two design points:
   - *nominal*: A_eff = 0.61 m², T_sys 70 K → SEFD ≈ 3.2e5 Jy (table above);
   - *worst case (design point)*: 0.61 m² physical × aperture efficiency ~0.5 →
     A_eff ≈ 0.3 m², T_sys 80 K → SEFD = 2·1.38e-23·80/0.3 ≈ **7.4e5 Jy** (2.3× nominal).
   t ∝ SEFD² → **5.3× longer**: J0437 ≈ **20 h per TOA** (SNR 7, dual pol),
   ρ ≈ −54 dB; normal pulsars 4 h (B0329+54) to hundreds of h (table above). σ_TOA per TOA unchanged (fixed SNR 7);
   per unit time √5.3 = 2.3× worse. The pipeline must work at the worst case; nominal is
   the bonus.
2. **Dual polarization on all 9 elements** (confirmed): n_pol = 2 in the radiometer
   equation (already in the table). Simulator: two independent noise streams per element,
   beamformed per polarization, total intensity |X|² + |Y|² after detection. Fast
   detected-power simulator: n_pol = 2 in the noise variance.
3. **Site-independent ("should work anywhere").** Station position (latitude, longitude,
   height) is an input parameter, needed anyway for the barycentric correction; no pulsar
   hard-coded; source list from visibility (elevation mask) at the site. J0437 (dec −47°)
   is never visible north of ~43° N → from Europe / most of North America **no feasible
   MSP**. Pipeline generic in P (ms to ~1 s).
4. **Normal pulsars are needed** (Jasper can advise the project): (a) geometry — one
   pulsar constrains position along one direction; 3-D position + clock needs ≥ 4 well-
   separated directions (3 with a known clock); (b) visibility — item 3. Cost: wide pulses
   → σ_TOA 0.3–0.7 ms at SNR 7 (100–220 km) vs J0437 ~12 µs (~3.6 km), Vela 146 µs;
   plus timing noise and glitches.

*Consequences of the worst-case design point.*
- A 20 h TOA is longer than one J0437 pass → **fold across passes with gaps** (TOA from
  several days of data). The fold phase must stay coherent over days → barycentric phase
  prediction (§10 agreed order item 4) is a hard requirement, not a refinement.
- Sub-ints within a pass (hours) are far below SNR 6–7 → combine folded profiles across
  sub-ints/passes with the phase model before the TOA fit; detection on the combined
  profile.
- Scintillation (J0437, B0950+08: low DM) dominates the per-pass SNR → weight passes by
  measured SNR.

---

## 12. Design rationale (why things are the way they are)

- **Complex data until detection**: coherent dedispersion needs phase; beamforming and
  RFI excision are linear and must precede squaring.
- **Square-law, then fold, then template match**: the pulsar signal is noise; information
  is in the variance → energy detection + profile weighting is optimal (NP).
- **Noise after dispersion**: physical order (sky → receiver); noise and RFI are never
  dispersed; dedispersion inverse-chirps them, as with real data.
- **Band-limited dispersion kernel**: models the receiver bandpass; exact because
  dispersion is phase-only and commutes with the bandpass (see 5.2).
- **fRef = fHigh**: top of the band keeps its time; bulk τ(fRef) recorded for later
  infinite-frequency conversion.
- **Linear phase assignment** in the fold: unbiased at sub-bin level; its neighbour
  covariance is tracked explicitly.
- **FFTFIT**: standard, continuous shifts, fast, ML for white noise; uncertainties from
  an observer-usable noise model incl. covariances.
- **Sub-integrations**: validation of error bars, drift detection (the navigation signal
  itself), outlier detection.
- **Ground-truth separation**: every processing stage must work unchanged on real data.
- **Parameters in one script** (`pipelineParams.m`): main and the MC cannot drift apart.
  Later split into simulation vs observation parameters for real data.
- **MC on separate files** (`data/mc/`): a later main run with skipped stages never
  picks up an MC realization (`checkConsistency` does not check the noise seed).

---

## 13. File list

| File | Role | Status |
|---|---|---|
| `main.m` | pipeline driver: run control, stages, `checkConsistency` | current |
| `pipelineParams.m` | all parameters + `ephem` (script, shared) | current |
| `runMonteCarlo.m` | noise-seed Monte Carlo, pooled validateTOA | validated (−5 dB) |
| `runSNRSweep.m` | phase D SNR sweep (pulsar + noise varied), summary + figure | validated (−25…+20 dB) |
| `loadInfo.m` | reload a stage's `_info.mat` | moved from main |
| `gaussianTemplate.m` | periodic Gaussian TOA template | moved from main |
| `noiseBandwidth.m` | (∫W²)²/∫W⁴ of the band tapers | moved from main |
| `generatePulsarSignal.m` | pulsar signal + ground truth | validated |
| `applyDispersionStream.m` | ISM dispersion (incl. `makeDispersionKernel`, `chooseBlockSize`) | validated |
| `addNoiseAndRFI.m` | receiver noise + RFI at RF, SNR predictions | noise validated (−5 dB, MC); RFI validated (5 Oct) |
| `rfiSource.m` | RFI source definitions | validated (5 Oct) |
| `applyIQmodulation.m` | downconversion to complex baseband | validated |
| `applyInverseDispersion.m` | coherent dedispersion | validated |
| `detectPower.m` | square-law detection, optional channels | validated |
| `foldProfile.m` | folding with phase model | validated |
| `estimateTOA.m` | FFTFIT TOAs + uncertainties | validated (noise-free; −5 dB MC) |
| `validateTOA.m` | TOA vs ground truth, MC pooling | validated |
| `detectPulsar.m` | NP detection (known / unknown phase), noise normalization, noise check | validated (6 Oct; P_FA via `runH0.m`) |
| `runH0.m` | long noise-only run: T0/Tmax statistics, exceedances, bin-bin correlation | validated (6 Oct) |
| `expectedPowerModel.m` | shared ground-truth power/variance model | written |
| `channelizeIQ.m` | oversampled polyphase filterbank, IQ → per-channel IQ files (channelized front end, unit 1) | validated (6 Oct) |
| `tests/testChannelizeIQ.m` | unit tests: tones, white noise, real-data Parseval | passes (Jasper's run, 6 Oct) |
| `dedisperseChannels.m` | coherent dedispersion per channel, common reference (channelized front end, unit 2) | validated (6 Oct) |
| `tests/testDedisperseChannels.m` | regression, commutation, end-to-end TOAs, time-bin noise | passes (Jasper's run, 6 Oct) |
| `tests/expBlankingVariance.m` | experiment 3a: detected-power variance after blanking (exact via mask convolutions) | passes (Jasper's run, 6 Oct) |
| `powerCovariance.m` | exact variance / lag covariances of detected time bins from the channel spectrum (unit 3b) | validated (6 Oct) |
| `detectChannels.m` | detectPower per channel → one [nChan × nBins] power file + noise stats (unit 3b) | validated (6 Oct) |
| `tests/testDetectChannels.m` | spectrum vs filter, layout, fold compatibility, measured V / X | passes (Jasper's run, 6 Oct) |
| `plotDispersionCheck.m`, `plotIQCheck.m`, `plotDetectedPower.m`, `plotFoldCheck.m` | checks | noise-free validated; noise versions written |
| `old/envelopeReconstruction.m`, `old/plotEnvelope.m` | legacy | to retire/update |
