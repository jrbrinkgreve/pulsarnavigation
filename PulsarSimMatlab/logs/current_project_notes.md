# Pulsar navigation – synthetic data & processing pipeline: project notes

Reference document for the MATLAB pipeline that simulates a pulsar signal through a
realistic receiver chain and processes it back into times of arrival (TOAs). It records
what exists, the conventions every stage follows, the physics and formulas behind each
step, what has been validated (with numbers), what went wrong along the way, and what is
still to do.

This file tracks the **current state** of the code (started 30 September 2026 as
`2026-09-30_project_notes.md`). Last updated: 1 October 2026, evening (after the
`pipelineParams` / Monte Carlo refactor; see `2026-10-01_session.md`). Status markers: **[validated]** = run in MATLAB and
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
VALIDATION
  9. validateTOA            (val struct + figure)                         [validated]
CHECKS (ground truth, optional via plots.*)
  plotDispersionCheck, plotIQCheck, plotDetectedPower, plotFoldCheck,
  expectedPowerModel (shared model)
```
\* validated on noise-free data. The receiver stage (3) and the noise-aware check
functions have been run once (−5 dB, L = 1 s, fLO 1.3 GHz, fs 800 MHz; confirmed by
Jasper, numbers not yet recorded in §7); the −20 dB reference test and the Monte Carlo
are still to do.

**Run control** (`main.m`): `runStage.sky` (generator + dispersion, the slow part),
`runStage.receiver` (noise/RFI + IQ; rerun alone to change SNR or RFI),
`runStage.process` (dedispersion + detection), `runStage.fold`, `runStage.toa`.
Skipped stages reload their `info` from disk; `checkConsistency` warns when files on disk
were made with different parameters (T, f_in, L, DM, band, SNR, receiver input, fLO, fs).
`plots.*` switch the check figures; `closeFigures` runs `close all`.

**Helpers** (in `functions/`, shared by main and the MC script): `loadInfo(dataFile)`
(reload `<name>_info.mat`), `gaussianTemplate(nBin, fwhmTurns)` (periodic Gaussian,
peak 1 at bin 1), `noiseBandwidth(fLow, fHigh, edgeWidths)` (= (∫W²)²/∫W⁴ for a
product of raised-cosine edge tapers). **Local function in main**: `checkConsistency`.

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

### 5.3 `addNoiseAndRFI(inFile, outFile, fs, ...)` + `rfiSource(type, ...)` [written]

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

(FFTFIT in noise-dominated data should reach these within ~1–2 %; B_eff ≈ 390 MHz, not 400.)

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

**Scenario in main** (illustrative L-band; `rfiOn` switch):

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

**Options:** `RefFreq`, `EdgeFrac`, `GuardTime`, `Nfft`, `MaxMemoryGB`,
`CheckLeakage`, `MaxLeakage`, `SaveInfo`, `Verbose`.

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
(on-pulse bins, model > 1 % of peak), tRef, fRef, turnRef; `toa.total` (whole fold):
phase, phaseErr, timeOffset (= phase/f0, relative to the ephemeris, **not** an absolute
TOA), timeOffsetErr, amp, ampErr, snr, redChi2.

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
| TOA | tRef + τ/fRef |
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
| Check-function centroids | per pulse mean +3.7 µs (noise ~6.2 µs each), total +3.68 ± 2.06 µs (1.8σ) | 0; watch in the MC |

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
is validated at this SNR. Mean predicted σ 2.907 µs ≈ 2 % above the radiometer optimum.
Caveat: the seeds share the pulsar realization; its pure self-noise share is ~4 % of the
per-bin variance at the peak here, so the TOAs are nearly independent. At high SNR the
MC must also vary the pulsar seed. The check-function centroid offset (+3.7 µs) is not
tested by the MC (no plot checks there).

**NumPy validations (algorithm ports):** dispersion kernel group delays exact to 1e-4 µs,
leakage ~1e-12, overlap-add = direct convolution (~1e-7); forward → IQ → inverse round
trip ~2e-6; IQ tones ~1.5e-6, burst centroid shift 0.04 ns, power ratio 2.000;
detection channel mapping + Parseval exact; fold indexing and nearest/linear bias;
FFTFIT MC: bias < 0.1 ps (noise-free), σ predicted 0.469 vs empirical 0.46–0.47 µs
(0.383 without covariance); red. χ² excess reproduced (1.065 median, tail-driven);
mixed signal+noise variance formula within 1 %.

**Runtimes (L = 1 s, M-series Mac, 24 GB):** dispersion 81.5 s, IQ 46 s, dedispersion
7.4 s, detection 1.2 s, fold 0.05 s; generator dominated by random-number generation.

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

---

## 9. Known limitations and modelling gaps

- **Processing assumptions not checked**: DM correctness (a wrong DM leaves a residual
  sweep), fRef convention vs TOA definition (constant offset), receiver sideband /
  spectral inversion (assumed f_RF = fLO + f_bb), true fs (trusted via info chaining).
- **Propagation**: no scattering (irreversible pulse broadening), no DM variations, no
  higher-order frequency terms, no scintillation.
- **Receiver**: ideal linear filters, infinite dynamic range; no amplifier compression,
  ADC quantization/clipping, intermodulation, bandpass ripple, gain drift.
- **Timing**: no clock jitter/drift, no relative motion (Doppler), single pulsar only.
- **Pulse**: stable Gaussian profile, no pulse-to-pulse jitter, no profile evolution with
  frequency, no polarisation. Gaussian-specific places to revisit for real data or other
  shapes: `ephem.profileFWHM` (→ template per pulsar), the generator and
  `expectedPowerModel` profile, the 1 %-of-peak on-pulse threshold in `estimateTOA`.
- **Estimator**: FFTFIT not optimal under self-noise; weighted (radiometer-model) fit
  not implemented.
- **Edges**: last ~6.3 ms of each file not fully supported; partial first/last turns
  skipped (`MinCoverage` = 1).
- **Scale**: single-FFT coherent dedispersion becomes memory-limited for large DM
  (e.g. DM ≈ 71 → ~90 ms sweep); channelized (filterbank) dedispersion needed.

---

## 10. Open items / immediate next steps

1. Run `main.m` once to confirm the refactor (`pipelineParams`, helpers in
   `functions/`): TOAs identical to before.
2. Record the −5 dB / L = 1 s run in §7.
3. ~~Run `runMonteCarlo.m`~~ (done, §7: ratio 0.950, χ²_red 0.91, 72 % within 1σ).
   Confirm seed 43 matches main (median SNR 93.4, σ 2.880 µs) from the Command Window.
   Check-function centroid offset (+3.7 µs, 1.8σ): run main with another `noiseSeed`.
4. −20 dB reference test (`subintPeriods = 10`, L = 1 s): off-pulse normalized residual
   std ≈ 1.00; 10 TOAs with SNR ≈ 12 and σ ≈ 25 µs; red. χ² ≈ 1.
5. Next phases (plan of 1 Oct): **D** SNR sweep in `runMonteCarlo.m` (rms TOA error vs
   radiometer optimum, −30…+10 dB); **E** noise normalization + NP detector
   (`detectPulsar.m`); **F** weighted radiometer-model fit in `estimateTOA`.
6. **Review (Jasper): bandpass filtering and noise bandwidth.** Not yet fully understood;
   it must be handled carefully. Points to go through together:
   - Which filters each component passes. Pulsar: dispersion taper W_f (8 MHz sin²
     edges) → IQ low-pass (flat to ±~355 MHz around fLO) → dedispersion taper W_i.
     Receiver noise: white at RF over 0–f_in/2 → IQ low-pass → W_i only. The IQ check
     figure shows it: before dedispersion the noise extends beyond 1.2–1.6 GHz; W_i
     removes it.
   - Why the radiometer bandwidth is (∫W²)²/∫W⁴ and not ∫W² or the nominal B: the
     variance of an averaged |z|² depends on how correlated the spectrum is.
   - Which B is used where: `estimateTOA` 'Bnoise' = W_i only (391.6 MHz, observer
     knowledge); the pulsar actually sees W_f·W_i (389.6 MHz); `predictPerformance` in
     `addNoiseAndRFI` uses a flat B = 400 MHz (~1–2 % optimistic SNR; part of the
     93.4 vs 97.9 gap). `expectedPowerModel` handles signal/noise filters separately
     (Js, Jn, css, csn, cnn).
   - Real data: the receiver bandpass is not ideal (ripple, slopes) and not known
     exactly → measure W² from the off-pulse spectrum and multiply it into W.
7. `rfiOn = true`: observe RFI effects (baseline rise, residual excess, TOA degradation).
8. Legacy `envelopeReconstruction` and `plotEnvelope` (now in `old/`): retire or update to
   `binTime0` and `'ieee-le'`.

Done: `MaxMemoryGB` default 16 GB in `applyDispersionStream`; generator progress print
per ~10 %; `*.asv` git-ignored.

---

## 11. Roadmap

**Block 1 – noise and detection (in progress)**
- Receiver noise + RFI [noise validated at −5 dB (MC); RFI not yet run]
- Noise normalization (baseline, variance → SNR units)
- Sub-integration length as a design parameter (SNR per sub-int ≳ 10–20)
- NP detector (matched filter on folded profile, threshold from false-alarm rate)
- Re-evaluate TOA estimator under noise; weighted radiometer-model fit
- Monte Carlo harness: over noise seeds [validated: `runMonteCarlo.m`]; TOA error vs SNR
  vs radiometer prediction [todo, phase D; vary the pulsar seed too at high SNR]

**Block 2 – processing robustness**
- DM check: sub-band TOAs (`NChan`), fit vs 1/f² → DM correction ± error
- Convention assertions: dedispersion takes `info_IQ` (fLO, fs, mapping, fRef)
- RFI excision before dedispersion (e.g. spectral kurtosis); test harmonic-PRF radar
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
- Scattering, DM variations; channelized dedispersion for high DM; polarisation
- Performance of generation for long simulations; optional fold-in-detection fusion for
  millisecond pulsars

**Open design decision:** a **reference scenario** (specific millisecond pulsar(s) with
P, DM, flux, profile; receiver: antenna/array gain, T_sys, bandwidth; ground vs
spacecraft; affordable observation time). It fixes realistic noise levels, sub-int
lengths, and turns Monte Carlo curves into "X ns per pulsar after Y minutes".

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
| `loadInfo.m` | reload a stage's `_info.mat` | moved from main |
| `gaussianTemplate.m` | periodic Gaussian TOA template | moved from main |
| `noiseBandwidth.m` | (∫W²)²/∫W⁴ of the band tapers | moved from main |
| `generatePulsarSignal.m` | pulsar signal + ground truth | validated |
| `applyDispersionStream.m` | ISM dispersion (incl. `makeDispersionKernel`, `chooseBlockSize`) | validated |
| `addNoiseAndRFI.m` | receiver noise + RFI at RF, SNR predictions | noise validated (−5 dB, MC); RFI not run |
| `rfiSource.m` | RFI source definitions | written |
| `applyIQmodulation.m` | downconversion to complex baseband | validated |
| `applyInverseDispersion.m` | coherent dedispersion | validated |
| `detectPower.m` | square-law detection, optional channels | validated |
| `foldProfile.m` | folding with phase model | validated |
| `estimateTOA.m` | FFTFIT TOAs + uncertainties | validated (noise-free; −5 dB MC) |
| `validateTOA.m` | TOA vs ground truth, MC pooling | validated |
| `expectedPowerModel.m` | shared ground-truth power/variance model | written |
| `plotDispersionCheck.m`, `plotIQCheck.m`, `plotDetectedPower.m`, `plotFoldCheck.m` | checks | noise-free validated; noise versions written |
| `old/envelopeReconstruction.m`, `old/plotEnvelope.m` | legacy | to retire/update |
