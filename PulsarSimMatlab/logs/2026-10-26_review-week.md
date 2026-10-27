# 2026-10-26 – Review week (E0), day 1

The review week was agreed on 26 Oct (roadmap.html#blockE, E0): five sessions, more if needed. Per part: the physics,
the code block by block, the tests and what they miss. Refinements are noted, not built (minor fixes excepted).

Jasper, 26 Oct: "lets also note all the smaller details from the other parts in the pipeline." So every smaller
detail goes onto one **running list** at the end of this log, numbered by day (D1-1, D1-2, …).

Earlier review material:
- the first review week (12–16 Oct, stages 1–11): logs `2026-10-12…16_review-week.md`, and notes §10 item 11;
- parts 2–3 of this review (23 Oct, the fast simulator and the RFI profiles): `2026-10-23_review-C-D.md`.

## The plan, re-made (Jasper, 26 Oct)

Jasper goes through the pipeline in the order of the pipeline overview page (docs/pipeline.html, "Stages and data
products"). He has done stages 1–3 (generator, dispersion, receiver noise and RFI) himself, and is now at stage 4,
the IQ generation.

| Day | Stages (page order) |
|---|---|
| 1 | 4 `applyIQmodulation` (IQ generation); 4a `makeClock`; 5 `applyInverseDispersion`, 6 `detectPower` (the full-band front end) |
| 2 | the channel path: 5a `channelizeIQ`, `detectRFI` (+ `periodicRFI`), `blankChannels`, 5b `dedisperseChannels`, 6b `detectChannels`, `blankingWeights` |
| 3 | 7a phase prediction: `makeTrajectory`, `delayModel`, `loadEphemeris`, `makePredictor`; 7 `foldProfile` (+ `'Predictor'`); 5–7 `processIQ` (reviewed below) |
| 4 | 8–9 TOAs: `foldToTOA` (reviewed below), `estimateTOA`, `detectPulsar`, `validateTOA`; 10 `navFit` in depth, `navFilter` briefly |
| 5 | the simulator rows: `trueResiduals` (`rfiProfile`, `simulateFold`, `adaptiveSubints`/`mergeSubints` were reviewed on 23 Oct); main.m run control and the drivers (`runScenario`, `runShip`); speed (profiling of 25 Oct); a main.m smoke run; one end-to-end walk; housekeeping; the choice of refinements before E1 |

## Reviewed early (26 Oct): the processing boundary (`processIQ`, `foldToTOA`)

**Physics.** What an observer does with the receiver's IQ samples, up to a TOA:
1. **Front end.** The input is complex baseband IQ: the 1.2–1.6 GHz band mixed down by the 1.3 GHz LO, 800 MS/s.
   It becomes detected power per time bin.
   - Full band: one coherent dedispersion of the whole 400 MHz band, then |y|² in bins of f_out.
   - Channel path: a polyphase filterbank (128 channels of 3.125 MHz, oversampled 4/3), then RFI excision, then
     coherent dedispersion per channel (aligned at refFreq), then |y|² per channel.
   - Excision comes before dedispersion: a radar pulse is a few samples wide there, while dedispersion would
     spread it over the filter's reach.
   - Dedispersion does spread every blank, so `blankingWeights` computes, per detected bin, the valid fraction W
     and the noise V, X.
2. **Fold.** The phase model φ(t) = Φ0 + f0(t − TRef) + ½F1(t − TRef)². Each time bin is split over two phase
   bins (linear assignment). Sums and weights are kept per sub-int.
3. **Inputs of the TOA step.**
   - B_noise, the noise-equivalent bandwidth of the detected power: full band from the observer's dedispersion
     taper, channel path that of one channel.
   - The weighting of the channels.
4. **TOAs (`foldToTOA`).** A template match gives the phase offset Δφ and TOA = tRef + Δφ/fRef. The error bar comes
   from the radiometer noise model, not from the scatter.
   - Detection is Neyman–Pearson at P_FA 1e-3: T0 at the known phase, and Tmax over all phases (Rice factor).
   - **Good TOA** = fitted, detected at unknown phase, and χ²_red < 2 (pChi2 is reported, not used).
- **The ground-truth rule holds by the signature:** IQ file, receiver settings, ephemeris and processing settings
  only.
- `processIQ` is also the MATLAB reference ("golden model") for F's FPGA or GPU port: the same stage boundaries and
  intermediate files.

**Tests.**
- No unit test calls `processIQ` directly.
- It is exercised by `main.m` (Jasper's run 23 Oct), the statistics drivers (`runMonteCarlo`, `runH0`,
  `runSNRSweep`, `runRFITest`: validated 19 Oct), `testRFIChain` (slow; channel path with excision; Jasper's run
  22 Oct), and the experiments.
- The `'Reuse'` branch and its checks (`checkProcessing`, `checkFold`) run only when main.m has
  `runStage.process = false`. No test covers them.
- `foldToTOA` runs in every fast-simulator test and in all of the above.

## Day 1, stage 4: the IQ generation (`applyIQmodulation`)

**Physics.**
- **What it does.** The receiver's real RF signal, sampled at 4 GHz, holds 0–2 GHz; the observed band is
  1.2–1.6 GHz. Three steps follow:
  1. mixing with a complex LO, e^(−i2π f_LO t) with f_LO 1.3 GHz, shifts the band to −100…+300 MHz (no inversion);
  2. a low-pass removes everything else, including the real signal's mirror, which lands at +1.1…+1.5 GHz (mod 4 GHz);
  3. keeping every 5th sample gives 800 MS/s complex, which covers ±400 MHz.
- **Why complex samples.** They keep the phase, which coherent dedispersion needs, since it is a phase filter. A real
  SDR also delivers IQ.
- **Gain convention: the complex envelope.** x(t) = Re{z(t) e^(i2π f_LO t)}, so z = LPF{2x e^(−i2π f_LO t)} and
  mean|z|² = 2·mean(x²) in band.
- **Timing.** The FIR's group delay, G = 1024 input samples = 256 ns, is removed exactly. Output sample k is input
  sample 5k, and the LO phase is referenced to input sample 0. Without this, every TOA would be 256 ns late.
- **Numbers, checked:**
  - Kaiser β = 7.857 for 80 dB;
  - transition width 9.79 MHz (flat to ±355 MHz, stopband from 365 MHz, below the output Nyquist of 400 MHz);
  - LO step per input sample r = 13/40 exactly, so the LO phase never drifts;
  - streamed by overlap-add (Nfft 2^22).
- **It is an ideal digital downconverter:** no I/Q imbalance, no DC or LO leakage, no phase noise, no quantization.

**Code** (applyIQmodulation.m):

| Lines | Block |
|---|---|
| 91–115 | Checks; the decimation D = round(f_in/fs) = 5 |
| 117–147 | The Kaiser low-pass (odd length, so the group delay is an integer); the alias and band-edge warnings |
| 149–205 | File size; block size from the memory budget |
| 224–264 | Per block: the LO phase from the absolute index (in double), mix, FFT overlap-add, carry; `writeDecimated` removes G and keeps every 5th sample on the global grid; flush the carry |
| 271–313 | info: mapping "f_baseband = f_RF − fLO (no inversion)", t0 = 0, the filter |

**Tests.**
- There is no MATLAB unit test (R-I1, noted 13 Oct).
- Checked by `plotIQCheck` (visual), by the NumPy ports (not in the repo) and end to end: main.m's TOAs are unbiased
  to 0.3σ (±0.7 µs), which bounds any timing error of this stage well below a µs.

## Day 1, stage 4a: the receiver clock (`makeClock`)

**Physics.**
- **What the clock does to a TOA.** Every time stamp comes from the receiver's oscillator, so its error ε(t) = clock
  reading − true time enters every TOA directly. 1 µs of clock error looks like 300 m of position.
- **The model:** ε(t) = Offset + Rate·(t − t1) + ½·Drift·(t − t1)² + x(t). x is the integral of the oscillator's
  random fractional frequency y(t), with the standard power-law spectrum S_y(f) = h0 + h₋₁/f + h₋₂/f².
  - **White FM** (h0): x is a random walk, var = h0·Δ/2 per step.
  - **Flicker FM** (h₋₁): no exact finite-state model exists, so it is synthesized by FFT over twice the span.
  - **Random-walk FM** (h₋₂): the frequency itself wanders, as an exact two-state process with q = 2π²·h₋₂.
- **The Allan deviation** follows: σ_y²(τ) = h0/(2τ) + 2 ln 2·h₋₁ + (2π²/3)·h₋₂·τ (Van Dierendonck et al. 1984).
- **The crystal preset** (that paper's 3e-11 oscillator), recomputed today:

  | τ | σ_y | τ·σ_y (time error scale) |
  |---|---|---|
  | 1 s | 4.3e-11 | 43 ps |
  | 1 h | 1.8e-10 | 0.66 µs |
  | 1 day | 8.8e-10 | **76 µs** |
  | 30 days | 4.8e-9 | 12.5 ms |

  Over days random-walk FM dominates. That is why the clock looks like a polar offset in the fit, and why the
  oscillator choice is step E5.
- **The SDR's ±2 ppm** enters only as `'Rate'` (deterministic); its noise coefficients are unknown (measure on the
  device).
- **Truth side only:** the observer never sees ε. navFit and navFilter estimate the clock; the observer folds with
  its own clock estimate (makePredictor `'Clock'`, day 3).

**Code** (makeClock.m):

| Lines | Block |
|---|---|
| 72–87 | presets, grid of `'Step'` (default 60 s), own random stream |
| 89–92 | white FM |
| 94–106 | flicker by FFT |
| 108–119 | random-walk FM, exact two-state steps |
| 121–134 | ε(t) as a function: deterministic part plus x linear between nodes; rate(t); sigmaY(τ) |

The left-out white-FM bridge between nodes is √(h0·Step/8) = 0.12 ns (checked).

**Tests.** testMakeClock, 5 checks (Jasper's runs, 22 and 23 Oct):
1. the deterministic algebra;
2. seeds and evaluation;
3. the Allan variance per noise type from 5 min to 17 h, within 3 % (40 clocks × 20 000 steps);
4. the Kalman relation from a known start, Q11 = 8.5 ns after 200 s (checked: 8.50 ns);
5. refusals.

Not covered: the flicker start transient (the paper's 11 ns at 200 s depends on it); `'Drift'` is only tested
deterministically.

## Stage 3 question: a carrier at the channel-48/49 edge (Jasper, 26 Oct)

Jasper: "there is a carrier at 1350 MHz being neglected, under the reason this is a 'special case'. But this is a
realistic case, so can you investigate whether the pipeline breaks at that point?"

**Background.**
- 1350 MHz = 1200 + 48 × 3.125 MHz, exactly the edge between channels 48 and 49.
- Each channel's dedispersion taper (sin², 62.5 kHz) is zero there.
- On 8 Oct the scenario carrier was moved to 1351.3 MHz so that it would not look harmless "for the wrong reason"
  (receiver.html: "an unrepresentative special case", i.e. too easy, not broken).

**Test** (Claude's run, 26 Oct, 14:33–14:47; a scratch copy of runRFITest, results in the scratchpad; the reference
`runRFITest_results.mat` untouched): main's settings (−5 dB pulsar, L 0.2 s, 15 sub-ints), channel path, without
and with excision:

| Case | Excision off: good, SNR, σ | Excision on: good, SNR, σ, blanked |
|---|---|---|
| noise only | 15/15, 92.4, 2.938 µs | 15/15, 92.0, 2.951 µs, 0.0114 % |
| cw 1351.300 MHz (−5 dB, the scenario's) | 15/15, 92.1, 2.950 | 15/15, 91.6, 2.963, 0.0113 % |
| **cw 1350.000 MHz (edge)** | 15/15, **92.4, 2.938** (= noise only) | 15/15, 92.0, 2.951 |
| cw 1350.030 MHz (in the taper) | 15/15, 92.1, 2.949 | 15/15, 91.6, 2.963 |
| cw 1349.970 MHz (in the taper) | 15/15, 92.0, 2.953 | 15/15, 91.7, 2.961 |
| cw 1350.100 MHz (past the taper) | 15/15, 92.1, 2.950 | 15/15, 91.6, 2.963 |
| cw 1351.300 MHz, +10 dB | 15/15, 92.1, 2.950 | 15/15, 91.6, 2.963 |
| **cw 1350.000 MHz, +10 dB** | 15/15, 92.4, 2.938 | 15/15, 92.0, 2.951 |

In every case: χ²_red 1.00, rms/pred 0.86–0.88, total-fold offset +0.15…+0.22 ± 0.71 µs, no extra flags, no channel
with r > 1.1.

**Conclusion: nothing breaks.**
- At the exact edge the carrier falls into the channel gap (the existing ~2 % gap of item 3) and costs nothing more.
- Inside the taper or past it, a carrier costs what any in-channel carrier costs: optimal weighting notches its
  channel, −0.4 % SNR (1 of 128 channels). That holds whatever its strength (−5 and +10 dB of the whole band, i.e.
  +16 and +31 dB in its channel).
- A constant carrier raises no flags.
- **Limits of the test:** one seed and 15 TOAs per case (effects below ~1 % are not resolved), and only constant
  carriers. A narrow *noise-like* signal straddling the edge would make both channels lopsided (the B5e/Inmarsat
  effect, χ² above 1). That is not tested here.

## Docs: the 4a and 7a pages merged into their stages (Jasper, 26 Oct)

Jasper: the separate "4a" in the docs sidebar "looks weird"; merge the pages, and "fold the folding page as well".
- `clock.html` became the section "4a · The receiver's clock" on `receiver.html` (anchor `#clock`, sub-anchors `clock-*`).
- `prediction.html` became the first part of `fold.html` ("7a · Phase prediction", anchors `#prediction`, `#prediction-why`,
  `#trajectory`, `#delay`, `#ephemeris`, `#makepredictor`). The fold's own sections now sit under "7 · Folding: `foldProfile`"
  (`#fold`), one heading level down; their anchors are unchanged.
- Text and figures unchanged (word diff against the old pages). Only the headings changed: the prediction page's section
  numbers 1, 2, 3, 3b, 4 became "7a", and the page titles were merged. Sidebar, home-page cards, contents and every link
  (about 25 files, figure scripts included) retargeted; a site-wide link check finds no broken page or anchor. Rendering checked.
- Committed `ddee4f4` (clock) and the prediction merge after it; changelog rows in `docs/maintenance.html`. No code change.
  The pending review-week edits (roadmap E0 row, notes, this log) were committed first as `bafc774`, on Jasper's go-ahead.

**Follow-up: a narrow noise-like signal straddling the edge** (Jasper: "also test a narrow noise-like signal (e.g.
200 kHz) straddling the 1350 MHz edge, as a scratch run").
- Claude's run, 14:49–17:03; slow, the machine was shared.
- Scratch copy; three passes per case: without excision, with excision, and excision + per-channel noise
  calibration (`noiseCal`).
- 200 kHz noise-like RFI, like one narrow Inmarsat carrier.

| Case | SNR (no exc. / exc. / cal) | χ² (no exc. / exc. / cal) | rms/pred (exc. / cal) | channels r > 1.1 (r, level) | error bars too small |
|---|---|---|---|---|---|
| noise only | 92.4 / 92.0 / 92.0 | 1.00 / 1.00 / 1.00 | 0.86 / 0.86 | – | 0.02 % |
| 1351.30 MHz, in ch 49, −21 dB | 92.2 / 91.8 / 91.7 | 1.02 / 1.02 / 1.00 | 0.87 / 0.87 | ch 49 (3.9, 2.0) | 0.28 % |
| **1350.00, straddling, −21 dB** | 92.0 / 91.7 / 91.4 | 1.02 / 1.02 / 1.00 | 0.93 / 0.89 | **ch 48 (2.0, 1.3), ch 49 (2.0, 1.3)** | 0.46 % |
| 1350.05, half in the taper, −21 dB | 92.1 / 91.8 / 91.7 | 1.02 / 1.02 / 1.00 | 0.91 / 0.88 | ch 49 (3.1, 1.6) | 0.34 % |
| 1351.30, in ch 49, −11 dB | 92.1 / 91.6 / 91.7 | 1.07 / 1.07 / 1.00 | 0.87 / 0.87 | ch 49 (10.3, 11.4; 23 flags) | 0.04 % |
| **1350.00, straddling, −11 dB** | 91.6 / 91.3 / 91.4 | **1.17 / 1.17 / 1.00** | 0.95 / 0.87 | **ch 48 (12.2, 4.2), ch 49 (12.1, 4.2)** | 0.51 % |

All cases: 15/15 good TOAs, offsets +0.09…+0.29 ± 0.72 µs.

**Reading: nothing breaks.** Straddling is the worse version of the known "lopsided channel" effect (B5e,
Inmarsat-narrow):
- Two channels are affected instead of one.
- The edge tapers cut the middle of the signal, so each channel's level stays low (4.2 instead of 11.4 at −11 dB).
  Optimal weighting therefore keeps these channels at a much higher weight (1/level²: 1/18 instead of 1/130), while
  their noise is ~12× the flat model.
- So χ² rises to 1.17 (1.07 for the same signal inside one channel).
- The TOA error bars are only ~0.5 % too small (2 of 128 channels).
- **The per-channel noise calibration brings χ² to 1.00 in every case.** This adds weight to the audit's open
  question 3 (should `noiseCal` be on by default?).

## Day 1, stages 5–6: the full-band front end (`applyInverseDispersion`, `detectPower`)

**Physics.**
- **Coherent dedispersion.** Interstellar dispersion is a pure phase filter: a frequency f arrives
  τ(f) = K·DM/f² late, with K = 4.148808e3 s MHz² pc⁻¹ cm³.
  - In the frequency domain the inverse is H = W(f)·exp(−i2π·K·DM·(f − fRef)²/(f·fRef²)). It removes the sweep
    exactly and keeps fRef (1.6 GHz, the top of the band) at its own time.
  - At DM 25 the sweep across 1.2–1.6 GHz is **31.5 ms** (checked). Frequencies below 1.6 GHz are advanced onto it.
  - The bulk delay of 1.6 GHz relative to infinite frequency, τ(fRef) = **40.5 ms**, stays in: TOAs refer to 1.6 GHz.
  - W is a sin² taper over 8 MHz at each band edge (2 % of the band). A sharp edge would ring (a sinc impulse
    response); the taper keeps the impulse response compact. The leakage beyond the overlap is measured (main.m:
    5.8e-13).
- **Overlap-save streaming.** Each output sample needs the input from nPast samples before to nFuture after. With
  fRef at the top: nPast = guard (2000 samples = 2.5 µs, 20 / edge width) and nFuture = sweep·fs + guard = **25.2 M
  samples** (checked).
  - The FFT size has the lowest total cost within 4 GB: 2^26 at DM 25.
  - The first nPast and the last nFuture samples (31.5 ms) are not "fully supported" and are excluded downstream.
- **Square-law detection.** |z|² averaged over time bins of 1/f_out = 1 µs (800 samples). The data rate collapses
  from 800 MS/s complex to 1 MS/s of power.
  - The relative variance per bin is 1/(B_noise·dt). B_noise = (∫W²)²/∫W⁴ = 391.6 MHz for the observer's taper
    (noiseBandwidth), about 392 independent samples per bin.
  - The bin time is its centroid, so a symmetric boxcar gives no timing bias.

**Code.**

applyInverseDispersion.m:

| Lines | Block |
|---|---|
| 107–134 | checks (band inside fLO ± fs/2, RefFreq inside the band unless allowed) |
| 136–148 | overlap sizes from τ(f) and the guard |
| 164–181 | FFT size by cost within the memory budget (`chooseNfft`) |
| 183–199 | H built in chunks, in double, only on in-band bins |
| 201–213 | leakage check |
| 242–271 | overlap-save loop on one preallocated buffer |
| 273–305 | info: bulkDelayRef, fullySupported, edgeWidth |

detectPower.m:

| Lines | Block |
|---|---|
| 89–96 | binLen = round(fs/f_out) |
| 98–133 | fast path (one channel) or sub-band FFT frames |
| 177–205 | stream: I² + Q² summed per bin in double |
| 207–238 | info: binTime0 (centroid), fullySupportedBins from the dedispersion's range |

**Tests.**
- The full band is the reference that `testDedisperseChannels` compares the channel path against.
- `detectPower` runs inside `testDetectChannels` and `testDedisperseChannels`.
- The 5 Oct full-band runs are validated (notes §7).
- Missing: an absolute round-trip test (dispersion → dedispersion gives back the band-limited input: R-V2, noted 16
  Oct), and the NumPy ports, which are not in the repo.

## Day 2, stage 5a: the filterbank (`channelizeIQ`)

**Physics.**
- **What it does.** It splits the 800 MS/s IQ stream into 128 channels of 3.125 MHz that tile 1.2–1.6 GHz. Each
  channel is a complex voltage stream at 4.17 MS/s.
  - It is exactly what 128 separate receivers would give (mix down to the channel centre, low-pass, keep every
    192nd sample), computed at ~1/256 of the cost by a polyphase filterbank.
  - Per output time: weight a frame of 4096 input samples with the prototype filter, fold it into 256 bins, and one
    256-point IFFT gives all channels at once, followed by a phase correction from the absolute sample index.
- **Why channels:**
  1. Narrowband RFI stays in one or two channels, where it can be excised or down-weighted.
  2. The dispersion sweep inside one channel is short: **0.38 ms at 1.2 GHz and 0.16 ms at 1.6 GHz at DM 25**
     (checked), against 31.5 ms over the whole band. Small FFTs per channel are enough, which solves the memory
     problem.
  3. Each channel gets its own noise model and weight (optimal weighting).
- **Why oversampled (4/3).** The channel rate (4.17 MHz) is above the channel width (3.125 MHz).
  - The prototype is flat to ±1.5625 MHz and stops from 2.604 MHz (80 dB), with 3853 taps, 16 per branch (checked).
  - Aliasing on decimation lands only outside ±1.5625 MHz, so every channel's useful band is alias-free. A
    critically sampled filterbank would alias into the channel edges.
  - The transition region is cut later, by the per-channel dedispersion's window.
- **The grid.** The channel centres are fLow + (j − ½)·3.125 MHz. The 1.3 GHz LO lies on a channel edge (between
  channels 32 and 33), so the half-slot offset is β = ½. The group delay (1926 input samples, 2.41 µs) is removed
  exactly.

**Code** (channelizeIQ.m):

| Lines | Block |
|---|---|
| 108–135 | the grid: K = 256, D = 192, the channel count, β from fLO |
| 137–157 | the Kaiser prototype (length from the transition width), the β-modulated taps, the coverage check |
| 207–241 | the stream: frames by index, ×taps, fold into K bins, IFFT, keep the 128 rows, exact phase correction, one file per channel |
| 243–281 | info: chanFreqs, the mapping (no inversion), fullySupported |

**Tests.** testChannelizeIQ (validated 6 Oct):
- tones against direct mix-filter-decimate receivers to 2.7e-7;
- passband ripple 1.2e-4;
- stopband −80.2 dB;
- noise power 1.0007;
- Parseval on real data 1.00007;
- 6 blocks.

Only today's grid is tested (R-C1: β = 0, oversampling 2, `BlockOut` bit for bit, the refusals).

## Day 2: RFI excision (`detectRFI`, `periodicRFI`, `blankChannels`)

**Physics.**
- **Why before dedispersion.** In the channel samples (0.24 µs each) a 2 µs radar pulse is a few samples in one or
  two channels. After dedispersion it would be smeared over the in-channel sweep (0.16–0.38 ms at DM 25).
- **Why blank instead of subtract.** At −54 dB one radar pulse holds the energy of ~80 000 pulsar pulses: even
  99.9 % subtraction would leave more RFI than pulsar. Blanking sets signal and noise to 0 together; the exact
  weights follow in blankingWeights.
- **The statistic.** Per channel, the power |x|² is averaged over windows of 1, 2, 4, 8 and 16 samples
  (half-overlapping).
  - Noise power per sample is exponential. A window's sum is Σλ_k·E_k, with λ_k the eigenvalues of the channelizer's
    sample correlation (ρ(1) = 0.05).
  - The exact survival function (matrix exponential) gives thresholds of 13.82 / 8.39 / 5.39 / 3.68 / 2.69 × the
    baseline at P_FA 1e-6 per window.
  - Expected false windows over 0.2 s, all channels: 106.7 / 106.7 / 53.3 / 26.7 / 13.3 (rechecked).
- **The baseline** is median(|x|²)/ln 2 per 16 ms block (robust; exact for exponential power), from fully supported
  samples. A second pass re-estimates it without the flagged samples.
- **Per-channel thresholds (since 21 Oct).** A narrow noise-like carrier is still Gaussian but correlated over ~20
  samples, so white thresholds flag its own peaks. From pass 2 a channel takes thresholds from its own measured
  autocorrelation (lags 0–15), but only if they move by > 1 %; clean channels stay bit for bit.
- **The guards:**
  - **time:** ±10 channel samples around each flagged window (the prototype's half-length / D);
  - **frequency:** a strong event (peak ≥ 100 × the baseline) is also blanked in ±7 neighbouring channels at the
    same time, against spectral sidelobes that a pulsar-locked radar would fold up coherently (B5b).
- **`blankChannels`** zeroes the masked samples in copies (the raw files stay).
- **`periodicRFI`** predicts and blanks every pulse of a periodic emitter, also the undetected weak ones. It is off
  by default since 9 Oct (the false-emitter bug; design parked).
- **Validity.** Up to a pulsar of about −5 dB. A brighter pulse flags itself (10 % of its 16-sample windows at 0 dB),
  which the weights do not model. The target (−25 … −54 dB) cannot trigger anything.

**Code.**

detectRFI.m:

| Lines | Block |
|---|---|
| 177–205 | memory check, guard, thresholds from the channelizer's ρ (`windowThresholds`, exact via `sumExpSurvival`) |
| 207–225 | baseline blocks |
| 227–308 | per channel and pass: own correlation and thresholds (pass ≥ 2), baseline per block, windows by cumulative sums, coverage with the guard, intervals, peaks → strong events |
| 309–318 | frequency guard; the mask is normalized |
| 320–366 | info |

blankChannels.m: block-wise copy with the masked samples set to 0 (exact, independent of the block size).

**Tests.**
- testDetectRFI, validated:
  1. noise only: baseline, flag counts = nWindows·P_FA at 1e-3 and 1e-4 (the exact distribution, not the Gamma
     approximation), blanked fraction;
  2. the −5 dB pulsar: flags on and off the pulse against the prediction;
  3. RFI: radar energy removed > 0.9999, impulses > 0.95, the carrier and GNSS not blanked; 3d the frequency guard;
  4. per-channel thresholds: three 200 kHz carriers; P_FA 1e-3 flags/expected 1.07 / 1.06 / 1.12 / 1.15 / 0.96
     against 1.8–63× with white thresholds.
- testBlankChannels: exact zeros, mask normalization, block size, empty mask, interface.
- testPeriodicRFI, testRadarEdges, testRfiGating, testRfiNoise; runRFITest (the B4/B5e tables), runLockedRadar.

## Day 2, stage 5b: coherent dedispersion per channel (`dedisperseChannels`)

**Physics.**
- **Same filter, smaller pieces.** Dispersion is a linear filter, so the full-band chirp
  H(f) = W(f)·exp(−2πi·K·DM·(f − f_ref)²/(f·f_ref²)) can be applied piece by piece. Each channel gets the same
  function of RF frequency, restricted to its 3.125 MHz. Filters commute, so channelize-then-dedisperse equals
  dedisperse-then-channelize inside each channel's flat band (test 2).
- **Why coherent inside a channel.** The smear inside one 3.125 MHz channel is 374 µs at the bottom and 159 µs at the
  top (DM 25). Only shifting the channels (incoherent dedispersion) would leave that smear.
- **No separate shift step.** f_ref is the band top (1.6 GHz) for every channel. Each channel's chirp then also holds
  its delay to 1.6 GHz: 31.5 ms for the bottom channel, 12.6 ms mid-band and 0.16 ms at the top. All channels come out
  aligned at 1.6 GHz on the same sample grid. The bulk delay τ(1.6 GHz) = 40.5 ms is not removed (the predictor
  carries it).
- **Channels tile the band.** The chirp is zero outside ±1.5625 MHz with sin² edges of 62.5 kHz (2 % of a channel).
  This cuts away the filterbank's transition region (oversampling 4/3).
  - Neighbouring channels do not overlap, so their noise is independent and the detected powers simply add.
  - Per channel ∫W² = 3.047 MHz, × 128 = 390.0 MHz: exactly the full-band path's 390.0 MHz (8 MHz edges).
  - B_noise per channel = 3.060 MHz, × 128 = 391.6 MHz.
  - The price: a 125 kHz notch at each channel boundary (D2-2, the edge-carrier test).
- **Why per channel: memory.**

  | DM | full band at 800 MS/s | bottom channel at 4.17 MS/s |
  |---|---|---|
  | 25 | FFT 2^26, 3.8 GB | 2^19, 29 MB |
  | 50 | 2^27, 7.5 GB | 2^20, 59 MB |
  | 100 | 2^28, 15 GB | 2^21, 117 MB |

  The FFT must exceed twice the overlap (sweep × rate + guards). With the 4 GB default the full-band path stops at
  DM ≈ 33.
- **The guard** is 20 / edge width = 320 µs (1334 channel samples) on each side, for the taper's ringing. The leakage
  self-check (one extra IFFT per channel, warning above 1e-6) is on by default.

**Code** (dedisperseChannels.m, 147 lines): a loop that calls the validated applyInverseDispersion per channel.
- The call at :84 passes the channel band fc ± dF/2, RefFreq 1.6 GHz and AllowRefOutsideBand. With the reference
  outside the channel, nPast is the guard only and nFuture is the delay to 1.6 GHz plus the guard
  (applyInverseDispersion.m:146–147).
- The common valid range at :102 is the channelizer's range minus the largest nPast/nFuture.
- B_noise at :103.
- Edge width, bulk delay and N come from channel 1: the same for every channel by construction.

**Tests** (testDedisperseChannels, validated):
1. the full-band path is bit-identical after the AllowRefOutsideBand change;
2. commutation in the interior channels: |ρ| > 0.9999, phase < 1e-3 rad, power within 0.2 %. This checks the chirp
   segments, the inter-channel delays and the phase reference at once;
3. end to end against the full band, same data: TOA differences against their predicted scatter (χ², mean), SNR
   against the taper prediction, fold noise ratio;
4. detected-power noise per time bin: in a 3 MHz channel the power is correlated over ~0.33 µs against 0.96 µs bins
   (bin rate 1.0417 MHz = channel rate / 4). The variance is below the radiometer value and the lag-1 correlation is
   positive; both match the exact prediction and nearly cancel within a phase bin.

## Day 2: missing unit tests, stage 2 (`testApplyDispersionStream`, R-D1)

Jasper, 27 Oct: filling in the missing unit tests is part of the review week ("let's consider filling in the
missing unit tests part of the review week"). Claude may run MATLAB itself; a test counts as validated only after
Jasper's run.

**The gaps (sweep of tests/ against functions/, earliest stages first).**
- Stage 1 has its test (R1, 13 Oct).
- Stages 2–4 had none: `applyDispersionStream` (R-D1), `addNoiseAndRFI` (R-N1: the noise level against the SNR
  definition; the older RFI types only as a regression), `applyIQmodulation` (R-I1).
- Stage 5 (`applyInverseDispersion`): a bit-for-bit regression and the comparison with the channel path only, no
  absolute round trip (R-V2). Stage 6 (`detectPower`): only indirect checks (D2-18). Stage 5a: only today's grid
  (R-C1).
- Helpers without a test of their own: `noiseBandwidth` (D2-17), `combineChannels` (exercised by the TOA tests),
  `processIQ`'s `'Reuse'` branch (D1-6). No test runs `main.m`.

**Physics checked.** From input file to output file the stage applies H_eff(f) = W(f)·exp(iφ(f)),
φ = 2π·K·DM·(f − fRef)²/(f·fRef²) (the kernel's s0 shift removed again). The group delay is exactly
τ(f) − τ(fRef). The first kernel had the wrong sign and its forward + inverse round trip still worked, so the test
compares every output with this formula, computed in the test in double precision; it never calls the kernel code.

**The test** (`tests/testApplyDispersionStream.m`, ~3 s, own files in tempdir; band 1.2–1.6 GHz at 4 GHz; DM 0.05
except check 3). All seven checks pass in Claude's runs (twice, also inside `runAllTests` with the generator test) and in
Jasper's run (27 Oct, the same numbers; with testGeneratePulsarSignal): validated. The two warnings in the
output are provoked on purpose (check 6 here, the generator's 25 % duty check).

| Check | Result | Limit |
|---|---|---|
| 1 format and info | 21 of 21 fields | exact |
| 2 transfer function from one impulse (RefFreq fHigh, 1.4 GHz, fLow; DM 0; EdgeFrac 0.01, 0.05) | relative error energy 0.9–4.1e-12; group delays within 7.4 ps | 1e-9 |
| 3 five tone bursts at DM 5 (1.25–1.55 GHz, σ 0.5 µs) | centroid − prediction ≤ 9e-4 samples; energy kept to 6e-7 | 0.01 sample; 1e-5 |
| 4 overlap-add vs one double FFT (automatic, M − 1 = 15 blocks, 300,007, 400,000 = 1-sample last block, one block) | relative rms 1.1e-6, file ends 8e-7 to 1.2e-6; a 5001-sample file (< s0, flush only) 3.1e-7 | 1e-5 |
| 5 power | kernelEnergy 0.19500000 = 2(B − 1.25·edgeW)/fs; white-noise power −0.33σ (σ 1.7e-3) | 1e-6; 4σ |
| 6 leakage | ≤ 4e-12 in 10 default runs; GuardTime 2/edgeW: 3.59e-8, warns, error energy / leakage = 1.00000 | 1e-10; 1 % |
| 7 input checks | 8 of 8 refused with their ids | — |

- **Prediction for the delays.** The burst's centroid moves by the group delay averaged over its own power spectrum
  (exact where |H| = 1). The plain τ(f0) − τ(fRef) is 2.2–5.2 samples off (curvature of τ across the burst's
  band). The plan said 4–10; the power spectrum is √2 narrower than the amplitude spectrum.
- **Beyond the plan:** from notes §10 D1's spec, the EdgeFrac cases, the file shorter than s0 and the 1e-10 leakage
  limit. The short-guard ratio limit was tightened from a factor 2 to 1 % after the first run gave 1.00000.
- **Not covered:** R-D2 (one refFreq for the simulation and the observer: a design question). DM 25 (main since
  20 Oct) runs the same code with a 5× longer kernel. The sky page's round-trip line (2e-6) still comes from the
  NumPy port; R-V2 replaces it.

**Also changed.** `runAllTests`: the test runs second, count 26. Docs: tests.html (table row, section `#disp`, the
check matrix regenerated by `gen_testmatrix.py`, the data-flow box "1 · sky stages", counts 26); sky.html
validation (the NumPy-port line replaced by the test's numbers); roadmap R-D1 done; figures README. Notes: §5.2
Tests paragraph, §10 D1 done, tests table row and the runAllTests row. `api.html` unchanged (tests are not in it).

## Day 2: missing unit tests, stage 3 (`testAddNoiseAndRFI`, R-N1)

Jasper: "lets go on to the next unit test"; plan given (physics, six checks), then "Go ahead".

**Physics checked.**
- **SNR.** `SNRdB` is the pulse-peak power spectral density over the noise's, in the band. For a white signal of std A
  and white noise of std σ_n at the same rate, ρ = A²/σ_n², so σ_n = A/√ρ.
- **The band's share.** White noise spreads over 0 … fs/2, so the band holds P_band = σ_n²·B/(fs/2) = 0.2 σ_n².
- **INR** = RFI power while on / P_band. Carrier amp √(2·INR·P_band) (power amp²/2); bpsk the same with ±1 chips;
  an impulse burst is white over the whole band, amp √INR·σ_n, so it puts the same share B/(fs/2) into the band.
- **Block size.** Phases from the absolute sample index, random parts from counter-based substreams.
- **Prediction.** SNR² = B·∫(ρp)²/(1 + ρp)² dt; weak pulses ρ·√(B·σ·√π) = 388·ρ for σ = 0.212 ms, B = 400 MHz.

**The test** (`tests/testAddNoiseAndRFI.m`, ~4 s, own files in tempdir). All six checks pass in Claude's runs (alone
and inside `runAllTests` with the dispersion test) and in Jasper's run (27 Oct): validated.

| Check | Result | Limit |
|---|---|---|
| 1 format and info | 17 of 17 (σ_n, noisePSD, P_band exact; no noise = input bit for bit) | exact |
| 2 SNR definition (white input of std A, −5 dB) | in-band ratio 0.31712 vs ρ 0.31623 (+1.26σ); variance −0.64σ; 8 sub-bands within 2.08σ; noise independent of the input | 4σ |
| 3 RFI vs INR | cw = formula to 1.5e-7 of amp, power 1.000000; bpsk chips ±1 to 5.9e-6, 1023 chips constant, power 1.00000; impulse 104 bursts, energy +0.53σ, in-band share 0.1991 (−0.43σ) | 1e-6, 1e-4 / 1e-4, 1e-3 / 4σ |
| 4 block size (noise + 6 sources, 2.5 ms) | 4e6, 1,000,003 and one block bit for bit | exact |
| 5 prediction | = quadrature to ~1e-16 (−54 dB, 0 dB, amplitude mode); weak forms −3.3e-6 / +2.2e-6 at −54 dB; 387.99·ρ | 1e-6; 1e-5 |
| 6 input checks | 7 of 7 refused with their ids | — |

- **Beyond the notes' N1 spec:** whiteness, input independence, the 0 dB (self-noise) and amplitude-mode predictions.
  `pulsed` power is not repeated: testRadarEdges 3 already compares the radar with a formula whose amplitude comes
  from INR. `noise` stays in testRfiNoise 3.
- **Not covered:** the receiver noise is white up to fs/2, not band-passed (known simplification, receiver page);
  no ADC (D1-11).

**Also changed.** `runAllTests`: third, count 27. Docs: tests.html (row, section `#rx`, check matrix regenerated, the
data-flow box "1 · stages 1–3", counts 27); receiver.html validation (the "no unit test yet" line replaced by the
test's numbers); roadmap R-N1 done; figures README; changelog. Notes: §5.3 Tests paragraph, §10 N1 done, tests
table row and the runAllTests row.

## Day 2: missing unit tests, stage 4 (`testApplyIQmodulation`, R-I1)

Jasper: "plan the next unit test (R-I1, applyIQmodulation)"; plan given, then "Go ahead".

**Physics checked.**
- **Complex envelope.** x = Re{z·e^(i2π fLO t)}; mixing gives 2x·e^(−i2π fLO t) = z + z*·e^(−i4π fLO t). The low-pass
  keeps z (the band at f − fLO) and removes the real signal's mirror at −(f + fLO). A tone A·cos(2πft + φ) becomes
  A·e^(iφ)·e^(i2π(f − fLO)t): amplitude A (factor 2), no inversion.
- **The low-pass.** Kaiser, 2049 taps, 80 dB: ripple and stopband both δ = 1e-4; flat to ±355.1 MHz, stopband from
  364.9 MHz; the band (−100 … +300 MHz) sees 9e-6 ripple. After ÷5 the Nyquist is ±400 MHz.
- **Timing.** Symmetric FIR: linear phase, delay exactly G = 1024 input samples, removed; with it removed the passband
  response is real, so a tone keeps its phase. LO phase from the absolute index (13/40 cycles per sample).

**The test** (`tests/testApplyIQmodulation.m`, ~2 s, own files in tempdir). All seven checks pass in Claude's runs
(alone and inside `runAllTests` with the stage 2–3 tests) and in Jasper's run (27 Oct): validated.

| Check | Result | Limit |
|---|---|---|
| 1 format and info | 22 of 22; 4,000,003 → 800,001 samples; no warning at main's settings | exact |
| 2 six tones (1.2 … 1.6 GHz, one at fLO → DC) | = the formula sample for sample within 0.9–5.3e-6 of A | 1e-4 |
| 3 burst timing | centroid on output sample k0, −8.5e-6 sample (204.8 late without removing G) | 1e-3 sample |
| 4 power ratio | 2.000003 | 2 ± 2e-4 |
| 5 suppression | mirror −134.3 dB, stopband tone −91.7 dB, aliasing tone −103.8 dB | −80 dB |
| 6 block length (2048 = 1954 blocks, 99,999, 1,000,003) | 2.8–3.0e-7 rms from one block | 1e-6 |
| 7 warnings and refusals | :rate, :band, :alias; 8 of 8 refused | — |

- **Predicted beforehand:** the suppression levels were computed from the filter's design (NumPy, before the run) and
  came out the same to 0.1 dB. The limits are the design limits, not the measured values.
- **Warnings silenced.** The three provoked warnings are switched off inside the check; `lastwarn` still records
  them (checked), so the output shows no warning text.
- **Not covered:** `'Gain', 'unity'` (R-I4: remove); the ideal receiver (R-I3); no ADC (D1-11).

**Also changed.** `runAllTests`: fourth, count 28. Docs: tests.html (row, section `#iq`, check matrix regenerated, the
data-flow box "1 · stages 1–4", counts 28); receiver.html IQ validation (the NumPy-port line replaced by the test);
roadmap R-I1 done; figures README; changelog. Notes: §5.4 Tests paragraph, §10 I1 done, tests table and the
runAllTests row.

## Day 2: missing unit tests, stage 5 (`testApplyInverseDispersion`, R-V2) and a leakage-check fix

Jasper: "plan the next unit test (R-V2, the dedispersion round trip)"; plan given (a new test instead of extending
the validated testDedisperseChannels), then "Go ahead".

**Physics checked.**
- **The inverse filter.** Per baseband bin, RF f = fLO + f_bb: H = W(f)·e^(−iφ(f)), the conjugate of the forward
  chirp; group delay −(τ(f) − τ(fRef)): the sweep collapses onto fRef; τ(fRef) (the bulk delay) stays in.
- **Overlap-save reach.** nFuture = ⌈(τ(fLow) − τ(fRef))·fs⌉ + guard, nPast = ⌈(τ(fRef) − τ(fHigh))·fs⌉ + guard, each
  ≥ guard; outputs nPast + 1 … N − nFuture fully supported.
- **Round trip.** Forward × inverse = W² (the chirps cancel): the IQ of the input, band-limited by W².
- **Channel path.** The same filter per channel with one common fRef; channels above fRef reach into the past.

**The test** (`tests/testApplyInverseDispersion.m`, ~2 s, own files in tempdir). All seven checks pass in Claude's
runs (also inside `runAllTests` with testDedisperseChannels and testDetectChannels) and in Jasper's run of the
same three (27 Oct): validated, the fix included.

| Check | Result | Limit |
|---|---|---|
| 1 format and info | 23 of 23 | exact |
| 2 impulse response (RefFreq fHigh, 1.4 GHz, fLow, 1.8 GHz above, 1.1 GHz below; DM 0) | error energy 4–23e-14; group delays within 13 ps | 1e-9 |
| 3 round trip RF → dispersion → IQ → dedispersion | input × W² to 1.93e-6 | 1e-4 |
| 4 channel path, RefFreq 1.4 GHz, DM 25 (4 fake channels) | each = formula (≤ 5.9e-14); nPast up to 53,013 above, nFuture up to 80,957 below; common range right; DM 0 taper only | 1e-9 |
| 5 overlap-save vs one double FFT (4 Nfft) | 2.4–2.8e-7 | 1e-5 |
| 6 leakage vs its formula | 1.170e-12 / 1.154e-12; 3.570e-8; 1.212e-6 (warns); noise error ≤ 0.065 × leakage | 1 % + 1e-13; ≤ 2× |
| 7 input checks | 10 of 10 | — |

**Found and fixed: the leakage check summed in single precision** (applyInverseDispersion.m:205–206, bug §8 #23).
- **Symptom (first run, check 6):** reported leakage 1.79e-7 at the default guard and 5.96e-7 at 0.5/edgeW (so the
  10× shorter guard did not warn), but 0 at RefFreq 1.4 GHz, DM 0 and DM 5. All are multiples of 2^-24 = 6e-8: near 1,
  1 − inside/E in single cannot resolve less. The true value (double sums) is 1.17e-12.
- **Fix (Jasper: "Fix it"):** `sum(…, 'double')` for E and inside. The |h|² values themselves are fine in single
  (1.165e-12 vs 1.154e-12 fully in double). No extra memory.
- **Validation:** the HEAD copy (renamed `_old`) and the new code on 7 settings (full band: default, RefFreq 1.4 GHz,
  Nfft 2^17, DM 0, RefFreq above the band, GuardTime 0.1/edgeW; one DM 25 channel): output files identical byte for
  byte, info identical except `leakage`; testDedisperseChannels' bit-for-bit regression on the frozen file passes;
  checkcode unchanged (two stale `%#ok` notes on the onCleanup lines, already in HEAD).
- **Plan corrected:** the plan said "error energy = leakage within 1 %, as in stage 2". That holds for a truncated
  kernel, not for overlap-save: a lag beyond the overlap disturbs only the start of each block (later outputs still
  find its true sample in the buffer). So check 6 compares `info.leakage` with the same quantity from the formula,
  and keeps the error as a bound (≤ 2 × leakage; measured 0.0004–0.065 ×).

**Also changed.** `runAllTests`: fifth, count 29. Docs: tests.html (row, section `#inv` with a History box, check
matrix regenerated, data-flow box "1 · stages 1–5", counts 29); fullband.html (leakage sentence ~1e-12, a History box
in the leakage design box, a new Validation list for stage 5); sky.html (the round-trip line from the NumPy port →
this test, 1.9e-6); roadmap R-V2 done; figures README; changelog. Notes: §5.5 (the fix and a Tests paragraph), §8 #23,
§10 V2 done, tests table and the runAllTests row.

## Running list of smaller details

| # | Part | Where | Detail | Weight |
|---|---|---|---|---|
| D1-1 | processIQ | processIQ.m:220 | With `noiseCal` and no excision, `foldArgs(1:end-2)` drops the last two elements, on the assumption that they are the `'NoiseCoeffs'` pair. This is positional and breaks silently if the argument order changes | minor robustness |
| D1-2 | processIQ | processIQ.m:198–223 | No `'Predictor'`: the fold uses the linear f0 / F1 model only. Real data (F) need the predictor (Doppler, delays); `foldProfile` supports it | for F |
| D1-3 | processIQ | processIQ.m:211–216 | The noise calibration's on-pulse window uses the f0 / F1 phase (P3-3); on real data the pulse drifts by ~1e-4·t | for F |
| D1-4 | processIQ | processIQ.m:114 | Excision runs only on the channel path; with `frontEnd 'fullband'` the `excision` setting is ignored without a message | minor, docs |
| D1-5 | processIQ | processIQ.m:302–323 | `checkFold` checks that a noise calibration was used, but not with which on-pulse window or coefficients: a changed `noiseCalOnPulse` keeps a stale fold without a warning | minor |
| D1-6 | processIQ | tests | The `'Reuse'` branch and its checks are untested; a reuse test belongs with the main.m smoke run (day 4) | test gap |
| D1-7 | foldToTOA | foldToTOA.m:67 | The good-TOA rule requires detection at UNKNOWN phase (threshold 4.15 at P_FA 1e-3) even when the predicted phase is known (3.09). This guards against wrong-phase TOAs. For navigation with a good predictor it discards TOAs detected only at the known phase; at worst-case flux the adaptive groups close at T0 ≥ 7, so both thresholds pass | design question (day 3) |
| D1-8 | foldToTOA | foldToTOA.m:61–66 | `estimateTOA` and `detectPulsar` each build the noise model of the same fold | efficiency, low |
| D1-9 | IQ stage | applyIQmodulation; deliverables/task2_4_signal_chain.tex | **The LO plan.** The simulator uses an in-band LO at 1.3 GHz, off-centre on purpose so that a mirrored frequency mapping would show (parameters.html). The 30 Sep set-up had 1.4 GHz. The task 2.4 write-up (another session, 21 Oct) *proposed* a high-side LO at 2.0 GHz with an inverted IF. Its log calls that a proposal, not a decision, and lists the alternative, the SDR's own mixer at 1.4 GHz. Jasper (26 Oct): "I thought we'd stick to 1.3/1.4GHz". So: keep the in-band LO. The relevant real-hardware gap is then D1-10 (I/Q imbalance; DC/LO leakage at the LO, in band, e.g. at 1.4 GHz in the protected 1400–1427 MHz band). The write-up's text presents the 2.0 GHz LO as the plan: Jasper may want to correct it there | LO decided in-band (Jasper); deliverable text to check |
| D1-10 | IQ stage | applyIQmodulation | An ideal DDC (R-I3): no I/Q imbalance (an image mirrored around the LO), no DC/LO leakage (at 0 Hz, which is inside our band with the 1.3 GHz LO), no LO phase noise | known gap |
| D1-11 | IQ stage | float32 output | No quantization. A real SDR delivers integer samples (8–16 bit). For noise-dominated data the SNR loss is small at 8 bits and more at 2–4, but a strong RFI line can clip a few-bit ADC. ADC bits and the RFI headroom belong in E4. Added 26 Oct (Q&A, Claude's estimates, not run): the weak pulsar makes quantization a fixed SNR factor with no TOA bias, and the channelizer re-Gaussianizes the samples, so V, X and the 6a thresholds stay valid (efficiency ~0.64 / 0.88 / 0.96 / 0.99 for 1 / 2 / 3 / 4 bits, Thompson, Moran & Swenson). The real issue is headroom: the ADC gain is set by the strongest RFI, not the noise. Keeping the realistic +40 dB radar (amplitude ~100x noise rms) unclipped leaves the noise at ~1 LSB with 8 bits (~10 % extra noise), ~4 LSB with 10 bits (<1 %), ~15 LSB with 12 bits (negligible). Alternative: let the radar clip; clipping spreads it over all channels but only during its 2 us pulses, so detectRFI's frequency guard should blank them (~1 % of the data at 973 Hz PRF with guards; untested). Model: an N-bit / full-scale / clip option at the output of applyIQmodulation (the SDR samples I/Q after the analog mixer); a fixed-point FPGA filterbank adds its own rounding later. Needs the SDR's bit count | for E4/F |
| D1-12 | IQ stage | LO and sample clock | One reference drives both the LO and the sampling. Its ±2 ppm offset shifts the band by ~2.6 kHz at 1.3 GHz, and it is the same error as the clock's rate in makeClock. The simulation keeps the LO exact and puts the clock error only into the time tags | for E4/F |
| D1-13 | IQ stage | applyIQmodulation 'Gain' | The legacy `'unity'` option is unused (R-I4) | cleanup |
| D1-14 | IQ stage | applyIQmodulation.m:134 | The alias warning (stopband below the output Nyquist) is stricter than our band needs: aliasing into −100…+300 MHz would need content at 700–1100 MHz, which is 80 dB down. Correct, just conservative | none |
| D1-15 | IQ stage | info.t0 = 0 | The IQ time origin is the RF file's start; real data get absolute time stamps from the SDR (UTC → TT; notes §10 item 13) | for F |
| D1-16 | clock | makeClock 'Drift' | Ageing is in the model but 0 in every run. A crystal ageing 1e-10/day gives ½·D·t² = 3.9 ms after 30 days, a quadratic the clock model must absorb. Realistic ageing per oscillator class belongs in E5 | for E5 |
| D1-17 | clock | not modelled | **Temperature.** A crystal's frequency follows its temperature. A daily temperature cycle (site enclosure, ship) gives a daily clock signature, the same period as the daily position signal (solar vs sidereal day: 4 min/day apart). It could bias the position, not just add noise. Not in makeClock, not in the limits | **E5, worth a number** |
| D1-18 | clock | makeClock presets | Only "perfect" and "crystal" (the 1984 paper's 3e-11 crystal, steadier than a stock SDR TCXO: 0.66 µs/h vs the log's ~3.6 µs/h estimate). TCXO, OCXO and Rb presets come with E5's numbers | for E5 |
| D1-19 | clock | navFit.m:41 | Flicker FM is left out of the fit's clock covariance (OK for the crystal at daily scales: ~3 µs against 76 µs; not for other oscillators) | for E5 |
| D1-20 | clock | scope | makeClock covers timescales ≥ Step (60 s). Sub-second jitter (ADC sampling) and LO phase noise are separate and not modelled (limits: negligible for sampling) | none |
| D1-21 | stage 3 / channel path | receiver.html (carrier note) | A carrier at a channel edge does not break the chain (test above). The docs' "unrepresentative special case" could say what happens there (removed by the edge gap) and point to this result. Untested: a narrow noise-like signal straddling an edge | docs text; optional test |
| D1-22 | main.m | main.m:15 (header) | The pipeline list still says "clock jitter, pulsar-Earth motion, polarisation?, 3x3 array element signals, multiple pulsars [todo]" under the receiver. The clock and the motion exist since block D (makeClock, makeTrajectory/delayModel/makePredictor); jitter is out of scope (D1-20) | minor text fix |
| D1-29 | channel path | noiseCal default | A narrow noise-like signal straddling a channel edge gives two lopsided channels with a high weight (the taper keeps their level low): χ² 1.17 at −11 dB, error bars ~0.5 % small. `noiseCal` fixes it (χ² 1.00); it is off by default. Decision for Jasper: noiseCal on by default (cost: an rfiProfile pass per observation). **Jasper (26 Oct): decide on day 5** | **decision, day 5** |
| D1-23 | full band | applyInverseDispersion, info.bulkDelayRef | TOAs refer to 1.6 GHz; τ(fRef) = 40.5 ms at DM 25 stays in. In the simulation truth and observer share it. For real data the DM's error and variations (2.1 µs per 1e-3 at 1.4 GHz) and the Doppler-shifted dispersion (±7 µs) enter through this term: notes item 13 / E2 | for E2/F |
| D1-24 | full band | fullySupported | 31.5 ms (15.7 % of a 0.2 s file) at the end is unusable. Real data are a continuous stream: the processing must carry the overlap across file boundaries (overlap-save across files), not drop the edges | for F |
| D1-25 | full band | memory | One FFT of 2^26 (3.8 GB) at DM 25; DM 100 needs 2^28 (~15 GB): notes item 5. The channel path avoids it | known |
| D1-26 | full band | detectPower 'NChan' > 1 | The sub-band mode (FFT frames) predates channelizeIQ. It looks unused by main and the tests: check, then keep or retire | cleanup check |
| D1-27 | full band | processIQ B_noise | The observer's B_noise (391.6 MHz) comes from its own taper only. The simulated pulsar also passed the generator's taper, so in-pulse the true value is ~0.5 % lower: high-SNR only, irrelevant for the weak signal | none |
| D1-28 | full band | no excision | Every RFI type passes the full band unblanked (D1-4): it is a reference path, not for RFI data | known |
| D2-1 | 5a filterbank | channelizeIQ.m:129–132 | β depends on the LO: any fLO on the 3.125 MHz grid gives β = ½ (1.3 or 1.4 GHz both do, D1-9). A fractional β (an LO off the grid) is supported by the code but never run or tested (R-C1) | test gap |
| D2-2 | 5a filterbank | header "tile the band without gaps" | True for the filterbank. But dedisperseChannels' edge taper (sin², 62.5 kHz at each channel edge) then removes ~2.5 % of each channel's W² bandwidth: the known channel gap (item 3; the edge-carrier test of 26 Oct). The header could say so | docs |
| D2-3 | 5a filterbank | channelizeIQ.m:311–330 | `kaiserLowpass` is copied from applyIQmodulation (on purpose: "so that validated stage stays untouched"). Two copies to keep in step | cleanup |
| D2-4 | 5a filterbank | output files | One file per channel: 128 files, 0.85 GB per 0.2 s. Fine for simulation; F streams in memory | for F |
| D2-5 | 5a filterbank | channelizeIQ.m:221–224 | The 4096-tap weighted sum and the fold are accumulated in single precision; the tests show 2.7e-7 relative error, far below the noise | none |
| D2-6 | drivers | runScenario.m:112, runShip.m:86 | Scripts share MATLAB's base workspace. runShip leaves `psr` as a 1×2 struct array; runScenario then fails at `psr.B0329 = struct(...)` ("Scalar structure required", Jasper's run 27 Oct). The same kind of leftover as the experiments' `restoreDir`. Fix options: `clear psr` at runScenario's top, or local names; for now run after `clear` (and cd back to the repo folder) | driver robustness (day 5) |
| D2-7 | excision | main.m reference run | main's −5 dB reference run sits exactly at the excision's validity limit: its 16-sample windows flag 109 against 13.3 expected (8×), on the pulse; SNR 92.4 → 92.0 (−0.4 %) with excision. Harmless for the weak target, but main's excision numbers are the pulse, not RFI | understanding; docs |
| D2-8 | excision | baseline median/ln 2 | Assumes exponential power. In a channel dominated by a constant-envelope signal (carrier, GNSS) the level is overestimated (up to 1/ln 2), so it flags less. Intended: such channels are notched by the optimal weighting | none |
| D2-9 | excision | detectRFI, periodicRFI, blankingWeights | Whole-file stages: one channel in memory (53 MB per 0.2 s; 267 MB per s). Real-time F needs streaming versions (notes §9) | for F |
| D2-10 | excision | periodicMask false | periodicRFI is off (the false-emitter bug, design parked 9 Oct), so weak periodic RFI (radar sidelobes at −35 dB: ~40 % of pulses caught) is only partly blanked. Doppler decorrelates topocentric RFI over long folds (block D), but that is unquantified (excision#locked) | parked |
| D2-11 | excision | per-channel thresholds | In the 26 Oct straddling test a 200 kHz carrier at −11 dB inside a channel still had 23 flagged windows against 2.4 expected (r 10.3). Residual data-dependent blanking at higher INR; small at 0.2 s | watch |
| D2-12 | excision | frequency guard ±7 | Chosen for the pulsar-locked radar (1 µs bias only below −51 dB); costs ~0.004 % of the data per guard channel. Tuned on one radar model | none |
| D2-13 | 5b dedispersion | full-band reference path | With MaxMemoryGB 4 the full-band path stops at DM ≈ 33 (FFT 2^26 = 3.76 GB at DM 25; 2^27 above DM 33). B0329+54 (26.8) fits; B1937+21 (71) or the Crab (56.8) need a larger budget for the full-band cross-check. Channelized: 29 MB at DM 25, 117 MB at DM 100 | understanding |
| D2-14 | 5b dedispersion | dedisperseChannels.m:102 | One common valid range for all channels: the bottom channel's 31.5 ms sweep (DM 25) cuts the end of every channel. The upper channels discard valid data, on average ~18 ms (~9 % of a 0.2 s file). Negligible in long runs; in streaming (F) the sweep is latency, not loss | none |
| D2-15 | 5b dedispersion | chirp DM | A DM error enters here. The chirp uses the ephemeris DM; a true DM off by δDM = 1e-3 shifts the TOA by 2.2 µs (band-mean of K·δDM/f²; notes §9) and smears the summed pulse by 1.3 µs across the band. The simulation generates and processes with the same DM, so no test can see it | block D (TOA corrections) |
| D2-16 | 5b dedispersion | time-bin noise | Power bins in the channelized path are correlated (lag 1 > 0, variance below the radiometer value; test 4). They nearly cancel inside a phase bin, so the radiometer noise model with B_noise 391.6 MHz holds for the fold. It would not hold for a much finer bin rate (bins ≲ 1/B_ch) | understanding |
| D2-17 | tests | noiseBandwidth.m | No test of its own. `processIQ`, `foldToTOA` and `dedisperseChannels` take B_noise from it, so an error would shift every σ and SNR prediction. Check against the exact integral of the tapers | test gap |
| D2-18 | tests | detectPower | Only indirect checks: testDetectChannels (rows = detectPower per channel) and testDedisperseChannels (end to end). No direct check of the detected power's mean and variance against their prediction | test gap |
| D2-19 | 5 dedispersion | applyInverseDispersion.m:143 | The default guard is 2001 samples, not 2000: ceil(20/8 MHz × 800 MHz) evaluates ceil(2000.0000000000002). Harmless (the tests use the same formula); round(…) or a small tolerance would give 2000 | none |
| D2-20 | 5 dedispersion | applyInverseDispersion.m:231, 241 | Two stale `%#ok<NASGU>` suppressions on the onCleanup lines (checkcode: "the message is no longer generated"); cosmetic, in HEAD before today | cleanup |
