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
| D1-11 | IQ stage | float32 output | No quantization. A real SDR delivers integer samples (8–16 bit). For noise-dominated data the SNR loss is small at 8 bits and more at 2–4, but a strong RFI line can clip a few-bit ADC. ADC bits and the RFI headroom belong in E4 | for E4/F |
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
