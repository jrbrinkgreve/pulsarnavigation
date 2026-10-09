# 2026-10-09 – handover: block B (RFI excision) closed for now — state, open items, ready designs

Start here when RFI work is picked up again. Day records: `2026-10-08_excision.md` (B1–B6, B5a–B5d),
`2026-10-09_B5e.md` (realistic scenario, today's findings). Full reference:
`current_project_notes.md` §5.3, §5.19–5.21, §7, §9, §10. Readable form: `docs/excision.html`,
`docs/limits.html`, `docs/results.html#b5e`.

## Why block B is closed now

Jasper (9 Oct, after B5e): "are we now not stuck in a very tiny hole about RFI?" Agreed:

- The essential part is done: strong pulsed RFI (radars) is detected and blanked before
  dedispersion; constant-envelope and noise-like RFI (GNSS, LTE, Inmarsat) only add noise, which
  the 'optimal' channel weighting handles; a radar locked to the pulsar period is handled by the
  frequency guard (and the periodic mask, when switched on).
- What remained were second-order effects, measured against RFI levels estimated to ±10 dB, at a
  pulsar 80,000× brighter than the target (−5 vs −54 dB) and over 0.1 s instead of hours. Tuning
  further has little value before real data from the hardware exist.
- The critical path to navigation is C (fast simulator, hours at −54 dB) and D (phase
  prediction: Earth's motion, clock). D also makes a locked radar unlikely: a ground radar keeps
  a fixed rhythm, the pulsar's apparent period drifts with Earth's motion.
- Decision: **periodic mask off by default** (`periodicMask = false`, 9 Oct, Jasper); the two
  open items below recorded with ready designs. **Do not start C yet** (Jasper, 9 Oct): ask first.

## What block B delivers

```
main.m, channel path (default):
channelizeIQ → detectRFI → [periodicRFI if periodicMask] → blankChannels → dedisperseChannels
→ detectChannels → blankingWeights → foldProfile('DataWeights') → estimateTOA / detectPulsar ('optimal')
```

| Step | What | Commit | Validated |
|---|---|---|---|
| B1 | `blankChannels`: mask → blanked copies of the channel files | 45d95bc | Jasper |
| B2 | `detectRFI`: multi-scale power threshold per channel, exact thresholds, guard 10 | 45d95bc | Jasper |
| B3 | excision in `main.m` (`excision = true`, default) | f4181f5 | Jasper |
| B4 | `runRFITest.m`: each RFI type with / without excision | f4181f5 | Jasper |
| B5a | radar `RiseTime` (0.1 µs default) | 724501d | Jasper |
| B5b / B5b-2 | `runLockedRadar.m`; `detectRFI` frequency guard (`FreqGuard` 7) | 724501d | Jasper |
| B5c | rotating antenna for any RFI type | 9b959f0 | Jasper |
| B6 | `periodicRFI` (periodic mask); **default off since 9 Oct** | 2082540 | Jasper (testPeriodicRFI) |
| B5d | `'noise'` RFI type (LTE) | 9750461 | Claude's run; **Jasper's run of `tests/testRfiNoise.m` pending** |
| B5e | realistic scenario (`rfiList`, `rfiRealistic`), `runRFITest` cases 14–15, per-channel noise check, periodic hook | f019a8a | Claude's runs at Jasper's request |
| close | `periodicMask` default false; this handover | (this commit) | – |

Defaults: `frontEnd 'channels'`, `weighting 'optimal'`, `excision true` (detectRFI windows 1–16
samples, PFA 1e-6, guard 10, FreqGuard 7), **`periodicMask false`**, `rfiOn false`,
`rfiList 'illustrative'`, radar `RiseTime` 0.1 µs.

Key results (−5 dB, seed 43; notes §7):
- No RFI: excision costs 0.0125 % of the data, TOAs unchanged.
- Radar: breaks the TOAs without excision; with excision fully restored. Realistic radar at
  +40 dB: all pulses blanked, no energy missed.
- Locked radar (worst case) at −54 dB: excision alone 147 µs bias → frequency guard 7: 2 µs →
  plus the periodic mask: 0.41 µs. Sidelobe-only locked radar: needs the periodic mask
  (1.2e-10 of the baseline with it).
- Realistic environment: −11 % SNR (~25 % longer observing), mostly LTE band 32.
- Lopsided channels: TOA error bars 0.39 % too small (realistic), 3.4 % with narrow
  Inmarsat-like carriers.

## Open item 1 — periodicRFI false emitters with dense detections (ready to fix)

**Symptom.** `runRFITest` case 14 / 15 with `periodicMask = true`: 29 % / 32 % of the data
blanked (SNR 69 / 66 instead of 82 / 80). Emitters reported: 351.000 Hz (radar 2, correct),
14108.5 Hz twice, 14467.6 Hz (false). Reproduce: set `periodicMask = true` in
`pipelineParams`, `runCases = [14 15]` in `runRFITest`.

**Diagnosis (confirmed 9 Oct with an instrumented copy; numbers in `2026-10-09_B5e.md`).**
1. *Pulse times on a grid.* detectRFI's windows step by 1, 1, 2, 4, 8 samples, so event centres
   lie on a 4-sample grid (time mod 8 = 0.5 or 4.5). bestPeriod's chance q = (2·TimeTol + 1)/P
   assumes continuous times; for P a multiple of 4 the true chance is (2·floor(TimeTol/g) + 1)·g/P
   = 12/P for g = 4. Inmarsat group: P = 288 = 72 × 4, best phase 138 of 2240 where q predicts 70
   → lp −14.2 (accepted). → the 14.5 kHz emitter.
2. *True period never a candidate.* Candidates = gaps to the next 3 pulses (÷ 1..MaxMissing).
   With a false flag every ~424 samples the next 3 pulses are flags, never the next radar pulse
   (4282 samples). An alias wins: radar 1's period = 14.5 × 295.33 samples, so every second radar
   pulse is on the 295.33 grid (68 of 707, lp −19.9; a real alias). → the 14.1 kHz emitters.
3. *No plausibility check.* Fraction of predicted pulses with a detection: false emitters
   4.7 % / 5.9 %; radar 2 (real) 54 %; validated B6 cases ~40 % (−35 dB sidelobes) to ~100 %.

**Fix design (agreed in outline, not started; one edit each, in this order):**
- **a. Grid-aware chance:** q = max over grids g = 1, 2, 4, …, max(stride) of
  (2·floor(TimeTol/g) + 1)·g / P (default strides → 12/P). Small sensitivity loss (radar 2 stays
  at lp ≈ −90 vs the −13.8 threshold).
- **c. Plausibility:** accept an emitter only if ≥ `MinHitFraction` (proposed 0.2) of its
  predicted pulses within the span of its fitted pulses have a detection. Rejects aliases and
  random rhythms; keeps all validated cases. With a + c the default path is safe.
- **b. Candidates from strong pulses:** also gaps between pulses with peak ≥ `StrongPeak`
  (proposed 30× the baseline; noise ~e^−30 per sample) and their next 3 strong pulses, so a strong
  radar's true period is found among dense flags (needed for B6's purpose: predict a strong
  radar's weak pulses).
- **Tests:** extend `tests/testPeriodicRFI.m` with synthetic events (no files needed:
  periodicRFI only uses info_chan's fs, N, nChan, t0 and info_rfi.events / guard / freqGuard;
  fix a would add info_rfi.stride): dense random flags on a 4-sample grid → no emitter; a radar + dense flags → its true
  PRF. Old tests must still pass. End-to-end: `runRFITest` 14–15 with the mask → ~0.3–0.5 %
  blanked, emitters radar 1 (973 Hz) and radar 2 (351 Hz) only. Then decide the default again.
- Effort ~1 session. Only needed if locked radars matter (after D) or to blank weak sidelobe
  pulses predicted from a main-beam passage.

## Open item 2 — per-channel noise calibration (lopsided channels; design sketch)

**Physics.** detectRFI's thresholds and the fold's noise model (V, X) assume a flat spectrum
within each channel. A channel partly covered by noise-like RFI (an LTE band edge, a narrow
carrier) has fewer independent samples per bin → false flags and noise above the model. Flat-
channel formula (RFI with PSD R × the noise over a fraction f of the channel): level
L = 1 + fR, noise ratio r = (1 − f + f(1 + R)²)/L². χ² counts channels equally; the TOA fit
weights them by F = 1/L², so the error bars are too small by √(Σ F r / Σ F). Per channel
F(r − 1) = f(1 − f)R²/(1 + fR)⁴, largest at R = 1/f: (1 − f)/(16 f) → worst case = narrow
noise-like carriers about as strong as their channel's noise.

**Measured (B5e):** LTE 0 dB 0.08 %; realistic scenario 0.39 % (formula exact); ten narrow
Inmarsat-like carriers χ² 1.22, error bars 3.4 % too small, 280–730 false-flag windows per
carrier channel (vs ~3).

**Design sketch (not agreed; design note first, it changes validated noise models):**
measure each channel's noise on the data (off-pulse bins via the ephemeris, or the variance of
the detected power relative to its mean², per channel and file) → scale V, X (and so the fit
weights and χ²) per channel; optionally detectRFI thresholds from the measured window
variance per channel. Ground-truth rule: fine (measured on the data). Test case: `runRFITest`
case 15 (χ² 1.22 → ~1.0, error-bar factor → ~1, false flags → ~3 per channel). Note: the r
check in `runRFITest` compares with the median channel and is invalid under heavy blanking.

## Other open items in block B (unchanged)

- Staggered-PRF radars in periodicRFI (v1: fixed PRF only); a frequency guard scaled with the
  event's strength; periodicRFI for long data (one row per pulse and channel).
- Spectral kurtosis / whole-channel flags: only if needed.
- **Receiver dynamic range** (not modelled): a +40 dB radar would clip a few-bit ADC and spread
  over the band — possibly as important as excision for real hardware.
- Ionosphere (delay ~3–70 ns at 1.4 GHz, Faraday rotation, scintillation; notes §9).
- The RFI levels of `rfiRealistic` are estimates (±10 dB); a measurement at the real site
  would replace them.

## Pending runs (Jasper)

- `run('tests/testRfiNoise.m')` (B5d, ~5 s).
- Optional: `runRFITest` with `runCases = [1 13 14 15]` (~10 min) to validate B5e yourself
  (now without the periodic mask by default — matches the B5e table in notes §7).

## How to run (RFI)

From `PulsarSimMatlab`: tests `testBlankChannels`, `testDetectRFI` (~33 s, ~2.5 GB in tempdir),
`testRadarEdges`, `testRfiGating`, `testPeriodicRFI`, `testRfiNoise`. Scripts: `runRFITest`
(cases 1–15; `runCases`; ~45 s per case, realistic 3–6 min; `data/rfi`), `runLockedRadar`
(cases 1–12, ~9 min; `data/locked`). Realistic scenario in main: `rfiOn = true`,
`rfiList = 'realistic'`.

## Working tree note (9 Oct)

Other sessions had uncommitted edits when this was written (this morning's consistency fixes
in the notes / `applyDispersionStream.m` help / `docs/maintenance.html`; the dispersed-domain
idea, `logs/2026-10-09_dispersed-domain.md`, `tests/demoBinningMatchedFilter.m`; f_out study
log, `docs/tests.html`). B5e and this close-out committed only their own lines.

## Next

Block B closed. Next block: C (fast simulator; design draft in notes §10 item 3) — **not
started; Jasper decides when** ("do not start part C yet", 9 Oct). Then D, E.
