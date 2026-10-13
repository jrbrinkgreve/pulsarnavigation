# 2026-10-13 – full read-through: things that look weird

Jasper (13 Oct, after the code review was committed): "go over the full project thoroughly
flagging anything weird" (he is going over the theory himself in parallel). Read-only:
nothing changed in code, notes or docs; no MATLAB runs. Items 1–33 of
`2026-10-13_code-review.md` are not repeated here.

Read in full: `main.m`, `pipelineParams.m`, every function in `functions/` (pipeline order,
plus the plot/check functions), the five run scripts, `tests/runAllTests.m`,
`testParams.m`, `makeTestData.m`, and the pass criteria of every test; notes §1–4 and the
docs' parameters, pipeline and scripts pages; a scan of all names cited in docs/notes
against the files; the `data/` folder (sizes, dates, which files nothing writes any more).

## Findings

### Can give wrong or misleading results without a warning

1. **`checkConsistency` (`main.m:285`) misses several sky/receiver parameters.** It compares
   T, f_in, L, DM, band, SNR, receiver input, RFI, fLO and fs with the reused files, but
   not `refFreq` (vs `info_disp.refFreq`), `dutycycle` / `genEnvMode` / `seed` / `A`
   (vs `info_gen`), `noiseSeed` (vs `info_rx.seed`) or `filterOrder` (vs
   `info_IQ.filterLen`). Examples with `runStage.sky = false`:
   - `refFreq` 1.6 → 1.4 GHz: the data were dispersed to 1.6 GHz and are dedispersed to
     1.4 GHz, so every TOA shifts by τ(1.6) − τ(1.4) = −2.48 ms. That is just under
     `validateTOA`'s T/4 warning, so nothing warns.
   - `dutycycle` changed: the template no longer matches the pulses on disk.
   - `noiseSeed` changed with `runStage.receiver = false`: silently not the seed-43
     realization.

   Same gap in the scripts. `runRFITest` and `runLockedRadar` do not check the sky files at
   all. `runMonteCarlo` checks T, L, f_in, seed, DM and band, but not refFreq or the duty
   cycle.
   Fix: a few more comparisons. Every field needed is already in the info structs.

2. **`runMonteCarlo`, `runSNRSweep` and `runH0` silently ignore `frontEnd`, `excision` and
   `weighting`.** They always run the full band and print nothing about it. They also pass
   `rfi` straight in (`runMonteCarlo.m:55`, `runH0.m:57`, `runSNRSweep.m:100,135`), so with
   `rfiOn = true` they run RFI through the full band with no excision.

   Their headers still say seed 1 / realization 1 "must match main.m" (`runMonteCarlo.m:10`,
   `runSNRSweep.m:28`). That is true only with `frontEnd = 'fullband'`: main's default
   channel path gives 8 TOAs at 2.956 µs, the full band 9 at 2.880 µs.
   `docs/scripts.html` says they are full-band only; the scripts do not. (Item 22 covers the
   missing channel-path validation; this is about the silent trap.)

### Modelling oddities (small effects, worth knowing)

3. **The band-edge taper shapes the pulsar but not the receiver noise.**
   `applyDispersionStream` applies the 8 MHz sin² edges W(f) to the sky signal. The receiver
   noise is added afterwards and is white. Notes §5.2 (line 287) says "W therefore models
   the receiver bandpass, not the sky". A receiver bandpass would shape the noise too, so the
   simulation treats it as a sky property after all.

   On the channel path, the fraction of the pulsar each edge channel keeps (mean W²) is:

   | Channels | Pulsar kept |
   |---|---|
   | 1 and 128 | 2.4 % |
   | 2 and 127 | 42 % |
   | 3 and 126 | 96 % |

   Their noise is normal. The 'optimal' fit assumes `ChannelGain` = 1 everywhere, so it
   weights these channels fully for a pulse that is almost absent. This costs about 0.5 % of
   SNR. It does not bias the TOA, because the pulse is symmetric. Options: noise shaped by the
   same W, or `ChannelGain` from the observer's bandpass.

4. **Nothing band-limits the receiver noise or the 'impulse' bursts before the IQ stage.**
   Both are white over 0–2 GHz. The IQ low-pass keeps ±360 MHz around 1.3 GHz, so
   `test_rx_IQ.dat` holds receiver noise from 0.94 to 1.66 GHz: 1.8× the in-band noise power.
   This is harmless today, because full-band dedispersion and the channelizer keep only
   1.2–1.6 GHz. But anything that uses the total IQ power, or an out-of-band RFI source,
   sees a receiver without an RF filter.

5. **The channel path loses the first turn of the file, which explains "8 TOAs vs 9".**
   The per-channel dedispersion guard is 20/edge width. Per channel the edge is
   0.02 × 3.125 MHz = 62.5 kHz, so the guard is 320 µs (1334 channel samples); the full
   band's is 2.5 µs. The channel path's first fully supported bin is therefore at 0.32 ms
   (power file: bin 338), against 3.5 µs for the full band (bin 4).

   The synthetic file starts exactly at phase −0.5 (`TRef` = T/2), so turn 0 is then
   incomplete and rejected (`MinCoverage` 1). This is why the docs show 8 TOAs for the
   channel path and 9 for the full band, without saying why.

   It is mostly an artefact of the file starting on a turn boundary; on long or real files
   it costs about 0.64 ms per file. Options: accept it, or a shorter per-channel
   `GuardTime` (the leakage check would show how short is safe).

### Stale or confusing text

6. **`docs/parameters.html:16` says "Values as committed in `9b959f0`".** The page content
   is current, but `9b959f0` (B5c) never touched `pipelineParams.m`; its latest change is
   `0499da7`.
7. **Notes §3 and §4 are stale.**
   - The pipeline map (lines 119, 122) gives the IQ and dedispersed files as "complex,
     500 MHz"; they are 800 MHz.
   - The "Reference configuration" table of §4 (L 1 s, fLO 1.4 GHz, fs 500 MHz, −20 dB,
     nBin 1024, 10 turns per sub-int) is the 30 Sep plan. Every validation since uses the
     seed-43 set, now frozen as `tests/testParams.m`, and §11 has a third "reference
     scenario" (−54 dB). Three different things are called "reference".
8. **`main.m` on the default channel path:**
   - `plots.detected` / `plots.fold` default to true, but these plots exist only for the
     full band. Every default run therefore prints "the iq / detected / fold check plots
     are for the full-band path only".
   - The best-achievable TOA error (`expectedPowerModel`) is printed only for the full
     band, so the default run has no optimum to compare against.

### Housekeeping

9. **Files in `data/` that nothing writes any more (~0.8 GB):**
   - `data/test_dispersed_IQ_modulated.dat` (381 MB, 30 Sep, an old file name);
   - `data/chan/test_dedispFB_chan_*` (407 MB, from the old tests' 12:31 run;
     `testDedisperseChannels` now writes it to `data/test/chan/`);
   - `data/chan/test_dedisp_chan_*`, also from the old tests. `main` writes these only with
     excision off, and `checkProcessing` would flag them as older than their input.

   `data/mc/` holds 14 GB of intermediate files (MC, sweep and H0 rx / IQ / dedispersed).
   Delete or keep: Jasper's choice.

## Checked and fine

No new correctness bug in the processing code. These were checked line by line against the
physics:

- **Dispersion:** chirp phase and group-delay sign, s0 bookkeeping, and output power
  2(B − 1.25·edge)/f_in. The inverse filter is the exact conjugate; the overlap-save index
  ranges and leakage check are right.
- **IQ:** baseband = RF − fLO with no inversion; the group delay is removed on a global
  decimation grid; no image noise aliases into ±360 MHz.
- **PFB:** the polyphase form re-derived from the per-channel mix-filter-decimate formula
  (β = ½ modulation, phase from the absolute index); the alias-free useful band was
  checked.
- **`detectRFI`:** eigenvalue (phase-type) tail of correlated exponentials, threshold
  bracket, median/ln 2 baseline.
- **`powerCovariance`:** complex-Gaussian |R|² form; V + 2ΣX → 1.
- **`blankingWeights`:** mean with the channelizer spectrum (2 Re over d > 0); |C|²
  covariances.
- **Fold:** the noise terms of the linear split, lag pairs across chunks and sub-int
  boundaries, the stable root for `tRef`.
- **TOA and detection:** FFTFIT sign and error propagation (C′ and C as linear forms of the
  bins), the optimal fit's sandwich covariance, and Rice's threshold with λ₂ from the
  template harmonics.
- **`expectedPowerModel`:** the factor 4 from the IQ gain of 2.
- **Radiometer prediction:** 388·ρ per pulse.
- **Text checks:** all function names cited in docs and notes exist (or are local
  functions or variables).

Already known and not repeated here: review items 10b–25 (e.g. item 20e: three copies of the
pulse envelope G(t), two with the old neighbour rule, giving identical results), the
truncation from data-dependent blanking (notes §9), R2 / R3, and excision options not covered
by the reuse checks.

## Items 1 and 2 done (Jasper: "do items 1 and 2")

Implemented as simple logic changes under the 13 Oct rule (implement, then validate each).
Validation: Claude's checks below, then **Jasper's run of `main.m` (13 Oct): no warnings**, so
item 1 in `main` is validated on the real `data/`. The script parts (stops in `runMonteCarlo`,
`runRFITest`, `runLockedRadar`; the warnings of item 2) rest on Claude's checks only.

**Item 1: reuse checks cover every sky and receiver parameter.**

- **New `functions/checkSkyFiles.m`.** It returns the sky parameters in which `info_gen` /
  `info_disp` differ from the current values ({} when all match):
  - generator: T, f_in, L, A, dutycycle, envelope mode, seed;
  - dispersion: DM, band, refFreq;
  - the dispersed file made from this generator file.
- **`main.m`.** `checkConsistency` uses it for the sky part, then checks the receiver part.
  New there: `noiseSeed` (only when the receiver file has noise; RFI needs noise anyway) and
  `filterOrder` (as the filter length `applyIQmodulation` makes: order + 1, made odd). The
  call passes the extra parameters, grouped as sky and receiver.
- **The run scripts.** `runMonteCarlo` replaces its own six-field check with it.
  `runRFITest` and `runLockedRadar`, which had no check, now stop with `<script>:stale` and
  the list of differences.

**Item 2: the full-band-only scripts say so.**

- **Warnings.** `runMonteCarlo`, `runSNRSweep` and `runH0` warn with `<script>:fullBandOnly`
  when `frontEnd` is not 'fullband' (excision and weighting are ignored). They warn with
  `<script>:noExcision` when `rfi` is not empty, because the full band has no excision.
  `runSNRSweep` warns only after its `replotOnly` early return.
- **Headers.** Each now says "full-band path only". The built-in "matches main" checks now
  say "with frontEnd = 'fullband'".

**Validation.** MATLAB R2026a, batch. Everything was read-only on `data/`, using scratch
copies in the session scratchpad.

- **`checkcode`.** It is clean on the new and changed files, except two messages that were
  already there at HEAD (checked against HEAD copies):
  - `runLockedRadar`: preallocation hint;
  - `runSNRSweep`: "cannot be reached" after `if replotOnly ... return`.
- **`checkSkyFiles` on the real info files.** The baseline gives {}. Each of the 11 inputs
  changed alone gives exactly its own name (fLow and fHigh give 'band'). Another
  dispersion input gives 'dispersion input'. A string vs char path gives {}.
- **`main.m`.** Scratch copies had every `runStage` false, so they only load infos and run
  the checks; one change was inserted before the `checkConsistency` call. `main:stale` was
  made an error so that the later `checkProcessing` warning could not mask it in `lastwarn`.
  A first try that read `lastwarn` showed that masking for the refFreq cases. All 16 cases
  are as expected:

  | Change | Warning |
  |---|---|
  | none (baseline) | none |
  | refFreq 1.4 GHz | refFreq |
  | dutycycle 10 | dutycycle |
  | genEnvMode 'amplitude' | genEnvMode |
  | seed 7 | seed |
  | A 2 | A |
  | noiseSeed 44 | noiseSeed |
  | filterOrder 1024 | filterOrder |
  | filterOrder 2047 (same 2049 taps) | none |
  | T, DM, fLO, snrDB (the old checks) | each its own name |
  | noise-free receiver info, other noise seed | none |
  | another dispersion input | dispersion input |
  | refFreq and noiseSeed together | refFreq, noiseSeed |

- **`runMonteCarlo`, `runRFITest`, `runLockedRadar`.** Scratch copies stopped right after
  the check. The baseline reaches the stop; refFreq 1.4 GHz and dutycycle 10 each stop with
  `<script>:stale` "(refFreq)" / "(dutycycle)".
- **`runMonteCarlo`, `runSNRSweep`, `runH0` warnings.** Scratch copies stopped right after
  the warnings, before any file is written. Each warning was turned into an error in turn,
  to see whether it is issued. All 12 cases are as expected:

  | Settings | fullBandOnly | noExcision |
  |---|---|---|
  | channels, no RFI (the default) | yes | no |
  | fullband, no RFI | no | no |
  | fullband, RFI on | no | yes |
  | channels, RFI on | yes | yes |

**Notes and docs updated.**

- **Notes:**
  - header;
  - §3 run control and helpers;
  - §5.11, §5.12, §5.14;
  - §10 item 12;
  - §12 (noise seed);
  - §13 (new row `checkSkyFiles.m`; rows `main.m`, `runMonteCarlo`, `runSNRSweep`, `runH0`,
    `runRFITest`, `runLockedRadar`).
- **Docs:**
  - conventions: the checked parameters, the info-chain diagram label (rendered: fits, clear
    of the arrow), the noise-seed sentence;
  - scripts: the front-end callout, the `runMonteCarlo` checks, the `runSNRSweep` built-in
    check, sky checks for `runRFITest` and `runLockedRadar`;
  - maintenance: a changelog row.
- `api.html` is to be regenerated after the commit (it reads git HEAD).

Committed on Jasper's request after his `main` run (`8a756ec`); `api.html` regenerated from
that commit (`d4e6344`).

## Item 5 questioned, then checked; quick fixes 6, 7, 8 (part)

**Item 5.** Jasper: "I don't think this is true; the final turn isn't fully completed +
dispersion, so it's flagged by FullySupported?" Checked on the saved folds (`data/test_fold.mat`,
`data/test_fold_chan.mat`): phase coverage per turn

| Turn | Full band | Channel path |
|---|---|---|
| 0 | 1.0000 | 0.9683 |
| 1–8 | 1.0000 | 1.0000 |
| 9 | 0.3701 | 0.3379 |

The final turn is flagged on both paths, as Jasper says: the full band is supported to
93.7 ms and the channel path to 93.4 ms, because the last 6.3 ms lack part of the sweep.
That is why both have fewer than 10. The one-TOA difference is turn 0:

- the channel path's support starts at 0.324 ms (power bin 338), the full band's at 3.5 µs
  (bin 4);
- the per-channel guard is 20/62.5 kHz = 320 µs (1334 channel samples, `nPast` in
  `test_dedispB_chan_info.mat`);
- the file starts exactly on a turn boundary.

So the finding stands, with the end turn added to the explanation.

**Quick fixes (Jasper: "further continue the quick fixes").** All text only.

- **Item 6.** `docs/parameters.html`: every value checked against `pipelineParams.m`
  (pulsar, simulation, receiver, noise/RFI, processing and its defaults, ephemeris, files,
  derived quantities, run control, overrides). The line now says "Values checked against
  `pipelineParams.m` at commit `8a756ec`". The docs checklist in `maintenance.html` now
  names that line, and a discrepancy row was added.
- **Item 7.** Notes §3: IQ and full-band dedispersed files "complex, 800 MHz" (were 500).
  Notes §4 rewritten:
  - the reference is now the seed-43 set (= `testParams.m` = `pipelineParams.m` on 13 Oct),
    as a full table with the dated changes;
  - the printed `main` results added;
  - the 30 Sep plan kept as one paragraph;
  - a pointer that §11's "reference scenario" is a different thing.
- **Item 8, part.** `main.m` comments now say which `plots.*` switches are full band only
  (iq, detected, fold) and which work on both paths (toa). The code is byte-identical apart
  from the comments; `checkcode` is clean. Defaults unchanged.
- **Item 8, other part, not done.** A best-achievable TOA error on the channel path would
  need a per-channel version of `expectedPowerModel` (per-channel tapers, optimal channel
  weighting). That is a new model, not a quick fix.
- **Item 5's fact documented.** Notes §4 (the L row) and the derived quantities on the
  parameters page now explain 9 vs 8 complete turns.

## Next

**Item 9 done** (Jasper: "commit this and delete the orphan files").

- **Moved to the Trash** (`~/.Trash/PulsarSimMatlab-orphans-2026-10-13/`, restorable until
  Jasper empties it), 789 MB:
  - `data/test_dispersed_IQ_modulated.dat` + info (30 Sep; no code reads it);
  - `data/chan/test_dedispFB_chan_*` (128 channel files + info, from the old tests' 12:31
    run; `testDedisperseChannels` now writes this set to `data/test/chan/`).

  131 files in all.
- **Kept, correcting the finding above:** `data/chan/test_dedisp_chan_*` is not an orphan.
  It is the source `tests/makeTestData.m` copies when `data/test/` is re-frozen, and
  `tests/makeFoldReference.m` reads its power file. `main` also writes it with excision off.
  Also kept: `data/mc/` (intermediates of the experiments).
- **Text fixed on the way.** The parameters page said the unblanked channel files are "the
  reference for the unit tests". Since 13 Oct they are the source of the re-freeze; the
  tests read `data/test/`.

The quick fixes and item 9 are committed together. Still open for Jasper to choose:

- 3, the edge-channel taper;
- 4, the out-of-band receiver noise;
- 5, whether to shorten the per-channel guard;
- 8, a channel-path optimum (needs a per-channel `expectedPowerModel`).

Items 3–5 are modelling choices and should be discussed before any change (review week:
note, don't implement).
