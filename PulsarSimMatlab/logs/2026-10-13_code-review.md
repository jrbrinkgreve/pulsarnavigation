# 2026-10-13 – code review of the whole codebase: list and the minor fixes

## State at start

Review week (handover `2026-10-12_handover-review-week.md`). `main` at `e7bb4be`, uncommitted:
yesterday's review log (stage 1), notes §10 item 11, roadmap, a whitespace edit in
`pipelineParams.m` line 25 (Jasper's; left alone). On 12 Oct Jasper asked Claude to go
through the whole codebase on its own (7 driver scripts, 30 functions, 19 test files,
~13.4k lines) and list everything worth changing; 13 Oct: "implement the minor changes,
validate each step", then "also small logic such as 16 and 17, and similarly simple ones".

## The list (33 candidates) and status

Status: **done** = implemented and validated on 13 Oct (Claude's checks below; Jasper's run of
`testRadarEdges`, `testRfiGating`, `testRfiNoise`, `testBlankingWeights`, `testChannelTOA` the
same day: ALL PASSED, same numbers);
**dropped** = measured, not worth it; **open** = for Jasper to choose.

Correctness and scaling
1. **measured 13 Oct, not changed** – `foldProfile` with `DataWeights`: `addTo` (local function)
   adds a whole-fold-size array ~60 times per chunk (~29k bins with weights), so the cost grows
   faster than the data. The review list called this quadratic, "hours at 10 s": **wrong** —
   memory bandwidth makes the full-array add cheap, and the 2 GB fold cap (`MaxMemoryGB`) bounds
   sub-ints × channels (128 channels, 1-turn sub-ints: ~1.4 s at most). Tried: accumulate only
   into the touched sub-ints (`addTo` with a cell range; bit-identical in 8 real cases and on
   synthetic data). Times HEAD → range version: 0.1 s / 128 ch 1.9 → 2.5 s (+35 %, the remap of
   every index costs more than the full add saves); 16 ch 1 s / 2 s 2.5 → 3.3 / 5.2 → 6.3 s;
   128 ch 1 s (100 sub-ints, at the cap) 34.3 → 26.8 s (−22 %; HEAD grows 18× for 10× the data,
   the range version 10.5×). Slower where it is used today, a small gain only at the cap:
   **reverted** (`foldProfile.m` = HEAD). If long channel-path folds with weights matter later
   (with item 32): accumulate per chunk into a local slice of the sub-ints it covers, built with
   local indices (no per-element remap), written back once per chunk (bit-identical if the
   slice starts from the global values).
2. **done 13 Oct, validated (Jasper's `main` run, only `runStage.toa`)** – reused processing files were never
   checked for staleness: `checkConsistency` (`main.m`) covered sky + receiver only. Now two more
   local functions in `main.m`: `checkProcessing` (each reused product made from the file before
   it — `inFile`, `rawFile` for the blanked channels, `powerFile` for the weights — and its info
   file at least as new as its input's; f_out as samples per bin, DM, refFreq, band, chanWidth;
   weight file N and Lmax vs the power file) and `checkFold` (power file, age, nBin,
   subintPeriods, ephemeris, data-weight file). Warnings `main:staleProcessing` /
   `main:staleFold`, only when a stage is reused. Not checked: excision / periodic-mask options
   (not stored in the infos). Validated with a scratch copy of `main.m` (all data paths
   redirected, L = 0.04 s): fresh build and unchanged reuse no warning; f_out 100 kHz →
   "f_out (power file: 4 samples per bin, this run: 42)"; subintPeriods 2 → staleFold; ephemeris
   DM 5.1 → "DM"; receiver rerun with processing reused → "test_rx_IQ_chan older than its
   input"; `checkcode` clean. Live case on the real `data/` (read only): "f_out (power file: 42
   samples per bin, this run: 4)". Before the fix that run went on to fold the 99 kHz data. **Live case found today:** main's channel-path
   products in `data/` (`test_dedispB_chan_power.dat` 9,920 bins, the weight file with
   Lmax 1, `test_fold_chan.mat`; all 9 Oct 10:49) are from the f_out study at ~99 kHz, not
   from today's f_out = 1 MHz. `main` with `runStage.process = false` would fold them
   without a warning. Fix: check the `inFile`/`sourceFile` links and key parameters.
3. **done 13 Oct, validated (Jasper's run: `makeTestData`, then `runAllTests` 14/14, 6.8 min)** – the unit tests read main's data files (`data/test_rx_IQ.dat`, `data/chan/test_*`),
   which `main` rewrites; e.g. after a main run with another SNR or RFI, test 1 of
   `testDedisperseChannels` (bit-identical vs `test_IQ_dedispersed.dat`, no longer rewritten
   by the channel path) fails for no real reason. Fix: frozen test inputs in `data/test/`.
   Done: `tests/testParams.m` (script) = the tests' own frozen settings, the seed-43
   reference (today's `pipelineParams` values), every path under `data/test/`; it errors
   `testParams:noData` when that folder is missing. `tests/makeTestData.m` (run once)
   first checks along the info chain that `data/` is that reference (pulsar, dispersion,
   noise seed 43 at −5 dB without RFI, IQ, full-band dedispersion and power, channels,
   channel dedispersion and power; each file made from the one before), then copies to
   `data/test/`: data + info of the IQ, full-band dedispersed and envelope files, info only of
   the generator / dispersion / receiver files (the tests read their ground truth, not their
   samples), the two channel sets and the channel power file, and the three frozen
   references of `data/mc/`; the path fields in the copied infos point into `data/test/`.
   Refuses if `data/test/` exists (`'Overwrite', true` to re-freeze on purpose). The 12 tests
   that read data run `testParams` instead of `pipelineParams` (one line each, plus header
   paths); `testChannelizeIQ` points at `data/test/`. `runAllTests` skips `testParams`.
   Found on the way: (a) `dutycycle` in a script is parsed as the Signal Processing
   Toolbox function when a function file uses it after `run`ning the script, so
   `makeTestData` reads the settings into a struct (local `readTestParams`); (b) `run()`
   changes into the script's folder, so the no-data check is relative to `testParams.m`
   itself, not to the current folder; (c) `testFoldWeights` compared the fold info's
   `sourceFile` (a path) with the reference: now by file name only (folds bit-identical).
   Validation (Claude's runs in an APFS copy of the repo): `makeTestData` 9 checks passed,
   272 files, 2.23 GB; every path field of the 10 copied infos points into `data/test/`;
   a second call stops with `makeTestData:exists`; with a changed source info (SNR) and
   `'Overwrite', true` it stops with `makeTestData:source` ("receiver noise") and leaves
   `data/test/` untouched (names, sizes, dates); without `data/test/` a test stops with
   `testParams:noData`. `runAllTests` 13/14 (only the `sourceFile` path in
   `testFoldWeights`), after the fix **14/14 with everything except `data/test/` moved out
   of `data/` and `pipelineParams` changed** (−20 dB, RFI on, fs 500 MHz, nBin 1024): 9.1 min.
   In the real repo: Jasper ran `makeTestData` and then `runAllTests`: 14/14 in 6.8 min
   (validates item 3, the 20 % generator test and the > 20 % warning).

Inaccurate text
4. **done** – `blankingWeights` help: true-spectrum correction "~4 %" → "~1 % (≤ 0.07 % for
   W ≥ 0.1)", as in notes §5.18 / §9.
5. **done** – warning `addNoiseAndRFI:weak` (below −60 dB, "pulse below single-precision
   resolution") removed and replaced by a comment: the noise dithers the float32 rounding,
   so the rounding error averages to zero, is uncorrelated with the pulse and adds white
   noise ~145 dB below the noise power; no floor at −60 dB or anywhere near the target.
6. **done** – outdated statements: `dedisperseChannels` ("Later, RFI excision ..."),
   `detectChannels` ("Blanking, step 2, will ..."), `estimateTOA` ("current noise-free
   data"; now with the phase-D numbers 1.5× / 2.4×), `pipelineParams` header (shared by all
   run scripts and the tests), `main` (`runStage.process` comment: also filterbank,
   excision, blankingWeights).
7. **done** – `rfiSource` help: `StartTime` is a phase (pulses at StartTime + k/PRF for every
   k, also before it); checked against the formula in `addNoiseAndRFI`.
8. **done** – `main` disk estimate per path: the channel path no longer counts the full-band
   dedispersed and power files (7.87 → 7.23 GB at 0.1 s with excision; full band unchanged,
   6.08 GB = the files on disk exactly).
9. **done** – `estimateTOA` local function `fitOne` re-indented (it was one level too deep).

Performance
10. (a) **done** – `'noise'` RFI: the filter spectrum `fft(h)` of a full block computed once per
    source instead of every block; bit-identical; stage −13 % with 3 of 7 sources 'noise';
    costs 67–134 MB per source (~1.7 GB in `runRFITest` case 15).
    (b) **open** – generate 'noise' RFI at complex baseband and upconvert (~100× less work;
    new realizations → B5d/B5e reruns).
11. **open** – `applyIQmodulation`: LO `exp` recomputed every block; fLO/f_in = 13/40, the LO
    repeats every 40 samples. Measure first (stage 46 s per second of data).
12. **open** – `applyIQmodulation`: filters at 4 GHz then keeps every 5th sample; decimation in
    the frequency domain (Nfft = 5·2^k) makes the inverse FFT 5× shorter. Validated stage.
13. **open** – `applyDispersionStream` FFT size: at L = 1 s 2^28 (11.8 GB) for ~5 % fewer flops
    than 2^27 (5.9 GB). Measure first.
14. **open** – `periodicRFI` period search O(n³) (each candidate period tries every pulse as
    the phase); sort + sliding window gives the same inliers in O(n log n). Low while the
    periodic mask is off.
15. **open** – per-channel loops (dedisperse, detect, blank, detectRFI, blankingWeights) with
    `parfor`; Parallel Computing Toolbox licensed (checked 13 Oct). Measure first.
16. **dropped** – larger `channelizeIQ` `BlockOut`: measured on the 0.1 s file, 1024 / 2048 /
    4096 / 8192 → 3.2 / 3.2 / 3.1 / 3.1 s (first run 3.6 s: warm-up), outputs byte-identical
    for all block sizes and identical to main's `data/chan` files. ~3 %: default kept.
17. **done** – generator neighbour rule: `max(1, ceil(6σ/T + 0.5))` → `floor(11.75σ/T + 0.5)`
    (every Gaussian above 1e-30 of the peak; 0 at 5 % duty). Bit-identical (5 % with A = 1
    and 1e20, 'amplitude' mode, 20 %, 30 %, 80 %); ~10 % faster at 5 % (RNG dominates).
18. **open** – per-channel dedispersion kernel also carries the delay to refFreq (126 ms at
    DM 100 for the bottom channel vs 1.5 ms own smear): integer shift + residual chirp.
19. **open** – `runH0` / `runSNRSweep` write a 1.6 GB zero file for an empty sky; an option
    for no input file in `addNoiseAndRFI` would avoid it.

Structure (duplicated code)
20. **open** – shared helpers instead of copies: (a) dispersion constant / τ(f) (6 files);
    (b) sin² band taper (5 functions); (c) `kaiserLowpass` (2 copies + `noiseBandpass`);
    (d) mask normalization (3 copies); (e) pulse envelope G(t) (3 copies); (f) cf32
    read/write (9 copies); (g) profile noise model (`estimateTOA` / `detectPulsar`).
21. **done** – `blankingWeights` temporary impulse files: unique names (`tempname`) instead of
    names from the output file; weight file bit-identical, no files left in tempdir.

Validation gaps
22. **open** – statistical validation (`runMonteCarlo`, `runH0`, `runSNRSweep`) only on the
    full-band path; main's default is channels + excision (with notes §10 item 7).
23. **open** – ground-truth check plots only for the full band (`plotFoldCheck`).
24. **open** – ephemeris = truth in f0 and DM too (extends R3): runs with deliberate small errors.
25. **open** – template not smoothed by the fold's time / phase binning.

Found while validating item 2
34. **done 13 Oct (Claude's validation; Jasper's run pending)** – `main.m` stopped with an error in
    `validateTOA` ("No valid TOAs to validate") when no sub-int gave a TOA (e.g. the 99 kHz fold:
    0 of 10 fitted, empty phase bins). Now `validateTOA` runs only if any TOA is valid, else
    "main: no sub-int gave a TOA (0 of N fitted; see the warnings above), nothing to validate";
    the channel-path summary prints "channels used none" instead of "zeros(1,0)". Validated:
    scratch set at f_out 100 kHz (0 of 4 fitted) runs to the end with both lines; a normal run
    validates as before (summary text unchanged, `mat2str`); `checkcode` clean. (Jasper reran
    `main` at 12:05 with `runStage.process = true`: `data/` is at 1 MHz again, 8 TOAs, 2.956 µs;
    reused, it gives no stale warning — item 2 without a false alarm on real data.)

Already recorded before the review
26. R1 generator unit test. 27. R2 remove `'amplitude'`. 28. R3 true pulse offset.
29. §10 item 7 one chain function. 30. "Run all tests" entry + unit tests for stages 1–4 and
`detectPower`. 31. `expectedPowerModel` still [written]. 32. Excision on long files.
33. Parked periodicRFI fix and per-channel noise calibration.

30 (the entry point) **done and validated 13 Oct (Jasper's run 13 Oct: the 13 tests and testGeneratePulsarSignal all passed)**: `tests/runAllTests.m`.
The 13 unit tests in one session in dependency order (testChannelizeIQ, testDedisperseChannels,
testDetectChannels write the shared channel files; then FoldWeights, ChannelTOA, OptimalTOA,
BlankingWeights, BlankChannels, DetectRFI, PeriodicRFI, RadarEdges, RfiGating, RfiNoise), each in
try/catch; table PASS / FAIL + time + first error line; error at the end if any failed; subset by
name (still in dependency order); demo*, exp*, make* not run; a new test*.m not in the list runs
last with a note. Validated in an isolated copy of the code and test data (APFS clones in
scratch, so `data/chan` was not touched): all 13 passed in 6.3 min (OptimalTOA 148 s,
BlankingWeights 127 s); a subset runs in dependency order; an unknown name is refused; a
deliberately failing extra test shows FAIL with its message, runs last with the note, and the
run ends with `runAllTests:failed`; `checkcode` clean. Unit tests for stages 1–4 and
`detectPower`: stage 1 done (R1, next paragraph); stages 2–4 and `detectPower` still open.

R1 (26) **done and validated 13 Oct (Jasper's run 13 Oct: the 13 tests and testGeneratePulsarSignal all passed)**: `tests/testGeneratePulsarSignal.m`,
~4 s, own files in tempdir, first in `runAllTests`. (1) format and ground truth (binary-exact T = 1/1024 s, f_in = 2^27 Hz: the 13th pulse centre on the last sample); (2) power profile, 5 %, 200 pulses at 100 MHz: folded mean x² in 1 µs bins vs A²·G from an independent exact periodic sum (21 images), exact per-bin σ: residuals mean +0.11 ± 0.08, std 0.97 ± 0.06; centre 0.500015 ± 0.000017 turns; width −0.3σ (FWHM 49.92 µs, nominal 50); (3) 'amplitude' mode, power FWHM / √2; (4) 20 % duty, the widest real pulses (Jasper: pulsars stay below ~20 %; was 40 %), all bins: the bins halfway between pulses check the neighbour sum (a generator copy without it fails test 4 only: residual std 5.3); (5) same seed identical, BlockSize 333,333 / 7,919 = one block bit for bit, A = 4 exactly 4 × A = 1 for normal floats (subnormals far off-pulse differ by ≤ 3e-45), another seed differs, global RNG untouched, clock seed reproduces; (6) bad mode / seeds refused; 25 % warns, 20 % does not (the warning was > 50 % until Jasper lowered it to > 20 %, 13 Oct). Found on the way: the
notes' open question (does single `randn` drawn in pieces equal one draw?) — yes, checked
directly (single and double, five splits), so the file is block-size independent. First run
failed only on the A-scaling check: for float32 subnormals (|x| < 1.2e-38, far off-pulse)
(4·env)·n and 4·(env·n) round differently (55,505 samples, ≤ 2.8e-45) — a test expectation,
not a generator bug; the check now requires exact scaling for normal floats. `checkcode` clean.

## How each change was validated (Claude's MATLAB runs, R2026a, 13 Oct)

- Text-only edits (4, 6, 7, 9): a script checks against git HEAD that every changed line is a
  comment (or, for 9, the same text with other indentation); `checkcode` messages identical
  to HEAD for every edited file. 9 also: `testChannelTOA` ALL PASSED (estimateTOA and
  detectPulsar regressions bit-identical to the pre-A2 references).
- 5 and 10a: `addNoiseAndRFI` vs the HEAD copy on a zero sky with all RFI types and three
  'noise' sources (13.5 MHz, 2 MHz, 200 kHz), −70 and −5 dB, block sizes 4e6 and 1e6:
  output bit-identical; HEAD warns at −70 dB, the new code does not; `info.rfi` identical
  (the cached spectrum is not saved).
- 8: estimate vs the files on disk (full band exact; channel files exact per type).
- 17: generator vs HEAD copy, byte-identical files and identical info in six cases.
- 21: weight file of 128 channels (four blanked) bit-identical to the HEAD copy.
- Unit tests rerun in one session: `testRadarEdges`, `testRfiGating`, `testRfiNoise`,
  `testBlankingWeights`: ALL PASSED, same numbers as recorded (e.g. blanked TOAs
  −0.251 ± 0.354 µs, noise ratio 0.9873, long-bin test 2.8e-8).

## Docs and notes

Notes: §3 (disk use per path), §5.1 (neighbour rule), §5.3 (StartTime, noise-type cache,
no −60 dB warning), §5.18 (temp files), §8 lesson 11 (test sessions), §10 item 12 (items 2,
30, R1, 34, 1, 3), §13 (runAllTests, testParams, makeTestData, generator test), header.
Docs for item 3: tests (how to run, shared files, data flow figure), conventions §8,
parameters §7, scripts (parallel runs), excision (blankChannels note), roadmap, maintenance.

**Diagrams (Jasper, 13 Oct):** 19 new figures on 14 pages, in the site's SVG style (light and
dark), checked by rendering each one (headless Firefox): flowcharts for run control,
ground-truth rule, info chain, receiver, overlap-save, bin stamps, filterbank, good TOAs,
test data, roadmap; computed plots (generators in `docs/tools/figures/`, numpy) for the
power profile and neighbour sum, dispersion delay, IQ spectra, excision of one channel,
folding, FFTFIT, T0/Tmax distributions (simulated: Rice factor 5.51 and the 4.15 threshold
reproduced, P_FA 1.0·10⁻³), the SNR and threshold sweeps, and the pulsar map. Doc errors
found while drawing, fixed: bin-start stamping makes TOAs *early*, not late (full band);
"others 5000+" in the scenario's figure of merit (B1937+21 is 300, two MSPs ~4000);
receiver options still mentioned the removed −60 dB warning; the sky page still called the
`MaxMemoryGB` help text out of date (fixed 9 Oct); roadmap status line (8 Oct) and
timeline. `.claude/launch.json` added (a local server for previewing the docs). Docs: `sky.html` (neighbour
rule), `receiver.html` (StartTime, noise-type cache), `roadmap.html`, `maintenance.html`
changelog. Not committed (Jasper decides); after the commit run
`python3 docs/tools/make_api.py` (help texts changed in 7 functions and 2 scripts).

## Next

Jasper ran the five tests (all passed) and asked for the medium items one at a time, in the
suggested order: 2 (stale-file check in `main`), 1, 30 + R1, 3. Item 2 done (above); Jasper's run
of `main.m` pending (with `runStage.process = false` it now warns about the 99 kHz products;
`runStage.process = true` remakes them at 1 MHz; done by Jasper 12:05). Item 34 done (above).
Item 1 measured and reverted (above). Items 30 and R1 done (below), validated by Jasper's run;
test 4 of the generator test then moved from 40 % to 20 % duty (Jasper: no pulsar is wider than
~20 %). Item 3 done (frozen test data, above; Claude's validation). Notes and docs updated,
diagrams added. Jasper then ran `makeTestData` and `runAllTests`: 14/14 (6.8 min), validating
item 3 and the generator changes. Jasper's `main` run with only `runStage.toa` (validates
item 2): `checkFold` silent; `checkProcessing` warned "test_rx_IQ_chan_blanked older than its
input". Cause: Jasper's 12:31 run of the old unit tests rewrote `data/chan/test_rx_IQ_chan_*`
after main's 12:05 blanking (before item 3 the tests wrote there). Contents unchanged: all
53,326,734 unblanked samples of the blanked files equal the rewritten raw files bit for bit
(Claude, numpy); 0.0125 % blanked as in B3; 8 good TOAs, median 2.956 µs, as at 12:05. A true
date warning on identical data, and it cannot recur (the tests now write only `data/test/`);
`runStage.process = true` once clears it. Item 34 (no valid TOA) stays Claude's check only.
Jasper then ran `main.m` with all stages, the sky included (13 Oct). Committed on Jasper's
request as `0499da7`; then `api.html` regenerated and the docs' commit line, changelog and
roadmap status pointed to `0499da7` (docs commit). Generator overlap warning
lowered from > 50 % to > 20 % duty (Jasper: yes, 13 Oct): help text, threshold, message;
test 6 now checks 25 % warns and 20 % does not (Claude's run and Jasper's `runAllTests`: passed, tests 1–5 the same
numbers; `checkcode`: test clean, generator only its pre-existing `onCleanup` suppression note).
