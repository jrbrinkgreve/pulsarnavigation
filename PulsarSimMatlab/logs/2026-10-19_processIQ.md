# 2026-10-19 – processIQ: design and roadmap

Continues `2026-10-16_review-week.md` (last section: the processIQ proposal).

## Choices (Jasper)

"Start at IQ for the boundary, and include both front ends."

## Read-only inventory of main's processing block (`main.m` 128–287)

- From `info_IQ`: `file`, `actualFsOut`, `fLO`, `t0` (receiver settings only).
- rx: `fLow`, `fHigh`. ephem (already a struct): f0, F1, TRef, DM, profileFWHM.
- proc: frontEnd, refFreq, f_out, nBin, subintPeriods, chanWidth, weighting, excision,
  excisionArgs, periodicMask, periodicArgs (chanWidth and refFreq in proc: filterbank and
  dedispersion/TOA reference choices).
- files: full band dedisp, envelope, fold; channels chanIQ, chanBlank, chanDedisp, chanPower,
  chanDedispB, chanPowerB, chanWeight, foldChan (an input: the run scripts use their own
  folders, and `makeTestData` copies `data/chan/test_dedisp_chan_*`).
- Stays in main (needs truth): bright-excision warning (snrDB), checkConsistency, check
  plots, validateTOA, expectedPowerModel. Moves in: checkProcessing, checkFold (observer only).
- Run scripts: runMonteCarlo / runSNRSweep / runH0 full band; runRFITest channelizes once and
  processes twice (with processIQ: +~4 s per case); runLockedRadar processes difference
  channels by superposition → stays as it is.

Proposed signature: `[toa, detection, out] = processIQ(info_IQ, rx, ephem, proc, files,
'Run', run, 'Verbose', v)`; `toa.good` inside (the good-TOA rule is method); `out` = the
intermediate infos for the plots and validation. Open (asked, not answered yet): split off
the fold → TOA part as its own callable piece (C's `simulateFold` output would use it).

## Where it sits in the roadmap (Jasper asked)

It is the open item "Factor the chain into one function" (§10 item 7), redesigned; it closes
A3b (A's leftover) and comes before C, whose validation plan needs repeatable multi-seed runs
of the default path. Not yet recorded before today: real IQ recordings and the FPGA/GPU port.

**Roadmap updated** (Jasper: "go ahead with the proposed roadmap"; text only):
`docs/roadmap.html` ("In short", A row, E row, open item → "The processing as one function:
processIQ" with the design and steps, block 4 line); notes §10 item 7 (redesign paragraph),
§11 block 4 (real IQ, FPGA/GPU port).

Next: step 1 (`functions/processIQ.m` + `main` + rx/proc/files structs in `pipelineParams`),
shown to Jasper before any run; validation = main bit-identical for both front ends.

## Step 1 written (Jasper: "go ahead and combine"; not yet run)

Asked what processIQ contains (answered: inputs, body steps 0–3, outputs, local functions,
what stays in main). Jasper: keep the TOA step as a separate function. First version had
processIQ end with the TOAs and main's runStage.process / fold / toa passed as `'Run'`;
Jasper asked why `runStage.fold` was still in main → simplified (Jasper: "include the
folding run with the main process_iq switch. Do keep runTOA as a separate function"; keep
the intermediate steps for checks before folding):
- `functions/processIQ.m` (new): `out = processIQ(info_IQ, rx, ephem, proc, files, 'Reuse',
  tf, 'Verbose', v)`, IQ → front end → fold; returns fold, info_fold, the intermediate infos
  (info_det, info_dedisp / info_chan, info_dc, info_w), Bnoise and weighting for the TOA
  step. Reuse = load the saved files and fold + checkProcessing / checkFold (moved from
  main; warnings now processIQ:staleProcessing / processIQ:staleFold, text "Set
  runStage.process = true"). Input field checks (processIQ:fields). Stage calls as in main,
  plus 'Verbose' passed through (default true = as before).
- `functions/foldToTOA.m` (new; name kept, says what goes in and out): template,
  estimateTOA, detectPulsar, toa.good, summary line (prefix now "foldToTOA:"); full band
  passes 'Weighting','equal' explicitly (= the default).
- `main.m`: `runStage.fold` removed (refolding alone saved ~2 s of ~45 s); processIQ call
  with `'Reuse', ~runStage.process`; foldToTOA under `runStage.toa`; check plots and
  validation unchanged (fold plot now also on reuse; plotDetectedPower after the fold);
  local functions checkProcessing, checkFold, infoTime removed (−190 lines); header.
- `pipelineParams.m`: structs rx, proc, files at the end (old variables kept).
- Docs: scripts.html (main's steps), pipeline.html (run table, one-switch note, design note
  simulate / process / evaluate), parameters.html (run table), detection.html (good TOAs in
  foldToTOA), conventions.html (where the reuse checks live), roadmap.html (IQ in, fold out);
  notes §4 run control and §10 item 7.

Validation plan (outputs must be bit-identical): (1) channel path: keep the 16 Oct
`data/test_fold_chan.mat`, run main, isequaln on the fold + the printed numbers of 16 Oct;
(2) full band: old main (HEAD) vs new with frontEnd 'fullband'; (3) reuse path
(runStage.receiver/process false): same TOAs, no stale warnings. Who runs it: asked.

## Validation (Claude's batch run at Jasper's request; Jasper: "this one time consider it validated by me")

One MATLAB batch session (driver and the HEAD copy of main in the session scratchpad, not in
the repo; plots off, sky and receiver files reused): old main (HEAD 97fb16d…a8aacf6) vs new
main, both front ends, plus the reuse path (`runStage.process = false`).

| Comparison | fold | toa | detection | data files (md5) | run time |
|---|---|---|---|---|---|
| old vs new, full band | identical | identical | identical | identical (dedispersed IQ, power) | 2.3 / 1.5 s |
| new vs new reuse, full band | identical | identical | identical | identical | 1.5 / 0.1 s |
| old vs new, channels | identical | identical | identical | identical (power, weights, ch 64) | 33.9 / 34.9 s |
| new vs new reuse, channels | identical | identical | identical | identical | 34.9 / 1.6 s |

Numbers: full band 9 TOAs (9 good), median σ 2.8799 µs, SNR 93.36, χ² 0.963, min pChi2
0.072, total −0.2069 ± 0.9587 µs (docs reference 2.880 µs, 93.4); channels 8 TOAs (8 good),
2.9559 µs, SNR 91.86, χ² 1.008, min pChi2 0.014, total +0.0564 ± 0.9826 µs (= 16 Oct).
The 14:09 fold file (Jasper's run) is identical to the new channel fold. No warnings in
any run (no stale warning on the reuse path). **Validated 19 Oct**; committed.

## Step 2 written: runMonteCarlo on processIQ (A3b) – Jasper validates

Plan explained (what the MC tests; fixed pulsar → ratio slightly below 1; 10 seeds ±0.08,
50 seeds ±0.035 on the ratio). Jasper: "write it", he runs it. `runMonteCarlo.m`: per seed
addNoiseAndRFI → applyIQmodulation (as before), then `processIQ` (own files in `data/mc/`,
channel files in `data/mc/chan/`) and `foldToTOA`; front end from pipelineParams (default
channels, excision, 'optimal'); per-seed line adds good TOAs, median χ², min pChi2, blanked
fraction; pooled: validateTOA (all, good), pChi2 (mean, fraction < 0.05, KS D), detection
rate, blanked fraction; expectedPowerModel full band only; warning fullBandOnly removed,
noExcision only for full band with RFI. The fold file is now saved per seed (was SaveFile
false; 110 MB on the channel path, overwritten). Docs `scripts.html` (front-end note, MC
section, timing row, "Planned").

Validation (Jasper's runs): (1) channels, nSeeds 10: seed 43 identical to main (8 TOAs,
2.9559 µs, +0.0564 ± 0.9826 µs); (2) `frontEnd = 'fullband'`: the old result (90 TOAs, ratio
0.950, χ² 0.91, 72 % within 1σ). Then the A3b run: nSeeds 50 on the channel path (~40 min).

### Run 1 (Jasper): channel path, 10 seeds, −5 dB (406 s)

Seed 43 = main (8 TOAs, 2.956 µs, SNR 91.9, χ² 1.0082, total +0.0564 ± 0.9826 µs) ✓.
Pooled 80 TOAs, all good, all detected: mean +0.075 ± 0.331 µs, rms 2.762 / predicted 2.957 µs,
**ratio 0.934** (±0.074; slightly below 1 as expected with the fixed pulsar), χ² 69.9 / 80,
|z| < 1 66.2 %, < 2 96.2 % → TOA error bars honest on the default path (first time over seeds).
Blanked 0.0124 % (0.0119–0.0131 %).
**pChi2 not uniform**: mean 0.291 (0.5 ± 0.032), fraction < 0.05 0.138, KS D 0.391 (> 0.218).
Per-seed median χ²_red 1.0047–1.0084 (mean ~1.0066; expected 1 ± 0.0009 for the mean of 80
sub-ints, σ per sub-int ~0.0082) → on-pulse residual variance ~0.65 % above the noise model on
the channel path. chi2Tail itself is right (testOptimalTOA test 4: uniform on model draws), so
the real chain has ~0.65 % more on-pulse variance than modelled. Candidates (not tested):
pulse-model mismatch (ChannelGain 1 in tapered edge channels, item 3: rough estimate
ε²·f·SNR²/dof ~ 0.003 at SNR 92), lags beyond Lmax 7 (captured 0.9984 → ~0.0016), white V, X
approximation of blanked bins (BW1). Proposed to Jasper before digging: full-band run
(planned) and a 10-seed channel run at −15 dB (pulse-model excess ∝ SNR², noise-model
excess not).

### Investigation of the χ² excess (Jasper: "go ahead and investigate"; Claude's runs, diagnostic)

**(a) On- vs off-pulse χ² on the two channel folds on disk** (seed 43 from main, seed 52 from
run 1; scratch script, model rebuilt from estimateTOA's fit: a_c, b, τ, band-limited template,
diagonal radiometer variance m²/(B·binDt)·W2/W²):
- off-pulse mean z² **0.9999 ± 0.0008** (3.65 M bins) → the noise model is right;
- on-pulse **1.0046 ± 0.0021** (0.54 M bins), rising with pulse height: template level
  0.01–0.3 → 0.997, 0.3–1 → 1.010–1.012 → pulse-related;
- per channel on-pulse: **channel 1 1.193, channel 128 1.167** (σ per channel ~0.023: 8σ, 7σ);
  edge channels 1–6, 123–128 mean 1.032, middle 1.0018 (±0.002, not significant); channels 1
  and 128 alone give ~0.0028 of the 0.0046 excess. Off-pulse edge 1.008, middle 0.999.
→ full-review item 3: the simulator tapers only the pulsar at the band edges (deep analysis
§2: pulse gain g = 0.023 in channels 1 and 128) while the fit assumes `ChannelGain` 1, so the
model puts a pulse where there is almost none. Excess ∝ (1−g)²·SNR_chan² → ∝ SNR².
(Channels 2/127, g 0.417 in the deep analysis, show no clear excess: 0.989 / 1.029.)

**(b) Channel path, 10 seeds at −15 dB** (scratch copy of runMonteCarlo with snrDB = −15;
372 s): 80 TOAs all good and detected, median SNR ~11, σ ~26 µs; ratio **0.971** (±0.077), χ²
77.3/80, |z| < 1 68.8 %, < 2 95.0 %; **pChi2 uniform**: mean 0.537 (0.5 ± 0.032), fraction <
0.05 0.013, KS D 0.090 (< 0.218); mean red. χ² 0.9987. → the excess scales away with SNR²
as item 3 predicts; at the weak-signal target (−50 dB) it is ~10⁻⁹ and irrelevant.
Side note: blanked fraction 0.0049 % at −15 dB vs 0.0124 % at −5 dB → at −5 dB detectRFI
still flags some pulse peaks (the known bright-excision limit; the warning starts above
−5 dB). Not pursued.

Conclusion: channel-path noise model and pChi2 are right; the −5 dB pChi2 non-uniformity
is item 3 (a simulator truth-vs-observer mismatch in two edge channels), already an open
decision for Jasper. No code change.

### Step 2 runs (Jasper: "let's do both"; Claude's batch runs from scratch copies of runMonteCarlo with only the overrides)

**Full band, 10 seeds, −5 dB** (167 s): 90 TOAs, all good and detected; mean −0.0003 ± 0.306 µs,
rms 2.762 / predicted 2.907 µs, **ratio 0.950**, χ² 81.8/90 (0.909), |z| < 1 72.2 %, < 2 95.6 %
= the old full-band Monte Carlo exactly (0.950, 0.91, 72 %) → runMonteCarlo on processIQ /
foldToTOA reproduces the old script. pChi2 uniform (mean 0.477, fraction < 0.05 0.067, KS D
0.059 < 0.206; mean red. χ² 1.0075): no edge channels, no item-3 excess. Best achievable
2.892 µs per sub-int (expectedPowerModel), measured median ~2.90 µs.

**A3b – channel path (default: excision, 'optimal'), 50 seeds, −5 dB** (3264 s; seeds 43–52
identical to Jasper's run 1): **400 TOAs, all good and detected; mean error −0.032 ± 0.148 µs
(unbiased); rms 2.821 / predicted 2.950 µs, ratio 0.956 (±0.034); χ² 366.2/400 (0.915,
P(≥) 0.89); |z| < 1 67.0 %, < 2 96.8 %.** Blanked 0.0126 % (0.0112–0.0150 % per seed).
pChi2 not uniform (mean 0.292, fraction < 0.05 0.205, KS 0.330; mean red. χ² 1.0063) = item 3
(see the investigation above). Both paths sit at the same ratio (0.950 / 0.956) slightly below
1, as expected with the fixed pulsar (its self-noise is not resampled) → **the default channel
path's TOA error bars are honest end to end over seeds; A3b done.**
