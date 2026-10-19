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
