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
