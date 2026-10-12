# 2026-10-12 – handover: the review week (understand and refine the whole pipeline)

Start here in the next session. Today's record: `2026-10-12_overview.md` (week-start overview,
then the fixes of the day). Full reference: `current_project_notes.md`. Readable form: `docs/`.

## Goal of the week (Jasper, 12 Oct)

"I don't fully understand blankingWeights, so let's keep this week to fully understanding and
refining the whole pipeline." Block C (fast simulator) is **not** started; Jasper decides when.
Aim: a codebase Jasper fully understands, one stage at a time (CLAUDE.md).

## State at hand-over

- Git: `main` clean. Five local commits since `849ccf5`: `f40df46`, `817cfe5`, `62e8d06`,
  `e947d4b`, plus the commit of this handover. Jasper pushes himself.
- No known bugs open. Fixed and validated today (Jasper's runs): `blankingWeights` long-bin
  crash (`fft(q, M, 1)`, test 4) and lag-array memory (bit-identical); B5b-2 (`testDetectRFI` 3d)
  and B5e (rerun of `runRFITest` 1/13/14/15, identical) validated; all stale text from the
  12 Oct consistency report fixed; legacy files stay in `old/`.
- Recorded limits, not for this week: excision on long files (`detectRFI`, `periodicRFI`,
  `blankingWeights` work on whole files; L = 1 s already exceeds `blankingWeights`' 4 GB default;
  design in notes §5.18 / §9); run time of `blankingWeights` ∝ bin length; the parked block-B
  items (periodicRFI false emitters, per-channel noise calibration; designs in
  `2026-10-09_handover-B-closed.md`).

## Proposed order (pipeline order; Jasper may change it)

| # | Stage | Code | Notes | Docs | Tests / checks |
|---|---|---|---|---|---|
| 1 | pulsar signal | `generatePulsarSignal` | §5.1 | sky | runs, §7 |
| 2 | dispersion | `applyDispersionStream` | §5.2 | sky | `plotDispersionCheck`, NumPy ports |
| 3 | receiver noise + RFI | `addNoiseAndRFI`, `rfiSource` | §5.3 | receiver | `testRadarEdges`, `testRfiGating`, `testRfiNoise` |
| 4 | IQ conversion | `applyIQmodulation` | §5.4 | receiver | `plotIQCheck` |
| 5 | channelizer | `channelizeIQ` | §5.15 | channels | `testChannelizeIQ` |
| 6 | RFI detection, blanking | `detectRFI`, `periodicRFI`, `blankChannels` | §5.19–5.21 | excision | `testDetectRFI`, `testPeriodicRFI`, `testBlankChannels` |
| 7 | per-channel dedispersion | `dedisperseChannels` (+ `applyInverseDispersion`) | §5.16, §5.5 | channels, fullband | `testDedisperseChannels` |
| 8 | detection | `detectChannels`, `powerCovariance` (+ `detectPower`) | §5.17, §5.6 | channels | `testDetectChannels` |
| 9 | data weights | `blankingWeights` | §5.18 | excision | `testBlankingWeights`, `expBlankingVariance` |
| 10 | fold | `foldProfile` | §5.7 | fold | `testFoldWeights` |
| 11 | TOA fit, detection | `estimateTOA`, `detectPulsar`, `combineChannels`, `validateTOA` | §5.8, §5.9, §5.13 | toa, detection | `testChannelTOA`, `testOptimalTOA`, `runH0` |
| – | drivers | `main.m`, `pipelineParams.m`, run scripts | §3 | pipeline, scripts | – |

Per stage: (1) the physics in plain words, (2) a walk through the code (what each block does
and why), (3) how it is tested and what the tests do not cover, (4) refinements Jasper wants
→ a list; each refinement is one agreed edit (explain first, wait for OK).

Refinement candidates already known:
- The same stage calls repeated in `main.m`, `runMonteCarlo`, `runSNRSweep`, `runH0`,
  `runRFITest`, `runLockedRadar` (notes §10 item 7): look at it when the drivers come up.
- Stages 1–4 and `detectPower` have no regression unit tests (validated by runs and NumPy ports
  only); there is no single "run all tests" entry point.
- `expectedPowerModel` is still marked [written] (notes §5.10).

## Working rules (CLAUDE.md and memory, unchanged)

- One stage at a time; explain the physics and the planned change, wait for OK; one logical
  change per edit; never edit files Jasper did not ask about.
- Jasper runs MATLAB himself; Claude runs it only when asked (Claude may run a unit's tests;
  "validated" only after Jasper's run).
- Ground-truth rule: processing stages use only ephemeris + receiver settings.
- Dated note in `logs/` after each session; docs and notes updated with the code; after a code
  commit rerun `python3 docs/tools/make_api.py` and commit `api.html`.
- Commit straight to main when asked; Jasper pushes.
