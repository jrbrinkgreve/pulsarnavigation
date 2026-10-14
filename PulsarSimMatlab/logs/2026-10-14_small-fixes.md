# 2026-10-14 – small fixes from the open-items list

Jasper asked for the list of things to fix (merged from notes §10 items 11–12,
`2026-10-13_code-review.md`, `2026-10-13_full-review.md`, `2026-10-13_deep-analysis.md`),
then: "fix the small things autonomously, and log the changes everywhere they need to be
logged". Done under the 13 Oct rule (small logic: implement, then validate each). A parallel
review-week session was editing `roadmap.html` and the notes at the same time; only
targeted hunks were added there. Its log (`2026-10-14_review-week.md`) was not touched.

## What counts as small

Two code items, both from the deep analysis:

1. its §9: `estimateTOA`'s convergence flag was dropped;
2. its §4: excision on a bright pulsar, options (a) a warning in `main` and (b) a stated limit
   in the help.

Checked and already done:

- I3 (ideal I/Q) is in notes §9 and on the limits page;
- the stale roadmap "commit pending" was fixed 13 Oct.

Left for Jasper, not small:

- design choices: the χ² threshold, read-through items 3/4/8, R2, R3, D2;
- new test functions: D1, N1, I1, C1 and `detectPower`;
- the code-review performance, structure and validation items.

## 1. `estimateTOA`: `toa.converged`

`fitOptimal` iterates fit → weights until τ and b change by < 1e-12, at most `MaxIterations`
(10). It set `r.converged`, but the caller never copied it. Now:

- `toa.converged` is kept per sub-int (false when invalid) and in `toa.total`;
- 'equal' (`fitOne`) does not iterate, so it is true;
- both fit result structs start with `converged = false`, so a failed total fit is safe;
- the summary prints "N sub-int(s) not converged in 10 rounds (fit -> weights; TOA kept)"
  only when that happens;
- the help text is updated (`MaxIterations`, outputs).

Validation (Claude's MATLAB runs, R2026a):

- New vs the HEAD copy (`estimateTOA_old`) on main's saved folds: channel path 'optimal',
  'equal', 'optimal' with `MaxIterations` 2, and the full band.
- In all four, every old field is bit-identical (`sameFields`, `isequaln` on info). The only
  new field is `converged`. No invalid sub-int is marked converged.
- Converged: all 8/8 valid on the channel path, 9/9 on the full band, and the total fold. With
  2 rounds: 0 of 8, and the summary line appears.
- `checkcode`: the same single pre-existing note (`mustBeInRange`) as HEAD.
- `testChannelTOA` and `testOptimalTOA` in one session: ALL PASSED (147 s).

## 2. Bright pulsars and excision

- `detectRFI` help: a "Limit" paragraph. It is valid for ρ ≲ 0.3 (`snrDB` ≤ −5). The peak flag
  rates from the deep analysis §4 are given (10 % at 0 dB, 63 % at +3 dB, almost all from
  +5 dB; the frequency guard from +10 dB; noise-free flags every on-pulse window). The mask
  depends on the pulse, so amplitude, SNR and χ² are wrong, while the TOA stays unbiased.
- `main.m`, right after `checkConsistency`: warning `main:brightExcision` when
  `~fullBand && excision && snrDB > -5` (Inf included). It is a driver-level check of a
  simulation setting, like the reuse checks; processing is unchanged. The default −5 dB stays
  silent.

Validation:

- `detectRFI`: one hunk of 9 added lines, all inside the help block, none removed;
  `checkcode` clean.
- `main.m`: the new block was extracted from the file by text and run over 8 cases. Full band
  never warns; excision off never warns; −5 and −54 dB do not warn; −4.9, −3, 0 and Inf do.
  `checkcode` is clean, as at HEAD.
- `main` itself was not run (it would rewrite `data/`).

## Docs and notes

- Notes:
  - §3: the `main:brightExcision` warning;
  - §5.8: outputs and the `converged` paragraph;
  - §5.20: the bright-pulsar limit;
  - §9: RFI excision limits (a);
  - §10 item 12: both findings marked done, the open list trimmed.
- Docs:
  - `toa.html`: weight iteration, per sub-int and total outputs;
  - `excision.html` §9: limits, a new item;
  - `limits.html`: the RFI excision summary;
  - `scripts.html`: the `main.m` steps;
  - `parameters.html`: the `excision` row;
  - `roadmap.html`: a new row "Small fixes from the deep analysis";
  - `maintenance.html`: a changelog row.
- After the commit: `api.html` regenerated, home page "checked against commit" line (done, see Status).

## Status

Committed on Jasper's request as `6ad5b4c` (only this session's hunks of the shared files `roadmap.html`,
`limits.html` and the notes; the review-week session's work stays in the working tree). Then
`api.html` regenerated from `6ad5b4c`, and the home page, roadmap status line and changelog
pointed to it (docs commit). Claude's validation only. Jasper's runs to call it validated:

- `runAllTests` (or `testChannelTOA` + `testOptimalTOA`) for item 1;
- a `main` run for item 2, as is (no warning at −5 dB), and optionally with `snrDB = 0` and
  only `runStage.receiver` / `process`, to see the warning.
