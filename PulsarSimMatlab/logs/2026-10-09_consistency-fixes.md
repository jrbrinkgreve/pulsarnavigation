# 2026-10-09 – fixes from the daily consistency check

## State at start
- Code at `d501ea7`, in sync with GitHub, working tree clean.
- Daily report `~/Documents/pulsarnavigation-reports/2026-10-09_consistency-check.md`: 5 real
  inconsistencies (stale text) plus minor notes. Jasper: "you can fix them now".

## Work done (text only; no MATLAB run, no behaviour change)
- `current_project_notes.md`:
  - B6 status: §5.21 heading, §10 status paragraph and §13 rows (`main.m`, `detectRFI`,
    `periodicRFI`) now say validated 8 Oct (Jasper's run of `testPeriodicRFI`), as in
    `2026-10-08_excision.md` and the handover.
  - §10: removed the leftover "Then B5d … B6 (spectral kurtosis …)" sentence; item 2 "Done /
    Open" updated (done: B1–B4, B5a–B5d, B6; open: B5e and/or per-channel noise calibration).
  - §3: map shows `periodicRFI` (between `detectRFI` and `blankChannels`, as in `main.m`), "next"
    line updated; `runLockedRadar.m` added to the script list.
  - §5.3 heading: B5a, B5c validated 8 Oct; B5d Claude's run.
  - §10: duplicated fragment "per ~10 %; `*.asv` git-ignored." removed; planned 3e test name
    noted as "became `tests/testBlankingWeights.m`".
  - §13: rows for the other `old/` files and for `docs/`.
- `functions/applyDispersionStream.m`: help text `MaxMemoryGB` default 8 → 16 (matches the
  `arguments` block). Comment only.
- `docs/maintenance.html`: that discrepancy marked fixed; changelog line.
- `2026-10-08_docs.md`: line for `5a6ae63` (landing page), which no log described.
- Memory (outside the repo): push item removed (done); docs-site memory and index line
  (B6 documented, style.css committed, help text fixed); status memory de-staled; §10 item
  number for the memory-efficient FFT (9 → 5).

## Open items / next steps
- After committing: `python3 docs/tools/make_api.py` so `api.html` shows the new help text.
- Unchanged from the handover: the decision B5e (measure first) vs per-channel noise
  calibration (fix first); Jasper's run of `tests/testRfiNoise.m`; the 7 Oct Q&A.

## Committed (later on 9 Oct, B5e session, at Jasper's request after a review)
- All fixes checked against the code / logs and committed, except the notes §13 row for the
  other `old/` files: commit 2f6eddc had meanwhile added an identical row, so this one was a
  duplicate and was dropped.
- Not part of this commit (other sessions' uncommitted work in the same files): the paragraph
  repairs in notes §5.7 / §5.8 / §10 (paragraphs split by the f_out commit's insertions), the
  `expSamplingRate` / `demoBinningMatchedFilter` rows, the dispersed-domain idea.
- `docs/maintenance.html`: the code-state cells of the two B5e changelog rows set to their commits.
- `api.html` regenerated afterwards (new `applyDispersionStream` help text).
