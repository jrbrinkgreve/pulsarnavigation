# 2026-10-09 – project review (started; finished 12 Oct in `2026-10-12_overview.md`)

Jasper asked for a full review of the project state. Stopped partway at Jasper's request;
continue next session. No code changed, no MATLAB run, nothing committed.

## Read so far
- Handover `2026-10-09_handover-B-closed.md`, logs README, today's logs (f_out study,
  dispersed domain, consistency fixes), notes §1–§5.7 and §9–§13, `main.m`, `pipelineParams.m`.
- Git: HEAD `849ccf5`, working tree clean, local main 1 commit ahead of origin on 9 Oct
  (pushed since; in sync on 12 Oct).

## State in short
- Block A (channelized front end, exact noise model, optimal weighting) complete and validated.
- Block B (RFI excision) closed 9 Oct; periodic mask off by default; two parked items with
  ready designs (periodicRFI false emitters; per-channel noise calibration).
- Next: C (fast simulator, design draft in notes §10 item 3), not started; Jasper decides when.
  D (barycentric phase prediction) is a hard requirement for the −54 dB design point.
- Known unfixed bug: `blankingWeights` line 203 `fft(q, M)` → row for scalar q (notes §10
  item 10; confirmed still in the code).

## Small stale-text findings (not fixed; need Jasper's OK)
1. Notes §3 map, line 112: `addNoiseAndRFI` still marked `[written]` (§5.3 says validated).
2. Notes §5.2, line 280: "`MaxMemoryGB` (default 8" (code and help text say 16).
3. Notes §9, line 1868: periodic mask "(default since B6-3)": now off by default.
4. Notes header (line 10): "Last updated: 5 October 2026" followed by a long run-on list.
5. `pipelineParams.m:102` comment "where RFI excision will go, block B" (done).
6. `main.m:33` header says the channel fold uses `NoiseCoeffs`; with excision (default) it
   uses `DataWeights`.

## To do next session
- Read notes §5.8–§5.21 and §7 (validation record); docs site spot-check (index, roadmap,
  results); skim the functions for consistency with the notes; check test coverage per
  function; then write the full review for Jasper.
