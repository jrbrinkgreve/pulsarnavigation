# logs

Dated development notes for the pulsar-navigation pipeline. They keep the context of the
project across sessions (for Jasper and for Claude).

## Convention

- One file per working session: `YYYY-MM-DD_<short-topic>.md` (ISO date, so files sort
  in time order).
- `2026-09-30_project_notes.md` is the full reference document (conventions, physics,
  stage reference, validation record, roadmap). It is the baseline; later logs record
  **what changed** relative to it, they do not repeat it.
- When the reference document gets too far out of date, write a new full snapshot
  (`YYYY-MM-DD_project_notes.md`) instead of editing the old one, so the history is kept.

## What a session log contains

1. **State at start**: what exists, what runs, open items carried over.
2. **Work done**: changes to code, with the reason for each.
3. **Results**: numbers from runs, compared to predictions.
4. **Decisions**: what was chosen and why.
5. **Open items / next steps**.

Status markers as in the reference: **[validated]** run and checked against ground
truth; **[written]** code exists, not yet run; **[todo]** not implemented.
