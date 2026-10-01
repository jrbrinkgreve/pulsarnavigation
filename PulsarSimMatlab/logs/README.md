# logs

Dated development notes for the pulsar-navigation pipeline. They keep the context of the
project across sessions (for Jasper and for Claude).

## Convention

- One file per working session: `YYYY-MM-DD_<short-topic>.md` (ISO date, so files sort
  in time order).
- `current_project_notes.md` is the full reference document (conventions, physics,
  stage reference, validation record, roadmap). It is kept **up to date with the
  code**: edit it whenever the code or a status changes. (Started 30 Sept 2026 as
  `2026-09-30_project_notes.md`; its history is in git.)
- Session logs record **what changed and why**, in time order; they do not repeat the
  reference document.

## What a session log contains

1. **State at start**: what exists, what runs, open items carried over.
2. **Work done**: changes to code, with the reason for each.
3. **Results**: numbers from runs, compared to predictions.
4. **Decisions**: what was chosen and why.
5. **Open items / next steps**.

Status markers as in the reference: **[validated]** run and checked against ground
truth; **[written]** code exists, not yet run; **[todo]** not implemented.
