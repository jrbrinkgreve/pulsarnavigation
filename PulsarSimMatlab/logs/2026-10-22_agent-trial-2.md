# 2026-10-22 – Second agent trial (evening): set-up, the first two jobs, linked sessions

Session note, written by the overview role (role 4): first as a subagent, then as a linked session. Nothing changed in
the code or the notes, and nothing was run in MATLAB. In the docs, `api.html` was regenerated, and the docs session is
fixing the pages (see "The docs job widened").

## State at start

- Block D closed on 22 Oct. The handover for block E is written and committed (65b1c95, `2026-10-22_handover-E.md`).
  E was planned to start on 23 Oct.
- Unit tests: `runAllTests` lists 25, and all pass in Jasper's runs.

## The second trial

Jasper (22 Oct, evening): "I'd like to try the agentic development again: can you spin up a main manager session, with
subagents that 1: do the coding, 2: do the testing, 3: do the docs writing, and 4: do the context checking and the
bigger picture reviewing as a second check?"

The first trial (20 Oct) stopped within the hour. Jasper said then: "nevermind, go back to the old development routine
... I liked the visibility of the development in the old form" (`2026-10-20_blockC.md`).

The roles, from `.claude/workflow.md`:

| Jasper's role | agent | owns |
|---|---|---|
| 1 coding | coder | the code files named in its task |
| 2 testing | evaluator | `.claude/runs/` and backups of `data/` results. It is the only agent that runs MATLAB. |
| 3 docs writing | docs | `docs/` |
| 4 context check and big picture | overview, with the reviewer for a diff's physics and code | `logs/` |
| on demand | researcher | nothing (read-only) |

The main session is the manager. It talks to Jasper, gets his OK, dispatches the tasks and makes every commit.

## Visibility rules added to workflow.md today

- One change is in flight at a time. Work runs in parallel only when a long run goes next to the next small task.
- Every agent writes its full report to `.claude/runs/YYYY-MM-DD_<role>_<task>.md` and returns the same text to main.
- Main relays each report in the session: what was done, the numbers, the reasoning and the open questions. It never
  reports only "done".
- Main keeps `.claude/board.md` current at every dispatch and every return.

Also new is flow step 1. Before a plan goes to Jasper, the overview agent checks it against the notes, the handovers and
the parked items. This is job 0 in `.claude/agents/overview.md`.

`workflow.md` and `overview.md` are changed but not committed yet. Main makes the commit.

## Board reset for E

`.claude/board.md` is local and not committed. The 20 Oct board (tasks of blocks C and D) closed with D. The new board
has four tasks:

1. Context check for E (overview): running.
2. Docs audit (docs): running.
3. Propose E1 with its physics to Jasper (main): starts after task 1.
4. E1 code (coder, then reviewer, then evaluator): starts after Jasper's OK.

## The first two jobs

1. **overview: context check for block E.** Inputs: handover-E, handover-D, notes §10 and the roadmap. Report:
   `.claude/runs/2026-10-22_overview_E-context.md`. The short version, which is the overview's view and has not been
   discussed with Jasper yet:
   - E1 as drafted (four northern pulsars at full scale) asks a question that the planning study `expNavPlan`
     already answered.
   - The ship has never been run at the worst-case design point.
   - The clock preset is steadier than a stock SDR TCXO.
   - Suggested first: Jasper decides route, flux and oscillator, and a worst-flux `runShip` run gives the baseline.
2. **docs: audit of all pages** for stale, missing or inconsistent text, from the brief of Jasper's 23 Oct reminder
   (scheduled task `remind-docs-check`). It started as report-only. Jasper widened the job to fixing the pages the
   same evening (see "The docs job widened"). Report: `.claude/runs/2026-10-22_docs_audit_partial.md`.

Both jobs ran as subagents inside the main session. Jasper could not open their conversations.

## The move to linked sessions

Jasper (22 Oct, evening): "can you move the operation to such a nested-session setup?"

Main proposed a docs session, and Jasper started it. Then Jasper: "i started the docs session. Can you also start
sessions related to the other funcrtions?"

Since then each role is a linked session that main starts, instructs and reports to. Jasper can open each one from the
sidebar, read it live and steer it. `workflow.md` describes the new set-up (still uncommitted).

| role | session's first task | output |
|---|---|---|
| docs | fix the pages from the audit, in four groups | `.claude/runs/2026-10-22_docs_update.md` |
| overview | bring this log up to date | this note |
| evaluator | prepare the worst-flux ship baseline without running it: command, run time, outputs, comparison numbers | `.claude/runs/2026-10-22_evaluator_ship-worst-prep.md` |
| coder | a read-only code map for the open E decisions (oscillator, profiles, metric, fusion) | `.claude/runs/2026-10-22_coder_E-code-map.md` |

Reviewer and researcher start on demand. All sessions share this checkout (no worktree). Only main commits. The session
ids are on `.claude/board.md`.

The evaluator's and the coder's first tasks change no code and run no MATLAB. They prepare for Jasper's decisions on E.

## The docs job widened

Jasper (22 Oct, evening): "i want the docs agent to review all documentation pages, update the visuals when needed or
lacking, and check missing or outdated information".

The audit:

- A subagent read all 30 hand-written pages against the code and the logs.
- It found 44 outdated, 23 missing, 14 inconsistent and 9 style items. Some rows cover several pages.
- All 873 internal links and anchors resolve.
- It did not check the rendering (themes, phone width) or the figures' contents. It judged the figures from their
  captions only. The style check is partial.
- It stopped at main's request when the job moved to its own session.
- Report: `.claude/runs/2026-10-22_docs_audit_partial.md`.

The subagent regenerated `docs/api.html` by accident:

- It ran `make_api.py --help`. The script does not parse `--help`, so it ran in full.
- The page now says "from commit 65b1c95" (it said 32026d1).
- The only change in content is runScenario's help header (tSub, sideTag). Commit 56ec1c8 changed that header after
  the last regeneration.
- Kept, because the page is now current.

The docs session now fixes the pages in four groups, group 1 first. Main commits after each group.

Correction: the first version of this note said the docs audit "reports only and makes no edits". The audit edited no
page, but it regenerated `api.html`, and the job now fixes the pages.

## Open

- Done: the memory index entry "Single visible session" is replaced by "Agent workflow, trial 2" (linked sessions, not
  subagents).
- The notes need a catch-up for blocks C and D. Their header still says "Last updated: 13 October", and block D is
  missing from them. This is the overview's job (board task 1d), and it waits for Jasper's OK.
- Jasper's decisions for E are open: route, flux and oscillator first (context check, §4). The worst-flux ship baseline
  (evaluator, 1b) and the oscillator and profile values (researcher, 1c) wait for them.
- (Committed later: see below.)

## Later the same evening: what happened, and the trial stopped

Written by the main session after the role sessions' work.

- **Commits.** `0717f20` holds the workflow, the overview role and this note. `aa4bd59` holds the docs review: 30 pages,
  `style.css`, 19 generators, the new `gen_clock.py` and `api.html`. Both are pushed. Main reviewed the docs only by
  spot checks (roadmap, limits, the clock figure in the preview).
- **The worst-flux ship baseline** was approved by Jasper ("Yes: once the evaluator's prep is in, start the worst-flux
  ship baseline (8 seeds x {self, calibrated})"). Claude's run, 19:48–21:49, two waves of four processes (peak 10.4 GB in
  one process). All 16 runs were clean. The median 2-D rms is 358 km self-calibrated and 298 km calibrated, against 372 km
  dead reckoning (nominal: 102 and 73 km). Details and conclusions: the update at the top of `2026-10-22_handover-E.md`.
  The docs got the result on navigation#ship-worst, limits, index, results and the roadmap.
- **The coder's code map** (`.claude/runs/2026-10-22_coder_E-code-map.md`) found four points for E:
  - `navFilter` and `navFit` leave the clock's flicker out. The reason (navFit :41) holds for the crystal only.
  - `runShip`'s file tag would overwrite D's `ship_14d_nominal*` files in any run with another clock, RFI or
    scintillation.
  - The real-time navigation error needs no filter change; runShip only has to keep what it computes.
  - It maps where ship sensors would enter `navFilter`.
- **The trial stopped.** Jasper asked whether the set-up worked well, then whether work quality was lost compared with
  one session. Main's answer named these risks:
  - Jasper's own understanding suffers when work arrives in large batches (52 files in one docs commit).
  - Code, tests and docs no longer come from one understanding.
  - Main reviews other sessions' work only by sampling.
  - The gain was an independent second check: the context check found three blind spots in D's handover.

  Jasper: "i think i can then better keep to 1 session for the potential loss of work quality."

  Development is back in one visible session from now on. `workflow.md` is marked "not in use".
- **Open after the trial:**
  - The ship batch's results.
  - Jasper's decisions for E: route, flux, oscillator.
  - The notes catch-up for C and D.
  - The docs session's 11 open questions (`.claude/runs/2026-10-22_docs_update.md`, end).
