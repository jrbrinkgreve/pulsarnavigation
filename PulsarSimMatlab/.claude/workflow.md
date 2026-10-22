# Agent workflow (first trial 20 Oct 2026, stopped for lack of visibility; second trial from 22 Oct)

One interactive **main session** (Jasper + Claude as manager) and one **linked session per role**,
which the main session starts and instructs (since 22 Oct evening; before that, in-process
subagents whose conversations Jasper could not open). Jasper talks mainly to the main session, but
he can open any role session from the sidebar, read it live, steer it, or put it in split view.
Role definitions: `.claude/agents/*.md`; each role session reads its file first.

Jasper's four roles (22 Oct): 1 coding = **coder**, 2 testing = **evaluator**, 3 docs = **docs**,
4 context check + big-picture second check = **overview** (with **reviewer** for the physics and
code of a diff). **researcher** is on demand.

## Roles

| role | owns (may edit) | does | never |
|---|---|---|---|
| **main** (coordinator + integrator) | `.claude/board.md`, git | explains physics/plans to Jasper, gets his OK, writes task specs, dispatches, relays results, makes every commit | edits code/docs/logs itself while a role owns them |
| **coder** | the code files named in its task (`functions/`, `tests/`, top-level scripts) | one logical change / one new function at a time, from an agreed spec | runs long MATLAB jobs, edits docs/logs, commits |
| **evaluator** | `.claude/runs/` (outputs), backups of `data/` result files | runs tests and experiments (MATLAB batch), reports numbers against expectations | edits code, fixes failures itself |
| **docs** | `docs/` (pages, figures, `api.html` via `make_api.py`) | keeps the site in step with each change, checks rendering | edits code/logs, runs MATLAB |
| **overview** | `logs/` (notes, dated session log) | records what happened, keeps the critical path, flags rabbit holes and inconsistencies | edits code/docs |
| **reviewer** | nothing (read-only) | physics/statistics/ground-truth review of plans and diffs | edits anything |
| **researcher** | nothing (read-only) | literature and catalogue values with sources | edits anything |

## Flow of one change

1. main ↔ Jasper: physics and planned change → Jasper's OK (CLAUDE.md). Before the plan goes to
   Jasper, overview checks it against the notes, handovers and parked items (context check).
2. reviewer on the plan when it is non-trivial (optional).
3. coder implements → reports files, backward-compatibility check, what to run and what to expect.
4. reviewer on the diff ∥ evaluator runs the tests / experiment.
5. main reports the numbers to Jasper ("Claude's run"; "validated" only after Jasper's own run).
6. docs ∥ overview update the site and the logs.
7. main commits (one commit per logical change), then docs regenerates `api.html` (second commit).

Independent tasks run in parallel (e.g. a long evaluator run while the coder does the next change).

## Visibility (the lesson of the first trial)

- Each role is its own linked session, so its whole conversation is visible.
- One change in flight at a time; parallel only for a long run next to the next small task.
- Every agent writes its full report to `.claude/runs/YYYY-MM-DD_<role>_<task>.md` (Jasper can
  open it) and returns the same text to main.
- Main relays each report in the session: what was done, the numbers, the agent's reasoning,
  open questions. Never only "done".
- Main keeps `.claude/board.md` current at every dispatch and return.

## Rules that keep parallel work safe

- **File ownership** as in the table: two roles never edit the same file. Shared facts go through
  the main session.
- **MATLAB**: only the evaluator runs it (the coder only when its task says "smoke run").
  One heavy job at a time (26 GB machine; a 1 s voltage-chain run needs 5–8 GB): check
  `pgrep -fl MATLAB` first. `matlab -batch ... > file 2>&1`, never a pipe.
- **data/**: result files are backed up before a run that would overwrite them, and restored.
- **No worktrees**: `data/` (GBs, untracked) is not in a worktree; everyone works in this checkout.
- **Git**: only main commits; straight to `main`, no branches/PRs.
- A role session is started when its role gets its first approved task, then continued with
  SendMessage, so it keeps its context. All role sessions share this checkout (no worktree).
- A role session that finishes a turn reports to main; main checks the result (diff, numbers)
  before relaying it as done.

## Board

`.claude/board.md` (local, not committed): tasks, owner, status, where the output is.
