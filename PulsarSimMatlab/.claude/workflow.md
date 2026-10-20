# Parallel workflow (trial, 20 Oct 2026)

One interactive **main session** (Jasper + Claude as coordinator) and role agents that the main
session instructs. Jasper talks only to the main session. Role definitions: `.claude/agents/*.md`.

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

1. main ↔ Jasper: physics and planned change → Jasper's OK (CLAUDE.md).
2. reviewer on the plan when it is non-trivial (optional).
3. coder implements → reports files, backward-compatibility check, what to run and what to expect.
4. reviewer on the diff ∥ evaluator runs the tests / experiment.
5. main reports the numbers to Jasper ("Claude's run"; "validated" only after Jasper's own run).
6. docs ∥ overview update the site and the logs.
7. main commits (one commit per logical change), then docs regenerates `api.html` (second commit).

Independent tasks run in parallel (e.g. a long evaluator run while the coder does the next change).

## Rules that keep parallel work safe

- **File ownership** as in the table: two roles never edit the same file. Shared facts go through
  the main session.
- **MATLAB**: only the evaluator runs it (the coder only when its task says "smoke run").
  One heavy job at a time (26 GB machine; a 1 s voltage-chain run needs 5–8 GB): check
  `pgrep -fl MATLAB` first. `matlab -batch ... > file 2>&1`, never a pipe.
- **data/**: result files are backed up before a run that would overwrite them, and restored.
- **No worktrees**: `data/` (GBs, untracked) is not in a worktree; everyone works in this checkout.
- **Git**: only main commits; straight to `main`, no branches/PRs.
- Agents are spawned once and then continued with SendMessage, so they keep their context.

## Board

`.claude/board.md` (local, not committed): tasks, owner, status, where the output is.
