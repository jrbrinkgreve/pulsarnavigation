# 19 Nov 2026 — ChannelGain, step 2: adaptiveSubints

Jasper: "commit this and continue" (after step 1 of the plan in
`logs/2026-11-17_budget-recheck-big30.md` §7). Step 1 committed as f1397ce.

## What changed

- **`adaptiveSubints` takes 'ChannelGain'** and passes it to `detectPulsar`, whose known-phase
  statistic T0 closes the groups: [] (default, gain 1 in every channel, results as before), a
  [1 x nChan] vector, or "measured" (each channel's mean power over the whole input fold, taken
  once). `sel.channelGain` records it. Give `foldToTOA` the same value, so that the groups close
  on the statistic the TOAs are then judged by.
- **`channelGainOf` is its own function** (`functions/channelGainOf.m`). It was a local function
  of `foldToTOA` since 17 Nov; both functions need it now. Same code, error identifier
  `channelGainOf:gain` (was `foldToTOA:gain`, two days old, used nowhere else).
- No caller passes the option, so every run script and result is unchanged.

## Tests (Claude's runs; not validated until Jasper's)

Log: `.claude/runs/2026-11-17_big30/tests_19nov.log`.

- `testFoldFaults`: ALL PASSED. New case 1c, on case 1b's tapered fold (300 sub-ints, groups
  closed at T0 >= 50): ChannelGain = ones gives the 30 groups and their T0 exactly as without the
  option; with the taper, every group's closing T0 equals `detectPulsar`'s T0 of the merged group
  with that gain (difference 0); "measured" gives `foldToTOA`'s gains exactly; T0 per sub-int
  ×1.0103 with the taper (predicted sqrt(n sum g²) / sum g = 1.0103). Cases 1, 1b and 2 as on
  17 Nov.
- `testSimulateFold`: ALL PASSED (16 cases; it calls `adaptiveSubints` and `foldToTOA` without
  the option, so this is the check that the default path is untouched).

## Docs

simulator.html (the option in `adaptiveSubints`' list), detection.html D5, tests.html (case 1c),
roadmap.html (FM-41 line), maintenance.html. `api.html` regenerated after the commits.

## Step 3 is not done: a decision for Jasper

Step 3 was "the run scripts give the taper to the observer". Three ways, with what each costs:

1. **Leave the run scripts as they are (Claude's recommendation).** The option exists for real
   data and for strong signals. At the weak-signal target the chi² excess is about 5·10⁻⁵ and the
   cost in sigma_TOA is the known +2.5 % of roadmap item 3. Nothing has to be rerun, and the
   critical path is the navigation redesign, which still waits for the D6 decision.
2. **An opt-in setting in `runScenario`** (default off): the taper passed to `adaptiveSubints`
   and `foldToTOA`. Small, and existing numbers stay. It gives the truth-side taper to the
   observer, which stretches the ground-truth rule, so it would be a named exception.
3. **Shape the simulated noise with the same band edge** (item 3's other option): then the
   taper is a receiver bandpass, "measured" sees it, and no truth leaks to the observer. It is
   the physically right picture of a real receiver, but it changes the truth side of both
   simulators and every frozen reference number.
