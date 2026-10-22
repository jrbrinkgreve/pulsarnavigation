# Handover for block E (written 22 Oct 2026, evening; start: 23 Oct)

Session note for 22 Oct and the starting point for block E. Jasper: "let's start with block E tomorrow".

## Update, 22 Oct late evening: read this first

Three things happened after this note was written. They change where E should start.

1. **A context check of this plan** (`.claude/runs/2026-10-22_overview_E-context.md`) found that E1 as drafted is off the
   critical path:
   - The planning study `expNavPlan` already answers what extra northern pulsars add.
   - The 'crystal' clock preset is the paper's 3e-11 oscillator, steadier than a stock SDR TCXO (0.65 against ~3.6 µs per
     hour, log estimate).
   - D's 73 km is a smoothed (hindsight) number.
   - runShip has no RFI or scintillation, and B0950+08 is a single Gaussian.
2. **The ship at worst-case flux** (Claude's run, 22 Oct 19:48–21:49; `.claude/runs/2026-10-22_evaluator_ship-worst.md`;
   `data/scenario/ship_14d_worst*.mat`):

   | 8 ships, 14 days | 2-D rms error, median | end error, median |
   |---|---|---|
   | dead reckoning | 372 km | 644 km |
   | self-calibrated, worst flux | 358 km | 622 km |
   | calibrated, worst flux | 298 km | 499 km |
   | calibrated, nominal flux | 73 km | 141 km |

   - The cause: B0329+54 gives 18–23 TOAs instead of ~115, and B0950+08 gives 1–3. The daily filter starts on day 6–9.
   - The pull spread is 1.6 east and 1.2–1.5 north, against ~1 for the twin, so the TOA chain loses most of it.
   - A start-up artefact: navFilter puts the port's 1 km prior at the first TOA (`navFilter.m:134`), which comes 9–14 h
     after departure at worst flux. Small fix, for Jasper to decide.
   - **Conclusion: at the design point, from 52° N, B0329+54 + B0950+08 do not navigate a ship.** The levers come before
     E1: collecting area or flux, the oscillator, the route and the pulsar set.
3. **The code map for the open decisions** (`.claude/runs/2026-10-22_coder_E-code-map.md`):
   - navFilter and navFit leave the clock's flicker out (navFit :41), which holds for the crystal only.
   - runShip's file tag would overwrite D's files in runs with another clock, RFI or scintillation.
   - The real-time error needs no filter change.

**Decisions for Jasper, in the order they block work:** route (52° N only, or south too), flux (the worst case is the
6 Oct rule), oscillator (TCXO, OCXO, Rb or CSAC). Then measured profiles, RFI and scintillation in runShip, the navigation
metric, and how much fusion goes into E.

Also open: the notes catch-up for C and D (header still 13 Oct), and the docs review's 11 questions
(`.claude/runs/2026-10-22_docs_update.md`, end).

## Where things stand

- Blocks A, B, C, D are closed. D closed 22 Oct evening; its handover with results and the ten parked items:
  `logs/2026-10-22_handover-D.md`. Full history of C and D: `logs/2026-10-20_blockC.md`. Project notes:
  `logs/current_project_notes.md` (open items in §10; another session added "Doppler-shifted dispersion delay, real data
  only" there and in limits.html on 22 Oct, `logs/2026-10-22_doppler-dispersion.md`).
- Roadmap (docs/roadmap.html) now has six labels: A-D closed, **E next**, **F later** (added by Jasper 22 Oct):
  - **E, the full navigation system in simulation:** several pulsars at full scale, the navigation solution (position +
    clock), the 3x3 array, the second polarization.
  - **F, the real implementation:** IQ recordings from antenna -> mixer -> SDR through `processIQ`; the processing ported
    to a GPU or FPGA and checked stage by stage against `processIQ`; in the end real time.
- Everything is committed and pushed (last commit 1c16951). No runs of this session are going. Untracked and not ours:
  `.claude/launch.json`, `deliverables/`, `logs/2026-10-21_task2.4-writeup.md`.
- Unit tests: `tests/runAllTests.m` lists 25; all pass in Jasper's runs (the D ones on 22 Oct).
- Docs: 31 pages in stages (7a prediction, 4a clock, 1-7 fast simulator, 10 navigation with six data figures);
  0 broken links. Reminder scheduled for 23 Oct 09:00: have an agent check all docs for missing or outdated text
  (scheduled task `remind-docs-check`; it only reminds, with a ready brief).

## What exists for E to build on

| Need in E | Exists | Where |
|---|---|---|
| Fast simulation of folds for hours to days | `simulateFold`, `adaptiveSubints`, `foldToTOA` | stage 1-7 page |
| Truth of a moving receiver and a clock | `trueResiduals`, `makeClock`, `makeTrajectory` | prediction, clock pages |
| Several pulsars on one timeline, shared clock | `runShip` (B0329+54 + B0950+08, simultaneous beams) | scripts page #ship |
| Navigation solution, fixed site / moving receiver | `navFit` (batch), `navFilter` (Kalman + smoother), both multi-pulsar, with `'Clock'`, `'Start'`, `'ClockPrior'`, `toa.tSub` | navigation page |
| Which pulsars help | `tests/expNavPlan.m` (navFit `'Predict'` on real pass schedules from 52 N) | navigation #plan |
| Exact Earth for real data | `loadEphemeris` (DE440), `delayModel` `'Ephemeris'`, `'DeltaT'` | prediction #ephemeris |
| Budgets per pulsar | `runScenario` (B0329+54 52 N, J0437 30 S; TOA rates, sigma) | scenario page |

The scenario budget already assumes the 3x3 array beamformed coherently (A_eff 0.61 m^2, SEFD 3.2e5 Jy nominal) and both
polarizations (n_pol 2) over 1.2-1.6 GHz. In E, "3x3 array" and "second polarization" therefore mean modelling them
explicitly (element signals, beam weights, calibration, the noise each brings), not a gain that is still missing.

## A first proposal for E (to agree with Jasper before any code)

Per CLAUDE.md: explain the physics and the planned change first, wait for the OK, one stage at a time.

1. **E1, several pulsars at full scale:** the ship (or a fixed site) with all pulsars visible from 52 N that the planning
   study names: B0329+54, B0950+08, B1642-03, B1929+10 (B1749-28 never rises above 20 deg there). Their profiles and
   fluxes from scenario.html's table, each in its own beam, one shared clock. Questions: ship or fixed site first;
   nominal or worst flux; how many days; whether J0437 (southern route) belongs in E.
2. **E2, the navigation solution:** position + clock from all of them with `navFilter`; the clock-synchronisation gain the
   planning study predicts (the polar component). Includes three parked D items that bear on it: the Earth's real axis
   in the fit's rows (for real data), the clock's wander inside a long TOA (an extra variance), and the self-calibrated
   ship's north gap against the ideal twin.
3. **E3, the 3x3 array:** what to model - element positions, steering weights for simultaneous digital beams, element
   phase calibration, the beam pattern (RFI from the side lobes), and noise correlated between elements. Needs Jasper's
   input on the hardware (element type, spacing, the SDR's 18 channels = 9 elements x 2 polarizations).
4. **E4, the second polarization:** the simulator already uses n_pol 2 as a sensitivity factor; the voltage chain and
   `processIQ` are single-polarization. Decide whether E models two IQ streams (and polarization-dependent RFI, Faraday
   rotation) or keeps the factor.

## Decisions Jasper may want to make first

- The first scenario of E: the ship from D (with more pulsars) or a fixed site.
- How far E goes into the array (a sensitivity model, or element-level voltages through `processIQ`).
- The oscillator: the planning study showed that a millisecond pulsar needs a better clock than the SDR's crystal; E is
  where an OCXO or a disciplined clock would enter the model.
- Which parked D items to take into E (handover-D lists ten).

## How to resume

Read this note, `logs/2026-10-22_handover-D.md`, the memory index, and docs/roadmap.html. Then propose the first step of E
to Jasper with its physics, and wait for the OK. Claude may run unit tests; "validated" only after Jasper's run. Keep
docs, logs and the notes updated with each step; commit straight to main.
