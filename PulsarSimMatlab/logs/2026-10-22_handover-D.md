# Block D handover (22 Oct 2026, prepared for the close-out)

Block D: motion, barycentric phase prediction, the clock, and the navigation fit. Started 21 Oct, built and
validated 21-22 Oct. Detailed log: logs/2026-10-20_blockC.md (D sections). Docs: prediction.html (stage 7a),
clock.html (4a), simulator.html (#truth, #motion), navigation.html (stage 10), localisation.html, scripts.html
(#scenario, #ship). This note collects the state, the results, and what is parked.

## State

| Step | What | Status |
|---|---|---|
| D1 | makeTrajectory, delayModel | validated (Jasper's runs, 22 Oct) |
| D2 | makePredictor, foldProfile 'Predictor' | validated |
| D3 | trueResiduals, makeClock, the simulator's truth in runScenario (motion) | validated; runScenario motion: Claude's runs |
| D4a, D4b | navFit (robust batch fit, the clock as correlated noise; 'Predict') | validated |
| D4c | planning study (expNavPlan); 'ClockPrior' for a calibrated start | validated |
| D4d | navFilter (Kalman filter + smoother, ship); runShip (closed loop) | navFilter validated; runShip: Claude's runs |
| D4 fix | TOAs as time-averages: navFit / navFilter toa.tSub (rows averaged over the sub-ints) | validated (Jasper's run, 22 Oct) |
| D5a-d | loadEphemeris (DE440), delayModel 'Ephemeris' (DE440, TDB - TT from DE440, precession + nutation), 'DeltaT' | validated |

Unit tests: runAllTests lists 25; all pass in Jasper's runs (the D ones on 22 Oct).

## Results (the numbers behind the docs)

- Fixed site, B0329+54, 30 days, perfect clock: 302 good TOAs (nominal), east +-31 km per 30 days; a 100 km offset at
  4.9 sigma in 100 days. Worst flux: 59 TOAs, +-71 km.
- SDR clock tracked day by day (runScenario motion.track): no TOAs lost; clock prediction error a few hundred us.
- Planning study (expNavPlan): extra pulsars add little to the 2-D position or the daily clock; with a real clock they
  separate the polar component from the clock's wander (B0329 alone 64 km -> 34 km with B0950+08, 24 km with three).
- Ship (runShip, 14 days, B0329+54 + B0950+08): first run (seed 1): dead reckoning 1103 km off at day 14, smoothed
  track 84 / 78 km rms. FINAL comparison self vs calibrated, 8 seeds, with the ideal twin: [to be filled in when the runs
  of 22 Oct evening finish].
- The analytic Earth orbit (default for simulations) is up to 3.2e5 km off the SSB (Uranus and Neptune left out); in a
  simulation that cancels. With 'Ephemeris': JPL Horizons to 4.5 cm (geocentre), 28-30 m (a site, 2026; UT1 - UTC left
  out); TDB - TT to 2.25 us rms of the USNO series.

## Findings worth remembering

- A TOA is a time-average over its sub-ints (often hours, sometimes two passes a day apart). One time tag damps the daily
  position signal (0.74 for 10 h of data) or, at the middle of two passes, flips its sign. Fixed by averaging the rows
  (toa.tSub); the run scripts tag each TOA at its data's mean hour angle (sideTag) for order and plots.
- A single ship run is far noisier than its pulls suggest: its errors are correlated along the track. Judge the filter
  by fixed-time pulls over many seeds, against the ideal twin on the same TOAs (runShip R.twin).
- A current drawn far out in the prior (seed 2: 0.76 m/s against 0.3 m/s) is shrunk by the prior, even with ideal data.
- The port clock must enter with its uncertainty ('ClockPrior'); taken as exact, the clock's error bars are 2.6x too
  small (the position barely suffers: the clock's wander absorbs the offset within a day).

## Parked (each with what it needs)

| Item | Why parked | What it needs |
|---|---|---|
| J0437's binary orbit | J0437 never rises at 52 N; needed for a southern route or a probe | an ELL1 binary model in makePredictor / trueResiduals (~3 s of delay per 5.74 days); test against a closed form |
| Real ship and probe trajectories | no data yet | makeTrajectory 'track' from logged positions (exists); a probe trajectory file reader |
| Earth orientation in navFit / navFilter rows | their rows use the Earth rotation angle only; with 'Ephemeris' data the axis turns 0.37 deg (~0.4 %, < ~6 km on the offsets they fit) | the same precession-nutation matrix in the rows (delayModel's earthOrientation as a shared function) |
| The clock's wander inside a long TOA | testNavFilter 7 north pulls 1.19: a split TOA spans a day, the crystal wanders a few hundred us | an extra variance per averaged TOA from the clock model over its span |
| Stored profiles shifted by the filter's correction | the fold is redone daily with the new predictor | shift stored sub-int profiles by the phase change instead of refolding (needed for MSPs, where the predictor matters inside a TOA) |
| UT1 - UTC and polar motion | 'DeltaT' takes TT - UT1 as one number; polar motion ~10 m | IERS Bulletin A values as a function of time |
| A probe's own clock | the code treats a probe's time stamps as TDB | the probe's proper time from its speed and potential (localisation.html) |
| The current as a Gauss-Markov state in navFilter | the filter models the wander as a random walk with the same daily change | a first-order Gauss-Markov velocity state (exact for runShip's truth) |
| Robust filter, one pass | navFilter weighs outliers once (navFit iterates) | an iterated smoother pass if outliers matter on real data |

## Open for Jasper

- Confirm the close of block D after the final ship comparison.
- Then block E (several pulsars at full scale, the 3x3 array, second polarisation, real IQ, the FPGA/GPU port).
