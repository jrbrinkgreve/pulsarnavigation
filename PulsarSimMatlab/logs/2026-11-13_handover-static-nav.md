# Handover: the static-site navigation redesign (written 13 Nov 2026)

Jasper: "let's write a handover and go to a next session". The next session starts from this note. It covers the redesign that
started on 4 Nov and what other sessions changed around it since.

## Start here
1. Nothing of the redesign is built. The proposal (D6a-e) and the geocentre predictor wait for Jasper's OK. Don't write navigation
   code before he gives it (CLAUDE.md: physics and plan first, then his OK).
2. Read `logs/2026-11-04_static-redesign.md`: the physics, the numbers, the steps, the 9 Nov correction under item 5, and the geocentre
   predictor section. The open decisions are below, in the order they block work.
3. The first question for Jasper: does the D6 plan stand, with the geocentre predictor in it? If yes, D6a and D6b come first, one
   function per step, his MATLAB run validating each.

## The objective
Locate a static antenna array anywhere on Earth from several pulsars, with digital beamforming and the full Earth motion. No dead
reckoning and no position prior. The field test will probably be in the Netherlands (~52 N). The ship comes after the field test.
Moving-platform code is developed later, if wanted.

## What changed since 4 Nov (other sessions)
- **5 Nov, stage 7 simplified, ship code archived** (`logs/2026-11-05_static-stage7.md`, commit c4aef71). `makeTrajectory` now has only
  'fixed' and 'custom'. `makePredictor` no longer splits segments: a segment that misses its 1 ns check is refused with
  `makePredictor:tol`. `navFilter`, `runShip`, `testNavFilter` and the old `makeTrajectory` and `makePredictor` are in
  `old/moving_platform/`, with a README on how to resume. `runAllTests` lists 30 tests and runs 33. Jasper's own run of 5 Nov: "all passed".
  This replaces two lines of the 4 Nov log: the bullet "Nothing is deleted: 'track', 'barycentric' ... stay" no longer holds, and step
  D6a may need no new code, since 'custom' with a zero position already gives the geocentre (the tests use it). `navFit` still asks for a
  'fixed' site (`navFit:site`, the fields ecef and enu), so D6b must accept the geocentre.
- **Kept on purpose on 5 Nov:** `makePredictor`'s 'Clock' option (about 30 lines, for the SDR oscillator). It can go only if the geocentre
  predictor is built, and that is Jasper's decision. testMakePredictor test 8 and testTrueResiduals depend on it.
- **9 Nov, correction to item 5** (`logs/2026-11-09_lo-and-clock.md`). A receiver clock's absolute offset moves only the first fit against
  the geocentre (286 m per second at 52 N), and a second pass against that fix removes it (0.03 m). An error in TT - UT1 turns the Earth
  and stays in every pass. The assumption is now: UT1 known to a fraction of a second. The clock's offset only needs a few seconds. An LO
  error acts through dedispersion only (0.2 us at 2 ppm, DM 25), not as a clock bias, unless the LO and the sample clock share a reference.
- **5 Nov, also:** docs Eq. (7.3) and the TOA optimisation problems, toa.html 8.2. A real-time estimate (about 350 GFLOP/s per beam
  without `blankingWeights`, 11 TFLOP/s with it as written). Jasper's rule from then: a working pipeline first, speed-ups later, each one
  tested.

## The proposal in short (details and numbers in the 4 Nov log)
- Physics: the site's delay is linear in its position, so the fold can use the geocentre (every site within 21.3 ms, under 0.1 turn for the
  normal pulsars). A 120 s sub-int stays sharp before the site is known, but a 10 h TOA does not, so a cold start searches the sub-ints.
  With the height known, the polar component follows from the equatorial one (north error = radial error / sin lat), with no absolute
  phases. Extra pulsars matter mainly against a daily clock term (FM-11).
- Expected 1 sigma, 52 N, east / north, height known, crystal clock: nominal 30 d about 30 / 47 km, worst 30 d about 67 / 107 km with all
  the data, against 100 / 194 km with SNR 7 TOAs. A cold start in a toy search: about 600 km after one day at worst flux, about 180 km
  after three. Covariance and a Python toy only; nothing simulated in MATLAB.
- Chain: fold against the geocentre; `siteSearch` (grid over X, Y and the clock rate); per-block residuals with no threshold; `navFit` with
  priors (height, per-pulsar constants); repeat from the new site. The TOA path stays as a check at nominal flux.
- Steps: D6a geocentre reference; D6b navFit priors and lat/lon output; D6c siteSearch; D6d the no-threshold residuals; D6e driver
  runStatic. One new function per step.

## Open decisions for Jasper (in the order they block work)
1. The D6 plan as a whole, and the four assumptions: on land with a known height; UT1 known to a fraction of a second (corrected 9 Nov);
   no position prior; the no-threshold fit as the main estimator.
2. The geocentre predictor (fold once, shift after the fold). He chose "Discuss it first" on 4 Nov; the four points put to him then were
   not answered:
   - Storage: the short sub-ints must be kept until the site is known. B0329+54, 30 days of 120 s sub-ints, 128 channels, 2048 bins:
     about 70 GB; channels combined, about 0.6 GB; also 512 bins, about 0.15 GB. Combining the channels after the fold loses the per-channel
     fit, such as a DM check. His main question to answer: is he happy to give that up?
   - Fold once is enough for the normal pulsars: the smear in a 120 s sub-int (site plus a 2 ppm clock) costs 0.03 % S/N for B0329+54 and
     0.1 % for B1642-03 from 52 N. J0437 from 30 S would lose about 36 %, so it would need 20 s sub-ints or a second fold. No millisecond
     pulsar is feasible from 52 N.
   - "Site-free" means a fixed reference: the geocentre in simulation (a true cold start), or a rough Dutch site in the field test, which
     also points the beam.
   - The fast simulator has no motion inside a sub-int (FM-10), so the smear figures are calculations. One short voltage-level test in
     `testFoldPredictor` would check them.
   Not settled then: the post-fold shift must use the exact delay at the clock-corrected time.
3. `makePredictor` 'Clock': keep or drop, which follows from decision 2.
4. Later, not blocking D6: the beam pointing for an unknown site (the 15 deg beam allows about +-330 km for under 0.5 dB, so "anywhere"
   needs a grid of beams, plus the array's heading and level; E4); the beam as a direction measurement (an idea, about 50 km after a day
   if level and heading are known; not simulated; E4); the target precision.

## On hold and waiting
- Tests: `testNavFitFaults` cases 2-4 (FM-9, FM-4, FM-26) and the FM-27 development wait for the new fit. Case 1 (FM-11) stays.
- Docs still describing `navFilter` or `runShip` as they stood (list in the 5 Nov log): roadmap, results, limits, glossary, lessons,
  conventions, formulas, receiver, simulator, history, the body of navigation.html and pipeline's Figure O.2. The roadmap still says the
  redesign comes "in the week of 9 November" (roadmap.html, near line 519). They wait for the redesign.
- `runScenario`'s motion block is a fixed site with an offset; its header still says "a ship that sailed south" for the 30 S site.

## Files
- Proposal and numbers: `logs/2026-11-04_static-redesign.md`. Stage 7: `logs/2026-11-05_static-stage7.md`. LO and clock:
  `logs/2026-11-09_lo-and-clock.md`. Archive: `old/moving_platform/README.md`.
- Python study of 4 Nov: `.claude/runs/2026-11-04_static-study/` (static_plan.py, cold_start.py and their outputs). That folder is
  gitignored, so it exists on this Mac only.
- The repository was clean on 13 Nov before this note (last commit 6810882, 10 Nov).
