# 5 Nov 2026: stage 7 simplified for a static site; moving-platform code archived

Jasper: "can we simplify the docs and the code for the 7a.5, where now ship movements etc are all elaborated", then "let's put it like
this: to prevent bloat, I want to work towards a static earth position, with code later developed if we want to move to a moving
platform", then "also, do this re-formatting for the entirety of "7: phase prediction and folding". you can archive the old code".

This replaces the 2 Nov rule "never delete or thin out ship material" for code: it goes to the archive, and moving-platform code is
developed later if wanted.

## Why this is safe for a static site
A fixed site's delay comes from the Earth's orbit and rotation, which are smooth, so one-hour Chebyshev segments of degree 12 pass the
1 ns check first time (two days: 48 segments, 1.86e-10 s). Segment splitting only ever fired at a kink, a ship's turn on a linear track.

## Code
- `makeTrajectory`: only `'fixed'` and `'custom'` (184 to 101 lines). Removed `'track'`, `'barycentric'`, the interpolant, the options
  Time/Pos/Interp and the output tRange, which nothing else used.
- `makePredictor`: no splitting, no MinSeg (190 to 181 lines). A segment that misses its check is refused with `makePredictor:tol`.
- `delayModel`: one help line. Its SSB branch stays, for the test probes.
- `testMakeTrajectory` (180 to 106 lines): fixed, custom, three input checks. `testDelayModel` test 5: a custom receiver at the site equals
  the fixed site, and the site's delay rate is v.n/c. `testMakePredictor` test 5: the kink now comes from a custom trajectory.
- `runAllTests`: testNavFilter out, so 33 tests. `runScenario`: a comment.
- Moved to `old/moving_platform/`, unchanged: `navFilter.m` (291 lines), `runShip.m` (342), `testNavFilter.m` (262). Copies from commit
  4433443 of the old makeTrajectory, makePredictor, their tests, testDelayModel, gen_predseg.py and the old stage page. README says how to
  resume. Nothing outside `old/` calls them any more.

## Kept on purpose
- `'custom'`: testDelayModel, testFoldPredictor, testTrueResiduals and testMakePredictor use it for the geocentre and for probes with a
  closed-form truth.
- `makePredictor`'s `'Clock'` option, about 30 lines: it is for the SDR's oscillator, not for a ship. It could go only if the geocentre
  predictor of `logs/2026-11-04_static-redesign.md` is built. Test 8 and testTrueResiduals depend on it. Jasper's decision.
- `runScenario`'s `motion` block and testSimulateFold's "moving receiver": despite the name they use a fixed site with an offset. The
  header of runScenario still says "a ship that sailed south" for the 30 S site. Not stage 7, left alone.

## Docs
- `fold.html`: 7a rewritten around the Earth carrying the receiver; Figure 7a.1 without the feedback loop; Figure 7a.3 removed; the options,
  outputs and test rows follow the code. Ship words remain only in the History boxes, which are records.
- `tests.html`: counts, Figure 2 regenerated (the generator first reproduced the old figure exactly), testNavFilter row removed, its
  section marked archived. Figures README, `delaymodel.py` (`track()` removed), `gen_testmatrix.py`; `gen_predseg.py` deleted.
- Short notes only in `navigation.html` (10.4), `scripts.html` (runShip), `pipeline.html`, `localisation.html`, `index.html` (one card).
- `maintenance.html` row; `api.html` regenerated, then once more from the commit so its header names it.
- Not rewritten, and still describing navFilter or runShip as they stood: roadmap, results, limits, glossary, lessons, conventions,
  formulas, receiver, simulator, history, the body of navigation and pipeline's Figure O.2. They wait for the navigation redesign.

## Checks
- MATLAB batch, my run (Jasper has not run these yet): testMakeTrajectory, testDelayModel, testMakePredictor, testFoldPredictor,
  testTrueResiduals, testMakeClock, testNavFit, testNavFitFaults: ALL PASSED. The other 25 tests and the scenario scripts were not run.
  Key numbers: 48 segments, 1.86e-10 s; delay rate vs v.n/c 3.0e-7 relative; the kink refused.
- Docs: no broken internal link or anchor on any page; fold.html and tests.html render without console errors. Other changed pages were
  checked for links only, not opened.

- Jasper's own `runAllTests`, after the push: "all passed" (his report; I did not see the output). The stage 7 badges hold, and the five
  stage 7 rows of tests.html now name this run.
