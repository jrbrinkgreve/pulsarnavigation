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

## Later on 5 Nov: the per-bin noise equation
Jasper asked how W, V and X couple to Section 7.4, then: "and yes, add the per-bin equation if it helps with clarity".

- `fold.html` 7.5: new Eq. (7.3), the form of Eq. (7.2) with per-bin V_k and X_k(L), and a sentence that W enters only through the
  denominator. The symbol list of 7.4 and the summary point to it. Eq. (7.2) and (7.3) are set on two lines: on one line (7.3) was cut off,
  and (7.2) already ran 55 px past the 643 px column, under its label. `maintenance.html` row. No code change.
- Found while answering, not added to the docs (he asked for the equation only): V_k + sum_{L != 0} X_k(L) = W_k, the per-bin form of
  V + 2 sum X = 1. Derivation: C = H diag(keep) H' with H a flat all-pass filter, so C^2 is proportional to C, and sum_b |C_ab|^2 is
  proportional to C_aa. Physically, a blank leaves mean W m and bandwidth W B, so (W m)^2 / (W B T) = W m^2 / (B T).
- Checked with Python, read only, on `data/rfi/rfi_dedisp_chan_weights.dat` (DM 25, 26.7 million bins). (V + sum X) / W is 0.998 at W = 1,
  0.994 at W 0.40-0.65, 0.986 at 0.20-0.40, 0.957 at 0.05-0.20 and 0.69 at 0.005-0.05. Not separated: how much of the low end is lags
  beyond Lmax and how much is the channel taper.
- Worked numbers from the same file, five time bins in one phase bin with a = 1. Clean (channel 1, bins 104166-170): weight 5.000,
  weight2 4.846, var(prof)/sigma^2 0.194, or 0.178 without X. Half-blanked (channel 2, bins 176090-94): weight 2.494, weight2 2.323,
  var(prof)/sigma^2 0.374, or 0.302 without X. So leaving X out costs 19 % instead of 8 %. With 363 time bins per phase bin, weight2/weight
  is 0.998 clean and 0.992 at mean W = 0.40.

## Later on 5 Nov: the optimisation problems behind the TOA fit
Jasper asked what optimisation problem FFTFIT solves, then: "can we couple the related optimisation problem as a written out link, and put it
in the docs?" No code change.

- `toa.html` 8.2: FFTFIT is least squares, Eq. (8.3). Parseval and the phase ramp give Eq. (8.4), the cross-correlation is Eq. (8.5), and the
  parabola in b gives Eq. (8.6). Only C(tau) depends on the shift, so the three-parameter problem becomes the highest peak of C(tau).
- `toa.html` 8.6: optimal weighting is weighted least squares, Eq. (8.10), with the weights held fixed each round. Eq. (8.11) gives
  chi2_min(tau) = const - N^2/Dn, which defines the N and Dn the page already named. With one channel and equal weights, N = (2/NBin) C(tau)
  and Dn does not depend on tau, so it becomes Eq. (8.6).
- Renumbered: old (8.3) to (8.5), (8.4) to (8.7), (8.5) to (8.8), (8.6) to (8.9), (8.7) to (8.12). Older logs quote the old numbers. `gen_toa.py`'s
  caption follows, and its output still equals the page's Figure 8.1 block. `formulas.html` got two rows; the summary box points to the equations.
- Stale number fixed: step 2 said 1/8192 turn = 1.22 us, left over from NBin 1024. At NBin 2048, Upsample 8 and a 10 ms period the coarse grid is
  1/16384 turn = 0.61 us.
- Checked in Python against direct least squares (2048 bins, a two-component template, 5 channels with blanks and unequal weights), not in MATLAB.
  FFTFIT vs direct: tau agrees to 3.3e-10 turns (the limit of the direct search), a and b to 1e-14. The Parseval form and chi2_min(tau) hold to
  4e-16 relative, the weighted chi2_min(tau) to 4e-16, and b_hat and a_c to 3e-14. One channel, w = 1: N / ((2/NBin) C) - 1 is 5e-15 and Dn is
  constant to 4e-16. With blanks and unequal weights Dn varied by a factor 1.22 with tau. C(tau) had 8 local maxima in that example.
