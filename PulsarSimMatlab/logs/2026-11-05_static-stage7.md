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
- Correction, later the same day: every bin with W below 0.9 in that file lies outside `fullySupportedBins` (338 to 175172), at the file end, where
  the zero padding acts as one very long blank. The fold never uses those bins. Inside the supported range the lowest W is 0.901, so "half-blanked"
  above describes the padding, not the RFI mask. The numbers stand as arithmetic on W, V and X. The check of V + sum X against W below W = 0.9 has
  only this padding behind it, not short blanks.

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

## Later on 5 Nov: is the pipeline real-time feasible?
Jasper: "is this computing pipeline real-time feasible on an FPGA/GPU? because it will need to process 1 second of data in .... 1 second", then
"later we might be able to do some computational improvements to the algoithm, but everything needs testing, and i want a "working pipeline" before a
"fast working pipeline"". So this is an estimate only. No code or docs change, and no speed-up started.

The work figures are counted, not measured on any FPGA or GPU. Counting rules: FFT 5 N log2 N flops, complex multiply 6. One beam is two polarisation
streams of 128 channels at 4.1667 MHz. The code runs one stream today, and the beamformer and second polarisation belong to block E.

| Stage | Work | How it was counted |
|---|---|---|
| `channelizeIQ`, all 18 streams | 1.9 TFLOP/s | 3853 taps, K 256, D 192: 107 GFLOP/s per stream; about 120 DSP slices per stream at 400 MHz |
| beamformer, per beam | 77 GFLOP/s | 18 complex multiplies per channel sample; about 100 DSP slices |
| `detectRFI`, per beam | 32 GFLOP/s | about 30 flops per sample: scales 1 to 16, two passes, a median per 16 ms block; thresholds in software |
| dedispersion, per beam | 217 GFLOP/s | overlap-save, 4096 points, 1560-sample filter at DM 25; the delay across the band as a delay line, 60 MB per stream |
| detection, per beam | 4 GFLOP/s | 4 flops per sample |
| fold with weights, per beam | 21 GFLOP/s | about 80 flops per channel-bin, 133 million channel-bins/s per stream |
| `blankingWeights` as written, per beam | 11 TFLOP/s | 86 FFTs of 2^20 points per blanked channel per 0.2 s, 125 channels |

- Without `blankingWeights` a beam is about 350 GFLOP/s.
- Data rates: 57.6 GB/s into the channelizer (18 x 800 MS/s x 4 B), which is 5 PB per day, so storage forces real time. After the beamformer, 4.3 GB/s per
  beam at 16 bits. After detection, 0.53 GB/s per beam if the polarisations are summed. A 500 MS/s SDR cuts the first figure to 36 GB/s; later stages keep the
  same channel rate.
- Only the stages up to the fold face the 1 s per 1 s limit. The TOA fit runs once per sub-int. Latency is about 50 ms of buffers plus the sub-int.
- Saved MATLAB times, one stream, 0.2 s, `data/rfi`: channelize 10.6 s, blank 0.9, dedisperse 2.6, detect 0.7, weights 62.6, fold 3.9. That is 81 s, about 400
  times real time, 77 % of it in weights. `data/rfi/profile`, 1 s of data, spent 285 s in `blankingWeights`; the older copy of that run took 1052 s.
  `detectRFI` and `estimateTOA` have no saved time. The implied rates, 18-29 GFLOP/s in weights, about 11 in dedispersion, 2-3 in the channelizer, say the
  counts are the right size. The machine's load during those runs is not known.
- Sparsity: the ten saved weight runs blank 0.005-0.07 % of samples on average, and 70-128 channels hold a blank. `blankingWeights` computes a whole
  channel if it holds one blank. In `data/rfi` and `data/chan`, read directly, 0.54 % of the channel-bins in the fully supported range have |W - 1| > 1e-6. The
  others equal the clean constants. If the work scales with those bins, a beam needs about 60 GFLOP/s. Not built, not tested, and none of the ten runs is a
  dense-RFI site. On a busier site the fraction grows, at worst back to the cost as written.
- Idea, untested, not for now: V + sum X is about W (above), so the fold's phase-bin variance might need W alone where a varies slowly across the lags.
- GPU: about 2 TFLOP/s sustained on FFT work is my estimate from memory bandwidth (1 TB/s for an RTX 4090; 4096 points, 2 read-write passes), not a
  benchmark. PCIe 4.0 x16 carries about 32 GB/s, so the 57.6 GB/s does not enter one GPU.
- Precedent: MeerKAT PTUSE (Bailes et al. 2020, arXiv 2005.14366, read through ar5iv): two GPUs per node, each half the band, one node per beam, doing coherent
  dedispersion and detection, with the beamformer on SKARAB boards. With the 856 MHz L band that is 428 MHz of two polarisations per GPU. Its real-time
  sentence (Section 2.5) I saw only through a summary of the page. Parkes Medusa (Hobbs et al. 2020) only through a search summary: its PDF could not be read here.
- The split this suggests: FPGA for the channelizer and beamformer; GPU from `detectRFI` to the fold; CPU for thresholds, predictor, TOA and `navFit`. No
  hardware chosen. An X410-class SDR has an RFSoC on board and 4 channels at 400 MHz, so 18 channels would take five units.
- Not speed questions: the port must give the same numbers (FM-37), and block processing (FM-20) is needed for any long recording, since `blankingWeights`
  refuses above 4 GB, under 1 s of data. FM-20 and FM-38 in the risk register are unchanged.
