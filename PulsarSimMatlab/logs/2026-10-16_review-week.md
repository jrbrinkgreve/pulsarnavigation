# 2026-10-16 – review week, continued: stage 6c fix, stage 7

Continues `2026-10-15_review-week.md` (stages 6a–6c).

## B-text fixed (Jasper: "fix it, and continue")

`functions/blankChannels.m` help, line 7: "~60 us per channel at DM 5" → "32-75 us per channel at
DM 5" (in-channel sweep 31.7 µs at 1.6 GHz, 74.7 µs at 1.2 GHz; notes §5.19 already had 130–310
samples). Comment only: `git diff` shows that one line, no code changed. `docs/api.html`
regenerated with `make_api.py --worktree` (only that line and the stamp "from working tree"
changed); rerun `python3 docs/tools/make_api.py` after the commit to stamp the commit hash.
Notes §10 item 11 "6c".

## Stage 7 – per-channel dedispersion (`dedisperseChannels` → `applyInverseDispersion`, §5.16, §5.5)

Walk-through given (code unchanged since `0499da7`; engine since 6 Oct `AllowRefOutsideBand`):
inverse chirp H = W·e^(−i2πK·DM(f − fRef)²/(f·fRef²)) per channel, band = centre ± 1.5625 MHz
(removes the oversampling transition, channels tile the band), common RefFreq 1.6 GHz outside
the channel → each filter also shifts its channel by τ(f) − τ(1.6 GHz) (up to 6.3 ms) →
channels aligned, no shift step. In-channel smear removed: 31.7 µs (132 samples) top channel,
74.7 µs (311 samples) bottom. Taper 2 % per channel (62.5 kHz) inside each channel, no overlap →
independent channels, ∫W² = 128 × 3.046875 = 390.0 MHz as the full band. Guard 20/edgeW =
320 µs = 1334 samples per channel (longer than the sweep, cheap). Overlap-save with
[past | step | future] buffer, leakage self-check, Nfft by cost. Alignment picture in the chat.
Tests: `testDedisperseChannels` (regression bit-identical, commutation |ρ| ≥ 0.999993, end-to-end
TOAs consistent, time-bin noise 0.884 vs 0.888 predicted).

Jasper: "fix V1, note V2". **V1 done** (text only, like I2; numbers checked against the saved
run `data/mc/mc_dedisp_info.mat`) notes §5.5 still has the 30 Sep set-up
(coverage 1.15–1.65 GHz at 500 MHz, nPast 1250, nFuture 3,152,472, reference run Nfft 2^26 in
8 blocks, fullySupported [1251, 496,847,529], run 7.4 s); current (fs 800 MHz, fLO 1.3 GHz):
coverage 0.9–1.7 GHz, nPast 2001, nFuture 5,043,956, Nfft 2^25 in 3 blocks at L = 0.1 s
(2^26 in 13 blocks, 3.76 GB at 1 s), fullySupported [2002, 74,956,045]; old ones kept as
history. **V2 noted** (notes §10 item 11 "Stage 7", roadmap row "Dedispersion test"):
(a) round trip forward → IQ → inverse = band-limited input (absolute check in MATLAB; replaces the
NumPy round trip 2e-6; could be part of D1); (b) a RefFreq inside the band on the channel path
(channels above it need past overlap; only the top-of-band reference has run); (c) DM = 0.

## Stage 8 – detection (`detectChannels` + `powerCovariance`, engine `detectPower`; §5.17, §5.6)

Walk-through given (code unchanged since `0499da7`): square-law |y|² per channel after all linear
steps; bins of round(fs/f_out) = 4 channel samples (1.0417 MHz, 0.96 µs; one warning);
channels interleaved into one [128 × nBins] float32 file for `foldProfile`; exact noise of
detected bins in a narrow channel: B·dt ≈ 3, power correlation time ~1/B ≈ 0.33 µs →
V 0.8877, X₁ 0.0458, X₂ 0.0051, X₃ 0.0020 (V + 2ΣX = 1, Lmax 7); a 5-bin phase bin then has
variance 0.969 of the independent-bin value (√ → error bars 1.6 % too large without the model);
the chirp is all-pass and drops out (not with blanking → stage 9). Tests: `testDetectChannels`
(formula vs the real impulse response 6e-10; interleave bit-identical; fold of the sum 1.2e-8;
measured V, X, 5-bin variance within errors).

Jasper: "fix S1". **S1 done** (text only, like V1; checked against the saved run
`data/mc/mc_envelope_info.mat`: binLen 800, binTime0 0.499375 µs, [4, 93,695]): notes §5.6 "binLen = round(fs/f_out) =
500" is the 500 MHz value; at 800 MHz it is 800 (still 1 µs bins; binTime0 0.4994 µs;
`fullySupportedBins` [4, 993,695] for 1 s unchanged because it is time-based; "Run: 1.2 s for 1 s"
was at 500 MHz). A sweep of §5 found no other 500 MHz numbers.

## Stage 9 – data weights (`blankingWeights`, §5.18)

Walk-through given (code unchanged since `0499da7`): a blanked input sample removes data from
every output the dedispersion filter spreads it into; for Gaussian input the covariance of the
blanked dedispersed samples is exactly C(a, a+τ) = Σ_l h(l)h*(l+τ)·keep(a−l) (the mask convolved
with q_τ, one FFT pair per lag τ = 0…31); W = mean of C(a,a)/m0 over the bin (with the true
channel spectrum R_x, ~20 extra lags), V and X(L) = (1/n²)Σ|C|² over the sample pairs (white
form); h of each channel by an impulse through `applyInverseDispersion` (same Nfft, reach
checked); unblanked channels W = 1, V0, X0; MinWeight 1e-6. Tests: `testBlankingWeights`
(direct covariance matrix 3e-8, R_x vs measured 7e-4, end to end unbiased incl. on-pulse
blanks, long bins) – thorough.

Computed for the explanation (Python toy with the same channel filter, scratch only, not in the
repo): bottom channel, a 20-sample (4.8 µs) blank → ~80 output bins (the 75 µs in-channel sweep)
each with W = 0.934 and V = 0.871 ≈ W² of the unblanked value; picture drawn in the chat.
Also checked: covariance beyond the stored Lmax = 7 lags near blanks: median 0.4–0.7 % of a
bin's long-sum variance (0.16 % unblanked), up to 1.3–1.6 % for bins with W ≥ 0.1 next to a
200-sample blank; error bars only, same order as the white V, X approximation (§9), consistent
with test 3's noise check.

**BW1 noted** (Jasper: "add BW1"): §9 "Blanking: channel-input spectrum", `docs/limits.html`
(blanking bullet), notes §10 item 11 "Stage 9". No code or test candidate.

## Stage 6b again, separate session ("Current state", 15–16 Oct): P3, P4 noted

A second session (started 15 Oct; it committed the 13–14 Oct review notes as `2f0d0e8`) gave its
own `periodicRFI` walk-through, not knowing about the 15 Oct one above. Two new points; Jasper:
"note P1 and the wording point". They were called P1 and P2 in that chat; renumbered **P3** and
**P4** because P1/P2 were already taken (15 Oct). Notes §10 item 11 "6b", §9 pointer, roadmap
parked-fix entry; nothing implemented.
- **P3** (fix item d) step 2 groups by strongest channel with single linkage (`ChanGroup` 2), so
  dense flags chain channels into one group; the 9 Oct 14.1 kHz group spanned channels 29–94
  (radar 1's sidelobes + LTE edge channels 81 / 94; `2026-10-09_B5e.md`), and an emitter blanks
  the union of its inliers' channels. Fixes a and c reject the false emitter but not the merging.
- **P4** (minor, text only) "mask does not depend on the noise values, so blankingWeights stays
  exact" (help step 4, §5.21, `docs/excision.html`): what holds is no threshold truncation (a
  predicted window does not depend on the noise inside it, so the kept noise is not selected,
  as `blankingWeights` assumes); the windows still depend on the data through the fitted t0, P.
  (The 15 Oct entry above repeats the original wording.)
Also checked in that walk-through, nothing new: q = (2·TimeTol + 1)/P instead of 2·TimeTol/P
(conservative, superseded by fix a); `FalseAlarm` is per group search (as the help says);
staggered PRF, long files and O(n³) cost already recorded.

## Sync with the other session (this session, after P3/P4)

Jasper relayed that the second session added P3 and P4 (above) and stopped. Re-read notes §10
item 11, the §9 pointer, the roadmap parked-fix entry and this log: both sessions' items are
present and consistent (E1, P1–P4, B-text, V1, V2, S1; BW1 still proposed). Correction to this
session's stage-9 explanation in the chat ("the mask is independent of the noise values, so the
weights are exact"): `blankingWeights` is exact for a given mask if the kept samples are
ordinary Gaussian noise. That holds for `periodicRFI`'s predicted windows (P4), not for
`detectRFI`'s mask, which removes the highest samples: the kept noise is slightly truncated
(§9 "RFI excision limits" (a): ~1e-4 of clean data, phase-uniform at weak SNR; at −5 dB the
pulse amplitude −0.4 %, symmetric). Nothing committed; the working tree holds both sessions' work.

## Stage 10 – fold (`foldProfile`, §5.7)

BW1 noted first (Jasper: "add BW1"): §9 "Blanking: channel-input spectrum", `docs/limits.html`,
notes §10 item 11 "Stage 9".

Walk-through given (code unchanged since `2f6eddc`, 9 Oct coverage warning): phase model
φ = Phi0 + f0·(t − TRef) + ½F1·(t − TRef)² from the ephemeris; turn = round(φ), so each turn runs
from −½ to +½ around its pulse and sub-int cuts fall at phase 0.5 (off-pulse); linear assignment
of each time bin to the two nearest phase-bin centres (no sub-bin bias; 'nearest' gave 0.5 ns);
kept per (phase bin, sub-int, channel): sum Σw·P, weight Σw·W, weight2 (variance), weightX
(covariance up to D = 3 phase bins: shared time bins + correlated time bins X(L)); empty
bins (W = 0) add no power; per-sub-int tRef / fRef for the TOAs. 2048 bins of 4.88 µs ≈ 5 time
bins; linear assignment variance (2/3)σ²/5.09 = σ²/7.6 per turn and phase bin. Picture of the
linear split drawn in the chat. Tests: `testFoldWeights` 1–7 (bit-identical references, algebra
AᵀCA to 6e-16, noise model on data, data weights, fake blanking unbiased, coverage warning).

Jasper: "fix FP1". **FP1 done** (text only; NBin 2048 checked in `data/test_fold_chan.mat`;
notes §5.7 points 4 and "Why weightX matters", `docs/fold.html` callout): §5.7 describes 1024 phase bins (9.77 µs; the
function default and the 30 Sep reference run) and gives "σ²/14.7 instead of σ²/9.77" for that
case; the pipeline uses 2048 bins (4.88 µs) since 7 Oct (A1a): add the current numbers, keep the
old example as history.

## Stage 11 – TOA fit and detection (`estimateTOA`, `detectPulsar`, `combineChannels`; §5.8, §5.13)

FP1 fixed first (see stage 10). Walk-through given (code unchanged since `6ad5b4c`, 14 Oct
`toa.converged`): model p = a + b·T(φ − τ) + noise with the exact Gaussian template (§12 design
choice: no template mismatch); TOA = tRef + τ/fRef at 1.6 GHz; error bar from the template
slope and the fold's full covariance (dᵀ·Cov·d / C″²); main's channel path = 'optimal' (one τ,
one b × ChannelGain, baseline per channel, weights 1/var iterated, sandwich error bars);
quality: red. χ² flag (> 2); detection: Neyman–Pearson T0 (known phase, η 3.09) and Tmax
(unknown phase, Rice η 4.15, exact by a 10⁶ Monte Carlo); good TOA = valid & detected
(unknown phase) & not χ²-flagged. Picture of the template fit and C(τ) drawn in the chat.
Tests: `testChannelTOA` 1–7, `testOptimalTOA` 1–10 (bit-identical references, explicit
covariance formulas to 1e-15, Monte Carlo pulls 1.006–1.04, H0 false alarms at design), plus
the validated runs (MC, SNR sweep, H0).

No new candidates. Open decisions from the deep analysis (13 Oct) that act in this stage, still
Jasper's: the χ² threshold (2 = 130σ on the channel path; dof-scaled option 1 + 5·√(2/dof));
items 3 / 4 / 8 (band taper on the pulsar only → `ChannelGain` 1 wrong in six edge channels,
the 2.2 % gap; receiver noise not band-limited; no best-achievable TOA error printed on
the channel path, item 8's open part). Links: R3, D2,
phase F (weighted radiometer fit, deferred), template exact by design (§12).

## Drivers – `pipelineParams.m`, `main.m`, run scripts (Jasper: "continue to the drivers")

Walk-through given (`pipelineParams` unchanged since `0499da7`, `main` since `6ad5b4c`; the
docs parameters page is current, stamped `8a756ec`, no change since): `pipelineParams` sets
truth (pulsar, simulation, noise/RFI) and observer values (receiver, processing, ephemeris) in
one script; five values feed both sides (T → f0/TRef, DM, dutycycle → template FWHM, refFreq,
fLow/fHigh), so today the observer's knowledge equals the truth exactly – each already recorded
(R3, D2, §12, item 3); picture drawn in the chat. `main`: runStage switches, disk estimate, sky
→ receiver (reuse checks `checkConsistency` / `checkSkyFiles`), bright-excision warning,
processing per front end (reuse check `checkProcessing` along the file chain), fold (reuse
check `checkFold`), TOA + detection + good TOAs + `validateTOA`; channel path prints no optimum
(item 8). Run scripts: `runMonteCarlo`, `runSNRSweep`, `runH0` full band only (warn since 13
Oct); `runRFITest`, `runLockedRadar` channel path, each with its own copy of the stage calls.

No new candidates. One connection proposed for Jasper's priority: §10 item 7 (one chain
function) + A3b (multi-seed Monte Carlo with blanking, "optional"): with one chain function the
statistical scripts get the default channel path (excision, 'optimal') for free; today the
default path's error-bar honesty and false-alarm rates are validated piecewise (noise model vs
data: `testFoldWeights` 3, `testBlankingWeights` 3; estimator vs drawn profiles:
`testOptimalTOA` 4, 8, 9) and by single-seed `main` runs, never end to end over seeds.

Noted (Jasper: "note the connection, then commit everything"): notes §10 item 7 (connection to
A3b) and item 11 (entries "Stage 11" and "Drivers"), roadmap row "Factor the chain into one
function". Then committed: the review notes of 14–16 Oct from both sessions and the text fixes
(B-text, V1, S1, FP1, BW1); `docs/api.html` regenerated from that commit afterwards.

## Next session (Jasper, 16 Oct: "let's continue that in a next session")

1. Reminder first: Jasper skipped parts of the walk-through and wants to read it fully
   (these logs 12–16 Oct, notes §10 item 11, docs stage pages; the chat pictures are only in
   the session – a compiled page on request). Ask which stages.
2. Open decisions: χ² threshold; full-review items 3, 4, 8; item 7 + A3b (priority vs C);
   parked periodic fix (P1–P4) and per-channel calibration.

## Open decisions, session 2 (16 Oct): χ² threshold → option C (`pChi2`)

Reminder given (read the walk-through fully; Jasper chose the decisions first). Order agreed:
χ² threshold now, then item 7 + A3b before C; items 3, 4, 8, P1–P4 and per-channel
calibration stay parked.

What is wrong with `MaxRedChi2` = 2 (explained): set on the full band (260 on-pulse dof,
σ 0.088 → 11σ); on the channel path 33,661 dof, σ 0.0077 → 130σ, so noise-model failures pass
(B5d 1.04 = 5σ, B5e 1.06 = 8σ, narrow carriers 1.22 = 28σ with error bars 3.4 % too small).
New point: √(2/dof) ignores the neighbouring-bin correlation of linear assignment (ρ ≈ ¼ →
scatter ~6 % larger). Design question underneath: a flag drops the TOA, so a sensitive
threshold would drop every B5e narrow-carrier TOA.

Jasper chose **C**: keep 2 as the drop rule, add a report. Done (code, test, docs; not yet
run):
- `functions/estimateTOA.m`: output `toa.pChi2` / `toa.total.pChi2` = P(χ² ≥ observed |
  template + modelled noise); local function `chi2Tail` (exact variance 2 Σ ρ² from the fit's
  lag covariances, moment-matched g·χ²_ν, `gammainc` upper tail); option `Chi2Alarm` (1e-3),
  summary line counting sub-ints below it. `redChi2`, `flagChi2`, good TOAs unchanged.
- `tests/testOptimalTOA.m` test 4: `pChi2` uniform on the 800 exact-covariance profiles per
  SNR level, 'equal' and 'optimal' (KS at 1e-3, D < 1.95/√n); header updated.
- Docs: `docs/toa.html` §4 new subsection "Noise-model check: pChi2"; notes §5.8 (new
  paragraph), §10 (decision marked in item 11 "Stage 11" and the deep-analysis list).
  `docs/api.html` to regenerate after the commit.

**Validated** (Jasper's runs, 16 Oct): `runAllTests("testOptimalTOA")` all passed (3.2 min); new
checks: pChi2 mean 0.498 / 0.525 / 0.500 / 0.504 (SNR 10 equal / optimal, SNR 40 equal /
optimal; 0.5 ± 0.010), fraction < 0.05 0.041–0.058, KS D 0.017–0.050 (< 0.069). SNR 10
optimal 2.5σ high, isolated (SNR 40 0.504, χ² 0.999): not pursued. `main` identical to the
8 / 13 Oct runs (8 TOAs, 2.956 µs, SNR 91.9, total +0.0564 ± 0.98 µs, ratio 0.825), no
P(χ²) alarm (χ² 1.008 ≈ 1σ). Committed.

## Item 7 + A3b: proposal revised (Jasper's answer pending, "I'll answer tomorrow")

First proposal (one chain function sky → TOAs) rejected by Jasper: the real IQ data enter
between generation and processing, so merging them is wrong; and the operational processing
will run on FPGA/GPU, not in MATLAB. Revised proposal: three roles, cut at `info_IQ`
(`main.m:112`; its fields are receiver settings only: file, fs, fLO, filter, t0):
1. simulate (truth; stand-in for antenna + ADC/DDC): generate → disperse → addNoiseAndRFI →
   applyIQmodulation, stays as calls in `main`;
2. process: new `processIQ(info_IQ, rx, ephem, proc)` → toa, detection, info (channelizeIQ …
   estimateTOA + detectPulsar, same intermediate files); the ground-truth rule enforced by the
   signature; the MATLAB golden model for the FPGA/GPU port; real IQ recordings can be fed in;
   removes the copies in runRFITest / runLockedRadar; multi-seed MC (A3b) = simulate(seed) →
   processIQ → evaluate;
3. evaluate (truth vs result): validateTOA, expectedPowerModel, check plots.
Open for Jasper: (1) boundary at IQ (recommended; a channel-data entry point later if the
FPGA does the channelizer) or at channels; (2) both front ends via `proc.frontEnd`
(recommended) or channel path only. Next step after his answer, read-only: the variables the
processing block of `main` (lines ~117–255) reads → proposed rx / ephem / proc fields.

## Docs: LaTeX equations (evening)

Jasper: "add LaTeX to the handbook; equations look cluttered and hard to tell apart; be complete
with the naming explanations near the equations". Plan agreed (download + physics.html first,
then the other pages one at a time; inline `<var>` math in running text left as is for now).
- KaTeX 0.16.11 (`npm pack`, files only: katex.min.js/.css, auto-render, woff2 fonts, LICENSE;
  600 KB) in `docs/vendor/katex/`; `nav.js` loads it only on pages that contain `\(` or `\[`
  (site still works offline from file://). Delimiters `\[ \]` and `\( \)` only, no `$`.
- `style.css`: KaTeX inside `.eq`; new `dl.where` (symbol | meaning, unit, code name).
- `physics.html`: all 12 display equations in LaTeX, each with a "where" list; the inline
  Isserlis covariance, ρ = A²/σn² and m(t) in LaTeX too. Content fixes on the way: the Isserlis
  sum now reads var(P̂)/P² (R is normalized, R(0) = 1, as in `powerCovariance`); the TOA list
  says f_ref is the *spin* frequency at t_ref (`fold.subint.fRef`), not the 1.6 GHz radio
  reference, and τ there is the template shift in turns, not the dispersion delay τ(f).
- `maintenance.html`: how to write an equation + changelog row.
Checked in the browser pane (local http server): 12 display + 91 inline renders, no KaTeX
errors, no sideways scroll at 375 px, light and dark. Not yet checked from file:// by Jasper.
Not committed yet (waiting for Jasper's look). Next: the other 12 pages with `.eq` blocks,
`formulas.html` last.

Committed (`c588d3d`). Then Jasper: "commit this and continue with all pages, double check the
things yourself". Done the same way on sky, receiver, fullband, channels, excision, fold, toa,
detection, validation, scenario (32 display equations; existing symbol bullet lists turned into
"where" lists) and `formulas.html` (all 43 rows LaTeX + a symbols line per row, names and
links unchanged, checked by script). Symbol meanings checked against the code help texts
(foldProfile, estimateTOA, channelizeIQ, detectRFI, periodicRFI, blankingWeights,
expectedPowerModel, addNoiseAndRFI, rfiSource). Clashes now stated where they occur: K (FFT
size vs dispersion constant), T (template vs period), Φ (normal CDF vs phase), L (frame length
vs lag), P (period vs probability), τ (shift vs lag), A/B/W in the sandwich (matrices).
Checks (browser pane, every page): 0 KaTeX errors, 0 unrendered delimiters, no equation wider
than its box at desktop width (two split over two lines: PFB formula, periodic-RFI test); at
375 px no page scrolls sideways (fix: KaTeX's hidden MathML escaped `.table-wrap`; `.eq` and
`.table-wrap` now position: relative), long equations scroll inside their box as before.
