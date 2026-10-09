# 2026-10-09 – Does a higher f_out improve detection? (meeting question + MATLAB sweep)

## Question (from Jasper's meeting)
"Increasing f_out improves detection: the matched filter / FFTFIT on the envelope averages
out more noise at a higher data rate while the signal stays the same."

## Answer
No, not in this pipeline, once the bins resolve the pulse. `detectPower` / `detectChannels`
average ALL samples into each bin (integrate-and-dump), so a 10× higher f_out gives 10× more
bins, each 10× noisier (variance m²/(B·dt)); the fold sums them back to the same noise. The
information is fixed by B·T (radiometer equation); f_out does not appear. Same statement in
the code: `powerCovariance`, V + 2ΣX = 1. Frequency view: the pulse lives below ~3 kHz,
the low-frequency noise level is set by B; higher f_out only adds noise-only frequencies
that FFTFIT does not weight.
Where it IS true: a detector that samples instead of integrating (e.g. a hardware power
detector + ADC without averaging): then S/N ∝ √f_out up to ~2× its video bandwidth.
Relevant for the small-antenna hardware: integrate before decimating. Also: bins wider
than ~W/5 smear the pulse (signal loss, not noise).
Python check (white band, weak pulse): integrating S/N 3.73 for f_out ≥ 10 kHz (the
radiometer limit); sample-only 0.37 at 10 kHz, 1.18 at 100 kHz, 3.73 at 1 MHz.

## MATLAB sweep (Claude ran it, Jasper asked)
Scratch script (not in the repo): `foutSweep.m` + `foutSummary.m` in the session
scratchpad; data files there too, `data/` only read (sky files, seed 42).
Channel path exactly as `main.m` (pipelineParams of today), but **excision off** (RFI is
off; and see the bug below). Per (SNR, noise seed 43–47) the front end ran once and every
f_out reused the SAME dedispersed channels → paired comparison. 5 seeds × 10 sub-ints.
nBin: 2048 and "matched" = 2^floor(log2(f_out·T)) (≥ 1 time bin per phase bin per turn).

−5 dB (pooled, matched nBin):

| f_out | bin | nBin | good | SNR | err-bar | rms err | ratio | red χ² | paired ΔTOA/σ |
|---|---|---|---|---|---|---|---|---|---|
| 2 kHz | 500 µs | 16 | 5/45 | 81.9 | 3.63 µs | 6.11 µs | 1.68 | 2.26 | 2.0 |
| 5 kHz | 200 µs | 32 | 45/45 | 91.3 | 2.32 | 4.72 | 2.04 | 1.16 | 1.2 |
| 10 kHz | 100 µs | 64 | 40/40 | 92.0 | 2.85 | 3.09 | 1.08 | 1.22 | 0.24 |
| 30 kHz | 33 µs | 256 | 40/40 | 91.5 | 2.96 | 3.04 | 1.03 | 1.06 | 0.13 |
| 100 kHz | 10 µs | 512 | 40/40 | 91.9 | 2.96 | 3.07 | 1.04 | 1.03 | 0.07 |
| 300 kHz | 3.4 µs | 2048 | 40/40 | 92.0 | 2.95 | 3.08 | 1.04 | 1.006 | 0.03 |
| 1 MHz | 0.96 µs | 2048 | 40/40 | 92.1 | 2.95 | 3.08 | 1.04 | 1.007 | ref |

−20 dB (weak signal; good TOAs = detected, χ² ok):

| f_out | nBin | SNR/sub-int | P_D known / unknown | good | rms err (good) | err-bar (good) |
|---|---|---|---|---|---|---|
| 2 kHz | 16 | 3.06 | 0.36 / 0.09 | 4 | 95 µs | 73 µs |
| 5 kHz | 32 | 3.49 | 0.60 / 0.18 | 8 | 73 | 50 |
| 10 kHz | 64 | 3.43 | 0.60 / 0.25 | 10 | 76 | 61 |
| 30 kHz | 256 | 3.47 | 0.57 / 0.25 | 10 | 75 | 63 |
| 100 kHz | 512 | 3.48 | 0.57 / 0.25 | 10 | 71 | 63 |
| 300 kHz | 2048 | 3.49 | 0.57 / 0.25 | 10 | 73 | 63 |
| 1 MHz | 2048 | 3.50 | 0.55 / 0.25 | 10 | 73 | 63 |

Conclusions:
- From ~30 kHz up nothing changes (SNR within 0.6 %, same TOAs to ≤ 0.13σ, same P_D):
  confirms the analysis. f_out is a resolution/cost knob, not a sensitivity knob.
- Below ~W/5 (W = 0.5 ms FWHM) it gets worse, mainly through the template: the Gaussian
  template is not smeared by the bin → red χ² rises and, at 200 µs bins, the error bars are
  2× too small (dishonest) although the SNR is still fine; at 500 µs bins SNR −11 % and
  most TOAs χ²-flagged.
- **Coverage trap: with nBin = 2048 and 1 turn per sub-int, every f_out ≤ 100 kHz gives
  0 valid TOAs** (MinCoverage 1 needs ≥ nBin time bins per sub-int: f_out ≥ nBin/(T·turns)
  ≈ 205 kHz here). `pipelineParams.m` in the working tree has f_out = 1e5 (uncommitted,
  not changed by Claude; committed value 1e6) → `main.m` would give no TOAs as it stands.
  Either back to 1e6, or nBin ≤ 512 with 1e5.

## Bug found (put on the agenda, not fixed)
`blankingWeights` crashes for long time bins (f_out 2 kHz with excision on): `fft(q, M)` of
a scalar q returns a row → `Kf .* Qf` is M×M. Fix `fft(q, M, 1)` (line 203). Added as §10
item 10 in `current_project_notes.md` and to the open-items table in `docs/roadmap.html`.

## Warning 1 done: `foldProfile:coverage` (Jasper: "start with warning 1")
Jasper asked for warnings on edge cases with critical performance impact. Two found:
(1) NBin too fine for the time bins → no TOAs at all (silent); (2) bins coarse vs the
pulse → error bars too small (silent). Agreed order: one at a time, (1) first.
- `functions/foldProfile.m`: check before folding, on the time grid only (sub-int 2, the
  first complete one, same assignment rule); warning with the limits (NBin ≤ time bins per
  turn, binDt ≤ 1/(f·NBin), or more turns). Help text of `'NBin'` extended. Fold outputs
  unchanged.
- `tests/testFoldWeights.m`: test 7 (5 cases: 10.08 µs / 2048 / 1 turn linear and nearest →
  warning, coverage 0.966 / 0.483; NBin 512, 0.96 µs, 10 turns → none). Claude's run:
  ALL PASSED (tests 1–6 still bit-identical). Jasper's run: passed → validated.
- Docs: `fold.html` (callout under "Choosing NBin"), `tests.html` (test 7, badge
  "Claude's run"), notes §5.7. `api.html`: `make_api.py` reads the committed code, so rerun
  it after committing (now only its commit hash changed).

## Follow-up: "a higher-rate matched filter IS less noisy" + professor's E/N0 point
Jasper: a matched filter / autocorrelation on a noisy signal sampled at 100 kHz is less noisy
than at 1 kHz; professor: the criterion is signal energy over noise PSD (detection theory),
and the square in the power should give a net gain beyond a linear system.
Experiment `tests/expSamplingRate.m` (added at Jasper's request; function with 4σ checks against
the predictions, saves `data/mc/expSamplingRate.mat`, one figure; ~20 s; Claude's run: ALL PASSED):
1. Known pulse (FWHM 5 ms, P 100 ms, 1 s) in white noise, three set-ups:
   A same σ per sample at every rate (the usual quick test) → SNR ∝ √fs (4.9 → 170): the
     noise PSD the filter sees, 2σ²/fs, drops 1000× from 1 kHz to 1 MHz;
   B fixed physical PSD, integrate-and-dump (= pipeline) → SNR 5.0 at every rate;
   C fixed PSD, keep one sample (aliasing) → ∝ √fs up to 5.0 at the analog rate.
   Measured = √(2E/N0) with the N0 the filter sees, in all three: E/N0 is the right
   criterion, and with the physical N0 fixed (T_sys) it does not depend on the rate.
   Same noise at 1 kHz and 100 kHz (case B): the two matched-filter outputs are the same
   curve.
2. Noise-like pulse (power 1 + a·g, B = 1 MHz), square law: i square → average (pipeline)
   5.2 at every f_out = radiometer bound a·√(B∫(g−ḡ)²dt); ii square → keep one per bin and
   iii average voltages → square both ∝ √f_out, reaching the bound only at f_out = B.
   Envelope-domain E/N0: E_env ∝ (signal power)², N_env = 2N²/B → d² = B∫(S/N)²dt: the
   square gives (S/N)² (the non-coherent loss vs a known waveform), and the number of
   independent cells is B·T, set by the pre-detection bandwidth, not f_out. The pipeline
   squares every IQ sample before binning, so it already sits at the bound.
   A higher rate only helps if it brings more RF bandwidth B (radiometer √B).

## Open items
- f_out: Jasper set it back to 1e6 (same day).
- Warning 2 done (below).
- Warning 2 (next, after Jasper's OK on 1): `estimateTOA:resolution` when
  max(binDt·f, 1/NBin) > FWHM/5 (template not smoothed by the bins → red. χ² up, error bars
  too small; sweep: 2× too small at FWHM/1.6). Later: a bin-smeared template.
- Unchanged: block B decision (B5e measure first vs per-channel noise calibration).

## Warning 2 done: `estimateTOA:resolution` (Jasper: "resume warning 2")
Calibration first (scratch sweep, f_out 1 MHz = fine time bins, nBin 24…2048, 5 seeds):
−5 dB rms error / error bar 1.14 (nBin 24, w/1.2), 1.16 (32), 1.06 (48–96), 1.05 (128),
1.04 (256, 2048); red. χ² 1.45 → 1.007; SNR unchanged (±0.5 %). −20 dB: no effect beyond
the scatter (10 good TOAs per point). Together with the f_out sweep (time bins w/2.5 → 2.04):
coarse TIME bins are the critical side; phase bins alone cost ≤ ~10 %. Rule kept as agreed:
warn when max(binDt·f0, 1/NBin) > w/5 (w = narrowest width at half maximum).
- `functions/estimateTOA.m`: check after the template normalization, local function
  `halfMaxWidth`, `info.binResolution` / `info.templateWidth`, help text (template input,
  info output).
- `tests/testChannelTOA.m` test 7 (six cases + template width, Gaussian and two-component).
  Claude's run: ALL PASSED (tests 1–6 unchanged, bit-identical). Jasper's run: passed → validated.
- Docs: `toa.html` (callout under "The template"), `tests.html` (row + item 7), notes §5.8.
  `api.html` after the commit (`make_api.py` reads committed code).
- Also added: `tests/expSamplingRate.m` (the E/N0 experiment above), registered in
  `tests.html` and notes §13.
