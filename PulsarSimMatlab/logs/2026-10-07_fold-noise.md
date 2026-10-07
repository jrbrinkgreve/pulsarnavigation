# 2026-10-07 – A1a: exact noise covariance in the fold

## State at start

- Everything through 5af516a committed (overview `2026-10-06_overview.md` written for this
  session). Channelized front end: `channelizeIQ`, `dedisperseChannels`, blanking
  experiment, `detectChannels` + `powerCovariance` validated. Next: A1 (fold with channels
  and exact noise); reference folds of the pre-A1 `foldProfile` saved in
  `data/mc/foldRef_pre3c.mat`.
- Jasper's local `main.m` change (`runStage.sky = true`) left alone, not committed. His
  run at 13:11 rewrote `test_rx_IQ`, `test_IQ_dedispersed` and `test_envelope`; same
  seeds → identical files (regression below bit-identical).
- Reminder given to read the overview first.

## Work done

- Physics explained and agreed: in a 3.125 MHz channel detected time bins (0.96 µs) are
  correlated (V 0.888, X(1) 0.046, …, Lmax 7). The covariance of phase-bin sums is
  cov(S_j, S_j+d) = σ²·ΣΣ a_kj a_k′,j+d c(k′ − k); with c = [1] it reduces to the old
  `weight2` (d = 0) and `weightX` (d = 1). Estimate before coding (10 ms / 2048 bins):
  phase-bin variance 0.987 of the old model, lag-1 correlation 0.258 vs 0.250, lag 2
  0.0009, lag 3 5e-6.
- A1 split (Jasper OK): **A1a** exact covariance (this session), **A1b** per-channel data
  weights (next).
- `functions/foldProfile.m`: option `NoiseCoeffs` = [V X(1..Lmax)] (default 1 = old
  behaviour); `weight2` = var/σ²; `weightX` [NBin × nSub × 1 × D], D = phase-bin lags
  (1 without lags → old shape); L = 0 terms times V; lag pairs L = 1..Lmax per chunk
  with a tail of Lmax time bins carried to the next chunk (local `addLagPairs`); pairs
  across sub-int boundaries skipped (documented); `info.noiseCoeffs`, `info.covLags`;
  memory estimate includes the D lag layers; the "too many lags" check only applies
  with lags (small NBin without lags must keep working).
- New `tests/testFoldWeights.m` (function): regression, algebra (Aᵀ·C·A), measured noise.
- Notes: §3 map, §5.7 (NoiseCoeffs + tests), §9 (narrow-channel item: fold part done;
  new item: covariance between sub-ints not stored), §10 labels/progress, §13 file list.

## Results (Claude's run, batch, 1.4 s)

- 1a–1d: the four reference folds reproduced bit for bit (fold and info).
- 2: weight2/weightX vs Aᵀ·C·A: 6.4e-16 (linear), 2.8e-16 (nearest), D = 4, nothing
  beyond lag D, zero lag coefficients = default bit for bit.
- 3: 128-channel seed-43 fold, 13,524 off-pulse bins × 128 channels, D = 3:
  variance 3.3402 ± 0.0038 vs new 3.3456 (−1.4σ) / old 3.3903 (−13.1σ);
  lag 1 0.8647 ± 0.0029 vs 0.8629 (+0.6σ) / 0.8477 (+5.9σ);
  lag 2 0.0028 ± 0.0029 vs 0.0030 (−0.1σ) / 0 (+1.0σ).
  New/old variance 0.9868, as estimated. Fold 0.26 s with lags vs 0.18 s without.

## Decisions

- Default `NoiseCoeffs = 1` and an explicit option (not read automatically from
  `info_det.noise`): the full-band path stays bit-identical and the A4 switch in
  `main.m` passes the coefficients explicitly.
- Lag dimension in dimension 4 of `weightX`, dimension 3 reserved for channels (A1b), so
  the layout does not change again.
- Covariance between consecutive sub-int profiles is not stored (only boundary bins at
  phase 0.5, negligible; notes §9).

Jasper ran `tests/testFoldWeights.m`: passed → **A1a validated**; committed (994ca4f).

## Later: A1b – data weights per channel

- Design agreed (Jasper: "go ahead"): blanking before dedispersion → per time bin and
  channel W (valid fraction, mean W·m), V, X(L) (exact from the mask, computed later in
  B). Fold: weight = Σ a·W (unbiased profile even for phase-correlated blanking),
  weight2/weightX from the per-bin V, X; per channel.
- `functions/foldProfile.m`: option `DataWeights` (info struct: file, nChan, N, Lmax,
  byteOrder; float32 [nChan × (2+Lmax) × N]); replaces `NoiseCoeffs` (error if both);
  weights [NBin × nSub × nW], weightX [NBin × nSub × nW × D], nW = nChan with weights
  else 1; one code path for scalars and per-bin arrays (helper `addTo`); `addLagPairs`
  takes per-bin X (tail carries it across chunks); chunks limited with a weight file;
  `info.dataWeights`, `nWeightChan`, `validFraction`.
- `tests/testFoldWeights.m`: tests 4–6 (constant stream = A1a fold bit for bit; random
  streams vs AᵀW and AᵀCA; fake blanking of whole detected bins in the 128-channel data:
  mean unbiased, noise as predicted). Temporary files (~0.5 GB) deleted at the end.
- Results (Claude's run, ~7 s, all PASS; numbers in notes §5.7): 4 bit-identical; 5
  errors ≤ 7e-16; 6 mean −0.0019 ± 0.0009 (window) / +0.0001 ± 0.0001 (rest) vs naive
  −0.574 / −0.150; noise all off-pulse −1.8σ / −0.0σ / −0.1σ (variance / lag 1 / lag 2),
  window −0.8σ / +0.6σ. Same noise realization as test 3 → the slightly negative variance
  offsets are correlated. Fold with weights 1.8 s vs 0.25 s.
- Noted (§5.7): `estimateTOA` / `detectPulsar` must not get a per-channel fold before A2
  (`fold.weight(:, s)` would take channel 1).

## Later: questions, A1b committed, A2a (estimateTOA combines channels)

- Jasper asked again what A1–A4 and A–E mean → answered in plain terms (A1 fold, A2
  detection + TOA fit with channels, A3 end-to-end blanking test, A4 switch in main;
  A front end, B excision, C fast simulator, D Earth's motion, E later). Asked whether B
  can be finished today: no — ~11–13 units left (A2 2, A3 2, A4 1, B ~6–8), about a week
  of sessions at today's pace; today realistic: A2.
- A1b committed (918cae3).
- A2 design agreed: p = Σ_c prof_c (equal weights); noise summed over independent
  channels, all lags d ≤ D, per-channel level m_c = a_c + b_c·template (same Fourier
  projection, τ fixed); exclusion rule; `Bnoise` = noise bandwidth of one channel of the
  fold; two units A2a `estimateTOA`, A2b `detectPulsar`. All current callers use
  single-channel folds (testDedisperseChannels sums the channel power first) → nothing
  existing changes.
- `tests/makeToaReference.m` (new, run before the edit): pre-A2 outputs of estimateTOA /
  detectPulsar on the frozen folds of `foldRef_pre3c.mat` → `data/mc/toaRef_preA2.mat`
  (A 2.880 µs as known). MinCoverage case set to 0.3 (the partial last sub-int has 37 %
  coverage; 0.5 did not exercise the gap filling).
- `functions/estimateTOA.m`: helper `combineChannels` (exclusion rule, sum, weights),
  `fitOne(ch, c)` with per-channel a_c, b_c, noise summed over channels and lags,
  'offpulse' lags 1..D, `nChanUsed`; header updated (Bnoise per channel).
- `tests/testChannelTOA.m` (new): Claude's run ALL PASSED, ~4 s (numbers in notes §5.8):
  full band bit-identical (5 cases); error bars = explicit covariance to 4e-16 (128
  channels with NoiseCoeffs; fake blanking with exclusion: 118 / 128 channels in odd /
  even sub-ints); per-channel vs summed power: TOAs equal to 5.7e-7 σ, error ratio 0.9996,
  red. χ² 0.983 vs 0.971; blanked vs unblanked: mean −0.20 ± 0.50 µs, error ratio 1.041.
  Observation: rms difference 1.33 µs > nested-estimator guess 0.8 µs (not exact for a
  non-optimal estimator; 8 TOAs) → A3.
- Tolerance choice: per-channel vs summed TOAs pass at < 1e-4 σ because the float32
  summed file alone moves TOAs by ~1e-6 σ (measured 5.7e-7).

## Later: A2b (detectPulsar combines channels)

- A2a committed (ce9babb). A2b note agreed: H0 level per channel = its baseline a_c;
  all lags in T0 and Tmax; `combineChannels` moved to its own file (shared).
- `functions/combineChannels.m` (new, moved out of `estimateTOA.m` unchanged; tests 1–3
  still pass identically). `functions/detectPulsar.m`: `testOne(ch, c)` with
  per-channel H0 noise, lag loops in T0 / Tmax (template products per lag precomputed),
  `Bnoise` per channel, `nChanUsed`; header updated.
- `tests/testChannelTOA.m` tests 4–6 (Claude's run, ALL PASSED, ~5 s; numbers in notes
  §5.13): detection regression bit-identical (5 cases); T0 / normProfile / Tmax =
  explicit H0 covariance to ≤ 4e-16; noise ratio channel path 0.9973 ± 0.0061 (blanked
  1.0049 ± 0.0067), summed path 0.9842, ratio 1.0133 (fold level 1.0134); same phase and
  detections.

## Later: A3a (blankingWeights) and a fold bug found with it

- A2b committed (d5afbab). A3 design agreed (interval mask, exact W/V/X from the mask
  and the channel filter, white-noise approximation stated; A3b multi-seed MC optional);
  Jasper: "go ahead. test thoroughly".
- New `functions/blankingWeights.m` (§5.18): h per channel via an impulse through
  `applyInverseDispersion` (same Nfft), one FFT convolution of the mask per sample lag,
  bin sums for V and X(1..Lmax), MinWeight 1e-6, unblanked channels from their own h.
  Two indexing errors fixed before the first run (3-D shape for the unblanked
  constants; empty-partner mask shifted the wrong way). Smoke run 2.3 s.
- New `tests/testBlankingWeights.m` (exactness vs direct covariance matrix, consistency,
  end to end with voltage blanking before dedispersion). First run stopped on a test bug
  (row vs column → 7×7 broadcast); then ALL PASSED, but two hints were followed up:
  - noise ratio 0.987 (blanked) vs 0.997 (unblanked): diagnosis per sub-int — the
    unblanked fold with the same 121 channels gives 0.991; blanking itself −0.003 ±
    0.0025 → statistics of this realization; test now has the paired comparison.
  - off-pulse mean +2.3σ: grouped by blanking type → channels 110–115 (100 µs blanks
    every 1.1 ms) +4.9e-3 ± 0.7e-3 (7σ), hidden in the all-channel average. Detected
    level: E[P] = W·m holds for W ≥ 0.1; nearly empty bins (W 1e-6…1e-2) ~4 % low (the
    white-noise approximation, as expected). Fold level: phase bins with valid fraction
    < 0.001 at +285 % → cause: bins set to W = 0 by MinWeight still added their leftover
    power (effective W 3.9e-7) to the fold sum. Confirmed by refolding with that power
    removed: −1.6e-4 ± 0.9e-4 (all channels −0.8e-5 ± 1.9e-5).
- Fix (Jasper OK): `foldProfile` with DataWeights: X(W == 0) = 0 before accumulating
  (header note); `testFoldWeights` test 5 checks sum = Aᵀ·(P where W > 0) (exact).
  `testBlankingWeights`: means per blanking type and per valid-fraction group,
  inverse-variance weighted; noise deviations as a fraction of the variance; paired
  noise ratio. All three tests in one session: ALL PASSED (53 s).
- Note: the TOAs of this test were not affected by the bug (channels 110–115 were
  excluded in every sub-int anyway); with longer sub-ints they would have been.
- Jasper asked to fix the other two findings too. They are model limits, not code
  bugs → explained with options (notes §9): (1) white-noise approximation — weak-signal
  requirement (phase-dependent mean errors ≪ 4e-6 of the baseline); exact mean with the
  true channel spectrum possible; (2) exclusion rule / equal channel weights → weighted
  multi-channel fit. Waiting for Jasper's choice.

## Later: the two A3a findings — (1) exact mean with the true channel spectrum

- Jasper: "go ahead with (1) and (2)". (1) first; (2) gets a design note first (core
  change to two validated functions).
- Channel spectrum from `info_chan.prototype` (|response|² aliased at the channel rate):
  flat ±1.6e-4 in the band, 0.5 at Nyquist; R_x significant to lag ~19 (< 1e-8 beyond).
- `blankingWeights`: new first input `info_chan`; local `channelAutocorr`; mean
  E|y|² = Σ_d R_x(d)(q_d ∗ keep·keep_d), normalized by the unblanked value (reuses the q
  FFTs of the lag loop, one extra FFT per d); option `InputSpectrum` 'channelizer' /
  'white'; `info.Rx`. V, X stay white (error bars only).
- `tests/testBlankingWeights.m`: new signature; test 1 W via u'·R·u (independent matrix
  form); new 2e (prototype autocorrelation vs measured on raw channel samples: 7e-4,
  white off by 0.05); (a0) per W group with the white-input W for contrast. ALL PASSED,
  66 s.
- Correction: the true spectrum changes W of nearly empty bins by ~1 % only; the
  earlier "~4 % deficit" was over-stated — with the honest channel-scatter σ it is
  −4 ± 4 % (bins inside a blank share their data; the bin-based σ of the diagnosis was
  too small). Notes §5.18 / §9 corrected.

## Later: (2) A2c-1 — weighted multi-channel TOA fit

- Jasper: "losing track … re-evaluate it once more, and if everything is right, go
  ahead". Re-evaluated: design holds (known gains s_c instead of free b_c per channel —
  weak signal; diagonal weights 1/var with exact sandwich error bars; 'equal' stays
  default); refinements: iterate weights to convergence; equivalence check via the
  'offpulse' model (uniform weights); bounded search + Newton polish.
- `functions/estimateTOA.m`: options `Weighting` ('equal' | 'optimal'), `ChannelGain`,
  `MaxIterations`; output `chanBaseline`; local `fitOptimal`, `wlsFit`, `objective`,
  `slope`, `chanNoise`, `quadForm`; 'equal' path untouched (bit-identical).
- First smoke run: 'offpulse' equivalence only to 1e-10 turns (search on values) → Newton
  polish on the analytic slope → 1.7e-16.
- New `tests/testOptimalTOA.m`: ALL PASSED (127 s); numbers in notes §5.8. Monte Carlo
  with the exact covariance of 64 fake-blanked channels: unbiased, pulls 1.006 (SNR 40) /
  1.04 (SNR 10, both estimators — known low-SNR effect), scatter optimal/equal 0.955–0.963
  as predicted by the error bars; on the blanked data all 128 channels used.
- Memory: MEMORY.md had been changed outside the session (overview reminder removed,
  report time 06:30) — taken as current.

## Later: reminder, Jasper's test runs, commit

- Jasper asked for a reminder tomorrow → one-time scheduled task
  `pulsarnav-continue-reminder` (8 Oct 09:00) + a memory note for the next session.
- Jasper ran `testFoldWeights`, `testChannelTOA`, `testBlankingWeights`, `testOptimalTOA`:
  all passed → **A3a (blankingWeights, exact mean), the fold rule (W = 0 adds no power) and
  A2c-1 (estimateTOA 'optimal') validated**; committed.

## Open items / next steps

1. A2c-2: `detectPulsar` with the weighted matched filter (same quantities N, Dn; H0
   noise with the exact covariance). Design note first.
2. A4 (switch in main.m: channel path, 'optimal'); A3b multi-seed MC optional; B, C, D.
