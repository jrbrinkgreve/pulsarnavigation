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

## Open items / next steps

1. (done) Jasper's runs of `tests/testFoldWeights.m` (A1a, A1b) and
   `tests/testChannelTOA.m` (A2a).
2. (done) Jasper ran `tests/testChannelTOA.m` (tests 1–6): passed → **A2b validated**; committed.
3. **A3**: exact W/V/X streams from a voltage-level blanking mask (promote the mask
   convolutions of `tests/expBlankingVariance.m` to a function), blank voltages before
   per-channel dedispersion, end to end over several noise seeds: TOAs unbiased, error
   bars honest. Design note first. Then A4; B, C, D.
