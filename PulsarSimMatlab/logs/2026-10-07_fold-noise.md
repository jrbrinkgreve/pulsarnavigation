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

## Open items / next steps

1. (done) Jasper's run of `tests/testFoldWeights.m` (A1a).
2. (done) Jasper ran `tests/testFoldWeights.m` (tests 1–6): passed → **A1b validated**.
3. **A2**: `estimateTOA` / `detectPulsar`: all lags d ≤ D, per-channel weights and
   variance, own baseline a_c per channel, exclusion of partially covered channels.
   Design note first. Then A3, A4; B, C, D.
