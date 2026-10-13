# 2026-10-13 – deeper analysis (while Jasper studies the theory)

Jasper (13 Oct, after the read-through fixes were committed in `1b65f5a`): "note the current
state down, and keep analysing further and deeper while i study".

The state at the start is in `2026-10-13_full-review.md`:
- items 1, 2, 6, 7, 9 done;
- 5 explained;
- 8 in part;
- open: 3, 4, 5 (guard), 8 (channel optimum).

Read-only throughout: no code, notes or docs changed. The numbers come from independent
Python/numpy re-computations (not the MATLAB code) and the saved info files in `data/`; no
MATLAB runs.

## Summary (most relevant first)

1. **Item 3 explains the whole 2.2 % channel-vs-full-band gap** (section 2). The pulsar-only
   band-edge taper costs 1.1 % through the channel tapers, plus 1.0 % from `ChannelGain` = 1
   in the six edge channels. Predicted +2.13 %, measured +2.19 %. With a realistic receiver
   it would be 0.37 %. `docs/limits.html` ("< 1 %") is wrong for the default path.
2. **Default excision blanks a bright pulsar** (section 4): on the pulse peak 10 % of
   windows at 0 dB, 94 % at +5 dB, plus the frequency guard from +10 dB; noise-free runs are
   the extreme case. Weak target unaffected. Options: a warning in main, a help-text limit.
3. **The χ² flag is nearly blind on the channel path** (section 12): threshold 2 = 130σ
   there; B5e's 3.4 %-too-small error bars (χ² 1.22, 28σ) pass as good.
4. **"−54 dB" carries a hidden ~3 dB margin** (section 8): J0437's worst case is −51 dB in
   the simulator's convention. The simulator at −54 dB needs ~62 h per TOA, not 20 h.
5. **Docs errors** (section 7): the Rice factor is read as "independent trials" on three
   pages (5.5 vs the true ~60 at η 4.15); LTE "(planned)"; 800 MHz "= the bandwidth"; the
   roadmap's "commit pending" (section 10).
6. **Item 5 settled**: no usable guard recovers turn 0 (section 3); keep it.
7. **Small code finding**: the optimal fit's convergence flag is dropped (section 9).
8. **Confirmed independently** (sections 1, 5, 11): every key constant, the formulas page,
   the ground-truth rule, no timing offset (−0.21 / −0.05 µs), full determinism. The design
   rule's 373·ρ is 360·ρ on the default path (~7 % more turns).

Interaction: fixing item 3/4 with a realistic receiver makes the six edge channels non-flat,
which is the parked per-channel calibration problem (section 6).

## 1. Key constants recomputed independently: all confirmed

- **`detectRFI` thresholds at P_FA 1e-6** (from the real channelizer prototype in
  `test_rx_IQ_chan_info.mat`; ρ(1) = 0.050, eigenvalues and phase-type tail via `expm`):

  | Window (samples) | η |
  |---|---|
  | 1 | 13.816 |
  | 2 | 8.391 |
  | 4 | 5.387 |
  | 8 | 3.681 |
  | 16 | 2.686 |

  Docs and help text: 13.82, 8.39, 5.39, 3.68, 2.69.
- **`powerCovariance`** (a continuous-spectrum integral with direct lag sums, not the code's
  FFT grid): V = 0.8877, X(1) = 0.0458, X(2) = 0.0051, Lmax = 7, captured 0.9984,
  Bnoise 3.0596 MHz per channel. Identical to the values saved in
  `test_dedispB_chan_power_info.mat`.
- **Detection thresholds** (template: nBin 2048, FWHM 0.05 turns): Rice factor 5.512;
  η 3.090 (known phase) and 4.155 (unknown phase).

  Side note, not a bug: the textbook continuous-Gaussian value 1/(2π√2 σ) is 5.30. The code
  is right for its own filter, which is zero-mean with harmonics k ≥ 1 only. Dropping the
  DC term raises the mean of (2πk)² (the spectrum is only ~5 harmonics wide), so 5.51 is
  correct.
- **Full-band optimum at −5 dB** (both tapers exact): SNR per pulse 96.53, σ_TOA 2.8920 µs
  (docs: 96.5, 2.892). The flat-band prediction (97.93, 2.8504 µs) is 1.45 % optimistic
  (docs: "~1.5 %").

## 2. The channel path's optimum, and why it is 2.2 % behind the full band (items 3 and 8)

The channel path's own optimum (item 8's missing number) is the Fisher information summed
over the 128 channels. Each channel has its own sin² edges of 62.5 kHz, the pulsar also has
the forward taper, and the channels are combined with optimal weights. At −5 dB:
**σ_TOA = 2.9236 µs, SNR per pulse 95.50**, 1.09 % above the full band's 2.8920 µs.

The fit also assumes `ChannelGain` = 1 in every channel. But the forward taper leaves the
pulse at these fractions of each edge channel's own noise level (the pulse gain g_c):

| Channels | g_c |
|---|---|
| 1 and 128 | 0.023 |
| 2 and 127 | 0.417 |
| 3 and 126 | 0.960 |
| all others | 1 |

So the fit uses wrong weights in six channels: SNR efficiency (Σg)/√(N·Σg²) = 0.9898, a
1.0 % loss. Expected channel-path error bar: 2.9236 / 0.9898 = 2.954 µs.

| Model | Full band | Channels | Channels with s_c = 1 | Gap to full band |
|---|---|---|---|---|
| as simulated (forward taper on the pulsar only) | 2.8920 µs | 2.9236 µs | 2.9536 µs | +2.13 % |
| no forward taper | 2.8807 µs | 2.8807 µs | 2.8807 µs | 0 |
| real receiver (band taper on pulsar **and** noise) | 2.8882 µs | 2.8989 µs | 2.8989 µs | +0.37 % |
| measured, seed 43 (`main`, median error bar) | 2.880 µs | 2.943 µs | | +2.19 % (with excision 2.956 µs, +2.64 %) |

**Finding.** The 2.2 % gap between the paths in every printed comparison (docs: "full band
2.880 µs, channel path 2.943 µs") is almost entirely item 3: the simulator tapers only the
pulsar at the band edges.
- In the full band the forward and inverse tapers coincide (the same 8 MHz edges).
- On the channel path the pulsar sees the band-edge taper and 256 channel edges, while the
  noise sees only the channel edges.
- The six edge channels then hold less pulse than the equal-gain fit assumes.

With a realistic receiver, where the bandpass shapes noise and pulsar alike, the channel path
would cost only 0.37 % (its 256 extra edges), and `ChannelGain` = 1 would be exactly right.
The model predicts the measured gap (+2.13 % vs +2.19 %), so the channel path itself is
working as well as its model allows.

Also: main's error bars sit within 1 % of the respective optima: full band 0.996, channels
1.007, channels with excision 1.011.

**What was known before.** Notes §5.16 (`testDedisperseChannels` test 3) already says that
in the simulation the pulsar also passed the forward 8 MHz taper, so the channel tapers cost
2.1 % of the signal, and that on real data both paths collect the same 390 MHz.

What is new here:
- the further ~1 % from `ChannelGain` = 1 in the six edge channels (the 'optimal' fit);
- that the two together predict main's measured 2.2 % gap;
- the real-receiver case (0.37 %);
- that `docs/limits.html` ("< 1 % in SNR") and the channels page's result section ("one TOA
  fewer because …", silent on the 2.2 % larger error bars) do not reflect this.

**Consequence for decisions.** Item 3 is not a cosmetic modelling note. It is what makes
the default path look worse than the reference path, and it biases any "channels vs full
band" comparison by about 2 %. The fix would be noise shaped by the same W (the real-receiver
row), or no forward taper on the pulsar. It is a sky/receiver modelling change, so its scope
and the bit-for-bit regressions (`data/test/`) need Jasper's decision first.

## 3. Item 5: can a shorter per-channel guard recover turn 0? No

The per-channel dedispersion filter of `applyInverseDispersion` was rebuilt in numpy for
channels 1, 64 and 128 (Nfft 2^20, 62.5 kHz sin² edges, chirp to 1.6 GHz). Leakage = filter
energy outside the lags [−nFuture, nPast] that overlap-save keeps:

| Guard (× 1/edge) | Time | Samples | Leakage, ch 1 | ch 64 | ch 128 |
|---|---|---|---|---|---|
| 20 (current) | 320 µs | 1334 | 1.3e-12 | 1.5e-12 | 3.3e-12 |
| 10 | 160 µs | 667 | 3.6e-11 | 4.0e-11 | 9.0e-11 |
| 5 | 80 µs | 334 | 9.7e-10 | 1.1e-9 | 2.4e-9 |
| 3 | 48 µs | 200 | 1.2e-8 | 1.3e-8 | 2.9e-8 |
| 2 | 32 µs | 134 | 9.0e-8 | 9.8e-8 | 2.2e-7 |
| 1 | 16 µs | 67 | 3.8e-6 | 4.4e-6 | 1.0e-5 |
| 0.5 | 8 µs | 34 | 1.3e-4 | 1.5e-4 | 3.3e-4 |
| 0.25 | 4 µs | 17 | 6.2e-4 | 7.0e-4 | 1.6e-3 |

- **Fall-off.** Leakage falls about as guard⁻⁵: the sin² edge has a continuous first
  derivative, so the ringing decays as 1/t³. The notes' "~1e-12" at 20/edge is confirmed.
- **What turn 0 needs.** Recovering it needs data from within about one phase bin
  (4.88 µs) of t = 0. The channelizer's own support takes 10 samples (2.4 µs), which leaves
  a guard of ≲ 2.5 µs (0.16/edge): leakage around 1e-3, far above `MaxLeakage` 1e-6. So no
  usable guard gets turn 0 back.
- **The smallest guard within 1e-6** is about 1.5/edge (24 µs). It would shorten the
  unsupported start and end from 0.32 ms to 0.024 ms each, which saves nothing that matters.
- **Conclusion for item 5.** Keep the guard; document the effect (done in `1b65f5a`). The
  lost turn is an artefact of the synthetic file starting exactly on a turn boundary
  (TRef = T/2, file start t = 0). With real data the first and last partial turns are lost
  on any path.

## 4. New: the default excision eats a bright pulsar (and a noise-free one)

`detectRFI` flags windows whose mean power exceeds η_n times the channel baseline. The pulsar
is noise-like, so on the pulse the sample power is exponential with mean 1 + ρp(t). Exact
flag probability per window at the pulse **peak** (same eigenvalue / phase-type method as
section 1; baseline 1):

| SNRdB | ρ | n = 1 | n = 2 | n = 4 | n = 8 | n = 16 | one sample > 100× (frequency guard) |
|---|---|---|---|---|---|---|---|
| −25 | 0.003 | 1.0e-6 | 1.1e-6 | 1.1e-6 | 1.1e-6 | 1.1e-6 | 5e-44 |
| −10 | 0.100 | 3.5e-6 | 4.1e-6 | 5.2e-6 | 7.3e-6 | 1.2e-5 | 3e-40 |
| −5 | 0.316 | 2.8e-5 | 4.2e-5 | 7.5e-5 | 1.7e-4 | 5.1e-4 | 1e-33 |
| −3 | 0.501 | 1.0e-4 | 1.8e-4 | 3.8e-4 | 1.1e-3 | 4.2e-3 | 1e-29 |
| 0 | 1.0 | 1.0e-3 | 2.2e-3 | 6.0e-3 | 2.2e-2 | 9.5e-2 | 2e-22 |
| +3 | 2.0 | 9.9e-3 | 2.4e-2 | 7.3e-2 | 0.24 | 0.63 | 3e-15 |
| +5 | 3.2 | 3.6e-2 | 8.9e-2 | 0.24 | 0.59 | 0.94 | 4e-11 |
| +10 | 10 | 0.28 | 0.55 | 0.86 | 0.99 | 1.0 | 1.1e-4 |
| +20 | 100 | 0.87 | 0.99 | 1.0 | 1.0 | 1.0 | 0.37 |

- **Weak pulsars.** At the target (≤ −25 dB) the pulsar cannot trigger anything, as
  `detectRFI`'s help text says. At −5 dB the known −0.4 % amplitude effect (notes §9)
  matches these rates.
- **Bright pulsars.** From about −3 dB the pulse blanks itself: 10 % of the 16-sample
  windows on the peak at 0 dB, 63 % at +3 dB, almost all from +5 dB. From +10 dB single
  on-pulse samples exceed 100× the baseline, so the frequency guard also blanks ±7
  channels around them.
- **Why the weights cannot repair it.** The blanking depends on the data (it removes the
  highest samples), but `blankingWeights` assumes a mask independent of the data. The
  profile peak is therefore suppressed and flattened. The TOA stays unbiased by symmetry,
  but the amplitude, SNR and χ² are wrong, and at high SNR the data disappear.
- **Noise-free (`snrDB = Inf`)** is the extreme case. The baseline is the median of
  pulsar-only power, ~0 off the pulse, so every on-pulse window is flagged. Nothing in main
  stops a noise-free or a bright run on the default channel path with `excision = true`.
  The full band has no excision, so it is not affected.
- **Where it would bite.** A noise-free or bright check of the default path; and porting
  `runSNRSweep` (−25…+20 dB) to the channel path (review item 22), where the points
  ≥ 0 dB would be destroyed by excision, not by the estimator.

**Options (not implemented, for Jasper):**
- (a) A warning in `main.m` when `excision` is on and `snrDB` > −5 dB or is Inf. This is a
  driver-level check that uses a simulation setting, like the reuse checks.
- (b) Say so in `detectRFI`'s help and on the excision page: valid for ρ ≲ 0.3.
- (c) Later, an ephemeris-aware detector, e.g. a baseline from off-pulse phases only, or
  higher thresholds on-pulse. That trades RFI sensitivity on the pulse.

(a) + (b) are cheap; (c) is design work and irrelevant for the weak-signal target.

## 5. Ground-truth rule and timing chain: clean

- **Ground-truth audit (grep of every processing function).** Searched for `info_gen`,
  `info_rx`, `info_disp`, `pulseCenters`, `snrDB`, `noiseStd`, `dutycycle`,
  `envelopeMode`, `kernelEnergy` and `SignalInfo` in `applyInverseDispersion`,
  `detectPower`, `channelizeIQ`, `detectRFI`, `periodicRFI`, `blankChannels`,
  `dedisperseChannels`, `detectChannels`, `powerCovariance`, `blankingWeights`,
  `foldProfile`, `estimateTOA`, `detectPulsar`, `combineChannels`, `gaussianTemplate`,
  `noiseBandwidth`. The only hit is a comment in `foldProfile`'s help. In `main.m`, truth
  reaches only the check functions (`plotIQCheck`, `plotDetectedPower`, `plotFoldCheck`,
  `expectedPowerModel`, `validateTOA`). The processing uses `ephem` and receiver settings
  (`info_IQ.fLO` / `.actualFsOut` / `.t0`, the channelizer prototype). The known exceptions
  are deliberate: R3 (`ephem.TRef` = the true epoch) and `ephem.profileFWHM` = the true duty
  cycle (catalogue knowledge).
- **Timing chain.** Every deliberate delay is removed by an integer sample count:
  - dispersion `s0`;
  - IQ group delay G = 1024 input samples;
  - channelizer G = 1926 input samples;
  - dedispersion: refFreq unshifted, chirp phase exact per bin;
  - detection: bin centroid time stamps;
  - all tapers real, i.e. zero phase.

  An independent numpy template fit (FFTFIT with a bounded search; my own code, not
  `estimateTOA`) of the saved total folds gives the following total-fold offset:

  | Fold | Offset | Rough σ (white-noise) |
  |---|---|---|
  | `test_fold.mat` (full band, 8 Oct) | −0.21 µs | ±0.6 µs |
  | `test_fold_chan.mat` (channels with excision, 13 Oct) | −0.05 µs | ±0.6 µs |

  Both are consistent with zero (`main`'s σ ~1 µs). The full-band value matches the notes'
  −0.21 µs mean error. No systematic offset at the sub-µs level, by an independent method.
- **Precision for long files** (no action). Absolute-index phases (`mod(n·r, 1)`) in double
  keep the LO and carrier phases to ≤ 2e-4 cycles up to L ~ 1000 s. The channelizer phase
  `mod(kb·c, K)` is exact (half-integer × integer < 2^53). The fold phase f0·t stays below
  1e-9 turns for a day of data. The simulator's other limits (disk, whole-file excision)
  come long before these.

## 6. Interaction to keep in mind: fixing items 3/4 touches the parked edge-channel problem

Today the receiver noise is white everywhere (item 4), also across the band edges and inside
the channelizer's transition regions. That keeps two models exactly self-consistent in
channels 1 and 128 too:
- `detectRFI`'s thresholds (white noise through the prototype, ρ(d));
- `blankingWeights`' channel spectrum R_x.

If the receiver bandpass (8 MHz sin² edges) shaped the noise, as a real receiver does (the
fix of items 3 and 4), channels 1–3 and 126–128 would no longer be flat. Channel 1's raw
spectrum would rise from 0 to ~0.11 across its 3.125 MHz. They would become "partly covered"
channels, the case of B5d/B5e (LTE edge channels: false flags, χ² 1.04, error bars 0.4–3 %
too small) that the parked per-channel noise calibration (item 33,
`2026-10-09_handover-B-closed.md`) was designed for.

So a realistic receiver bandpass should come together with per-channel calibration, or the
six edge channels should be excluded or given their measured spectrum. That is one reason not
to treat items 3/4 as a quick fix.

## 7. Docs physics pages: two errors found

- **The "Rice factor" is misread as a number of trials (3 pages).** `docs/detection.html:43`
  ("roughly the number of independent phase trials per turn"), `docs/glossary.html:56`
  (same) and `docs/physics.html:191` ("η = 4.15 … about 5.5 effective independent trials
  per turn") all say this.
  - **What 5.51 is.** √λ₂/2π = 5.51 is the mean number of up-crossings of the *mean level*
    per turn.
  - **The effective trial count.** At a threshold η it is N ≈ P_FA/Q(η) ≈ 5.51·η·√(2π).
    At η = 4.155: Q = 1.63e-5, the Rice term is 9.8e-4, so **N ≈ 61**. Check: 5.5
    independent trials would need η = 3.565 at P_FA 1e-3; 60 trials give η = 4.149.
  - **Scope.** The code and the threshold are right (`runH0` validated them); only the
    interpretation in the text is wrong. Text fix: "the Rice factor 5.51 is the up-crossing
    rate of the mean level per turn; at η = 4.15 the phase search acts like ~60
    independent trials".
- **`docs/physics.html:204`** lists the noise-like RFI type (LTE / OFDM) as "(planned)". It
  exists since B5d, 8 Oct (`rfiSource('noise', ...)`). The table should list five types,
  matching the page's own "five types".
- **Small inaccuracy, `docs/physics.html` §4.** "Complex samples keep the phase … at a sample
  rate equal to the bandwidth of interest (800 MHz here …)". The bandwidth of interest is
  400 MHz. 800 MHz is needed only because the LO sits off-centre (band at −100…+300 MHz) and
  the IQ low-pass needs a margin (cutoff 360 MHz).

**Everything else on the physics page checked and right:**
- delays 8.103 / 14.406 ms and the 6.30 ms smear at DM 5;
- B_n 391.6 MHz;
- 388ρ;
- σ_TOA ≈ √2·σ_t/SNR;
- 1 µs ≈ 300 m;
- the Earth-motion drift ~100 µs/s (v/c = 1e-4);
- η = 3.09 / 4.15;
- the "80,000 pulsar pulses" per radar pulse (computed: ~94,000 for a 2 µs pulse at
  INR +20 dB vs a −54 dB pulse of 0.5 ms FWHM; same order).

**`docs/formulas.html`** (all entries re-derived): no errors. Checked:
- the fold's 1 : 2/3 : 1/6 weights and correlation 0.25;
- the baseline loss 0.962;
- the radar energy PW − RiseTime/4;
- B_n 391.6 / 389.6, mean-power integrals 390.0 / 388.4;
- the ionospheric 0.69 ns per TECU at 1.4 GHz.
- **`docs/limits.html`, "Signal and receiver"** says "Receiver noise does not pass the
  dispersion taper in the simulation: < 1 % in SNR". That is true for the full band:
  +0.39 % in σ_TOA against a pulsar without the forward taper (section 2: 2.8920 vs
  2.8807 µs). It is not true for the default channel path, where the same effect costs
  **+2.5 %** (2.9536 vs 2.8807 µs, with `ChannelGain` = 1). That is item 3, which the
  page's sentence should point to.
- **`docs/limits.html`, "Edges and scale"**: "Disk: ~61 GB per second of simulated data"
  holds for the full band (60.8 GB/s). The default channel path writes 63.4 GB/s, or
  72.3 GB/s with excision (blanked copies and the weight file; `main`'s own estimate gives
  the same). Small; the page could give both.
- Also checked on the limits page, right: the f^−1.6 pulsar spectrum costs 0.88 % with
  equal weights (page: "~1 %"); 1e-9 clock rate over 20 h = 72 µs; the 6.3 ms unsupported
  end; the B5d/B5e numbers.
- Glossary: only the Rice-factor entry is wrong (above). The rest checked: B_n values, 388ρ,
  the good-TOA definition, the TOA convention, NBin 2048, LO 1.3 GHz, frequency guard
  ≥ 100× and ±7.

## 8. The "−54 dB" target is in the scenario's ρ convention, not the simulator's

`docs/scenario.html` defines ρ = S_peak/SEFD with S_peak = S_mean·P/W_eq and
W_eq = 2·W50, which it calls a "margin for wings and components". The simulator's `snrDB` is
the peak of a Gaussian of FWHM W50: S_peak = S_mean·P/(σ√(2π)) = S_mean·P/(1.0645·W50). For
the same pulsar the scenario's ρ is therefore 1.88× (2.7 dB) lower. Recomputed for J0437
(S1400 150 mJy, P 5.76 ms, W50 0.141 ms, dual pol, 400 MHz, SNR 7):

| Design point | SEFD | Scenario ρ, time | Simulator-convention ρ | Matched-filter time | With the 0.962 baseline loss |
|---|---|---|---|---|---|
| nominal | 3.17e5 Jy | −50.1 dB, 3.9 h | −47.4 dB | 2.8 h | 3.0 h |
| worst | 7.24e5 Jy | −53.7 dB, 20.4 h | **−51.0 dB** | 14.6 h | 15.8 h |

The scenario's own numbers reproduce exactly. In the simulator's convention, J0437's worst
case is −51 dB, so "the target −54 dB" that notes and docs quote as the simulator's operating
point carries a hidden ~3 dB margin.

The simulator's own pulsar at −54 dB (T 10 ms, 5 % duty, one polarization) needs
SNR/turn = 1.49e-3 → **~62 h per TOA at SNR 7**, not the scenario's 20 h. The shorter
period, wider duty cycle and single polarization of the simulated pulsar add to the 2.7 dB.

Nothing is wrong. Every choice is on the conservative side, and the "≪ 4e-6 of the baseline"
bias requirement is then conservative too (J0437 worst case: ρ = 7.9e-6 in the simulator's
convention). But the two conventions are never stated side by side, so the numbers invite
being mixed. Suggestion (text only): one sentence on the scenario page and in notes §11
saying that the ρ there includes W_eq = 2·W50, and that in `snrDB` terms (Gaussian peak)
J0437's worst case is about −51 dB.

## 9. Small code finding: the optimal fit's convergence flag is dropped

`estimateTOA`'s `fitOptimal` iterates fit → weights → fit and stops when τ and b change by
< 1e-12, at most `MaxIterations` (10). It sets `r.converged` (`estimateTOA.m:495`), but the
caller never copies it into `toa` and never prints it. A sub-int whose weights had not
converged is reported as a normal valid TOA. That is probably rare and harmless: at weak SNR
the weights hardly depend on b. But it is invisible. A 1e-12 relative criterion on b can also
fail at rounding level, so it might run all 10 rounds without anyone knowing.

Option (small logic, Jasper's call): a `toa.converged` field, and in the summary line the
number of sub-ints that did not converge.

## 10. Planning notes: block C draft, item 22, stale roadmap row

- **Block C (fast simulator) draft** (`docs/roadmap.html`) is consistent with what is
  validated, with two details for when it starts:
  - its time-bin noise σ² = m²/(n_pol·B·binDt) should use the noise bandwidth B_n
    (391.6 MHz full band, 3.0596 MHz per channel), not B;
  - its "asymptotic weights" W : W₂ : W_X = n : 2n/3 : n/6 hold for independent time bins
    (the full band). The default channel path has correlated bins (V 0.888, X(1) 0.046,
    Lmax 7); there they come from `foldProfile` with `NoiseCoeffs` or `DataWeights`.

  Working in S_mean/SEFD, as the draft does, also avoids the ρ-convention mix-up of
  section 8. If C is calibrated on the channel path, the 2.2 % of section 2 (the
  pulsar-only taper) would be copied into it unless item 3 is decided first.
- **Item 22, sizing.** On the default channel path, end-to-end error-bar honesty rests on
  seed 43 alone (8 TOAs: rms/pred 0.84, known only to ±25 %). The unit-level Monte Carlo
  (`testOptimalTOA`, exact covariance draws) covers the estimator, not the voltage chain
  (channelizer, per-channel tapers, data-dependent excision). For a check to ±5 / ±3 / ±2 %
  (1σ) it needs 200 / 556 / 1250 TOAs, i.e. 25 / 70 / 157 noise seeds of 8 TOAs. At ~1 min
  per seed on the channel path with excision, the ±3 % check is about an hour. A
  channel-path `runMonteCarlo` would have to stay below ~−3 dB because of section 4.
- **Stale: `docs/roadmap.html:151`.** The code-review row still says "commit pending"; it
  was committed as `0499da7`. The roadmap and home page also have no entry yet for the
  read-through commits `8a756ec` / `1b65f5a`: reuse checks, full-band-only warnings,
  orphan files. The maintenance changelog has them.

## 11. Determinism, and the design rule's "SNR per turn ≈ 373·ρ"

- **The whole chain is deterministic, confirmed on the real files.** Jasper's `main` runs
  regenerated `data/test_rx_IQ.dat`: the full chain with the sky at 13:57, the receiver
  again at 15:28. It is bit-identical to the copy frozen in `data/test/` at 12:05, and so
  are the channelizer outputs (channels 1, 64, 128). Same seeds → same bits through
  generator, dispersion, noise, IQ and channelizer. The frozen test data therefore stays
  reproducible from `pipelineParams` (re-freezing would give the same bits).
- **373·ρ is the flat-band value.** It is 388·ρ (flat band, weak signal) × 0.962 (baseline
  loss). Weak-signal SNR per pulse with the real tapers, effective bandwidth Js²/Inn, × 0.962:

  | Case | B_eff | SNR per turn |
  |---|---|---|
  | flat band (the formula) | 400.0 MHz | 373.2·ρ |
  | full band, as simulated | 388.4 MHz | 367.8·ρ |
  | full band, realistic receiver | 389.6 MHz | 368.5·ρ |
  | **default channel path, as simulated** (`ChannelGain` 1) | 372.3 MHz | **360.1·ρ** |
  | channel path, realistic receiver | 386.7 MHz | 367.0·ρ |
  | either path without the forward taper | 391.6 MHz | 369.3·ρ |

  Turns per sub-int go as 1/SNR². The default path therefore needs ~7 % more turns than
  N = (7/(373ρ))² predicts, and the full band ~3 % more. That is small against the
  scintillation and flux uncertainties of the scenario, but it is one more place where
  item 3 (the pulsar-only taper) shows up. With a realistic receiver both paths are at
  ~367–369·ρ.

## 12. The χ² quality flag is nearly blind on the channel path

`estimateTOA` flags a sub-int when its reduced χ² over the on-pulse bins (template above 1 %
of its peak, 263 bins per channel) exceeds `MaxRedChi2` = 2. That fixed threshold was set on
the full band. The statistical scatter of χ²_red is √(2/dof):

| Path | dof per sub-int | σ(χ²_red) | Threshold 2 = | A 5σ threshold would be |
|---|---|---|---|---|
| full band | 260 | 0.088 | 11σ | 1.44 |
| channel path, 128 channels | 33,661 | 0.0077 | 130σ | 1.04 |

The measured channel-path values, in units of that σ:
- seed 43 with 'optimal': 1.008, 1σ, right;
- B5d LTE edge channels: 1.04, 5σ;
- B5e wide blocks: 1.06, 8σ;
- B5e narrow Inmarsat carriers: 1.22, 28σ, with error bars 3.4 % too small.

None reaches 2, so none is flagged: on the default path the flag catches only gross failures
(a radar spike the fit locked onto), not the noise-model failures it could easily see.

A design choice for Jasper: keep 2 ("profile is not template + noise" in the gross sense),
or scale with dof, e.g. flag if χ²_red > 1 + k·√(2/dof) with k ~ 5. The latter would mark
the B5e/LTE sub-ints, whose error bars are a few % too small. The parked per-channel noise
calibration (item 33) would then turn such flags back into good TOAs.

## Docs errors fixed (Jasper: "fix the docs errors and commit")

All text only. Fixed:
- **The Rice factor**: detection, glossary and physics pages now say up-crossing rate per
  turn and ~60 trials at η 4.15.
- **Physics page**: the noise-like RFI type exists (since 8 Oct); 800 MHz vs the 400 MHz band.
- **Limits page**: the pulsar-only taper costs +0.4 % (full band) and +2.5 % (channel path),
  the whole 2.2 % gap; disk per path.
- **Channels page**: both reasons for the channel path's result (turn 0; the 2.2 %).
- **Scenario page**: a callout on the two ρ conventions (J0437's worst case is −51 dB in
  `snrDB` terms; the simulator at −54 dB needs ~62 h).
- **Formulas, fold, results**: 373·ρ is the flat-band value (368·ρ / 360·ρ with the tapers).
- **Maintenance**: discrepancy and changelog rows.

**Not touched:** `docs/roadmap.html`'s "commit pending" (section 10). Another session is
editing the roadmap and the notes in parallel (Jasper's review-week window: stage 2/3 notes
D1, D2, N1 and `logs/2026-10-13_review-week.md`, uncommitted). Its changes are left
untouched and uncommitted; only this session's hunks are committed.

## 13. Detection thresholds: exact Monte Carlo confirms Rice (wider round)

10⁶ noise-only profiles in numpy: 2048 phase bins, Gaussian noise with the fold's
linear-split correlation (0.25 between neighbours, a 2-tap filter), the zero-mean 5 %
template, the same discrete search over the 2048 bin shifts as `detectPulsar`:

| Test | η | Measured P_FA | Nominal |
|---|---|---|---|
| known phase | 3.090 | 1.06e-3 | 1e-3 |
| unknown phase, Rice for 1e-3 | 4.155 | **1.00e-3 ± 0.03e-3** (MC threshold 4.156) | 1e-3 |
| unknown phase, Rice for 1e-4 | 4.676 | 0.95e-4 ± 0.10e-4 (MC threshold 4.664) | 1e-4 |
| white noise for comparison | 4.155 | 1.07e-3 | 1e-3 |

- The bin correlation smooths T(τ) slightly, so there are fewer up-crossings.
- Rice is exact at 1e-3 and slightly conservative at 1e-4, and becomes more accurate further
  into the tail. The limits page's "P_FA 1e-6 by extrapolating Rice tails" is therefore on
  firm ground for Gaussian fold noise. This complements `runH0` (full band, 10³ profiles) and
  `testOptimalTOA` test 8.

## 14. Further: what block D's phase predictor must model (J0437, a 20 h fold)

The fold's phase-model error budget for J0437 is W50/10 ≈ 14 µs over the fold
(W50 = 0.141 ms). Delays of the effects (sinusoidal approximations):

| Effect | Delay amplitude | Max v/c | Change in 1 h | Change in 20 h |
|---|---|---|---|---|
| Earth orbit (Roemer) | 499 s | 9.9e-5 | 0.36 s | 7.2 s |
| Earth rotation (R_E/c) | 21 ms | 1.5e-6 | 5.6 ms | 42 ms |
| **J0437 binary orbit** (x = 3.367 lt-s, P_b = 5.741 d) | 3.37 s | 4.3e-5 | 0.15 s | **3.0 s** |
| J2145−0750 binary (x 10.16 lt-s, P_b 6.84 d) | 10.2 s | 1.1e-4 | 0.39 s | 7.6 s |
| J1713+0747 binary (x 32.3 lt-s, P_b 67.8 d) | 32.3 s | 3.5e-5 | 0.13 s | 2.5 s |
| clock rate 1e-9 / 1e-10 | – | – | – | 72 / 7.2 µs |
| station position error 1 / 4 / 40 km | ≤ 3.3 / 13 / 133 µs, varying daily with Earth rotation | | | |

(Binary parameters are approximate catalogue values.)

**New for the plan.** The pulsar's own binary motion is in no note, docs page or roadmap.
Block D lists Earth's orbit and rotation, station position and the clock. But the primary
reference pulsar J0437−4715 is a binary: its orbit moves the arrival time by ~3 s within one
20 h fold, about 2·10⁵ times the budget. Two more pulsars in the scenario table are binaries
too (J2145−0750, J1713+0747; the others are isolated). So:
- the D phase predictor needs a binary timing model (circular-orbit ELL1 is enough for
  J0437 at the µs level; the catalogue parameters are precise enough);
- end-to-end tests of D need binary motion in the simulator (generator or block C), the same
  machinery as Earth's Doppler, with a different period and amplitude.

**Position prior (bootstrapping).** A station position error δr makes a daily-varying delay
error of up to δr/c, which smears a long fold. For J0437 the fold needs the position to
≲ 4 km (W50/10). Normal pulsars, with W50 of a few ms, tolerate ~100 km. A navigation
solution from scratch therefore has to start from the wide-pulse pulsars (or a coarse prior)
and refine with J0437, or search over position like an acceleration search.

**Clock.** A stability of ~1e-10 over a fold keeps the drift within the budget (7 µs in
20 h); 1e-9 needs a fitted rate (the limits page says so).

## 15. Further: a first 4-pulsar fix from the scenario numbers (block E expectations)

Weighted least squares for 3-D position and the clock. Assumptions:
- one TOA per pulsar with the scenario's σ_TOA and time per TOA, observed one after the
  other;
- directions from the J2000 names;
- no timing noise or glitches (normal pulsars and Vela have them), visibility ignored.

| Pulsars | 3-D position error | PDOP | Time (worst / nominal) |
|---|---|---|---|
| J0437 + Vela + B0950+08 + B1749−28 | 311 km | 3.3 | 184 h / 34 h |
| J0437 + Vela + B0950+08 + B1642−03 | 313 km | 3.3 | 366 h / 68 h |
| J0437 + Vela + B1642−03 + B1749−28 | 541 km | 4.6 | 395 h / 74 h |
| northern only: B0329+54 + B0950+08 + B1642−03 + B1929+10 | 1550 km | 10.8 | 920 h / 171 h |

Reading:
- The first fix is limited by the normal pulsars' 0.3–0.7 ms TOAs, not by J0437 (12 µs) or
  Vela. Repeated fixes improve as 1/√N, so ~30 km would need ~100 fixes: about 2 years in the
  worst case (100 × 184 h), ~5 months nominal (100 × 34 h).
- The faint MSPs (B1937+21 at 5 µs, J1713+0747, J2145−0750) would fix the geometry
  precisely but take months to years per TOA (scenario page).

This is a rough expectation, not a design. It suggests that block E's achievable precision
with this hardware is set by the normal pulsars and by integration time. That is worth
knowing before D and E are scoped: e.g. averaging many J0437 + Vela TOAs to constrain two
directions very well, with the clock and the third direction from the slow pulsars.

## 16. Wider round: docs numbers, links, figures, statuses

- **Numbers recomputed from the code's formulas, all right:**
  - forward kernel 25,229,772 taps; at L = 1 s: Nfft 2^28, block 243,205,685 samples
    (90.6 %), 17 blocks, 11.8 GB;
  - full-band dedispersion nPast 2001 / nFuture 5,043,956 (= the info file);
  - τ(1.6 GHz) 8.103 ms, τ(1.2 GHz) 14.406 ms;
  - IQ transition 9.8 MHz; channelizer prototype 3853 taps;
  - `detectRFI` n_eff 15.8, and the Gamma approximation's false-flag excess 1.02 / 1.04 /
    1.08–1.10 at P_FA 1e-3 / 1e-4 / 1e-6 (help text and excision page).
- **Option defaults:** all 109 options with defaults (`arguments` blocks) compared with every
  docs mention "Option (value)": no mismatch, only notation (16 ms = 16e-3, 2²², 10⁻⁶) or
  context (the glossary's NBin 2048 is main's value).
- **Links:** all 407 internal links and `#anchors`, and `nav.js`, resolve.
- **Figures:** `gen_results.py` plots exactly the saved phase D/E sweep
  (`data/mc/snrSweep_phaseDE.mat`: rms, robust, predicted, optimum per SNR). `gen_scen.py`
  uses the scenario table.
- **Docs error found and fixed** (`docs/scripts.html`): it said the threshold sweep's results
  are in `snrSweep_phaseDE.mat`, but that file holds the phase D/E sweep and
  `snrSweep.mat` the threshold sweep (checked in the files). Maintenance row added.
- **Probably stale status, Jasper to decide:** notes §5.10 / §13 and `docs/validation.html`
  label `plotDetectedPower` and `plotFoldCheck` "noise version written". Notes §7 records
  their noise versions running at −5 dB on 1 Oct with correct statistics: detected-power
  residual std 0.9995 (73,604 bins), fold std 0.987 / 0.990. Whether that run counts as
  validation (rule: Jasper's run) is for Jasper.
- **Cross-reference to the parallel review-week notes** (stage 2, D2: the truth convention
  is coupled to `refFreq`). Since `8a756ec` a `refFreq` that differs between the sky files and
  `pipelineParams` is at least caught by `checkSkyFiles` (warning in `main`, stop in the
  scripts). The coupling itself remains as D2 describes.
