# 17 Nov 2026 — budget recheck against a paper's claim, and a 30 m² scenario run

Jasper is reading a paper that says TOA integration times "aren't that bad": 30 m² of aperture at
1.2–1.6 GHz gets 0.3 km from J0437−4715 after 1 h. He asked for a recheck of the noise and power
budgets, also for a cascade of LNA and SDR, then "can you run it in a scratchpad?" and "Yes, log it".
The paper itself was not named or read. No repository code or docs changed. The MATLAB runs below are
Claude's (asked for by Jasper), not his validation.

## 1. The budget recomputed (Python, by formula)

| Quantity | Recomputed | Notes §11 |
|---|---|---|
| SEFD nominal (70 K, 0.61 m²) | 3.17·10⁵ Jy | 3.2·10⁵ |
| SEFD worst (80 K, 0.3 m²) | 7.36·10⁵ Jy | 7.4·10⁵ |
| A_eff of a 22 dBi beam at 1.4 GHz | 0.58 m² | 0.61 m² |
| B0329+54, SNR 7, box rule, nominal / worst | 0.78 / 4.2 h | 0.8 / 4.2 h |
| B1937+21, box rule, nominal / worst | 606 / 3270 h | 670 / 3600 h |

Power per polarisation in 400 MHz: receiver noise −94.1 dBm; B0329+54 mean −156.1 dBm (−61.9 dB),
peak −44.6 dB (box rule's rho). The arithmetic holds; the long times are not an error.

Levers on B0329+54's 0.78 h (one change each): SNR 5 → ×0.51; W = W50 → ×0.50; flux ×2 → ×0.25;
T_sys 30 K → ×0.18; 400 MHz band with S 1.5 Jy, B 100 MHz, T_sys 110 K (assumed) → ×0.18;
A_eff 10 m² → ×0.004.

## 2. Receiver cascade (Friis; T_ant 17 K; SDR NF 8 dB assumed)

| Layout | T_rx | T_sys | time vs 70 K |
|---|---|---|---|
| A: 0.2 dB feed loss, LNA 0.35 dB/19 dB, LNA 0.54 dB/22.8 dB, filter 2 dB, cable 3 dB, SDR | 40 K | 57 K | ×0.66 |
| B: one LNA 0.5 dB/20 dB at the element, cable 3 dB, SDR | 86 K | 103 K | ×2.2 |
| C: 1 dB of cable before one LNA | 162 K | 179 K | ×6.5 |
| D: 1.5 dB filter before two LNAs | 155 K | 172 K | ×6.0 |
| E: SDR alone on 3 dB of cable | 3360 K | 3380 K | ×2300 |

One 20 dB stage leaves the SDR adding 32 K; loss before the first LNA costs 75 K per dB. Layout A
puts the noise at −58 dBm at the SDR input and a +40 dB radar at −18 dBm. Jasper's actual chain
(components, losses, SDR noise figure) is not known yet.

## 3. Why the radiometer rule is optimistic: the box

The rule's factor sqrt((P − W)/W) comes from the on/off derivation: all flux in a box of width W
(S_peak = S_mean P/W), every instant in it weighted equally. For any shape, with the profile as the
weight, the factor is sqrt((<s²> − <s>²)/<s>²). B0329+54's five components (block C log): area/peak
width 14.15 ms but ∫p² dt = 6.71 ms, 47 % of a box (a Gaussian: 71 %). Shape factor 23.0 against
the box's 53.1 at W = 13.2 ms.

| Step | Time for SNR 7 |
|---|---|
| Box rule | 0.78 h |
| Five-component shape | 1.73 h |
| With the baseline fitted | 1.81 h |
| Simulated (20 Oct) | 1.93 h |

The last 7 % is not explained (guess: the stopping rule closes on whole 120 s sub-ints).
Correction to the chat: an earlier answer blamed the lower peak; the rule's W = 2 W50 already had
the peak about right. The cause is the efficiency of the shape.

## 4. The paper's claim against our numbers

Scaling the simulated J0437 nominal case (13 µs at ~15.75 h, 52 µs·sqrt(h)) by the area ratio
30/0.61 = 49 gives 1.05 µs after 1 h = 0.31 km. The claim agrees with our budget at T_sys ≈ 70 K;
the paper has 49× the aperture (about 440 dual-pol elements, or a 6–7 m dish).

## 5. The 30 m² run (fast simulator) — a throwaway test

Jasper, later the same day: "this was just a test, forget those scratch numbers". The numbers in
this section are not a design point and not a project result. They do not go into the docs, the
roadmap or the scenario tables. What stays is the chi² finding and its cause (section 6).

A scratch copy of `runScenario.m` (repository script untouched) with: `sefd.big30` = 2k·70/30 =
6443 Jy; the rho loop over all fields of `sefd`; `TargetSNR` and `MaxSubints` settable; results to
its own folder. rho: J0437 −34.22 dB, B0329+54 −27.98 dB. Variant clean, motion off, 70 K, NPol 2.

| Case | TOAs | sigma | c·sigma | std z | chi² mean | pChi2 < 0.05 | sigma·sqrt(h) | run |
|---|---|---|---|---|---|---|---|---|
| J0437, 30 d, fixed 12 sub-ints (1 h) | 367 | 1.074 µs | 0.32 km | 1.079 | 1.0018 | 9.8 % | 1.07 | 2.1 min |
| J0437, 30 d, TargetSNR 7 (one 300 s sub-int each) | 4410 | 3.714 µs | 1.11 km | 1.035 | 0.9986 | 3.3 % | 1.07 | 5.9 min |
| B0329+54, 10 d, fixed 30 sub-ints (1 h) | 195 | 21.32 µs | 6.39 km | 1.032 | 1.0853 | 100 % | 21.3 | 2.5 min |

Predicted by scaling: 1.05 µs, 3.6 µs, 20.5 µs; the runs are 3–4 % above. No |z| > 5 anywhere.

**Findings.**
- The simulator reproduces the paper's 0.3 km after 1 h for J0437 at 30 m².
- B0329+54 at 50× the signal: every TOA has pChi2 < 0.05 (mean red. chi² 1.085), while the TOAs
  stay honest (std z 1.03 ± 0.05). Diagnosed in section 6: the channel gains, not the bins.
- J0437's error bars are a few % small: std z 1.035 ± 0.011 (4410 TOAs), 1.08 ± 0.04 (367).

**Not covered.** A position (the site is known; this is a range along one line of sight), RFI and
scintillation (clean only), the clock (motion off; at 1 µs a crystal would dominate).

**Files** (gitignored, this Mac only): `.claude/runs/2026-11-17_big30/` — `runScenario30.m`,
`test1d.m`, `full30.m`, the two logs and `out/` with five result files; the chi² test:
`runScenario30b.m`, `chi2test.m`, `chi2test.log`, `out2/`; section 7: `runScenario30c.m`,
`chi2test_gain.m`, `chi2test_gain.log`, `out3/`, `testFoldFaults_17nov.log`.

## 6. Why chi² flagged every B0329+54 TOA (Jasper: "can you look into why")

**Cause: the fit assumes the same pulse amplitude in every channel, and the simulated truth does
not have it.** Correction (same day): this was already an open item, roadmap item 3 ("Band-edge
taper on the pulsar only", read-through of 13 Oct, `logs/2026-10-13_full-review.md`), which also
says the pChi2 excess grows with SNR² and is irrelevant at the weak-signal target. The chat first
called it FM-41. FM-41 is the related case of a receiver bandpass, which scales pulse and noise
alike; the taper scales the pulse only, so a channel's mean power does not show it. The runs below
confirm item 3 with numbers; they found nothing new.

- Truth: `runScenario`'s `chanGain` gives each channel the fraction of the pulse the forward
  band-edge taper leaves (8 MHz at each band edge). Three channels at each edge are affected:
  g = 0.023, 0.417, 0.960; the other 122 have g = 1. Mean 0.975, sum of (g − mean)² = 2.51.
- Fit: `foldToTOA` calls `estimateTOA` 'optimal' without 'ChannelGain', so s_c = 1 for all channels
  and one amplitude b. The fitted b is the mean, and each edge channel keeps a residual
  b (g_c − mean) × template. Its chi² grows with the signal squared:
  excess ≈ sum (g_c − mean)² × SNR² / (nChan × mean² × dof), SNR the TOA's.
- Sizes (dof = on-pulse bins × 128 − 3; chi² uses diagonal variances, the reported SNR the full
  covariance, so this estimate is a lower bound by up to ~1.5):

  | Case | SNR per TOA | dof | predicted excess | measured chi² |
  |---|---|---|---|---|
  | B0329+54, 30 m² test, 1 h | 247 | 21117 | ≥ 0.060 | 1.085 |
  | J0437, 30 m² test, 1 h | 86.5 | 64765 | ≥ 0.0024 | 1.0018 |
  | J0437, 30 m² test, 300 s | 25.0 | 64765 | ≥ 0.0002 | 0.9986 |
  | B0329+54, real array, SNR 7 | 7 | 21117 | 0.00005 | — |

  The 128 fitted baselines are not taken off the dof, which puts the clean value near 0.999.
- Paired test (Claude's run, 3 days, 58 TOAs each, same seeds; `chi2test.m`):

  | Truth | red. chi² mean (min–max) | pChi2 < 0.05 | median SNR | std z |
  |---|---|---|---|---|
  | band-edge taper (as `runScenario`) | 1.0856 (1.048–1.109) | 100 % | 247.0 | 1.048 |
  | flat, g = 1 | 0.9994 (0.966–1.024) | 3.4 % | 253.4 | 1.048 |

  With a flat truth the excess is gone. The bins and the template are therefore not the cause.

**What it means.**
- It is a model mismatch in the fit, not a wrong TOA: std z is the same in both rows, and the SNR
  loss is 2.5 % (the mean gain 0.975).
- At the real array's signal it is invisible (excess 5·10⁻⁵), so no existing result changes.
- It would matter for any strong signal: a bright pulsar on a large aperture, a calibration source,
  or real data with a sloped bandpass, where the gains differ in all channels and not only six.
  pChi2 is the alarm the docs say to watch on real data; with a real bandpass it would ring for
  this reason first.
- FM-41's listed fix (the gain from each channel's mean power) would not remove this flag, because
  the taper leaves the mean power alone. The fit needs the gains as known values (section 7).

## 7. Code change: `foldToTOA` 'ChannelGain' (Jasper: "fix the chi-squared bug", then "Yes, step 1")

Plan put to Jasper first (three steps); he approved step 1.

- **Step 1, done:** `foldToTOA` takes 'ChannelGain' and passes it to `estimateTOA` and
  `detectPulsar`: [] (default, nothing passed, results bit for bit as before), a [1 x nChan]
  vector of known gains, or "measured" (each channel's mean power over `fold.profTotal`,
  normalized to mean 1; FM-41). `info.channelGain` records it. "measured" is not the default: it
  would shift the frozen reference numbers slightly, and leftover RFI raises a channel's mean power.
- **Test:** `testFoldFaults` case 1b. Claude's run of the whole test: ALL PASSED
  (`.claude/runs/2026-11-17_big30/testFoldFaults_17nov.log`). Not validated until Jasper's run.
  - default = gain 1 exactly (TOAs, error bars, T0); = G exactly; "measured" = the direct call;
  - taper on the pulse alone, total fold at SNR 276: red. chi² 1.2934 with gain 1 (excess 0.293,
    predicted ≥ 0.190; pChi2 2.5·10⁻⁶¹); 1.0317 with the taper as ChannelGain (pChi2 0.027, one
    draw, about 2 sigma); SNR ×1.0102; TOA shift 0.01 sigma.
- **Check on the scratch scenario** (3 days, 58 one-hour TOAs, same seeds, tapered truth, the
  taper given to the fit; `chi2test_gain.m`): red. chi² mean 0.9994 (0.966–1.024), pChi2 < 0.05
  for 3.4 %, median SNR 249.5, std z 1.042. Against 1.0856 and 100 % with gain 1.
- **Docs:** detection.html D5 (paragraph), toa.html (s_c entry), tests.html (case 1b),
  roadmap.html (FM-41 line, item 3), maintenance.html (changelog). `api.html` reads the committed
  code: regenerate with `docs/tools/make_api.py` after the commit.
- **Not done (steps 2 and 3, each needs Jasper's OK):** `adaptiveSubints` passing the gains to
  `detectPulsar`; `runScenario` (and the experiments) giving the taper to the observer, which
  touches the ground-truth rule.

## Open

- Which paper it is, and its T_sys, threshold and pulse model.
- Jasper's real receiver chain, for the cascade.
- Steps 2 and 3 of the ChannelGain plan; "measured" as a default (FM-41) after a look at RFI.
- J0437's error bars a few % small in the throwaway test (std z 1.035 ± 0.011): not looked into.
