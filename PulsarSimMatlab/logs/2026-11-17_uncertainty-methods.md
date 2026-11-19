# 17 Nov 2026 — note: three ways the project gets a TOA uncertainty

Jasper: "can you note the difference in uncertainty estimation methods between direct analytical
calculation via the radiometer equation, following the SNR expression, and the empirical fit
uncertainty?", then "and can you reason and comment on it". No code or docs change; nothing was run.
Read for this note: `addNoiseAndRFI.m` (predictPerformance), `estimateTOA.m`, `docs/results.html`
(sections 3 and 4), the scenario rule in `current_project_notes.md` §11.

## The three

**A. Analytic, before observing (radiometer equation → SNR → sigma).**
- Scenario rule (notes §11, `docs/scenario.html` R.1):
  SNR = (S_mean / SEFD) sqrt(n_pol B t) sqrt((P − W_eq) / W_eq), sigma_TOA = sqrt(2) (W50 / 2.355) / SNR.
  The sqrt(2) sigma_p / SNR is the Cramér–Rao bound for a Gaussian pulse of width sigma_p in white
  noise (the integral of p'^2 is the integral of p^2 over 2 sigma_p^2).
- The exact form of the same idea: `predictPerformance` in `addNoiseAndRFI.m` and
  `expectedPowerModel`: Fisher information B ∫ (rho p')^2 / (1 + rho p)^2 dt per pulse, with
  self-noise in the denominator. This is the "opt" column of the SNR sweep.
- Inputs: assumed flux, SEFD, bandwidth, time, pulse shape. It is the best any estimator can do.

**B. The fit's own error bar, from the data (`estimateTOA`, toaErr; the "pred" column).**
- var(tau) = d' Cov(p) d / C''(tau)^2: the fit's sensitivity to each profile bin, times the
  covariance of the bins.
- Cov(p) comes from a noise model, not from the scatter of the data: 'radiometer' (default),
  var = m^2 / (Bnoise binDt) per time bin with m = a + b template the FITTED power model, carried
  through the fold weights (neighbouring-bin covariance from linear assignment, blanking). Options
  that are more empirical: 'offpulse' (variance and lags from off-pulse residuals), and `noiseCal`
  (rfiProfile measures each channel's V and X on this observation's off-pulse bins).
- Inputs: the measured baseline and amplitude, the real data weights after excision, the real
  channel noise bandwidths. One value per TOA.

**C. Measured scatter against the truth (Monte Carlo; simulation only).**
- rms of (TOA − truth) over seeds, compared with B ("rms/pred"), plus red. chi^2 and pulls.
  It is the only one of the three that tests the others, and it does not exist for real data.

## What the record says (docs/results.html)

- B against C, −5 dB reference, 90 TOAs: rms / predicted = 2.762 / 2.907 us = 0.950 (expected
  1 ± 0.075). SNR sweep, SNR per sub-int ≳ 10: rms/pred 1.03–1.06, no outliers.
- B against A: mean predicted sigma 0.5 % above the exact optimum (2.907 against 2.892 us);
  pred/opt 1.00–1.01 up to −5 dB. At high SNR they part: 1.5 at +10 dB, 2.4 at +20 dB, because
  FFTFIT weights for white noise and self-noise is not white (not the project's regime).
- Below SNR ≈ 10 both A and B fail against C: the fit locks onto noise peaks. Outliers 5.6 % at
  SNR 3.8 and 57 % at 1.2, rms/pred 11–12, while the core stays right (robust 83.9 against a
  predicted 82.2 us).
- A's simple rule against the fast simulator with measured profiles: J0437 3.9 h → 15.5–16 h per
  TOA, B0329+54 0.8 h → 1.93 h. The single Gaussian with W_eq = 2 W50 missed J0437's wings and
  overstated B0329+54's peak by 3 dB.
- Long runs: the core of B0329+54's pulls has width 1.01, but 1.9 % of its TOAs lie beyond 3 sigma.

## Differences, in short

| | A analytic | B fit error bar | C scatter |
|---|---|---|---|
| When | before observing | per TOA, from the data | after many repeats, with truth |
| Knows the real flux, T_sys, blanking, bandpass | no (assumed) | yes (fitted, measured) | yes |
| Estimator | the optimal one (a bound) | the actual one (FFTFIT) | the actual one |
| Sees outliers at low SNR | no | no | yes |
| Sees a wrong noise model | no | only through chi^2 | yes |
| Sees systematics (template, DM, clock) | no | no | only those in the simulation |
| Available for real data | yes | yes | no |

## Comments (Claude's reasoning, not measured unless a number is quoted above)

1. They answer different questions. A: how long must I observe? B: how much weight does this TOA
   get in the navigation fit? C: are B's error bars honest? The project already uses them this way.
2. A and B agreeing to 1 % in the weak-signal regime is expected, not a coincidence: with receiver
   noise dominating, the noise is white and FFTFIT is the optimal estimator, so its variance is the
   bound. The agreement is a check that the noise model in B is right.
3. B's error bar is itself noisy. It scales as 1 / b, and the fitted amplitude b has a relative
   error of about 1 / SNR: about ±14 % at SNR 7. A is noise-free but only as right as its assumed
   flux, and catalogue fluxes and scintillation move that by a factor 2.
4. Selection makes B slightly optimistic near the threshold: TOAs pass the detection test more
   often when noise pushed b up, and a larger b gives a smaller error bar. Not quantified here; the
   existing threshold sweep (SNR ≥ 6–7) is where it would show.
5. Neither A nor B covers the non-linear regime. Both are second-derivative (Gaussian) statements
   around the true peak. That matters for the redesign: a 120 s sub-int has SNR below 1, so no
   per-sub-int error bar from B means anything. The uncertainty must then come from the joint fit
   over all sub-ints (its curvature at the solution), and from a Monte Carlo of that fit.
6. For real data C is gone, so B needs stand-ins: red. chi^2 and pChi2 per TOA (already there),
   the chi^2 of the navigation fit's residuals, repeated TOAs of the same pulsar, and pulsar-free
   recordings (FM-34). A measured noise model (`noiseCal`, 'offpulse') moves B from "the model
   says" toward "the data say" and protects it against RFI residue and gain noise (FM-42).
7. What none of them contains: template error, mode changes (FM-3), DM error, the clock, the
   timing model. Those enter as extra terms in the navigation fit (FM-27: a per-pulsar error scale
   and an outlier fraction), not in the TOA error bar.
