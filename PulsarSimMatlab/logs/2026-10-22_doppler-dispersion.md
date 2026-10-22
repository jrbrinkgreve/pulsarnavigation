# 2026-10-22 – Q&A on block C and Earth motion; open item: Doppler-shifted dispersion delay

A question-and-answer session (no code changed, nothing run).

## State at start

Block C closed (22 Oct); block D: D1–D5 validated, close-out waits for the ship comparison
(`2026-10-22_handover-D.md`). Uncommitted at the start: `runScenario.m`, `tests/expNavFit.m`
(other work, not touched here).

## Questions answered

1. **How block C is constructed.** `simulateFold` step by step (schedule, phase offset, mean
   profile in Fourier space, asymptotic weights, circulant noise by FFT, RFI profile with the
   observer's clean model, scintillation field, assembly into foldProfile's fields and a
   separate truth), the helpers around it (rfiProfile, mergeSubints, adaptiveSubints,
   foldToTOA, runScenario), and the build order C1 → C2 → C4a → C4b → C3 with its validation.
2. **How the Earth's motion against the pulsar enters** (Jasper: "in the end, we need to
   compensate this while folding"). It is compensated in the fold, in block D: `makePredictor`
   (the assumed trajectory and clock estimate) and `foldProfile 'Predictor'` (phase per time
   bin). `simulateFold` draws the fold after that compensation, in the predictor's frame. Only
   the residual (`trueResiduals`: true minus assumed trajectory, unestimated clock) enters, as
   `'PhaseOffset'` per sub-int (`runScenario` motion.on). Two gaps named: the voltage generator
   has a fixed period (predictor folding is checked on synthetic power series only), and the
   Doppler-shifted dispersion delay below.

## Finding: Doppler-shifted dispersion delay (real data only)

The dispersion delay holds at the barycentric radio frequency, which the Earth's orbit shifts
by up to 1e-4. `delayModel` has no dispersion term; `makePredictor` treats the DM delay at
refFreq as a constant. Missing at the observer: ≈ 2(v·n̂/c)·K·DM/ν_ref², annual.

| Pulsar | DM | K·DM/ν² at 1.6 GHz | v·n̂/c max | missing term | in-band sweep change |
|---|---|---|---|---|---|
| B0329+54 | 26.8 | 43.4 ms | 8.2e-5 | ±7.1 µs | ±2.8 µs |
| J0437−4715 | 2.64 | 4.3 ms | 3.7e-5 | ±0.3 µs | ±0.1 µs |

Earth rotation adds ~0.05 µs (daily) for B0329+54. No effect in the simulator: truth and
prediction both leave the term out.

## Decisions

Jasper: record it as an open item for real data (notes §10 "Open items", next to the
ionosphere item); no code change now. Then: "add a row for the Doppler-shifted dispersion delay
to docs/limits.html" → done: a row in the Propagation table after "DM variations" (size,
matters for, plan), checked in the rendered page. The size figure at the top of the page and
the "In short" box are unchanged.

## Open items

- Doppler-shifted dispersion delay: a DM term in `delayModel` at the barycentric frequency
  when real data come (notes §10).
