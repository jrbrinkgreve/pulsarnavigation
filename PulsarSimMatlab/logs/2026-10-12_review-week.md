# 2026-10-12 – review week: stage by stage

Goal (Jasper): understand the whole pipeline, one stage at a time; this week is for
understanding. Refinements found along the way are **noted, not implemented** (Jasper, stage 1:
"stop for now and note that it should change in the future, just doing a review for my own
understanding now"). Stage order and method: `2026-10-12_handover-review-week.md`.

## Stage 1 – pulsar signal (`generatePulsarSignal`, notes §5.1)

Walk-through given:
- Physics: incoherent emission → Gaussian noise whose variance is modulated once per turn,
  x = A·√G·n, E[x²] = A²·G. Reference: T = 10 ms, 5 % → power FWHM 0.5 ms (σ 0.21 ms), white
  0–2 GHz at f_in = 4 GHz (band 1.2–1.6 GHz = 1/5), L = 0.1 s → 4e8 + 1 samples, 10 pulses at
  5, 15, …, 95 ms.
- Code: FWHM → σ ('power' mode: √G amplitude so the power profile has the FWHM); periodic sum
  over the nearest pulse ± `nNeighbors`, normalized by the exact peak (neighbour terms e⁻¹¹¹²
  at 5 %, irrelevant for us); N = round(L·f_in) + 1; private `mt19937ar` stream; blocks of
  4e6 samples, double time, single samples, float32 little-endian; info + `_info.mat`.
- Tests: no unit test (validated by end-to-end runs and the plot checks).

Refinements noted (notes §10 item 11, docs roadmap), not started:
- R1 unit test `testGeneratePulsarSignal` (incl. same seed, different `BlockSize` → identical).
- R2 remove the legacy `'amplitude'` mode: larger than first said — six downstream readers branch
  on `info.envelopeMode` (`addNoiseAndRFI`, `expectedPowerModel`, `plotFoldCheck` 2×,
  `plotIQCheck`, `plotDetectedPower`) plus `pipelineParams` and three drivers.
- R3 a true pulse offset unknown to the observer: every end-to-end run so far has true offset 0
  (`ephem.TRef` = exact true epoch). `validateTOA` matches TOAs to the nearest true centre (fine)
  but reports the whole-fold offset as if the truth were 0 → needs `info_gen.timeOffset`.
- Physics gaps (white spectrum, one stable Gaussian, no jitter/scintillation): already in §9,
  for D/E.

No code changed. Next: stage 2 (dispersion, `applyDispersionStream`).
