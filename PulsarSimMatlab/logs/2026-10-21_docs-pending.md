# Docs changes on hold (21 Oct 2026)

Jasper is rewriting the docs site for readability in a separate session. Claude's doc changes wait here
and go into the pages afterwards, in the tone of the rewrite. Each item: the page, where, and the
content (facts and numbers; wording to be adapted). Remove an item once it is in the page.

## 1. limits.html, "Estimation and detection", after the error-bar entry (was commit 008cfed, reverted)

**Rare outliers at SNR ~7 (20–21 Oct).** With adaptive groups about 0.5 % of B0329+54's TOAs land one core
width off (0.005–0.009 turns, 5–9σ), at SNR 7.3–7.9; χ² cannot flag them. The template's autocorrelation
has no second peak there: the narrow core sits on a broad component, and noise now and then moves the best
match onto that shoulder, which the error bar (from the curvature at the peak) does not cover. Fixed groups
reach SNR ~5.5 and add true jumps to the side peak at ±0.025 turns (19σ). None at SNR 10 (0 of 600).

**Decision (21 Oct, Jasper): keep SNR 7 for more TOAs, and use redundancy.** Block D's navigation filter
checks each TOA against the prediction from the other TOAs and the ship's dead reckoning (B0329+54 nominal:
TOA error 220 km, prediction ~70 km; outliers 1100–2000 km = 4.8–8.6σ from the prediction; a 4σ gate
rejects good TOAs about once in 4–5 years), with an error model that has a small outlier part (rate and
width from the simulator). No fixed groups for multi-component pulsars. Unhandled, one 6σ outlier among a
day's 10 TOAs would bias that day by ~130 km.

## 2. roadmap.html, "Block D: motion and barycentric phase", after its first paragraph (was 008cfed)

**Agreed for D (21 Oct): an outlier-tolerant fit.** The adaptive SNR-7 groups give the most TOAs, but
~0.5 % of B0329+54's land one core width off (link to limits#estimation). The navigation fit uses
redundancy: each TOA is checked against the prediction from the others plus the ship's dead reckoning
(innovation gating at ~4σ), with an error model that includes a small outlier part. D's first deliverable
— a 1-D range from B0329+54 recovering a 100 km site offset — is tested on simulated TOAs that contain
these outliers.

## Still to come (results of the runs of 21–22 Oct; added here as they land)

- Validation run (Claude's run of the frozen 171a716): tests.html / results.html status.
- expSimVsChain at DM 25 / L 0.2 s (after main on the full band): tests.html.
- New RFI profiles (per-channel detectRFI thresholds, 30 lags): simulator.html.
- testRFIChain on the new code (narrow Inmarsat): tests.html, simulator.html assumptions b and d.
- 100-day runs: scenario.html (final budgets), roadmap status.
