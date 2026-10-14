# 2026-10-14 – docs: the remaining pages in the new format

Continues `2026-10-13_docs-rationale.md` (the stage pages, committed in `4c8e87a`). Jasper: "continue with all the other
pages in that format, i like it!" No code change, no MATLAB run.

## Done

Every page except the generated `api.html` now follows the agreed format where it applies: "In short" first, plain
steps, exact formulas folded into `<details>`, "Design choices" boxes (the choice, the alternative, why).

| Page | What changed | New figure(s) (`docs/tools/figures/`) |
|---|---|---|
| `physics.html` | rewritten: "In short" per section, plain steps, design-choice boxes (square first; coherent dedispersion; ρ as a density ratio; detect → fold → template; FFTFIT not the centroid; good-TOA rule; blank before dedispersion) | `gen_navgeom` (lines of position), `gen_whysquare` (replaces the hand-drawn voltage sketch), `gen_radiometer`, `gen_foldgrow`, `gen_flanks` |
| `conventions.html` | "In short: the ten rules"; each "Why" box as choice / alternative / why; plain steps | `gen_timeprec` (float32 vs double, a summed clock drifts) |
| `pipeline.html` | "In short", how to read the map, design-choice boxes, anchors `#map #stages #why #run #sizes` | `gen_datarate` (bytes per stage vs 24 GB) |
| `parameters.html` | "In short", new section "The numbers at a glance" (`#plan`, `#scales`) | `gen_freqplan`, `gen_timescales` |
| `scenario.html` | "In short", contents, how the time per TOA is computed, new section "Who can see what" | `gen_visibility` (hours above 20° vs latitude; matches deep-analysis §17) |
| `scripts.html` | "In short", design-choice boxes (own files, skipping, H0 = zero sky, common random numbers, leftover in baseline units) | `gen_experiments` (what each script varies) |
| `tests.html` | "In short", design-choice boxes (functions; pass criteria from the predicted σ) | `gen_testmatrix` (regression / algebra / physics per test) |
| `results.html` | "In short", section "At a glance" (`#glance`) | `gen_honesty` (every honesty check on one axis) |
| `limits.html` | "In short", contents, section anchors, `#sizes` | `gen_effects` (size of each left-out effect as a timing error) |
| `lessons.html` | "In short", section anchors | `gen_neighbour` (lesson #8: shared noise, up to 1.5×) |
| `history.html` | "In short"; 9 Oct studies (f_out, dispersed domain); the review week 12–13 Oct; decision rows 9–13 Oct; 12 missing logs in the log table | `gen_timeline` |
| `roadmap.html` | "In short"; the two "why this order" boxes with their alternatives (only these hunks are mine; the D1/D2, N1, I1–I3 rows belong to the review-week session) | – |
| `formulas.html` | "In short"; an "In words" column for all 43 formulas | – |
| `glossary.html` | "In short"; seven terms (blocks, common random numbers, line of position, phase D/E, PRF, TECU, W50) | – |
| `index.html`, `maintenance.html` | stage-card text and status line; the page format documented; changelog and a discrepancy row | – |

Site check (all 26 pages): tags balanced, ids unique, all internal links and anchors resolve, all CSS classes defined, all
SVGs labelled (only the home hero has no caption, by design). Defaults check: one gap, fixed: the sky page did not state
`SaveInfo` (true) and `Verbose` (false) of `generatePulsarSignal`. Figures rendered in dark and light mode before insertion;
their numbers come from the parameters, the code's formulas or the recorded results (cross-checked: visibility hours equal
deep-analysis §17; data volumes 60.8 / 63.5 / 72.5 GB per second equal the limits page).

## Reasons to confirm (reconstructed, not yet in notes §12)

Most boxes restate reasons already recorded (notes §12, the decision log, the logs). These are new or extended:

1. **A template fit, not the highest bin**: the highest bin is limited to the bin width and jumps with the noise (the
   centroid part, ~2× noisier, is recorded).
2. **Pass criteria from the predicted σ**: a hand-picked tolerance can be loose enough to hide a real error or tight enough
   to fail a correct code by chance (recorded: "guessed tolerances prove nothing").
3. **A private stream per random process, not the global generator**: with one global stream every draw would depend on
   every other draw.
4. **Headerless float32, not a self-describing format**: the description of each file lives in its `_info.mat`.
5. **Detect → fold → template is also the cheapest order**: detection cuts the data rate ~800× on the full band, and
   folding turns hours of data into one profile before the fit.

## Confirmed and committed

Jasper: "the reasons are correct, add them to notes and commit". The five reasons are in notes §12 ("Added 14 Oct"). Committed
with only this session's hunks of `roadmap.html` and the notes; the review-week session's uncommitted work (notes §10, roadmap
rows, its log) stays in the working tree for that session.
