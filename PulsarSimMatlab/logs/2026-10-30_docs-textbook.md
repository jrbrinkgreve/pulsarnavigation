# Docs rewrite as a textbook (30 Oct 2026 –)

Jasper, 30 Oct: "I like this style much better! Can you back up the docs, and rewrite all pages in this style? Textbook,
linked, yet relatively short format. Make sure to keep track of consistency and accuracy." Then: "I already approved of
the style, so you can go ahead and do the whole thing."

The style is the one of stage 5a on channels.html (commit 77bf36d): numbered sections that each open from the one
before and close on the next question (one or two sentences); numbered equations and figures referenced from the
text; derivations with each step said in words; a short worked example where it helps; a short Summary box; captions
that say only what the figure shows. The 21 Oct sentence rules stay (one fact per sentence, no dates or review codes in
the running text, History in collapsed boxes, tables instead of slash lists). Stage pages keep their DM 5 figures and
worked numbers, with DM 25 in note boxes.

**Backup.** Git tag `docs-before-textbook` on 77bf36d; a copy in
`~/Documents/pulsarnavigation-backups/docs-before-textbook_2026-10-30_77bf36d.tar.gz`.

**Scope.**

| Group | Pages | Treatment |
|---|---|---|
| A, explanatory | channels (rest), sky, receiver, fullband, excision, fold, simulator, toa, detection, validation, navigation, scenario, physics, localisation, conventions, pipeline | full rewrite |
| B, reference | index, parameters, scripts, tests, results, lessons, limits, roadmap, formulas, glossary | light restyle: purpose and links up front, duplication trimmed, tables kept |
| C, records and generated | history, maintenance, api | unchanged, apart from links |

**Checks per page** (scratchpad tool `tb.py`, against the tag):
1. every number, `code` name, link and id of the old page that is missing from the new one is listed, and each is
   accounted for below (moved, merged, or dropped as a repeat);
2. current settings are checked against the code (`pipelineParams.m`, the function headers), not copied;
3. ids that other pages link to stay; all links of the docs resolve (baseline 30 Oct: 1041 links, none broken);
4. rendered in the browser pane: KaTeX errors, console errors, figure numbering.

Figures keep their SVG. A caption that changes is listed here, so that the `gen_*.py` scripts can be brought in line
in the final pass (they write the caption into the snippet).

## Progress

| Page | Words before → after | Numbers / names not carried over, and why | Checks | Notes |
|---|---|---|---|---|
| channels | 5221 → 5383 (5a already done in 77bf36d; the rest restructured, two 5a design boxes that repeated 5a.1 and 5a.5 removed) | none missing | facts ✓, links 1041/0 ✓, render ✓ (75 equations) | 5b.1–5b.4 and 6b.1–6b.3 numbered; `powerCovariance` moved inside 6b (id `powcov` kept on 6b.2); the DM-memory figure moved to 5b.4 as Figure 5b.1; captions unchanged apart from the label (gen_dmmem, gen_powcov: add the label in the final pass) |
| sky | 4293 → 3846; visible 3810 → 3269; captions 583 → 446 | none missing; `channels.html#why` → `channels.html#dedisp-mem` (the memory figure moved; id added on channels 5b.4) | facts ✓, links ✓, render ✓ | 1.1–1.4, 2.1–2.4; three design boxes folded into their sections (random noise → 1.1, exact dispersion → 2.1, band-limited filter → 2.2); "History: why the current form" now a collapsed box. **Captions shortened** for Figures 1.1–1.3 and 2.1–2.3 (gen_voltage, gen_profile, gen_rfspec, gen_waterfall, gen_disp, gen_kernel) and 2.4 (hand-written overlap-add) |
| receiver | 6768 → 6477; visible 6063 → 5772; captions 904 → 738 | none missing; id `clock-why` dropped (only the old table of contents linked it); new ids `snr`, `scenarios` | facts ✓, links ✓, render ✓ (4 equations, 9 figures); **settings checked against the code**: both RFI lists in `pipelineParams.m` (every source, frequency, rate and INR), `noiseSeed` = 43, `rfiOn` false by default, `RiseTime` 0.1 µs, `SidelobeDB` −30, `ScanPeriod` Inf, `BlockSize` 4e6, `Seed` 1, makeClock presets (h₋₁ = 6.492e-22, shown as 6.5e-22), `Step` 60 s | 3.1–3.5, 4.1–4.3, 4a.1–4a.3; eight design boxes folded into their sections (RF after dispersion → 3.1, SNR as PSD ratio → 3.2, own streams → 3.3, INR vs all noise → 3.4, two lists → 3.5, analog receiver / 800 MHz / low-pass / odd length → 4.2, factor 2 → 4.1, flicker → 4a.2); the realistic scenario is now a paragraph with id `realistic` inside 3.5. **Captions shortened** for Figures 3.1 (hand-written), 3.2 (gen_rxorder), 3.3 (gen_snrscale), 3.4 (gen_rfiscene), 3.5 (gen_radar), 4.1 (gen_iq), 4.2 (gen_iqfilter), 4a.1 (gen_clock), 4a.2 (gen_clockdraw) |
| fullband | 2884 → 2916; visible 2559 → 2591; captions 431 → 375 (already compact; the gain is structure) | none missing | facts ✓, links ✓, render ✓; outline checked (box placement); settings checked: `MaxMemoryGB` 4, `MaxLeakage` 1e-6, `CheckLeakage` true, `EdgeFrac` 0.02, `AllowRefOutsideBand` false, `BlockSamples` 2²², `FineBinsPerChan` 8, `NBin` 2048 | 5.1–5.4, 6.1–6.3; **fixed**: the legend said `ephem.DM` = 5, now "5 in the worked numbers, 25 by default"; design boxes folded in (round trip → 5.1, overlap-save / leakage check / fully supported → 5.2, coherent → 5.4, square-fold-match → 6.1, integrate-and-dump / centroid → 6.2, 1 µs → 6.3). **Captions shortened**: 5.1 (hand-written overlap-save), 5.2 (gen_support), 5.3 (gen_smear), 6.1 (hand-written bin stamps), 6.2 (gen_fout) |
| excision | 7798 → 7431; visible 7196 → 6512; captions 599 → 527 | none missing | facts ✓, links ✓, render ✓ (4 equations, 7 figures), outline ✓ (tables and boxes in place); settings checked: detectRFI `Scales` [1 2 4 8 16], `PFA` 1e-6, `BaselineTime` 16 ms, `Passes` 2, `FreqGuard` 7, `FreqGuardMin` 100, `ChannelThresholds` true, `MaxMemoryGB` 4; blankingWeights `MinWeight` 1e-6, `InputSpectrum` 'channelizer'; periodicRFI `TimeTol` 4, `MinPulses` 6, `FalseAlarm` 1e-6; `periodicMask` false | Already the 21 Oct model page; now E1–E9 with subsections E2.1–E2.5, E4.1–E4.3, E8.1–E8.2, equations (E.1)–(E.4) and Figures E.1–E.7 (letter prefix: not a stage). Callouts folded into the text; blankChannels' outputs and tests and blankingWeights' implementation list moved into collapsed boxes. **Captions shortened**: E.1 (hand-written chain), E.2 (gen_blanksub), E.3 (gen_excision), E.4 (gen_windows), E.5 (gen_blankvar), E.6 (gen_locked), E.7 (gen_periodic) |
| fold | 7197 → 7076; visible 6233 → 5706; captions 567 → 498 | none missing; id `why` dropped (its design boxes now sit in the sections they explain; nothing outside the page linked it) | facts ✓, links ✓, render ✓ (2 equations, 7 figures); settings checked: foldProfile `NBin` 1024 (main 2048), `SubintPeriods` 1, `Assign` 'linear', `NoiseCoeffs` 1, `MaxMemoryGB` 2, `UseSupportedOnly` true; makePredictor `SegLength` 3600, `Degree` 12, `Tol` 1e-9, `MinSeg` 1; delayModel `DeltaT` 69.184, `Moon`/`Planets` true; makeTrajectory `Height` 0, `Interp` 'linear', `Frame` 'ECEF', `VelStep` 1 | 7a.1–7a.5 (each function a numbered section instead of its own h2), 7.1–7.7. **Captions shortened**: 7a.1 (hand-written chain), 7a.2 (gen_dopsmear), 7a.3 (gen_predseg), 7.1 (gen_fold), 7.3 (gen_nearest), 7.4 (gen_subint); 7.2 (hand-written) unchanged |

**Convention added 30 Oct, applied to every page:** the "Options and outputs" lists (and channels 5b's "Outputs") are now a
collapsed box with that title, built automatically by `tb.py build`. They repeat the function reference (api.html), so they
stay one click away but out of the running text. Visible words after this change: channels 5004, sky 3169, receiver 5544,
fullband 2498, excision 6445, fold 5706.
| simulator | 6566 → 6590; visible 5717 → 5539; captions 348 → 275 | none missing | facts ✓, links ✓, render ✓ (3 figures), outline ✓ (13 boxes placed; "option a" history moved with its design choice into F4); settings checked: `scintOn` false, `scint` 740e6 / 2290, simulateFold `FWHM` 0.05, `NBin` 2048, `NPol` 1, `ObserverNoise` 'model', `SignalGain` 1; adaptiveSubints `TargetSNR` 7, `Weighting` 'equal', `MaxSubints` Inf, `AcrossGaps` false; rfiProfile `MinValid` 1−1e-6, `Lmax` [] (= the clean model's), `ChunkBins` 2^18, `OnPulse` 0 | Letter prefix F (stages 1–7 in one step): F1–F12, F4.1–F4.2, Figures F.1–F.3. Design choices "RFI from a measured profile" and "the fold carries the observer's model" folded into F4; id `limits` now on the paragraph in F12. **Captions shortened**: F.1 (gen_simdraw), F.2 (gen_adaptive); F.3 (gen_simvschain) unchanged apart from the label |
| toa | 3435 → 3538; visible 3335 → 3338; captions 259 → 221 (already compact; the gain is structure) | none missing; id `outputs` now on the paragraph after the collapsed options | facts ✓, links ✓, render ✓ (79 KaTeX, 3 figures); settings checked: estimateTOA `NoiseModel` 'radiometer', `Weighting` 'equal' (main 'optimal'), `MaxIterations` 10, `MinCoverage` 1, `Upsample` 8, `OnPulseFrac` 0.01, `MaxRedChi2` 2, `Chi2Alarm` 1e-3; `noiseCal` false | 8.1–8.7, equations (8.1)–(8.7), Figures 8.1–8.3; "In short" box turned into a short introduction. **Captions shortened**: 8.1 (gen_toa), 8.2 (gen_chi2), 8.3 (gen_chanweights) |
| detection | 2064 → 2270; visible 2021 → 2120 (+5 %: a short page that gained an introduction and a summary) | none missing; new ids `question`, `noisecheck`, `chans` | facts ✓, links ✓, render ✓ (48 KaTeX, 4 figures); settings checked: detectPulsar `PFA` 1e-3, `Phase` 0, `Weighting` 'equal', `MinCoverage` 1, `OnPulseFrac` 0.01; foldToTOA `Weighting` 'equal', `Verbose` true | Letter prefix D: D1–D7, equations (D.1)–(D.5), Figures D.1–D.4. **Captions shortened**: D.1 (gen_rice), D.2 (gen_det), D.3 (hand-written gate), D.4 (gen_pd) |
| validation | 1887 → 2022; visible 1882 → 2017 (+7 %: introduction and summary on a short page) | none missing; new ids `rule`, `plots` | facts ✓, links ✓, render ✓ (2 figures) | 9.1–9.5, equations (9.1)–(9.2), Figures 9.1–9.2; the two "why" boxes and the "In short" box folded into 9.1 and the introduction. **Captions shortened**: 9.1 (gen_pulls), 9.2 (gen_qq) |
| navigation | 6392 → 6391; visible 5486 → 5191; captions 526 → 502 | none missing | facts ✓, links ✓, render ✓ (6 figures), outline ✓ (9 boxes, 10 tables in place); settings checked: navFit `Start` "self", `Robust` true, `OutlierFraction` 0.02, `OutlierScale` 4, `MaxIter` 50, `Predict` false; navFilter `ClockPrior` [0 0], `Height` true, `PositionPrior` 1e7, `VelocityPrior` 0.3, `VelocityNoise` 1.2e-7 | 10.1–10.5; six design boxes folded into their sections (outliers by redundancy → 10.1, mixture and clock-as-noise → 10.2, averaged row → "A TOA is a time-average", ship frame / port clock / same as batch → 10.4). **Captions shortened**: 10.1 (gen_navsignal), 10.2 (gen_navavg), 10.3 (gen_navfit), 10.4 (gen_navplan), 10.5–10.6 (gen_navship) |

**Stage pages done (30 Oct):** channels, sky, receiver, fullband, excision, fold, simulator, toa, detection, validation, navigation.
Every one: no number, code name, link or linked id lost; all links resolve; rendered without KaTeX or console errors.
| physics | 5220 → 5201; visible 4980 → 4960 (dense background chapter; structure gained, length kept) | none missing | facts ✓, links ✓, render ✓ (103 KaTeX, 7 figures); DM 25 numbers checked against the sky page (40.52 / 72.03 ms, sweep 31.51 ms) | Letter prefix P: P1–P10, equations (P.1)–(P.11) incl. (P.7) inside its collapsed box, Figures P.1–P.7. Ten "In short" boxes and seven design boxes folded into the running text, with a bridge at the end of each section. **Captions shortened**: P.1 (gen_navgeom), P.2 (gen_whysquare), P.3 (hand-written sweep), P.4 (gen_radiometer), P.5 (gen_foldgrow), P.6 (gen_flanks), P.7 (gen_scint) |
| localisation | 3071 → 3226; visible 2866 → 3021 (+5 %: summary and bridges) | none missing | facts ✓, links ✓, render ✓ (3 figures) | Letter prefix L: L1–L6, equation (L.1) (the residual, now a numbered display equation), Figures L.1–L.3; "In short" box → introduction; the frames design box folded into L2; the consistency box folded into the end of L5 as the bridge to L6. **Captions shortened**: L.2 (gen_timesys), L.3 (gen_delays); L.1 (hand-written chains) unchanged apart from the label |
| conventions | 2187 → 2059; visible 2153 → 2025 (−6 %) | none missing | facts ✓, links ✓, render ✓ (3 figures); seed rules checked against runScenario (500 + pulsar, 600 + pulsar) and runShip (10 + seed, 1000 + seed); `tests/make*Reference.m` exist (fold, rfi, toa) | Letter prefix C: C1–C10, Figures C.1–C.3; nine design boxes folded into their rules. **Captions shortened**: C.1 (hand-written truth diagram), C.2 (hand-written info chain), C.3 (gen_timeprec) |
| pipeline | 2979 → 3009; visible 2883 → 2809 | "72" → **72.5 GB**: gen_datarate.py gives 60.8 / 63.5 / 72.5 GB for full band / channels / channels + excision at L = 1 s, and the page's own "14.5 GB at 0.2 s" is 72.5/5; "about 72 GB" was truncated. *Final pass: the same "about 72 GB" on conventions and simulator* | facts ✓, links ✓, render ✓ (4 figures, 30 clickable boxes) | Letter prefix O: O1–O6, Figures O.1–O.4; six "why this order" boxes → paragraphs in O4; the processIQ design box → O5 text. **Captions shortened**: O.1–O.3 (hand-written diagrams), O.4 (gen_datarate) |
| parameters | 2606 → 2645 (light restyle, reference page) | none missing | facts ✓, links ✓, 10 tables in place; **every value re-checked 30 Oct** against `pipelineParams.m` (unchanged since 95a93bc) and `main.m`'s switches | "In short" box → introduction; the excision-files design box → a sentence under Files; figures numbered 1–2 (reference pages use plain numbers); captions unchanged apart from the label (gen_freqplan, gen_timescales) |
| scenario | 4328 → 4357; visible 4174 → 4203 (a results page: 7 tables carry most of it) | none missing | facts ✓, links ✓, render ✓ (2 figures), 7 tables and 3 boxes in place | Letter prefix R: R1–R9 with R7.1–R7.3, equation (R.1), Figures R.1–R.2; the two "real profile" warning boxes became subsections R7.2 and R7.3 (they are results, not warnings); three "In short" boxes and the worst-case design box folded into the text. Captions unchanged apart from the label (gen_scen, gen_visibility). Checked: every page's closing "next page" sentence follows nav.js's order |
| index | 1854 → 1874 (landing page: targeted edits only) | "27" → 30 (status date); "72" → **72.5 GB** and 14 → **14.5 GB** in the quick start (as pipeline and parameters) | facts ✓, links ✓ | "Start here" describes the new chapter structure of the stage pages; E0 row mentions the docs rewrite; "checked against" now 77bf36d (code unchanged since 7274ced apart from a comment in runAllTests.m); the hero figure has no caption and stays unnumbered |
| scripts, tests, results, lessons, limits, roadmap, formulas, glossary | ±0 (light restyle, `light.py`) | none missing on any of them | facts ✓, links ✓; rendered: formulas (307 KaTeX, no errors), tests (2 figures) | Figures numbered 1… per page (captions unchanged apart from the label: gen_experiments, gen_testmatrix, gen_honesty, gen_neighbour, gen_effects, gen_roadmap and the hand-written ones); every "In short" box unwrapped into plain text; "Design choice:" dropped from box titles. Tables and wording unchanged |
| history, maintenance, api | records and generated: only the 30 Oct entries (final pass); api unchanged | — | links ✓, rendered | maintenance: page format and writing style now describe the textbook form; changelog row. history: session entry, decision row, log row |

## Final pass (30 Oct)

1. **Reference pages rendered** (scripts, results, lessons, limits, roadmap, glossary): no KaTeX or console errors, figures
   numbered, no "In short" left. Together with the earlier checks, every page of the site has now been rendered.
2. **Disk per second made consistent** with `gen_datarate.py` (60.8 / 63.5 / 72.5 GB): "about 72 GB" → 72.5 on conventions,
   simulator (twice), roadmap; limits "~63 / ~72" → 63.5 / 72.5. The full band stays "61" (60.8 rounded). The 0.2 s values on
   parameters (12.2 / 12.7 / 14.5 GB) were already right.
3. **Figure scripts brought in line.** All 72 `gen_*.py` now write the page's caption word for word, 71 of them with its label (the history timeline is
   unnumbered). 67 scripts changed; the four 5a scripts already matched. Captions that compute numbers keep their f-string fields: the field values were recovered
   from the pre-rewrite page and put back where the shortened caption still states them (navfit's `err[99]` and
   `est/err` and roadmap's date by hand; roadmap's caption also took the page's hand edit "the review week E0 (after the
   review of C and D)"). **Proof:** every script was run from the scratchpad and its output compared with the page: 73 of 73
   captions and 73 of 73 SVG bodies identical. The figures README now states the label convention, and its table names each
   figure's current label (it still had the old positions such as "physics.html (7)").
4. **Cross-page numbers** collected by script (taps, sweeps, channel numbers and widths, bins, phase bins, tests, daily sine,
   leakage, η, memory limits, TOA hours). All agree or differ for a stated reason (DM 5 worked numbers; the 1024-bin
   historical runs; the 25 tests of 23 Oct). **One fix:** conventions said "dispersion kernels of 25 million taps" (the DM 5
   kernel, on a page that is not a stage page); it now says 126 million taps at the default DM 25.
5. **maintenance.html** describes the textbook form (page format, writing style rules 2 and 7, the reference example) and has
   the changelog row; **history.html** has the 30 Oct session entry, the decision and this log.
6. **Links:** 1052 checked, none broken.

**Open, for Jasper (not changed):**
- *"In short" headings.* The 5a model page has an `<h3>In short</h3>` opening each stage, and so do sky, receiver, fullband,
  the rest of channels and excision. Fold, toa, detection, validation and navigation open each stage with a plain paragraph
  instead. Both read fine; one form everywhere would be more consistent. Recommendation: keep the plain opening paragraph
  (it reads more like a textbook) and drop the heading, keeping its id on the paragraph.
- *Roadmap figure.* Figure 1 on the roadmap page is the roadmap "on 23 October" (`TODAY` in `gen_roadmap.py`), while the
  page's status is 27 Oct (E0). Regenerate it with updated lists at the next roadmap change.

Nothing is committed. The rewrite is in the working tree; the tag `docs-before-textbook` is the restore point.


## Later the same day: paragraph flow (draft, not committed)

The rewrite above is committed as `7fccde4`. Jasper then asked for more flow between paragraphs ("the RFI excision page still
reads like a collection of paragraphs, not a complete reasoning"). Drafted on two pages, uncommitted in the working tree:

- `excision.html` (Jasper's example) and `toa.html`: every paragraph opens from the one before; bold-label openers ("Blank,
  do not subtract.", "Caveat.", "Test:") and the excision "In short" headings became sentences; the E6 finding list became the
  argument it is; every section and subsection ends on the question the next one answers, and the next opens by answering it
  (second pass, after Jasper: "right direction, not quite there yet … also the links between E4.2 and E4.3").
- Check: every number, code name, link and id kept (token diff against `7fccde4`); both pages render without KaTeX or
  console errors. Words: excision 7431 → 7993, toa 3538 → 3778.
- Open, for Jasper: his verdict on the second pass, then the rollout to the other 26 pages ("ask me next time"). Not changed:
  the "5 Oct 2026" date in toa's error-bar box (running text should have no dates).
