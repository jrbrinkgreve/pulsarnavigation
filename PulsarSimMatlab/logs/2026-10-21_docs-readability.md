# 2026-10-21 — Docs readability rewrite

Docs only; no code changes. Parallel to the validation runs on the frozen copy (see 2026-10-20_blockC.md).

## State at start

Jasper: the docs site had become "compact and comprehensive, but very hard to comprehend". The 13 Oct page format
(In short → physics → plain steps → design choices → tests) was still in place; the prose inside it had drifted, mostly
because every code change since 13 Oct appended a dated paragraph instead of being worked into the text.

Measured on the 25 hand-written pages: sentences were short (12–16 words), but each carried several facts —
5–9 parentheses per 100 words, stage codes in the running text (simulator 21, excision 17), dates (tests 34, scripts 28,
excision 20), slash lists ("1.02 / 1.04 / 1.08 at 10⁻³ / 10⁻⁴ / 10⁻⁶"), duplicated step lists, and sections that started
with the detail instead of the problem.

## Work done

- **Restore point:** git tag `docs-before-readability` on 305af83 (all 27 pages before the rewrite; local only).
- **Writing style** agreed (plan approved, then the excision page as the pilot: Jasper "much better") and written into
  `docs/maintenance.html#style`: one fact per sentence; each section = problem → what goes wrong → what the code does →
  what it buys; no dates, stage codes or "decided by" in the explanation (collapsed `History` boxes instead); define
  terms before use; every number gets a sentence saying what it means, tables instead of slash lists; full sentences;
  say each thing once; tests = one conclusion sentence + table.
- **Pages rewritten** (one commit each, 21 Oct): excision (pilot, full rewrite), index, physics, sky, receiver,
  channels, fold, TOA, detection, full band, simulator (full rewrite), scenario, limits, results, conventions, pipeline,
  parameters, maintenance.
- **Check per page:** a scratch script compares every number and code name in the page text before and after; every
  item that disappeared was a deliberate change (a date, "pre-8-Oct", "MSP" spelled out, "at 800 MHz"). Rendering
  checked on the docs server (KaTeX errors 0, no raw TeX, all in-page and cross-page anchors resolve, phone width
  without sideways scroll on the excision page).

## Content corrections found while rewriting

- index: status block was still "13 October / review week" → brought to 21 Oct (block C in progress, processIQ and
  A3b rows, simulator card). Disk use: the default path (channels + excision) writes ~72 GB per simulated second
  (main.m's own formula), not the 61 GB of the full band; 0.2 s needs ~14 GB.
- physics: section 10 said scintillation is "not yet simulated" → it is in the fast simulator (not the voltage chain).
- limits: adaptive sub-ints were still "next" (done); J0437's wings "not modelled yet" (8-component fit in the fast
  simulator since 21 Oct); the flat-spectrum box listed per-channel thresholds as a possible remedy (exist since 21 Oct;
  the realistic scenario has not been rerun with them yet — stated so on the excision and limits pages).
- fold: "an adaptive choice of sub-int length is on the roadmap" → `adaptiveSubints` exists (observer side).
- scenario: three budget tables in a row with no indication which is current → a note says so (B0329+54 five
  components, J0437 eight components; the single-Gaussian table is the reference).

## For Jasper to check

- Two statements are new in wording (facts from the logs): excision/limits "measured before the per-channel
  thresholds existed; the realistic scenario has not been rerun with them"; excision History "tested in Claude's run,
  not yet in Jasper's own validation run".
- tests.html: test counts differ between sections ("all sixteen", "All 14 tests pass … 13 Oct", "the fourteen unit
  tests") — a record of different dates, not changed.

## Second round (Jasper: "continue working like this for the full site")

- tests.html: every test starts with its conclusion; sub-tests as numbered lists (testSimulateFold 1–15 under
  subheadings), paired numbers as tables (slash lists 67 → 1). testRFIChain run time ~2.5 h (two cases of 72 and 78 min,
  log of 20 Oct), not ~1 h; the overview figure is the 13 Oct set of 14 tests, now said in its caption.
- scripts.html: dense paragraphs split, H0 false alarms and the weak-regime sweep as tables, history into boxes;
  `runScenario` variants were still "follow once testRFIChain is in" → all variants done (scenario page).
- roadmap.html: status brought to 21 Oct — `processIQ` row had badge "next" though done 19 Oct; block C row said C3a/C3c
  "written" while the status paragraph said validated → rewritten from the tests table (Jasper's runs: testSimulateFold
  1–15, testRFIModel 1–5; Claude's runs: expSimVsChain, expSimRFI, testRFIChain, scenario); status list for 21 Oct;
  21 Oct timeline row; In short no longer lists the J0437 profile and detectRFI thresholds as open decisions. Figure
  caption (said thresholds, lags and J0437 adoption were open, while the figure's own items are marked done) fixed in the
  page and in `docs/tools/figures/gen_roadmap.py`.
- validation, lessons, glossary, history: light pass. history's In short said the review week "goes through" the stages
  (present tense); it now says the page covers sessions up to 13 Oct and points to the roadmap timeline and logs.
- Whole-site check: 592 links resolve, KaTeX errors 0 on every page, balanced details/table tags.

## Open items

- `tests/runAllTests.m` help text says testRFIChain takes ~1 h; the full run is ~2.5 h (code comment, not edited).
- history.html stops at 12–13 Oct (sessions and decision log); 14–21 Oct are only in the logs and the roadmap timeline.
- lessons.html "Findings that changed the design" has nothing after block B (e.g. the lag cut-off, detectRFI false flags on
  narrow carriers, bin resolution at low SNR).
- After each code change, weave the change into the page in this style rather than appending a dated paragraph
  (memory updated).
