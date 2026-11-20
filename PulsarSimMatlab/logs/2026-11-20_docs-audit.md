# 20 Nov 2026: read-only audit of the docs site (no file in docs/ or the code changed)

Jasper: "can you write a prompt for yourself to audit the docs while im away for a bit", then "Yes, run it".
The audit compares `docs/` with the code at `cb05d21`, the notes and the logs. Nothing was fixed, MATLAB was not run,
nothing was committed. Every finding quotes the page with its line in the HTML file and names the source it was checked
against. "Checked" means compared with the source or recomputed by hand. "Suspected" means not traced to a source.

## Summary

| Severity | Count | Meaning |
|---|---|---|
| wrong | 1 | the page states something the code contradicts |
| stale | 20 | true when written, overtaken by a later decision or commit |
| inconsistent | 4 | two places in the docs disagree, or a list is incomplete |
| broken | 0 | links, anchors, nav, figure and equation numbering are all clean |
| minor | 12 | wording, rounding, a date in running text |
| code comments | 3 | the fix is in a `.m` header, not in docs; two of them show in `api.html` |
| already known (log of 13 Nov) | 12 | waiting for the navigation redesign; listed once at the end |

Counts are per bullet below. The stale ones cluster: twelve of the twenty are the ship and site wording on four pages.

The five that matter most:

1. **`parameters.html:226` gives the wrong DM limit.** It says "about DM 30 is this machine's limit". The sizing rule of
   `applyDispersionStream` allows the 2^28 FFT (11.8 GB) only up to DM 26.6. B0329+54's own DM, 26.8, needs 2^29 =
   23.6 GB, above the 16 GB default budget and at the edge of the 24 GB machine.
2. **The static site at an unknown place is on two pages only.** `index.html` and `roadmap.html` name the objective of
   2 Nov. `scenario.html` (lede and R6), `localisation.html` and the end of `navigation.html` still say that the array
   is on a ship, that the site moves, and that the endpoint is a fixed position in the Netherlands. Only part of this
   is in the known list of 13 Nov.
3. **The pages read as if millisecond pulsars were optional.** They are part of the objective, and have been from the
   start. `scenario.html:202` calls J0437 "the only realistic millisecond pulsar", and six places park J0437's binary
   orbit "because J0437 never rises at 52° N". With the site anywhere, that reason no longer holds.
4. **Two test rows and the home page disagree on who ran what.** `tests.html:447-448` gives `testNavFitFaults` and
   `testFrontEndFaults` as Claude's runs of 2 and 3 Nov. `index.html:122` and `roadmap.html:84` say all 33 tests passed
   in Jasper's run of 5 Nov. One of the two is wrong, and "validated" depends on it. This needs Jasper's word.
5. **`history.html` stops at 4 Nov.** Sessions, decisions and the log list miss 5 to 20 Nov, including the decisions to
   go static-only and archive the ship code, and the 4×4 array. The decision
   log also lacks the objective change of 2 Nov itself.

What held up: all 1078 links and anchors, the 29 nav entries, figure and equation numbering, every options list but
one against the `arguments` blocks, `api.html` against the 43 functions, the test tables against the test files, and
every worked number I recomputed on the sky, receiver, full-band, channels, excision, fold, TOA, detection, scenario
and formulas pages. No DM-5 number without its DM was found on the pages read in full.

## Findings per page

Severity in bold. "Fix" is one line and is a suggestion only.

### index.html

- **stale, checked.** Line 119: "with all 18 data streams beamformed". Line 148 of the same page says "4×4 array" and
  `roadmap.html:31` says 32. The update of 20 Nov (`logs/2026-11-20_array-4x4.md`) changed the block E row but not this
  sentence. Fix: 32.
- **stale, checked.** Line 94: heading "Status · 30 October 2026". The text under it gives the status of 13 Nov
  (line 124), 5 Nov and 2 Nov. Fix: date the heading as the text.
- **stale, checked.** Line 218: "Docs checked against commit 77bf36d (30 Oct 2026; the code is unchanged since
  7274ced)". `functions/` changed in `c4aef71` (5 Nov), `f1397ce` (17 Nov) and `982d045` (19 Nov). The maintenance
  checklist names this line (item 8). Fix: update after the next check.

### parameters.html

- **wrong, checked by hand.** Line 226, DM row: "The forward dispersion still fits the 2^28 FFT (~12 GB) that the DM-5
  runs used; about DM 30 is this machine's limit." The rule is `nMin = 2^nextpow2(2M - 2)`
  (`functions/applyDispersionStream.m:328`, design grid `:390`), 44 bytes per FFT sample (`:156`), budget 16 GB
  (`:101`), error above it (`:357`). By hand: DM 25 gives 126,068,856 taps and 2^28; DM 26.6 still 2^28; DM 26.8 gives
  135,144,373 taps and 2^29 = 23.6 GB; DM 30 the same. `sky.html:772` has it right ("At DM 30 this stage needs about
  24 GB"). Fix: "DM 26.6 is the limit of the 11.8 GB step; up to DM 53 the stage needs 23.6 GB".
- **stale, checked.** Line 15: "as of commit 95a93bc (21 Oct 2026, its last change; checked 22 Oct and again 30 Oct)".
  `pipelineParams.m` last changed in `44f5679` (20 Nov, comment lines 64-66 only; all values are the same). Fix: name
  the commit and say the values did not change.
- **minor, checked.** Lines 354-364, "Where main.m overrides function defaults". Since 19 Oct `main.m` calls none of
  these functions but `addNoiseAndRFI`. `processIQ` passes `NBin` and `TRef` (`functions/processIQ.m:198-199`), `T0`
  (`:143`) and the bin rate (`:160`); `main.m:158` passes `Weighting` and `Bnoise` to `foldToTOA`. The values are
  right. Fix: "main.m and processIQ pass".
- **stale, checked.** Lines 219-220 and 226: "The scenario runs use B0329+54 from 52° N" and "the pulsar for the
  Netherlands scenario". True of the runs, but the objective of 2 Nov is not mentioned. Part of top finding 2.

### pipeline.html

- **stale, checked.** Lines 245-246: "makeTrajectory describes where the receiver is: a site on the ground, a ship or
  a probe." The function has `'fixed'` and `'custom'` only (`functions/makeTrajectory.m:24-27`). Not in the known list.
- **stale, checked.** Lines 347-349: "runScenario, which draws a whole observing campaign, and runShip, which closes
  the navigation loop on a ship". `runShip` is in `old/moving_platform/`. Line 257 of the same page says so; this
  sentence does not. Not in the known list.

### sky.html

No findings. Recomputed: 25,229,772 and 126,068,856 taps, the design grids and their memory, the block-size table,
Σh² = 0.195, 10^17 multiply-adds, the group delay, the smear in the lowest sub-band. Options match the code.

### receiver.html

- **minor, checked.** Line 586: "3×3 when the list was made, 4×4 since 20 November 2026". A date in the running text
  of an explanatory page; writing-style rule 3 puts dates in a History box.
- Numbers recomputed and right: 388·ρ, −124 dBW, the 9.8 MHz transition, 256 ns, the clock's 0.66 µs per hour, 76 µs
  per day and 13 ms per 30 days. The scenario tables equal `pipelineParams.m:55-95` source by source.
  `h₋₁ = 6.5·10⁻²²` is the code's 6.492e-22 rounded (`functions/makeClock.m:77`).

### fullband.html

No findings. Recomputed: 2^26 and 3.8 GB at DM 25, `nFuture`, the 2^25 case, 374 µs / 75 µs / 1.5 ms in-channel
smear, 0.4994 µs. Options match the code.

### channels.html

- **inconsistent, checked.** Line 945, caption of Figure 5b.1: "a channel's filter stays near 100 MB". Line 889 says
  "about 120 MB" at DM 100 (2^21 × 56 bytes = 117 MB). The changelog of 13 Oct corrected the text, not the caption.
- Recomputed and right: 3853, 8989 and 1285 taps, the channel centres and DFT rows, 27 against 2700 multiplications,
  416,667 samples, DM 33 for the 4 GB step. The two numbers 0.18σ (line 965) and 0.26σ (line 1170) are two different
  runs, both in the notes (§5.16 and §10); not a finding.

### excision.html

- **minor, suspected.** Line 228, caption of Figure E.3: "At the default DM 25 each output mixes about 940 input
  samples." The definition used for "660-1560" (lines 90 and 162), the sweep across the full channel width, gives 984
  at 1.4 GHz. 940 matches the tapered width of 3.0 MHz. Two definitions on one page.
- Options of `detectRFI`, `periodicRFI`, `blankChannels` and `blankingWeights` match the code. Recomputed: the guard of
  10 samples, 80,000, 120×, 450×, 41 ms, 4.8 GB, 59 and 295 µs per channel.

### fold.html

- **stale, checked.** Lines 263-264: "J0437's binary orbit was parked at the close of block D, because J0437 never
  rises at 52° N". The site is not known beforehand since 2 Nov, and millisecond pulsars are part of the objective;
  most of them are binaries. The same reason stands in `navigation.html:857`, `roadmap.html:413`,
  `limits.html:17` and `:163`, and `localisation.html:507`. Fix: one sentence in each place that the reason has gone and
  the item is open.
- **minor, checked.** Line 24: "±8 ms over a day", and line 254 "±7.6 ms per day for B0329+54 from 52° N". Both hold at
  52° N; at the equator the daily term reaches 21 ms.
- **minor, checked.** Line 311, History of `makeTrajectory`: names a `'spline'` option. The function has none now
  (`functions/makeTrajectory.m:37-45`); it went with the `'track'` type.
- Stage 7a reads correctly for a static site. Options of the four functions match the code, and the test tables match
  the test files (3, 11, 8 and 6 numbered checks).

### simulator.html

- **minor, checked.** Line 506: "since 19 Nov" in an options list.
- Not read line by line (see the end). Headline numbers agree with the other pages: 16 h and 82-95 h per J0437 TOA,
  72.5 GB per second, 0.05 s for 20 hours. The options of `adaptiveSubints` match the code.

### toa.html

No findings. Options match `estimateTOA`. Recomputed: the 0.61 µs coarse grid, 1/1700, ±0.008 and 130σ, w/102.

### detection.html

- **inconsistent, checked.** Lines 251-252, options of `foldToTOA`: "Weighting ('equal'; …), Verbose (true)". The
  function also has `ChannelGain` (`functions/foldToTOA.m:70`), which the paragraph at line 241 describes. Fix: add it
  to the list.
- **minor, checked.** Lines 241 and 246: "Since 17 Nov" and "Since 19 Nov" in the running text.
- Recomputed and right: 3.09, 4.15, 3.57, 4.90, the expected false alarms 90 / 9 / 0.9, 0.962, 1700 km.

### validation.html

No new findings. Recomputed: ±0.075 for 90 TOAs, ±25 % for 8.

### navigation.html

- **stale, checked.** Lines 881-882, "Next": "Block E (revised 23 Oct) leads to a fixed position in the Netherlands from
  several pulsars with the joint navFit. The navigation filter and the ship runs are parked". The objective changed on
  2 Nov, the redesign D6a-e waits for approval (`logs/2026-11-13_handover-static-nav.md`), and the ship code is
  archived, not parked.
- **stale, checked.** Line 191: "the prediction from the other TOAs and from the ship's dead reckoning is good to about
  70 km". Section 10.1 is the live text, not the archived Section 10.4.
- The page has no word on the objective of 2 Nov, the D6 proposal, or millisecond pulsars. The 13 Nov log already says
  the body waits for the redesign. The options of `navFit` match the code.

### localisation.html

None of these is in the known list of 13 Nov. All **stale, checked** against `fold.html` (already static) and
`logs/2026-11-13_handover-static-nav.md` (no dead reckoning, no position prior).

- Line 28: "The receiver sits on a ship on a rotating …".
- Line 116: "A millisecond pulsar therefore needs a prior position, for example from dead reckoning, that is good to a
  few hundred km." With no prior, and millisecond pulsars part of the objective, this is an open design question
  (`logs/2026-11-20_array-4x4.md`, last paragraph).
- Lines 121 and 138: "ship's position is given in latitude and longitude", "A ship and a ground site".
- Line 232, table row: "Roemer, the ship | The ship's position on the Earth: what navigation wants to find".
  `fold.html:255` has the same row as "Roemer, the site".

### scenario.html

- **stale, checked.** Lines 30-31: "The array will go on a ship, and J0437 never rises at 52° N, so the scenario runs
  time B0329+54 from the Netherlands". Lines 34-35: "on a ship (R6)". Lines 222-226, R6: "the array will go on a ship.
  The site therefore moves, and a fixed flat installation is not an option … For now the location is the Netherlands".
  The array is static at a site not known beforehand; the ship comes after the static field test. Only "a ship that
  sailed south" (line 226) is in the known list.
- **stale, checked.** Lines 199-207: the faint millisecond pulsars "take them months to years", "J0437 is the only
  realistic millisecond pulsar", "there is no feasible MSP". The numbers are right. The wording reads as if they were out; they are part of the
  objective, and long integration times are accepted.
- **minor, checked by recomputation.** Line 78: the column "ρ [dB]" has no label for the design point. Its values are
  the nominal ones (J0437 −50). The text and the lede speak of −54, the worst case in the same convention.
- **minor, checked.** Line 34: "3×3 until 20 November 2026" in the running text.
- Recomputed and right: both SEFDs, 3.9 and 21 h, every σ_TOA, the figures of merit, every visibility limit, 62 h, and
  the conversions between the two ρ conventions.

### roadmap.html

- **stale, checked.** Line 490: "X and Y of one element share its position, so the 9 element positions are the only
  spatial information." 16 since 20 Nov; missed by that update.
- **stale, checked.** Line 13: "Status as of 13 November 2026", and the dated table ends with "9–13 Nov" (line 85).
  Other parts of the page were edited on 16, 17, 19 and 20 Nov. No row for the direct-sampling note, `ChannelGain`
  (validated 19 Nov), or the 4×4 array. Millisecond pulsars as part of the objective are nowhere on the page.
- Lines 413-414: the binary-orbit reason, see `fold.html`.

### history.html

- **stale, checked** against `git ls-files logs` and `git log`. The sessions end at "4 November" (line 378). The
  decision log ends with two rows of 2 Nov (lines 504-506) and has no row for: the objective change of 2 Nov, static-only
  and the archive (5 Nov), "the run scripts pass no gains" (19 Nov), the 4×4
  array (20 Nov). The log list ends at `2026-11-04_docs-number-audit.md` (line 569); twelve tracked logs are missing:
  `11-04_static-redesign`, `11-05_static-stage7`, `11-09_lo-and-clock`, `11-13_docs-status`,
  `11-13_handover-static-nav`, `11-13_tsys-questions`, `11-16_direct-sampling`, `11-17_budget-recheck-big30`,
  `11-17_uncertainty-methods`, `11-19_channel-gain-step2`, `11-19_idea-earth-motion-baseline`, `11-20_array-4x4`.
  The three `task2.4` logs are ignored by git and should stay out.

### tests.html

- **inconsistent, checked; needs Jasper.** Line 447: `testNavFitFaults` "passes (Claude's run, 2 Nov)". Line 448:
  `testFrontEndFaults` "passes (Claude's run, 3 Nov)". `index.html:122`: "33 since 5 Nov, all passed in Jasper's run of
  that day". `roadmap.html:84`: "33 unit tests, all passed in Jasper's run of 5 Nov". Commit `aa55f7d` re-dated only the
  five stage 7 rows. If the run of 5 Nov was the full `runAllTests`, these two rows are validated and should say so; if
  not, the home page and the roadmap overstate. Asked on 20 Nov, Jasper answered "Not sure". The record does not settle
  it: `logs/2026-11-05_static-stage7.md:51` has "Jasper's own runAllTests, after the push: 'all passed' (his report; I
  did not see the output)". A plain `runAllTests` on that day did include both tests (they were committed on 2 and
  3 Nov, and the runner picks up every `test*.m`), but a subset call would not. Open until he reruns the two.
- **minor, checked.** Line 138: "The table lists every test in the order runAllTests runs it". The order list in
  `tests/runAllTests.m` differs in three places: `testOptimalTOA` before `testBlankingWeights`, `testPeriodicRFI` before
  `testRadarEdges`, `testMakeClock` before `testNavFit`.
- Every file in `tests/` is on the page, and the counts are right: 30 in the list, three fault tests, 33 in the default
  run, `testRFIChain` as the thirty-fourth. `testNavFilter` is marked archived.

### conventions.html

- **inconsistent, checked.** Line 135: for the sky files the check "compares T, f_in, L, A, dutycycle, envelope mode,
  seed, DM, band and refFreq". It also compares the profile components (`functions/checkSkyFiles.m:1`, `main.m:123`),
  as `parameters.html:225` says.

### formulas.html

- **minor, checked.** Line 95: "For a moving receiver, the TOAs must be referred to the barycentre". That holds for
  every receiver on the Earth.
- Recomputed and right: 0.69 ns per TECU, 0.962, the half-turn limits 107,000 km and 863 km, ±31 km, the Allan
  variance.

### glossary.html

- **stale, checked.** Lines 35 and 58: "dead reckoning" and "Kalman filter, smoother" describe `navFilter` in the
  present tense, with no word that it is archived. Same family as line 108, which is in the known list.

### lessons.html

- **minor, suspected.** The table ends at row 19 (27 Oct). The χ² alarm found on 17 Nov that led to `'ChannelGain'`
  (`logs/2026-11-17_budget-recheck-big30.md`) has no row. Whether it counts as a bug is Jasper's call.

### scripts.html, results.html, limits.html, physics.html, maintenance.html, api.html

- `scripts.html`: the table equals the eight scripts in the folder, and `runShip` is marked archived. No findings.
- `limits.html`: lines 17 and 163, the binary-orbit reason (see `fold.html`).
- `api.html`: stamped `44f5679`, 43 functions, the same set as `functions/*.m`. Current; `cb05d21` only regenerated it.
- `results.html`, `physics.html`, `maintenance.html`: swept only, see the end.

## Fixed the same day (Jasper: "Yes, quick fixes")

Text only, on five pages and one figure script. Checked afterwards: 1079 links, none broken; tags balanced on every edited
page. The pages were not opened in a browser by me, and the figure was not regenerated (only its caption changed, in the
page and in the script alike).

- `parameters.html`, DM row: the 2^28 FFT holds up to DM 26.6; from there to DM 53 the stage needs 23.6 GB, so B0329+54's
  own DM does not fit. Commit line: `44f5679`, a comment only, values those of `95a93bc`, checked again 20 Nov.
- `index.html`: 32 data streams; heading "Status · 20 November 2026" and "(status 20 Nov)"; the last line now reads
  "Docs last audited against commit cb05d21", with a pointer to this log.
- `roadmap.html:490`: 16 element positions.
- `detection.html`: `ChannelGain` in the options list of `foldToTOA`.
- `channels.html`, caption of Figure 5b.1, and `docs/tools/figures/gen_dmmem.py`: "a channel's filter needs about 120 MB
  there" (at DM 100).
- `maintenance.html`: one changelog row.

Still open: everything else above. The ship and site wording, the pages' wording on millisecond pulsars, and the history
page need new text, so they get a draft first.

## Fixed later the same day: the static site in the wording

Jasper, on the draft: "the end goal right now is just a static site". So the draft's sentence "A ship is the goal after the
static field test" went, and so did the same clause on the roadmap. Applied to seven pages and one figure script; 1080
links, none broken; tags balanced; not opened in a browser by me.

- `scenario.html`: lede; "at a static site" in R1; R6 is now "The site" (anchor `#ship` kept), with 52° N and 30° S as
  the two latitudes of the runs and the ship in the History box. The reason for "constant gain over a pass" (the ship)
  is gone; the sentence now only says what the simulations assume. R4: "J0437 is by far the fastest millisecond pulsar
  to time", and for the north "B1937+21, the fastest there, needs 670 h of data at nominal flux" (the page's own table).
- `localisation.html`: lines 28, 121 and 232 ("Roemer, the site"); line 116, a millisecond pulsar's first position is an
  open point of the redesign, no dead reckoning; "site, 100 km" in Figure L.3 and in `gen_delays.py`; the binary-orbit row.
- `navigation.html`: "Next", and the binary-orbit row of the parked table.
- `pipeline.html`: `makeTrajectory`'s receivers; the scripts sentence without `runShip`.
- `fold.html`, `limits.html` (twice), `roadmap.html`: the binary orbit's reason. `roadmap.html` also: "the ship comes
  after the static field test" removed, and "return for the ship" is now "return if a moving platform is developed later".
- `logs/current_project_notes.md` §11: two dated updates next to the old statements, which stay.

Left as they are, on purpose:
- Figure L.1 on `localisation.html` (box "ship or probe", caption "A ship and a ground site"). It draws the general
  trajectory interface and has a "site" box of its own, so a relabel would change what it shows. It needs a decision.
- `localisation.html:101-103` ("A moving ship adds its velocity … a Kalman filter") and `navigation.html:191`: how
  navigation works, which waits for the redesign.
- Not asked, not done: whether the binary orbit moves from "parked" to an open roadmap item.

Still open from the list above: `history.html` (stops at 4 Nov), the two test rows, the minor items, and the known list.

## Correction, the same day: millisecond pulsars

Jasper: "MSP's were in the objective from the first moment; thats not a recent decision. Remove notions like 'on ... it
was decided to include MSPs'". The first version of this report called it a decision of 19 Nov, after the wording of the
notes. Reworded here, in `logs/current_project_notes.md` §1 and in one sentence of `logs/2026-11-20_array-4x4.md`. No docs
page had that wording. The commit message of `eb772a0` still has it and stays as it is.

## Code comments that `api.html` shows (the fix is in the code, so not touched)

- `main.m:15-16` and `:39-40`: "clock jitter, pulsar-Earth motion … [todo]", "barycentric / timing corrections,
  residuals [todo]", "navigation solution (multi-pulsar) [todo]". `makeClock`, `delayModel`, `makePredictor` and
  `navFit` exist.
- `tests/runAllTests.m:5-6`: "About 10 min for all 31". The default run has 33 tests and takes about 20 min
  (`tests.html:30`). Not in `api.html`, but the same kind of comment.
- `runScenario.m:6-7`: "J0437-4715 from 30 S (a ship that sailed south)". Already in the known list.

Outside the docs: `logs/current_project_notes.md`, heading of §5.29, says "testMakeTrajectory (7 tests)". The file has
three since 5 Nov.

## Already known: the list of `logs/2026-11-13_docs-status.md`, still open

These wait for the navigation redesign and are not counted as new. Checked today: all are unchanged.

`scenario.html:226`; `pipeline.html` stage table row 10 (line 296) and Figure O.2; the body of `navigation.html`;
`glossary.html:108`; `limits.html:179`; `formulas.html:121`; `conventions.html:291`; `receiver.html:1136`;
`simulator.html:248`; `validation.html:48`; `results.html:656` (needs a line that the code is archived); the roadmap
figure and the block E table (lines 239-247 and 279), which show the plan of 23 Oct.

`logs/2026-10-21_docs-pending.md` holds no open item.

## Not covered

- Not read line by line, only swept by keyword and for recurring numbers: `physics.html`, `results.html`,
  `limits.html`, `lessons.html`, `glossary.html`, the body of `simulator.html`, the per-test sections of `tests.html`
  below its table, the body of `history.html` and of `roadmap.html` below the block table, `scripts.html` below
  `main.m`, and Sections L2 onward of `localisation.html`.
- Check 5 of the prompt is half done. The test table was compared with the test files. The claims of `results.html`
  and `validation.html` were not traced to their logs one by one.
- Test results quoted on the stage pages (for example "1.2·10⁻¹²" on `sky.html:758`) were not compared with the tests'
  printed output; that needs a MATLAB run.
- Figures were not rendered or compared with their generators. The pages were not opened in a browser.
- `api.html` was compared by function set and commit stamp, not help text by help text.
