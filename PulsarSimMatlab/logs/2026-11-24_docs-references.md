# 24 Nov 2026: references for the docs, stage 1 (inventory and verification; no docs page changed)

Jasper: "can you go through the documentation, find, and double-check those found references for the information on the
pages? examples are rice' formula, the FFTFIT paper, everything. You can add the citation as [1] etc, with all citations
used in the docs in a new page. Also keep track yourself of these references". Plan approved the same day, in three
stages: (1) inventory and verification, (2) the references page and two pilot pages, (3) the rollout. His two choices:
named results and literature numbers get a citation, not general textbook statements; a mismatch between a source and a
page goes on a list and the page stays as it is.

This log covers stage 1. No file in `docs/` or the code was changed, MATLAB was not run, nothing is committed.

## What was done

- **Inventory.** All pages except `api.html`, `maintenance.html` and `history.html` were searched for named methods,
  theorems, data sources and numbers from the literature, and the code headers for the sources they already name. The
  pages now cite about a dozen sources as loose text; about 70 are needed.
- **Verification.** 69 references are in the new register, `logs/references.md`, each with a state:

  | State | Count | Meaning |
  |---|---|---|
  | A | 17 | text opened, the statement the docs rely on found in it |
  | B | 9 | record confirmed, content confirmed through a second source |
  | C | 30 | record confirmed, cited only as the origin of a named method; text not opened |
  | D | 13 | record or content not found yet; not to be cited as checked |

- **How.** Records from Crossref (by DOI) and the arXiv listing. Texts read from the arXiv HTML versions, the Internet
  Archive (Rice 1945), JPL's and NAIF's own pages, and one scanned report (Van Dierendonck et al. 1984, rendered page by
  page). Three page fetches by a helper; two of them failed and were replaced by the routes above. Where a number on a
  page is the source's number rescaled, the arithmetic was redone by hand and is marked "own arithmetic" in the register.

## Mismatches and open points (nothing changed for these)

Each line: where, what the page or code says, what the source says.

1. **Rice's formula is in the 1945 instalment.** `detection.html:73` says "Rice (1944)". The level-crossing result is in
   Part III, Section 3.3, eq. (3.3-14), printed in Bell Syst. Tech. J. 24(1), January 1945. The 1944 instalment holds
   Parts I and II. The reference will name both instalments.
2. **"The two largest terms" of TDB − TT.** `localisation.html:492`, `formulas.html:102` and `fold.html:339` give
   1.657 ms·sin g + 14 µs·sin 2g as the two largest terms. In Kaplan 2005, eq. (2.6), the second largest is a 22 µs
   term with another argument (575.3385 T); the 14 µs term is third. The page's "up to 60 µs with the two terms" was not
   checked.
3. **`L_C` is attributed to the wrong resolution, in the code.** `functions/loadEphemeris.m:33`: "L_C = 1.48082686741e-8
   the IAU's mean rate (Resolution B3, 2006)". The value is that of IAU 2000 Resolution B1.5 (Irwin & Fukushima 1999),
   as quoted in Kaplan 2005. Resolution B3 of 2006 defines TDB through L_B = 1.550519768·10⁻⁸. The number is right; the
   comment shows in `api.html`.
4. **The dispersion constant has two conventions.** The pages use K = 4.148808·10³ s MHz² pc⁻¹ cm³, the Handbook's
   value. Kulkarni 2020 reports that pulsar-timing work fixes it at 1/(2.41·10⁻⁴) = 4.149378·10³ instead, which is
   1.4·10⁻⁴ larger. A catalogue DM made with that convention then gives a sweep that differs by about 4.6 µs across
   1.2–1.6 GHz at DM 26.8 (own arithmetic). It cancels in the simulator and matters for real data. Which convention
   the ATNF values follow was not checked.
5. **B0329+54's scintles.** `scenario.html:307`: "about 7 MHz × 13 min (ATA 2025, Wang et al. 2005)". Wang et al. give
   14 MHz and 16.9 min at 1540 MHz. Brown et al. 2025 (the "ATA" paper) give 2.6 MHz at 1140 MHz and 14.7 MHz at
   1668 MHz, so 7 MHz at 1.4 GHz is an interpolation (notes §11), and the 13 min (780 s) was not found in what was
   read. The sentence reads as if both papers gave 7 MHz × 13 min.
6. **"HWHM" and "1/e" for J0437's scintles.** `parameters.html:262` and `pipelineParams.m:107` call 740 MHz the HWHM
   and 2290 s the 1/e time. Neither Sosa Fiscella et al. 2021 nor Keith et al. 2013 states these definitions in the text
   read. Wang et al. and Gwinn et al. do use them for their own values.
7. **DE440's accuracy for the Earth.** `fold.html:351` ("accurate to centimetres"), `fold.html:408` ("well below a
   metre") and `localisation.html:491` ("a few centimetres"). In the part of Park et al. 2021 that was read, the numbers
   are: DE440 and DE441 differ by less than 1 m, the inner orbits are tied to ICRF3 to about 0.2 mas, and ranging
   residuals are about 0.7 m for Mercury. No accuracy for the Earth was found. Not confirmed, not refuted.
8. **Two citations are informal.** "IAR 2020 (arXiv:2010.00010)" is Sosa Fiscella et al. 2021, ApJ 908, 158. "ATA
   2025" is Brown et al. 2025, ApJ 995, 120.
9. **The Handbook's year.** Its copyright page says 2005; it is often cited as 2004. The register uses 2005 and says so.
10. **The catalogue version.** `scenario.html:16` says "ATNF catalogue v2.6.5". Not checked; the research of 21 Nov
    used v2.7.0.
11. **The GNSS levels.** `receiver.html:589` says they are "the published minimum received power + ~3 dB". Only GPS
    L1 C/A (−158.5 dBW) was confirmed, through a secondary source. The other five signals, the Inmarsat band and the
    radar characteristics are in state D.
12. **Already known from the research of 21 Nov** (`reports/Pulsar navigation project research.md`, section 2), not
    repeated here: J1713+0747's W50 of 0.30 ms on the scenario page against 0.112 ms measured, and three catalogue rows
    that lean on one paper.

Confirmed exactly, for the record: B0329+54's five components in `runScenario.m:118` equal Table 2 of Shang et al. 2017
number for number, and the code is right to treat the widths as standard deviations in degrees. The clock preset
`h = [1.8e-21 6.492e-22 1.37e-24]`, "eqs. 6–8", "eq. 60", the 8.50 ns and the 11 ns at 200 s are all as in Van Dierendonck
et al. 1984. J0437's 740 MHz and 2290 s are a verbatim quote. "Stated good to about 10 µs" is Kaplan's sentence.

## Stage 2, the same day (Jasper: "Yes, stage 2", and "let's fix the 13 open ones afterwards")

- **`docs/references.html`** (new): 50 numbered references in the order in which the pages first cite them along the
  sidebar, each with its link and one line that says how far it was checked ("Source read" 17, "Through a second
  source" 8, "Record only" 25) and what was confirmed. A last section lists what has no confirmed source yet, without
  numbers. These are the site's first links to the internet; they are on this page only.
- **Numbers** are fixed from today (register, table "Numbers"). The 13 open references get [51] onward when they are
  confirmed. Six helper sources are not numbered; two of them are named inside the entry they support.
- **Sidebar:** one line in `docs/nav.js`, after Glossary. **Style:** four short rules in `docs/style.css` (`a.cite`,
  `.ref`, `.refnum`, `.refnote`); the plan said one rule, the entries needed a hanging indent as well.
- **Pilot citations**, six in all. `toa.html`: FFTFIT [7], the Kolmogorov–Smirnov test [39], the sandwich estimator
  [40], and "(arXiv:1703.03582, Table 2)" is now "(Shang et al. 2017 [38], Table 2)". `detection.html`: the
  Neyman–Pearson test [6] and "Rice (1944) [8]". The year on the page stays as it was (mismatch 1); the entry names both
  instalments and says that the formula is in the 1945 one.
- **Details trimmed from the entries because they were not confirmed:** the title of Fairhead & Bretagnon 1990, the
  end pages of White 1980 and Manchester et al. 2001, a publisher and a web address for IS-GPS-200.

Checks: 1222 links, none broken; 50 entries numbered 1 to 50 without a gap; every citation's number equals its entry's;
the text of both pilot pages differs from the last commit only by the six insertions. In the browser (local server,
dark theme): the sidebar shows References between Glossary and Docs maintenance, previous and next links are right, the
entries hang under their numbers, [7] on the TOA page lands on its entry, no console or KaTeX errors, no sideways
scroll at phone width. The light theme was not looked at; the new rules use the site's colour tokens only.

Not done yet: the other 44 citations (stage 3), the "cited in" back-links of those entries, the paragraph on references
in `maintenance.html`, and the 13 open references.

## Stage 3, the same day (Jasper: "Yes, roll out")

- **116 citations on 19 pages**: physics, localisation, conventions, parameters, sky, receiver, fullband, channels,
  excision, fold, simulator, toa, detection, navigation, scripts, limits, scenario, formulas, glossary. All 50 numbered
  references are cited at least once.
- **Where a citation sits.** In the running text or a table cell, right after the name or number it supports. On the
  formulas page in the first cell of the row; in the glossary in the term. Not in headings and not in figure captions,
  so no figure script had to change.
- **No sentence was changed.** The informal citations keep their wording and get a number: "IAR 2020 (arXiv:2010.00010)
  [10]", "ATA 2025 [47]", "Rice (1944) [8]". The one rewording stays the pilot's: "(Shang et al. 2017 [38], Table 2)".
- **Numbers.** The first three were swapped (Lorimer & Kramer [1], Hankins & Rickett [2], Dicke [3]), because the first
  numbering counted a table-of-contents line as the first appearance. Nothing had been committed. The numbers now follow
  the true order of first appearance on all pages, checked by script.
- **`references.html`**: every entry lists the pages that cite it, each link going to the nearest section anchor.
- **`maintenance.html`**: a "References" item under "How the site is built" (the markup, the numbering rule, the register,
  where a citation may sit), and the changelog row of today extended.
- Two planned citations were dropped for lack of a place in the text: `price2016` (no fitting sentence on the channels
  page) and a second one for pulsar spectra on the limits page. [36] (Fairhead & Bretagnon) is cited next to [13] on the
  fold page, where the USNO series is named.

Checks: 1443 links, none broken; 50 entries numbered 1 to 50; every citation's number equals its entry's; 109 text
changes against the last commit on the cited pages, each one an added citation (script); numbers in order of first
appearance (script). In the browser, with the local docs server: formulas (329 KaTeX elements, no error), glossary and
physics render with their citations, no console errors. Not looked at one by one: the other 14 pages.

Still open: the 13 unconfirmed sources (register, state D), and the twelve mismatches above, which are Jasper's to decide.

## The open sources, the same day (Jasper: "let's fix the 13 open ones afterwards")

Nine got a number, [51] to [59], appended after the fifty: Rao 1945 and Cramér 1946 (the Cramér–Rao bound), Kay 1993 (its
textbook treatment), the WGS84 standard, Meeus 1998 (nutation), Giorgini et al. 1996 (JPL Horizons), ITU Resolution 222
(the Inmarsat band), Kaiser 1974 (next to [29]) and Huber 1967 (next to [40]). 21 more citations on 9 pages; 137 in all.

How far each was checked is weaker than for the first fifty, and the entries say so. Six are "Record only", and for five
of those the record itself comes from a secondary reference list, not from Crossref (Rao, Cramér, Kaiser, Huber, Meeus).
Three are "Through a second source": Kay (the publisher's table of contents), WGS84 (the ellipsoid constants in a copy of
the specification, equal to `makeTrajectory.m:86`), ITU Resolution 222 (a search summary; its full title was not
confirmed). None of the nine was opened in the original.

Settled without a number: the catalogue version (2.6.5 exists, 26 Aug 2025; noted inside [45]), Irwin & Fukushima 1999
and the IERS Conventions 2010 (records confirmed, no page names them), and two textbooks that turned out not to be
needed.

Still without a confirmed source (register, last table; References page, section 3): the GNSS levels other than GPS
L1 C/A, with four carrier frequencies found in secondary sources only; the radar characteristics; the troposphere's
delay, the sky temperature and the constant 40.3; a textbook source for the matched filter.

Mismatch 10 is closed: "ATNF catalogue v2.6.5" is a real version. Version 2.7.0 came out on 5 Oct 2025, so the page did
not use the newest one.

Checks after this step: 1488 links, none broken; 59 entries numbered 1 to 59; 137 citations, each number equal to its
entry's; every numbered entry cited and listing its pages.

## Not done in stage 1

- The 13 references in state D (register, last table).
- The content of the 30 references in state C was not opened. They are cited for a name only.
- Texts behind a paywall were not read: Taylor 1992, the Handbook, Brown et al. 2025 beyond its abstract.
- No page was read line by line a second time. A statement that needs a source but uses none of the searched names
  can still be missing from the inventory.

## Next

Stage 2 after Jasper's look at the register: `docs/references.html`, the sidebar entry, and the citations on
`toa.html` and `detection.html` as a pilot. Numbers `[n]` are given then.

## Files

- New: `logs/references.md` (the register), this log.
- Not touched: `logs/2026-11-02_risk-register.md` and `logs/2026-11-24_detection-literature.md`, which carry another
  session's uncommitted work.
