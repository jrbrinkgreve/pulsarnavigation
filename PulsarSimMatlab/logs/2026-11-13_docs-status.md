# 13 Nov 2026: progress bullets and a docs staleness check (no code change)

Jasper: "can you first give my a list of bullet points for my progress update? Also review the docs whether anything it out of date
there", then "Yes, fix those" (the stale status lines on the home page and the roadmap). No MATLAB run, no code change.

## The check
A keyword sweep of all docs pages (navFilter, runShip, testNavFilter, dead reckoning, 'track', 'barycentric', dates, test counts)
against the code and the November logs. It was not a full read of every page. `api.html` is current: no code changed since c4aef71.

## Fixed today
- `roadmap.html`: the status line (was "as of 27 October", "now at day 2"); a paragraph on the 2 Nov objective and the archived ship
  code; four rows in the dated table (30 Oct to 13 Nov); a note under the block E endpoint; the sentence that said the redesign comes
  "in the week of 9 November".
- `index.html`: the block E paragraph names the 2 Nov objective and the waiting redesign; the ship code is said to be archived;
  31 tests then, 33 since 5 Nov.
- `maintenance.html`: one changelog row.
- Not claimed anywhere: that E0 is closed. The pages say days 1 and 2 are done; whether the review week is finished is Jasper's call.

## Found, not fixed
- `scenario.html:226` and the header of `runScenario`: "a ship that sailed south" for the 30 S site.
- Still describing the archived ship code as live, waiting for the redesign: `pipeline.html` (stage table row 10, Figure O.2), the body
  of `navigation.html`, `glossary.html:108`, `limits.html:179`, `formulas.html:121`, `conventions.html:291`, `receiver.html:1136`,
  `simulator.html:248`, `validation.html:48`; `results.html:656` needs a line that the code is archived.
- The roadmap figure (the calendar SVG) and the block E decision table still show the 23 Oct plan.
- In the logs only, not in the docs: the static redesign proposal, the real-time estimate of 5 Nov, the LO and clock findings of 9 Nov.

## Checks
Internal links and anchors of the two edited pages checked by script; the pages were not opened in a browser.
