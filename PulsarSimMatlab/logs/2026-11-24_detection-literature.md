# 24 Nov 2026 — literature research on pulsar detection and how it performs (no code or docs change)

Jasper: "for more deep research, can you look for papers on pulsar detection and how they perform? we are in quite a
weird situation, where we do academic research, but are also building a real thing at the same time". No code, docs,
notes or earlier logs were changed. Nothing was run in MATLAB. Nothing is committed.

## The session

- Brief for the researchers (`00_brief.md`): the request, the hardware as of 23–24 Nov, what the research of 21 Nov
  already covers, and the rules (read-only, at most 15 searches each, no downloads, no Browser pane, no project
  details in queries).
- Round 1: six read-only researchers in parallel, one per angle: detection theory and algorithms; predicted against
  achieved sensitivity; transit instruments; first light on new instruments; verification and test practice;
  research while building an instrument.
- Round 2: two verifiers re-opened the primary sources for 28 claims the report rests on. Reason: every source is read
  through a fetch tool that returns a model-written summary. No number was found invented; the errors were in what a
  number was taken to mean.
- My own part: the brief, `coordinator_checks.md` (arithmetic and geometry checks, one conflict between notes
  settled), the brief for the report writer with the corrections listed, and a read of the whole report with its
  arithmetic rechecked (thresholds in SEFD, the derating compositions, the detection probabilities, the single-pulse
  figure, the Sun and Cas A levels). No change to the report was needed.

## Where the files are (git-ignored, this Mac only)

- Report: `reports/Pulsar detection performance literature.md` (about 12,000 words, eight sections). Section 6 is the
  project side, section 7 publishing, section 8 what could not be verified.
- Notes: `research_notes/Pulsar detection performance literature/` (six topic notes, two `gap_verification_*.md`,
  `coordinator_checks.md`). The verification files and the coordinator's checks win over the first-round notes.

## Findings that bear on the project

Labels: published = in a source and verified in round 2 unless marked; estimate = arithmetic in the notes or mine;
simulated = the project's runs (Claude's runs, not yet Jasper's).

| # | Finding | Kind |
|---|---|---|
| 1 | The known-phase test in `detectPulsar` is the most powerful test for a known profile in white noise, with detection probability Φ(S − η) (Morello et al. 2020). Other statistics need more signal on the B0329+54 profile: best boxcar 5.3, H-test 6.8, against 4.5 for the matched filter with the phase searched (50 % detection, false-alarm 10⁻³). | published; the per-statistic numbers are estimates |
| 2 | The simulated S/N 11.5 per pass equals the radiometer equation with the 10 % width W10 (11.3). With W50 the same pass reads 25. A published S/N in the W50 or single-pulse convention is about 2.2 times higher for the same data. | estimate, checked by me |
| 3 | Dishes meet their SEFD within about 20 %. New arrays with a measured figure came in at 1.2 to 3 times the ideal SEFD (Apertif 85 Jy against 70 Jy on calibrators and 130–200 Jy from pulsars; BURSTT 11–12 kJy against 6 kJy; CHIME Pathfinder a few times the expectation). | published |
| 4 | Composition of the loss terms that apply to a template fold: 0.72 (0.37–0.98) of the nominal sensitivity if the 3 m² and 70 K were not measured for the whole beamformed panel on the sky, 0.90 (0.74–0.98) if they were. In S/N per pass: 8.3 (4.3–11.3) or 10.4 (8.6–11.3). Before red noise and before about ±50 % from the pulsar per pass. A range, not a prediction. | estimate |
| 5 | In SEFD: S/N 7 per pass needs at most 1.05·10⁵ Jy (0.026 m²/K); 90 % detection with the phase searched needs at most 1.3·10⁵ Jy; the 50 % point 1.65·10⁵ Jy (0.017 m²/K). Beam centre, mean flux. | estimate from the simulated 11.5 |
| 6 | Red noise: the published loss factors (PALFA 1.1–2 at 0.1–2 s, CHIME stacking an order of magnitude at 1 s and longer) are for blind Fourier searches. No source gives a measured loss for a fold with a known ephemeris. The simulator has no slow gain variations, so the size for this receiver is unknown until an off-source power spectrum is recorded. | published; the gap is the finding |
| 7 | Daily-transit timing is routine (Molonglo, CHIME, Pushchino, Nançay), at 30 to 1400 times the panel's sensitivity per transit. UTMOST-NS, a fixed array with no moving parts, accepts a TOA above S/N 7 and gets one from 53 % of its observations. No transit instrument with a stated SEFD above 10⁴ Jy reports a TOA. | published; the scaling is an estimate |
| 8 | No transit instrument read weights sub-integrations by the beam when forming a TOA. An unweighted fold between the half-power points loses 6.6 % against the weighted full pass; a pass is worth 0.753 of the half-power crossing time, which reproduces the 0.55 h. | estimate |
| 9 | √N for adding many passes is not established in the stacking literature (measured on quiet sky up to about 49 days; the Pushchino gains are planning values). Adding folded profiles in phase gives √N by construction in white noise, but needs the clock to hold the phase to about 1 ms over the span. | published; second part reasoning |
| 10 | Confirmation in the literature rests on null tests, repetition and independent software. For B0329+54 at 1300–1450 MHz the DM test is weak (the sweep is 12.9 ms; a fold at DM 0 keeps about half the S/N). Strong: sidereal recurrence, off-source hours (about 43 slices a day), the empty second beam, and the period against the Doppler prediction (±82 ppm over the year; one pass measures about 0.6 ppm). | published practice; strengths are estimates |
| 11 | A 20 ppm crystal moves the pulse 40 ms across one pass. A constant rate offset is a period offset, so a period search recovers the detection at a higher threshold, but the absolute period test is then lost. The RFSoC 4x2 paper (Effelsberg, A&A 2025) does not state its clock reference. | estimate; published gap |
| 12 | Time-stamp faults were found by timing, never by a detection: Effelsberg 97 ms and 409 ms, uGMRT 0.67 or 1.34 s, Apertif over 1 ms after power cycles, UTMOST-NS 4–18 µs. | published |
| 13 | For a hardwired analogue sum, path errors cannot be repaired afterwards. EMBRACE found up to 15 cm, the MWA's Engineering Development Array up to 20 cm. At 1.4 GHz an rms of 1 cm costs about 4 % and 2 cm about 16 % of S/N (exp(−σ²/2)). BURSTT found its array 3.0° off in azimuth after deployment; 2.6° costs the panel 11.5 → 7.2 (simulated). | published; the percentages are estimates |
| 14 | An on-sky SEFD on a continuum source needs the panel tilted: at 52° N the Sun stays 28.6° from the zenith, Cas A passes 6.8° north (beam response 0.04), Cyg A 11.3° south. | estimate, checked by me |
| 15 | The both-beams result of 23 Nov assumes the panel turned about its normal so that the second beam lies on B0329+54's track: panel diagonal about 69° east or west of north, panel horizontal. With the diagonal at azimuth 45° the pulsar misses the second beam by about 8°. Already in the log of 23 Nov; now an explicit mounting requirement. | estimate, checked by me |
| 16 | Publishing: the recurring sequence is a forecast paper before the hardware works, then a system paper built around a verification section (BURSTT 2022 / 2025 / Nov 2026; SEXTANT 2015 / 2016 / 2018). The Journal of Navigation judges simulation-only papers less favourably. No published radio work reports TOAs from a purpose-built aperture of square metres against their predicted error bar. | published; the gap from a limited search |

## What I told Jasper during the run and then corrected

- Red noise: I passed on "real sensitivity 1.1–2 times worse at 0.1–2 s, and B0329+54 sits in that range". The
  verification showed these are blind-search figures at DM above 150; for a known-ephemeris fold there is no number.
- "No transit instrument above 10⁴ Jy has published a TOA": holds only for instruments that state an SEFD; the
  Woodchester 2.4 m dish times Vela daily and states none.

## Candidate items, none decided and none implemented

For Jasper to pick from; each is in section 6 of the report.

1. Ask the antenna team whether 3 m² and 70 K were measured for the whole beamformed panel on the sky, and whether the
   per-tile path lengths were measured on the bench. This single answer moves the range in row 4.
2. Ask for the planned mounting: tilt (2.6° for B0329+54 at the first beam's centre) and rotation (row 15), and whether
   the panel can be tilted once to a continuum source for an SEFD measurement.
3. Clock for the first detection: external reference, or a period search with a trial-adjusted threshold.
4. Detector side, when hardware data exist: an empirical null from off-source and wrong-ephemeris folds, injection of
   the simulator's pulsar into recorded noise, a count of the trials, the minimum flux at a named detection probability.
5. Simulator side: slow gain variations (red noise) and the refractive flux modulation of B0329+54 are not modelled.
6. A record of measurements from now on for a later instrument paper (table in section 6 of the report).
7. Risk register: rows for on-sky SEFD shortfall of an analogue array, red noise in the fold, panel orientation, and
   the weak DM test are not added yet.

## Later on 24 Nov: the risk register updated

Jasper, asked whether to add the findings to the register: "Yes, update it". Item 7 above is done.
`logs/2026-11-02_risk-register.md` has a new section "24 Nov 2026: the register against the detection literature":

- Six rows annotated with a paragraph marked **24 Nov**: FM-1, FM-9, FM-10, FM-14, FM-42, FM-48.
- One rating changed: FM-14 (the analogue beams deliver less than their budget), 3 × 3 = 9 → 4 × 3 = 12. A draft for
  Jasper's review, like the other ratings.
- Two new scenarios: FM-49, the panel mounted a few degrees off its tilt or rotation (9); FM-50, a first detection that
  is not one, or a real one that cannot be shown to be one (9). The register now has 50 rows.
- Actions: R12 changed, R16–R21 new (the question to the team, the off-source spectrum, the confirmation set, the clock
  for first light, the survey of the mounting, the S/N convention). None is decided or built; not in notes §10 or the
  roadmap.
- A correction to R12 of 23 Nov: it asked for drift scans of the Sun through both beams, but at 52° N the Sun never
  comes within 28.6° of the zenith and crosses neither fixed beam. The measurement needs the panel tilted.
- Two estimates made for the register, mine: a period search over ±25 ppm is about 60 trials (steps of 0.8 ppm, a
  quarter of the pulse width of drift per pass) and raises the unknown-phase threshold from 4.5 to about 5.3; the
  beam gain off centre is 0.93 at 1°, 0.76 at 2°, 0.62 at 2.6°, 0.53 at 3° (Gaussian, 6.3° wide).
- Checked after the edit: 50 rows, ten cells each, P × I equal to the score in every row, highest score first.
- Docs: not changed. The pages name the register as the record of 2 Nov (39 scenarios) and were not changed for the
  updates of 21 and 23 Nov either.

## Not verified

The sources were read through the fetch tool's summaries in both rounds; nothing was read from a PDF. Six sources are
dated after mid-2026 and rest on the researchers' readings only. EMBRACE publishes no system temperature or effective
area in the texts read. Section 8 of the report has the full list.
