# 21 Nov 2026 — deep research on the whole project (literature only; no code or docs change)

Jasper: "can you do a general deep research on the project? there's a lot of usage left, you can continue until
about 80 fable or 80% overall%." No repository code, docs, notes or earlier logs were changed. Nothing was run in
MATLAB. Nothing is committed.

## The session

- Read first: notes §1 and §11, the risk register (2 Nov), the static-nav handover (13 Nov), the budget recheck
  (17 Nov), and the logs of 13–20 Nov.
- Round 1: fourteen read-only literature researchers in parallel, one per part of the project. This is a one-off
  exception to the single-session rule, taken because "deep research" was asked for. They shared one web-search
  rate limit and exhausted it, so every first-round notes file ends its sections with a Gaps list.
- Round 2: six researchers on the gaps that decide the most. Then one writer for the report.
- My own part: the project brief for the researchers, three checks against the code, one combined calculation, and a
  list of first-round statements that the second round corrected (`coordinator_checks.md`).

## Where the files are (all untracked)

- Report: `reports/Pulsar navigation project research.md` (about 25,000 words, 223 source links, eight sections; section 7
  lists the options per roadmap step, section 8 what could not be verified). I read it in full against the researchers'
  reports and changed four places: the heading on the "30 m²" paper (it said the number was "the project's own, not a
  paper's"; the finding is "not found as worded"), the researcher count, and two statements of "a factor of six" that
  now give the actual factors (×6.7 for the clock line; ×3.8 east and ×6.1 north for the zenith panel).
- Notes, one file per topic, about 190,000 words: `research_notes/Pulsar navigation project research/`. The six
  files named `gap_*.md` are the second round and win where they disagree with a first-round file.
  `coordinator_checks.md` lists the corrections.
- Front-end and hardware findings that refer to the Task 2.4 deliverable: `logs/2026-11-21_task2.4-research.md`
  (ignored by git, like the other Task 2.4 logs).
- **The repository is public.** The notes and the report quote the Task 2.4 writeup in places. Jasper, asked whether
  to add `research_notes/` and `reports/` to `.gitignore`: "Yes, ignore them". Done the same day (two lines in the
  root `.gitignore`, not committed). Both folders are therefore on this Mac only.

## What kind of numbers these are

Every number below is from a paper, a data sheet, a catalogue, or a researcher's own estimate, as marked. None is a
project result, and none was simulated in the pipeline. The fetch tool passes pages through a summarising model, so
quoted values are paraphrase-level unless a researcher cross-checked them.

## Findings that touch project numbers

| # | Topic | The project says | The research says | Kind |
|---|---|---|---|---|
| 1 | Array area against direction | 0.61 m² for every pulsar, 24 h a day (`scenario.html`, notes §11) | Near broadside a 4×4 lattice of 7–9 dBi patches realises 0.38–0.49 m² (band mean, 90 % efficiency); at 44° from the axis 0.14–0.20 m². Day-averaged SEFD for a zenith panel at 52° N: 2.7× nominal for B0329+54, 5–8× for the pulsars at declination 8–22°. A pointed panel: 1.5–2.4×. | researcher's pattern model, ±1 to ±1.5 dB; checked against the MWA tile's quoted area to 5–7 % |
| 2 | Position with that array | 30 / 47 km east / north after 30 d, nominal | A one-day fit weighted by the zenith panel's area: east ×3.8, north ×6.1, so roughly 110 / 290 km. A panel that keeps the pulsar near its axis: ×1.45. | my own scale estimate (white noise, constant T_sys); `coordinator_checks.md` §3 |
| 3 | Hour-angle coverage | not stated | B0329+54, 30 d: full day 28.6 / 44.8 km; ±6 h 42 / 113 km; ±3 h 118 / 558 km; ±1 h no fix. Zhou et al. 2026: with a constant zenith angle the position signal is essentially zero. | researcher's white-noise model (reproduces the 4 Nov table within 5 %); Zhou read in full |
| 4 | SDR reference | "±2 ppm TCXO class" (`makeClock.m:40`, `receiver.html` 4a, FM-10) | The RFSoC 4x2 has no TCXO: a 48 MHz crystal, ±15 ppm, ±25 ppm over temperature. 54–144 ms drift per hour, not 7.2; FM-10's rate search about 500 trials, not 40. All three candidate boards take an external 10 MHz. | reference manual rev. A5 (first-hand); ppm from distributor listings |
| 5 | FM-11, oscillator slopes | TCXO 10⁻⁸–10⁻⁷ per K; OCXO "about 10⁻⁹ per K" | TCXO confirmed. OCXO 2·10⁻¹¹ to 1.1·10⁻¹⁰ per K, so the register is conservative by 10–50×. Daily term on B0329+54 for ±1 K: TCXO 70–703 km, single-oven OCXO 0.14–0.8 km. | data sheets; slopes are range totals divided by the range, no room-temperature slope is published |
| 6 | FM-11, cost of a free solar-day line | "99 % of the term looks like position" | The cost is 1/√(1−ρ²) on the position error: ×6.7 (one pulsar, 30 d), ×2.1 (100 d). Two equal pulsars 6 h apart in right ascension ×1.40; 12 h apart ×1.00. A faint second pulsar helps little. | researcher's derivation and numerical model |
| 7 | Clock rate in `navFit` | one common `rate` column (`navFit.m:158`) | Stale spin frequencies (10–150 µs/day) and old position epochs differ per pulsar. A phase, frequency and derivative per pulsar costs the position under 0.1 % beyond three days. | three researchers independently; cost from a numerical model |
| 8 | Time-scale units | not handled | ATNF lists B0329+54 in TDB and B0950+08, B1133+16, B1642−03, B1749−28, B1929+10, B1937+21 in TCB. Unconverted: 1.340 ms/day between the groups. Not a live bug: truth and observer share the values today. It becomes one with real data or a PINT check. | catalogue; my grep of the code |
| 9 | UT1 | `delayModel` default UT1 = UTC; FM-31 "up to 0.9 s" | CGPM Resolution 5 (Oct 2026): continuous UTC from 20 May 2027, \|UT1 − UTC\| up to 3600 s. UT1 must come from IERS Bulletin A (good to 1.4 ms at 10 d, 4 ms at 40 d). UT1 − UTC is −0.129 s now: 37 m at 52° N. | BIPM text read as a close paraphrase; Bulletin A |
| 10 | FM-5, position epoch | 47 µs per year of annual term for B0329+54 | 37.9 µs per year for its geometry. Proper motion left out today: about 557 µs (167 km) annual amplitude; B1133+16 about 8 ms, B1929+10 about 2.5 ms. | researcher's calculation from ATNF v2.7.0 |
| 11 | FM-4 | "B0329+54 has glitched" | No entry in the Jodrell Bank or ATNF glitch tables. Timing noise under 0.2 µs over 30 d after a ν, ν̇ fit. | glitch tables; Hobbs et al. 2010 scaled |
| 12 | FM-3, mode changes | abnormal 15 %, "spells of hours", offset unknown | Occupancy 13–18 % (measured). Abnormal spells average 26 min, longest seen 2.5 h. Offset 0.15 ± 0.04 ms from a 1400 MHz component table (0.8 ms as the pessimistic bound). The abnormal mode is about 1.7× brighter. | papers; the offset is a researcher's template fit, table signs to be checked against the PDF |
| 13 | Scintillation, FM-7 / S2-5 | diffractive only; B0329+54 7 MHz × 13 min | 7 MHz confirmed. Refractive modulation index 0.31–0.47 at 1540 MHz is missing. 5th-percentile flux about 0.45 × mean for 2 h, 0.75 × for 30 d. | measured index; percentiles from a researcher's Monte Carlo |
| 14 | Pulsar table (notes §11) | J1713+0747 W50 0.30 ms; B1929+10 29 mJy; B1642−03 26 mJy; J2145−0750 5.5 mJy; B1937+21 W50 ~0.05 ms; "uncertainty of t ~×2" | 0.112 ms (measured); 41 ± 6 up to ~90 mJy; 11–55 mJy; 8–10 mJy; 37–39 µs; a factor 2 in flux is a factor 4 in time | papers and EPN profiles; ATNF v2.7.0 has the same S1400 as v2.6.5 |
| 15 | RFI scenario (`pipelineParams.m:75–95`) | 16 sources; radars with 2 µs pulses; mobile downlink 15/10/15 MHz; DME "just below 1.2 GHz" | Solid-state radars send 17.5–410 µs chirps at 4–10 % duty; `detectRFI`'s longest window is 3.8 µs. Licensed downlink layout is 15/15/10 MHz. DME replies reach 1213 MHz (Pampus on 1212 MHz). Missing: GLONASS L1 and L3, Galileo E5b, wind profiler (Cabauw, 1290 MHz), radar satellites, the L1 aggregate's real width. | ITU-R M.1463-3, Staatscourant 2020 41318, interface documents; levels are estimates |
| 16 | FM-15, polarisation | "the size is unknown" | For a 10 % differential gain: B0329+54 11–13 µs band-averaged, B0950+08 23–28 µs, B1929+10 31–50 µs, J0437 1.65 µs, B1937+21 0.4 µs. | computed from published Stokes profiles; MSP values within 2–2.6× of van Straten 2006 |
| 17 | FM-2, the Sun | 4·10⁵ Jy at minimum, "perhaps twice" at maximum; sidelobe level unknown | 45–55 sfu at minimum, 110–150 sfu at maximum (yearly medians, 1415 MHz). First sidelobe of a uniform 4×4: −11.3 dB. Main-beam SNR factor 0.17–0.22 at maximum. | NOAA daily values 1979–2026 |
| 18 | Doppler-shifted dispersion | "±7.1 µs" for B0329+54 | ±9.3 µs at 1.4 GHz; 7.1 µs is the 1.6 GHz value | researcher's arithmetic |
| 19 | Polar component | 19 Nov idea note: "or pulsars at different declinations" as a separate route | Declination spread alone does nothing while each pulsar has its own free constant. The 4 Nov redesign log, item 3, words it correctly. | Thompson, Moran & Swenson ch. 12; a 1980 VLBI analysis |
| 20 | Stream counts | FM-9, FM-14, FM-38 and the 5 Nov log say 18 | 16 elements since 20 Nov (32 streams if dual-polarised) | bookkeeping |

## Confirmed by the literature

- `delayModel`'s terms agree in form and sign with TEMPO2 and PINT. Expected difference against PINT today: 0.1–0.2 µs
  daily, mostly UT1 (estimate).
- `estimateTOA` 'equal' is Taylor's FFTFIT; `detectPulsar` is the optimal known-template test; the √2·σ_p/SNR rule is
  Cordes and Shannon's bound; the box rule being optimistic by 2.2 in time is consistent with the literature.
- "5–10 % small error bars, 1.9 % beyond 3σ at SNR 7" is consistent with the literature's breakdown at S ≈ 8–10.
- T_sys 60–80 K is realistic near broadside (measured arrays: 45–70 K at half-wavelength pitch).
- The daily-sine position sensitivity: Zhou et al. 2026 give cos δ / (c√2) per equatorial component; applied to
  B0329+54 with the project's TOAs that is 27.8 km, against about 30 km east.
- `siteSearch` and the no-threshold fit have published relatives: SEPO (Zheng et al. 2019), Han et al. 2020 (a grid of
  trial sites with χ² and an elevation model, Parkes), and gamma-ray timing without TOAs (Ajello et al. 2022).
- The sky is 3.2–6.3 K in a 15° beam for all pulsars except B1749−28 (about 11 K).
- Zero-fill plus zero weight for lost samples (FM-9) is what DSPSR and other backends do.
- B0329+54 is the only northern pulsar with TOAs within a pass; the next ones need 13–87 h on source.

## Not found, or not verified

- The paper behind "30 m², 1.2–1.6 GHz, J0437−4715, 0.3 km after 1 h" (17 Nov log). The closest is Kaur, Hobbs et al.
  2025 (arXiv:2506.08056), Table 1: 0.27 km in 1 h, but for a 256-element array at 700 MHz with 64 MHz. Its method at
  30 m² and L band would give about 0.08 km. The 17 Nov log's "the simulator reproduces the paper's 0.3 km" holds only
  if the paper states 30 m² at L band. **Needs the title or authors from Jasper.**
- No ground position fix from radio pulsars with a small aperture exists in what was read. No pulsar detection with a
  patch array or any array of about 1 m² or less was found; the smallest confirmed at 1.4 GHz use 3 m dishes.
- No published treatment of a daily oscillator term aliasing into position, or of a search over site and clock rate
  together.
- No measured board-to-board phase repeatability after a restart for several RFSoC 4x2 boards.
- No measured efficiency or T_sys against scan angle for an L-band array near one wavelength pitch.
- No multi-epoch 1.4 GHz flux for B0329+54 newer than the Urumqi series; no recent public ephemeris for B0329+54
  (catalogue epoch 1986).
- LNA and filter data-sheet values at 1.2–1.6 GHz (vendor sites blocked automated reads).
- Sala et al. 2004, Buist et al. 2011 and Chen 2018 were not read; two summaries of Sala disagree.
- ESA NAVISP-EL1-079 "CAPPA" (TU Eindhoven with TU Delft, from 21 April 2026, 12 months): nothing published beyond
  the project page.

## Open items (none is decided; each needs Jasper's OK)

1. ~~`research_notes/` and `reports/`: keep them out of the public repository.~~ Settled: ignored by git (above).
   This log itself is not ignored and quotes no deliverable text; the Task 2.4 detail is in the ignored log.
2. The title of the 17 Nov paper.
3. The array's orientation and usable scan range decide the position budget more than T_sys does (rows 1–3). The fast
   simulator has no A_eff(direction, frequency); the scenario runs assume full gain for 24 h. Candidates for step E4
   or earlier: a gain-against-hour-angle factor in `runScenario`, and the choice between a zenith panel, a pointed
   panel and 4:1 analogue sub-arrays.
4. `navFit`: a phase and rate per pulsar (row 7); the solar-day line and its harmonics as optional columns (row 6).
5. E1/E2 inputs: proper motion and position epoch, TCB → TDB, UT1 from Bulletin A (rows 8–10); the PINT check needs
   a Linux container or Rosetta on this Mac (native ARM Python lacks the extended precision).
6. ~~Risk-register updates.~~ Done the same day at Jasper's "update the risk register with the findings"
   (`2026-11-02_risk-register.md`, section "21 Nov 2026: the register against the literature"): 31 rows annotated with
   a paragraph marked 21 Nov, six ratings changed (FM-4 12 → 15, FM-18 12 → 15, FM-7 6 → 10, FM-12 6 → 8, FM-3 10 → 5,
   FM-32 8 → 4), four scenarios added (FM-43 the array's area against direction, 20; FM-44 mixed TCB/TDB units, 16;
   FM-45 the reference not locked, 12; FM-46 ADC calibration and the 5G alias, 9), the table re-sorted, and eleven
   actions R1–R11. The ratings stay a draft. The actions are not yet in notes §10 or `roadmap.html`.
7. The notes §11 pulsar table (row 14) and the refractive term in the scintillation runs (row 13).
8. The RFI scenario (row 15) and a site survey early (FM-18); long radar pulses against `detectRFI`'s windows.
9. Block F before hardware: public baseband that downloads today without a login (B1133+16 from Green Bank at
   1125–1875 MHz, 17 GB per 22 s file; B1642−03 from Dwingeloo at 438 MHz, 12 GB, with the observatory's own DSPSR
   fold as a reference). Downloading needs Jasper's go-ahead.
10. Docs: the unknown-phase threshold 4.15 is the single Gaussian's; the five-component template's is about 4.5 (the
    code computes it from the template, so only the docs' number is template-specific).

## Process notes

- Two researchers wrote helper scripts in the session scratchpad (outside the repository); their models are written
  out in the notes.
- The fetch tool saved copies of PDFs it could not read into its own cache outside the repository. In round 1 some
  researchers extracted text from public PDFs (data sheets, ITU recommendations, GNSS interface documents). In round 2
  the instruction was not to download; one researcher still read three PDFs from that cache, and once opened the
  Browser pane on a vendor site, where a data-sheet link raised a save dialog. It closed the pane and did not retry.
- Usage: 14 % overall and 25 % Fable at the start; 30 % and 57 % at the end. The 80 % allowance was not used up: the
  second round closed the gaps that could be closed without reading PDFs by hand.
