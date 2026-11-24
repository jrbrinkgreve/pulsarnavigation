# 24 Nov 2026 — the panel is fixed for now; two observations from the reruns explained (no code change)

Jasper answered the open points of 23 Nov, and asked: "5) can you investigate". No code changed. The analysis below
reads the result files of 23 Nov; no new simulation. The docs change for points 1–3 is described at the end.

## Jasper's answers (his words)

1. The beam width: "is 3 m^2 in a square array implies 6 deg, thats okay as a working assumption for now."
2. The mount: "its fixed for now, but we can still show tracking as a what-if scenario."
3. Splitter or switch: "splitter is an okay working assumption."
4. His own validation runs: "sure, i can run some time."
5. The two unexplained results: "can you investigate."
6. The overwritten roadmap note: "tell me about that roadmap in 2 hours!" (a reminder is set in the session for 10:54).
7. The numbers not yet replaced in the docs: "what do we need to update the docs?"

What follows from 1–3: the fixed panel is the case, with a beam of about 6° and both beams live all the time. Its
figure of merit is the SNR per pass (11.5 for B0329+54, TOA error 459 µs; scratch runs). The runs in which the pulsar
is followed (0.23 h per TOA, ±9.4 km east in 30 days) are a what-if for a tracking mount. A fixed panel with one beam
on the pulsar gives no position; with both beams on B0329+54, 108 km along one axis after 30 days.

## A. J0437−4715: the error bars look 23–33 % small, but the core is right

The spread of (TOA − truth)/σ was 1.33 (30 d, 207 TOAs) and 1.23 (100 d, 709 TOAs) in the clean runs.

| File (good TOAs) | spread of z | robust spread (1.48 × MAD) | spread of \|z\| < 3 | beyond 3σ | beyond 4σ | beyond 5σ | largest |
|---|---|---|---|---|---|---|---|
| J0437, design point, clean, 100 d (709) | 1.233 | 1.013 | 1.060 | 2.54 % | 0.99 % | 0.56 % | 5.9 |
| J0437, design point, clean, 30 d (207) | 1.332 | 1.138 | 1.067 | 4.35 % | 1.45 % | 0.48 % | 5.7 |
| J0437, design point, rfiCalScint, 100 d (1009) | 1.133 | 1.020 | 1.027 | 1.59 % | 0.40 % | 0.20 % | 5.2 |
| J0437, old nominal, clean, 100 d (76) | 1.129 | 0.952 | 1.073 | 1.32 % | 0 | 0 | 3.0 |
| B0329+54, design point, clean, 100 d (8649) | 1.124 | 1.013 | 1.026 | 1.34 % | 0.36 % | 0.14 % | 13.2 |

- The core is honest: the robust spread is 1.01 and the spread inside 3σ is 1.06. The total spread is raised by a
  tail: 18 of 709 TOAs lie beyond 3σ.
- Those 18 errors are not spread out. They sit at ±0.0065 to ±0.013 turns, around ±0.009 turns (±52 µs), on both
  sides. The template's two sharp components are at −0.0008 and +0.0091 turns (FWHM 0.0100 and 0.0085, amplitudes
  0.46 and 0.31), 0.0099 turns (57 µs) apart. So these are **wrong-peak fits**: the noise lines the template up one
  sharp peak over, at a secondary maximum of the correlation. About 2.5 % of the TOAs at SNR 7.3.
- It is not new with the array: the SNR per TOA is the same as before (7.1–7.3). The old runs had 22 and 76 TOAs,
  where 2.5 % is one or two, so it could not be seen. The rerun has ten times the TOAs.
- With scintillation the tail is smaller (1.6 %), because many TOAs then come from bright states with a higher SNR.
- What it means: J0437−4715's stated error bar is right for 97.5 % of its TOAs, and the rest are off by about 57 µs
  (4–6σ). `navFit`'s mixture (2 % at 4σ) is close to this by accident. Options, none built: a higher target SNR for
  pulsars with several sharp peaks; a wrong-peak term in the TOA's likelihood (the literature's phase-likelihood
  shape, 21 Nov research); a component of ±57 µs in the fit's outlier model.

Seen on the way, in every file old and new: the spread of z rises with the TOA's reported SNR (lowest third 0.75–0.81,
highest third 1.28–1.36). My reading, not tested: the error bar goes as one over the fitted amplitude, so a TOA whose
amplitude came out high by noise gets too small an error bar. This is the known 5–10 % at SNR 7, sorted by SNR.

## B. B0329+54 with scintillation: 10–11 % of the TOAs have pChi2 < 0.05 because the channels differ in amplitude

| File | mean reduced χ² | its spread | pChi2 < 0.05 | slope of χ² against SNR² |
|---|---|---|---|---|
| B0329+54, design point, clean, 100 d | 0.9989 | 0.0169 | 4.3 % | (−0.6 ± 3.7)·10⁻⁵ |
| B0329+54, design point, rfiCal, 30 d | 0.9986 | 0.0166 | 4.3 % | (8.5 ± 6.8)·10⁻⁵ |
| B0329+54, design point, rfiCalScint, 100 d | 1.0058 | 0.0174 | 10.4 % | (12.7 ± 2.5)·10⁻⁵ |
| J0437−4715, design point, rfiCalScint, 100 d | 0.9987 | 0.0094 | 3.7 % | (−0.3 ± 1.8)·10⁻⁵ |
| B0329+54, old nominal, rfiCalScint, 100 d | 0.9993 | 0.0103 | 4.6 % | (30 ± 15)·10⁻⁵ |

- Cause: scintillation gives every channel its own pulse amplitude (B0329+54's scintles are 7 MHz wide, about two
  channels), while the 'optimal' fit assumes the same amplitude in every channel (`ChannelGain` 1). The misfit adds to
  χ² in proportion to SNR². It is the mechanism of the band-edge taper of 17 Nov, now in every channel.
- Four checks agree. (1) The mean χ² rises by 0.0070, which is 0.40 of its spread; a shift of 0.40σ turns 5 % into
  P(z > 1.245) = 10.7 %, the observed fraction. (2) χ² grows with SNR² only in the scintillation runs (5σ). (3) The
  slope, 1.3·10⁻⁴ per SNR², times the degrees of freedom (about 6600, from the spread of χ²) gives a variance of the
  channel amplitudes of about 0.8–1.0, the value of fully modulated scintillation. (4) J0437−4715, whose scintles
  (740 MHz) are wider than the band, shows nothing.
- Why now: a TOA takes 10–12 minutes at the design point, less than the scintillation time of 13 minutes, so the
  amplitudes do not average out. At the old design point a TOA took 2.5 hours and averaged about eleven of them.
  Within the new runs the short TOAs show more (11.9 %) than the long ones (9.0 %).
- What it means: the TOAs are not worse (the spread of z is 1.11 with and without scintillation), and no TOA is
  dropped, since the quality flag does not cut at 0.05. But pChi2 is the alarm for real data, and with scintillation it
  will ring twice as often as stated for this reason alone. Options, none built: expect it and widen the alarm for
  scintillating pulsars; or estimate the channel amplitudes over a scintillation time and give them to the fit
  (`'ChannelGain'`), which would also add sensitivity in principle (up to √2 for fully modulated scintles) if the
  amplitudes could be measured at this signal level, which is doubtful.

## C. The tail fractions for the docs (read from the same files)

The pages still quote October's 1.9 % beyond 3σ and 0.3–0.7 % beyond 4σ. At the design point, B0329+54: 1.34 % beyond
3σ, 0.36 % beyond 4σ, 0.14 % beyond 5σ, spread inside 3σ 1.026 (clean, 100 d, 8649 TOAs); with RFI and scintillation
1.41 %, 0.30 %, 0.13 %. The largest single error is 13σ (one TOA in 8649). `navFit`'s outlier fraction of 2 % at 4σ
stays a fair description for B0329+54.

## What the remaining docs numbers need (Jasper's point 7)

| Not replaced | What it needs |
|---|---|
| Tail fractions beyond 3σ and 4σ | nothing more: section C above; a docs edit |
| Locked-radar biases (excision tables, limits, results) | `runLockedRadar.m` has its levels typed in (−5 … −60 dB): a small code change to add −38.0 and −44.2 dB, then its run (voltage chain, DM 5, 0.1 s per level) |
| The RFI list's link budget (receiver page) and the RFI variants built on it | the antenna's gain towards the horizon and the place of the clean 150 MHz from the team, or a site survey; then `rfiRealistic` in `pipelineParams.m`, `runRFIProfile` (25–50 min per case) and the RFI reruns (about 1.5 h) |
| The single-pass numbers as repository results (they are scratch runs) | the drift through a fixed beam in the simulator (risk register R1): a code change with its own plan |
| The archived ship results | nothing, unless the ship comes back |
| The 6-day SDR tracking table (scripts page) | a 6-day run of the SDR case at the design point, about 2 min |
| `validated` badges | Jasper's own runs of `runScenario`, `tests/expNavFit.m`, `tests/expNavPlan.m` |

## The docs change of 24 Nov

One docs subagent, briefed from this log. The pages now lead with the fixed panel and name every result with the pulsar
followed a what-if for a tracking mount; the 15° alternative is gone; the tail fractions of section C replaced
October's; the two "not looked into" sentences became the explanations of sections A and B. 19 pages and the captions
of six figure generators; no figure regenerated. My checks: all 29 pages parse and their links resolve; I reworded one
sentence on the TOA page that the subagent flagged (it said a wrong peak is dropped by the χ² threshold, which finding A
contradicts) and wrote the changelog entry in `maintenance.html`. Also recorded: notes §11 and the risk register
(R14 answered for now; R1 becomes the first development).

Left by the subagent as questions: `Chi2Alarm` defaults to 10⁻³ while the scenario tables count pChi2 < 0.05; the
23 Nov history entry still mentions that a wider beam was run and dropped.

## The roadmap's dated note, restored (later on 24 Nov)

The docs update of 23 Nov had replaced the numbers inside the roadmap's "Physics notes from the revision", a note of
23 Oct. At 10:54 (the reminder Jasper asked for) I showed him the old and the new numbers and two options. Jasper:
"option 2, restore it as a dated record". Done: the three notes are back word for word from commit `4a28aa2`, under the
heading "Physics notes from the revision (23 Oct 2026)" with one line saying which receiver they were written for, and
one paragraph "Since then (23–24 Nov 2026)" after them with what changed (analogue beams; the planning study at the
design point: ±9.5 / ±5.2 km east, polar 28 km, 28 → 14 km with the three extra pulsars, 36 and 24 TOAs), labelled as
part of the what-if of a tracking mount. Checked in the preview: renders, the link resolves, no console errors.
The changelog entry of 24 Nov in `maintenance.html` has a sentence on it.

For a next docs job: tell a subagent explicitly that a dated note keeps its numbers and gets the new ones beside it.

## Unit tests and single-pass checks (later on 24 Nov; Claude's runs)

Asked "is there any other compute you /we can do?", I listed four runs; Jasper: "Tests + pass checks".

**`runAllTests`: 33 of 33 passed** (15.3 min, with the other job running beside it). It is the first full run since
`runScenario.m`, `tests/expNavFit.m` and `tests/expNavPlan.m` changed on 23 Nov. The tests do not call those three, so
this shows that nothing else moved; it is not a test of the new design point. Not Jasper's run.

**Single-pass checks, fixed panel** (scratch copy `runScenario3m2.m`; B0329+54, a 6.3° beam, design point, clean; 240
passes per case, new seeds; a gain g on the pulsar is run as SEFD / g; "detected" is the pipeline's good-TOA flag):

| Case | Gain | SNR per pass, median (10–90 %; lowest) | Detected | SNR ≥ 7 | TOA error | spread of error/σ |
|---|---|---|---|---|---|---|
| beam centre | 1 | 11.40 (10.02–12.78; 8.44) | 100 % | 100 % | 462 µs | 1.06 |
| second beam, 20° off the normal (cos 20°) | 0.94 | 10.72 (9.34–12.10; 7.76) | 100 % | 100 % | 492 µs | 1.07 |
| off centre: zenith-facing, no tilt (2.6° off) | 0.62 | 7.17 (5.81–8.54; 4.17) | 99.2 % | 57 % | 735 µs | 1.15 |
| low-flux day (5th percentile, literature) | 0.45 | 5.23 (3.88–6.63; 2.41) | 74.6 % | 7 % | 927 µs | 1.70 |
| off centre on a low-flux day | 0.28 | 3.54 (2.47–4.79; 1.85) | 17.9 % | 0 % | 1092 µs | 2.80 |

The medians follow the gain, as the scaling of 23 Nov said (7.1 and about 5 predicted). New: without the 2.6° tilt only
57 % of the passes reach SNR 7; on a low-flux day a quarter of the passes is lost and the errors are 1.7 times the
stated ones; the second beam is nearly as good as the first, so the two-beam position of 108 km becomes about 112 km
(scaled from the two TOA errors, not simulated). Docs: scenario R7.4 (the table), limits (one sentence), the changelog.

Not run, offered: both beams over 100 days with an SDR clock; the 6-day SDR tracking table.

Files: `.claude/runs/2026-11-24_passchecks/` (`passchecks.m`, `passchecks.log`, `runAllTests.log`).

## Files

`.claude/runs/2026-11-24_investigate/investigate.m` and `investigate.log` (gitignored, this Mac only).
