# 23 Nov 2026 — the array is about 3 m² with two hardwired analogue beams (no code or docs change)

Jasper, after talking with the antenna engineering team, then "record it". Recorded here and in notes §11
("Array, update 23 Nov 2026"). No code, docs or risk-register change. Nothing was run in MATLAB: every budget number
below is a scaling of the simulated 1.93 h per SNR 7 TOA of B0329+54 at 3.2·10⁵ Jy (block C), not a new simulation.

## What Jasper said

- "we have about a 3m^2 effective area array!"
- "we get a 4x4 'patches' bowtie array, but beamforming is done differently: each patch contains 4x4 mini-antennas,
  dual pol. these are sum-and-delay analogly beamformed to 0,0, and theta,phi = 20,45 deg using a power splitter. These
  dual beams x dual pol output 4 data streams per patch, and between the 16 patches the same thing again happens So we
  get 4 analog outputs, already hardwired to analog beamform to 0,0 and 20,45 degrees."
- "the engineering team prefers to talk in SEFD en sensitivity compared to SNR, so update me on the terminology".
- The mount: "that is still open to decide, but ideally, they want to detect a pulsar within a single pass through the
  beam".
- LNA placement, and whether the 3 m² includes losses: not known, "but measurements say there is about 70K of system
  temperature at the total output".
- The second beam: "idk what they made the second beam for and why its fixed". "phi is in the xy plane, so its
  'diagonally' across the array patches /elements."
- The delays: "yes fixed for good".
- The band: still 1.2–1.6 GHz, "but there's quite some RFI, so we get about 150 MHz of 'free area' with Tsys actually
  at 70K".
- The two beams: "The beam switch is called 'simultaneous', and implemented with a power switch per element that
  'doesnt impact SNR', but im afraid how it exactly impacts the integration time." Later: "I very briefly saw a diagram
  with 'SPDT' in it."

## The hardware as I read it

- 16 tiles in a 4×4 grid; each tile has 4×4 dual-polarized bowtie elements: 256 elements, 512 element outputs.
- Analogue delay-and-sum in two stages, inside each tile and then across the 16 tiles, each stage hardwired to the same
  two directions: (θ, φ) = (0°, 0°), the panel's normal, and (20°, 45°), tilted 20° along the diagonal of the grid.
- Outputs of the whole array: 2 beams × 2 polarizations = 4 analogue signals. One RFSoC 4x2 has four ADC inputs.
- No digital beamforming and no electronic steering. Pointing anywhere else needs the panel to move.

My arithmetic from 3 m² and 256 elements (not from the team): cell 0.0117 m², pitch 0.108 m (0.43–0.58 λ over
1.2–1.6 GHz, half a wavelength at 1.4 GHz), tile side 0.43 m, array side 1.73 m. Final beam, 0.886 λ/D: 7.3° / 6.3° /
5.5° at 1.2 / 1.4 / 1.6 GHz. Tile beam: about 22–29°. Directivity for 3 m² at 1.4 GHz: 29.2 dBi. At this pitch there is
no grating lobe for these two beams, so the one-wavelength-pitch problem of FM-43 does not apply.

## Terminology: SEFD and sensitivity (the team's language)

These describe the instrument alone. SNR also needs a pulsar, a bandwidth and a time.

| Quantity | Definition | 0.61 m², 70 K (the notes until now) | 3 m², 70 K |
|---|---|---|---|
| Sensitivity | A_eff / T_sys | 0.0087 m²/K | 0.043 m²/K |
| SEFD | 2 k T_sys / A_eff = 2760 Jy / (A_eff / T_sys in m²/K) | 3.2·10⁵ Jy | 6.4·10⁴ Jy |
| G/T at 1.4 GHz | 10 log(4π A_eff / λ² / T_sys) | 3.8 dB/K | 10.7 dB/K |

- SEFD is the flux density of a source that doubles the system noise; lower is better.
- From SEFD to SNR: SNR = (S / SEFD) √(n_pol B t) × a pulse-shape factor. The simulator's ρ is S_peak / SEFD.
- Time per TOA goes as SEFD² / B; a position error as SEFD / √B.
- From now on: receiver budgets in SEFD and A_eff / T_sys first, SNR derived for a named pulsar, bandwidth and time.

## A single pass through a fixed beam (B0329+54)

The relations, so each number can be checked:

- SNR² per hour at the beam peak = (49 / 1.93 h) × (3.17·10⁵ Jy / SEFD)² × (B / 400 MHz).
- The pulsar drifts at 15°/h × cos δ = 8.7°/h. With a Gaussian beam of FWHM 6.3° the pass lasts 0.72 h between
  half-power points. SNR² adds with the beam gain squared, so the pass is worth 0.753 × 0.72 h = 0.55 h at the peak.

| Case (SEFD 6.4·10⁴ Jy, beam centre unless said) | SNR per pass | TOA per pass |
|---|---|---|
| 400 MHz, average flux | 18 | 285 µs |
| 150 MHz, average flux | 11 | 460 µs |
| 150 MHz, zenith-facing panel without a tilt (B0329+54 culminates 2.6° off the beam centre) | 7 | 730 µs |
| 150 MHz, low-flux day (5th percentile, 0.45 of the mean; research of 21 Nov) | 5 | — |

SEFD needed for SNR 7 in one pass with 150 MHz at beam centre: ≤ 1.0·10⁵ Jy (≥ 0.027 m²/K; T_sys ≤ about 110 K with
3 m²). For a bare detection at unknown phase (4.5): ≤ 1.6·10⁵ Jy. Real small systems obtain 0.5–1.0 of the ideal SNR
(research of 21 Nov), so a margin of 1.6 at average flux is thin. If the other 250 MHz is usable at a higher T_sys, the
channel weighting recovers part of the loss; not estimated.

Other northern pulsars: the next-brightest needed 13–87 h on source at 3.2·10⁵ Jy and 400 MHz (research of 21 Nov),
which is 1.4–9.6 h at 6.4·10⁴ Jy and 150 MHz, against 0.55 h per pass. None of them detects in one pass; their passes
have to be added.

## The position with fixed beams

A pulsar seen every day at one hour angle gives one TOA per day and no position: the site's delay at that hour angle is
a constant, which merges with the pulsar's phase constant (notes §11; Zhou et al. 2026 in the research report). The
30 / 47 km budget follows the pulsar through the day. Of that pair only the ±30 km per equatorial component is
simulated (block D, `expNavFit`); the 47 km north with the height known is the covariance study of 4 Nov, not a
simulation (`2026-11-13_handover-static-nav.md`). On 21–23 Nov I called the pair "the simulated budget" in chat; that
was too loose for the north component.

| Configuration, 150 MHz, 30 days, B0329+54 | East / north, 1σ (scale estimate) | How obtained |
|---|---|---|
| Fixed panel, one beam on the pulsar | no position | — |
| Fixed panel, both beams on the pulsar | about 120 km, east–west only | below |
| Panel on a tracking mount | about 10 / 16 km | 30 / 47 km × (6.4·10⁴ / 3.2·10⁵) × √(400 / 150) |

Both beams on one pulsar: B0329+54 is 20° from the zenith at an hour angle of ±33.5° (2.2 h before and after
culmination), so the panel can be turned about its normal to put the second beam on its track. The difference of the
two daily TOAs measures −0.17 a + 0.55 b of the daily sine (a radial, b east): 1.07 µs per km east. With 460 µs and
about 510 µs per pass the difference is 690 µs per day, 640 km per day, 120 km after 30 days. It needs the delay
difference between the two beam paths, a constant of the hardware, measured once; and the oscillator over 2.2 h, which
an oven oscillator holds.

A mount for B0329+54 can be one motor on a polar-aligned axis, since the pulsar stays at one declination.

## Splitter or switch: the working assumption is a splitter

**Working assumption (Jasper, later on 23 Nov):** "i think assuming a decent power splitter is okay, beams are on
continuously without switching". So both beams are live all the time, each with the full 3 m² and 70 K, and the
integration time per beam is not reduced. Every estimate in this log already assumes that. It is an assumption, not a
confirmation by the team; what the SPDT in the diagram does stays unknown. The three cases, for the record:

- A power splitter after an LNA copies each element's amplified signal to both delay networks. Signal and noise divide
  together, so the SNR is unchanged and both beams are on all the time: no effect on the integration time. The measured
  70 K at the output includes whatever the splitter adds.
- A switch that alternates between the two networks gives each beam half the time: √2 in SNR per beam, twice the
  integration time. It also chops the signal, which puts sidebands of strong RFI at the switching rate into the clean
  150 MHz. A known switching pattern can be treated as blanking (zero weight, as tested for sample gaps on 3 Nov), but
  the pattern must then be recorded with the data.
- A switch set once and left there gives one beam at a time, at full sensitivity.

"SPDT" in a diagram is a single-pole double-throw switch; it is never a splitter. What it switches is not known. The
common roles in such a chain: the beam selection itself; a one-bit delay line (two SPDTs that choose between two line
lengths per element), which would make the two directions selectable, not simultaneous; or a calibration switch at the
LNA input (antenna against a load or a noise source), which has nothing to do with the beams. Four outputs that are live
at the same time need a splitter somewhere, whatever the SPDT does.

## The second beam and polarization (my geometry, ideal crossed dipoles along the grid's x and y)

- In the diagonal plane (φ = 45°) the two polarizations are foreshortened equally at 20° (power × 0.94 each), so their
  gains stay equal. In a principal plane they would differ by about 12 % in power.
- They are no longer orthogonal on the sky: 93.6° apart (cos α = −sin²θ sin φ cos φ / (1 − ½ sin²θ)). That leaks about
  6 % of the linearly polarized flux into the summed power: a few µs for B0329+54, against 460 µs per pass (FM-15).
- Real bowties inside the array will differ from this picture.

## What it changes in the project (nothing changed yet)

- The objective of 2 Nov says "locking onto multiple pulsars using digital beamforming". This hardware has two fixed
  analogue beams. The objective's wording is Jasper's to change.
- Step E4 (the digital beamformer, element calibration, simultaneous beams on chosen pulsars) does not fit this
  hardware. Spatial RFI nulling is not possible. One board suffices, so the multi-board synchronisation of FM-14 and
  action R8 fall away.
- Every budget in the notes and the docs is at SEFD 3.2·10⁵ Jy (worst case 7.4·10⁵ Jy), 400 MHz and a pulsar followed
  through the day. A reference scenario at 6.4·10⁴ Jy, 150 MHz and fixed beams needs new runs.
- Risk register: updated later on 23 Nov at Jasper's "update the risk register for the new array"
  (`2026-11-02_risk-register.md`, section "23 Nov 2026: the register for the 3 m² array"): 12 rows annotated, three
  ratings changed (FM-9 16 → 12, FM-38 12 → 6, FM-15 9 → 4), two scenarios added (FM-48 a drift scan is not what the
  pipeline assumes, 10; FM-47 the beams alternate instead of being split, 6), action R1 changed, R8 dropped, R12–R15 new.
  As listed before the update: FM-43 changes character (hour-angle coverage instead of scan loss); FM-14
  shrinks; FM-15 has a size for the second beam; FM-18 gets the measured 150 MHz; FM-1 and FM-7 matter more, because a
  single pass on a low-flux day is at SNR 5; FM-46 (ADC calibration, alias band) and FM-45 (reference lock) stay.
- The fast simulator has no drift through a beam: a gain that changes within a pass, and one pass per day.

## Scratch runs for the new array (later on 23 Nov; Claude's runs, not Jasper's validation)

Jasper: "can you re-run the numbers using this new default case?" Run in a scratch copy of `runScenario.m`
(`runScenario3m2.m`, derived from the 17 Nov scratch copy); the repository's script and its defaults are untouched, and
the docs still say "not simulated". Receiver: SEFD 6.44·10⁴ Jy (3 m², 70 K); 150 MHz as 48 of the 128 channels (33–80,
1300–1450 MHz, my assumption for where the clean part sits), flat gain; two polarizations; B0329+54 from 52° N; receiver
noise only (variant "clean": no RFI, no scintillation), a perfect clock. ρ = −37.98 dB (was −44.95).

| Case | Setting | Simulated | Scaled estimate given before |
|---|---|---|---|
| A | pulsar followed through the day, SNR 7 TOAs, 30 d | 0.233 h of data per TOA (7 sub-ints of 120 s), 2596 TOAs (87 per day), σ 705 µs, spread of error / σ 1.10 | 0.21 h per TOA |
| B | one pass per day through a fixed beam, 120 d | SNR per pass 11.5 (10–90 %: 10.4–12.6; lowest 9.0), all 120 detected, σ 459 µs, spread of error / σ 1.12 | SNR 11, 460 µs |
| C | position, followed, 30 d, true site 100 km east | error bars 9.4 km east, 11.9 km radial (15 km north with the height known) | 10 / 16 km |
| D | position, one fixed beam, 30 d | the fit is singular (RCOND 10⁻¹⁷): no position | no position |
| E | position, both fixed beams (the second 2.23 h after culmination), 30 d | 108 km along [east 0.958, radial 0.287]; the other axis undetermined | 120 km, east–west only |

How a pass was simulated: the fast simulator has no gain that changes within a pass (FM-48), so a pass is a window of
0.546 h at full gain, which carries the same SNR² as a drift through a 6.3° beam (0.753 × 0.725 h); one TOA per window
(16 sub-ints). The second beam has full gain in case E (about 6 % lower in reality). No worst case exists for the new
array, so none was run.

**The east estimate of case C, and a lesson on seeds.** Case C returned east 73.9 ± 9.4 km for a truth of 100 km
(−2.8σ). My first hypothesis, a bias from `adaptiveSubints`' stopping rule (a group closes when T0 at the predicted
phase reaches the target), was wrong: the same truth with fixed groups of 7 sub-ints (case F) gave 67.3 ± 9.4 km. But C
and F share their seeds (`1e5 ip + 1000 d + c0`), so they share their noise, and one low draw shows in both. The
October runs did the same at the old array's level (36.5 ± 29.6 km for 100 km and 947.2 ± 29.8 km for 1000 km "with the
same seeds", notes §7.4). Eight independent draws of the same truth (600 s sub-ints in fixed groups of two, SNR 9 per
TOA, σ 586 µs): east 101.4, 94.9, 104.0, 85.8, 92.3, 75.3, 88.3, 103.6 km: **mean 93.2 km, spread 10.0 km** (stated
error 9.4–9.5 km), error of the mean 3.5 km; radial mean +0.9 km, spread 12.5 km (stated 12.0 km). So the error bars
are right, and the 26–33 km shortfall of C and F was mostly their noise draw. The mean is still 1.9σ below 100 km
(−6.8 ± 3.5 km; the slope of TOA error on true residual −0.059 ± 0.029): a shrink of a few per cent towards the assumed
site is not excluded and not established. Not pursued further. For any later comparison of two settings: runs that
share seeds are not independent.

Files (gitignored, this Mac only): `.claude/runs/2026-11-23_array3m2/` — `runScenario3m2.m`, `run3m2.m` and
`run3m2.log` (cases A–E), `check3m2.m` and `check3m2.log` (case F, the regressions), `seeds3m2.m` and `seeds3m2.log`
(eight draws), `out/`. Run times: A and C about 6 min each, the rest seconds; the eight draws 24 min.

## Code change: `runScenario.m` gets the new array as its default design point

Asked "Make 6.4·10⁴ Jy with 150 MHz the default design point in runScenario.m?", Jasper: "Yes, make default", and
while I worked: "you can still include RFI to make sure that chain remains tested, however".

- A third design point, `"array3m2"`, is the default of `fluxes`: SEFD `2 k · 70 K / 3 m²` = 6.44·10⁴ Jy on channels
  33–80 of the 128 (1300–1450 MHz). Where the clean 150 MHz sits in reality is not known; the placement is my
  assumption, and the code comment says so.
- `"nominal"` (3.2·10⁵ Jy) and `"worst"` (7.4·10⁵ Jy) stay under their names with all 128 channels. I did not
  redefine `"nominal"`: `expNavFit` and the saved result files go by that name, and the old numbers must stay
  reproducible.
- Per design point the script now takes its channels (`chanSel`), the band-edge taper values of those channels, the
  RFI profile cut to them (new local function `pickChannels`: `chanFreqs`, `level`, `noiseCoeffs`, `validFraction`,
  `nValidBins`, `nChan`), and the scintillation band at their frequencies. So the variants `rfi`, `rfiCal`,
  `rfiCalScint` and `rfiCalScintFixed` run at the new default. Inside 1300–1450 MHz the model has the radar at
  1332.9 MHz, the carrier at 1388.9 MHz and the broadband impulses.
- The result file also stores `sefd` and `chanIdx`. The header describes the three design points and says that the
  script still follows the pulsar through its pass, which for the new hardware is the array on a tracking mount. A
  drift through a fixed beam is not in this change (risk register R1).

Validation (Claude's runs, `.claude/runs/2026-11-23_array3m2/val/`; not Jasper's validation):

| Check | Result |
|---|---|
| Committed script against edited script, `nominal` clean, 3 d | 30 and 30 TOAs, identical in σ, error, SNR, hours, T0, χ², merged count |
| The same, `worst` `rfiCal`, 2 d | 3 and 3 TOAs, identical |
| The same, `worst` `rfiCalScint`, 2 d | 3 and 3 TOAs, identical |
| Edited script at its default (`array3m2` clean, 30 d) against scratch case A | 2596 and 2596 TOAs, identical |

The RFI chain at the new default (B0329+54, followed through the day; 10 days each, 30 for clean):

| Variant | Good TOAs per day | Hours of data per TOA, median (10–90 %) | σ median | spread of error / σ | mean reduced χ² | pChi2 < 0.05 |
|---|---|---|---|---|---|---|
| clean | 86.5 | 0.233 (0.167–0.300) | 705 µs | 1.10 | 0.999 | 3.9 % |
| rfi | 84.3 | 0.233 (0.167–0.300) | 705 µs | 1.11 | 0.985 | 0.6 % |
| rfiCal | 84.7 | 0.233 (0.167–0.300) | 706 µs | 1.13 | 0.998 | 5.1 % |
| rfiCalScint | 89.8 | 0.167 (0.100–0.400) | 703 µs | 1.08 | 1.006 | 12.2 % |

The modelled RFI inside the clean band costs about 2.5 % of the TOAs. With scintillation 12 % of the TOAs have pChi2
below 0.05 (5 % expected); not looked into. Run times: 6.6 min for 30 d clean, 6.7 min for the three variants of 10 d.

## The reruns at the new default, for the docs (Claude's runs, evening of 23 Nov)

Jasper: "the documentation has many numbers still referencing the old situation with 0.61 m^2 area etc. update those
as well in every page". A read-only inventory by a docs subagent found 215 such numbers and 14 questions
(`.claude/runs/2026-11-23_array3m2/docs_inventory.md`). Jasper's answers to the four that block the work:

1. One design point: no worst case for the new array; the old nominal / worst pair stays only as a dated record.
2. "we can get the tracking numbers for integration time, plus the 15-deg beam pass numbers. be mindful about how to
   report these two different figures of merit". So: time per TOA from the runs that follow the pulsar, SNR per pass
   from the fixed-beam passes, in separate tables under their own names. He wrote 15°; 3 m² implies about 6°. Both
   were run; the beam width is to be confirmed with the team.
3. The stage pages' "−54 dB" target is replaced by the new levels (J0437−4715 −44.2 dB, B0329+54 −38.0 dB).
4. The research's other outdated numbers (±2 ppm, ±7.1 µs, 4.15) are not part of this update.

My defaults for the smaller questions, told to him: J0437−4715 from 30° S stays and is rerun; the archived ship
results stay a dated record; unit tests that carry 734 µs as a setting stay; superseded scenario tables stay as
history; block E's endpoint and steps stay annotated.

Code edits the reruns needed (both keep the old design points selectable): `tests/expNavFit.m` takes the design point
from a `fluxes` setting; `tests/expNavPlan.m` has a third typed-in input set for "array3m2" (B0329+54 0.233 h per TOA
at 705 µs from the new run; B0950+08 1.30 h, B1642−03 5.6 h, B1929+10 11.5 h per TOA by formula × 0.108).

**TOA budget, pulsar followed through its pass** (`runScenario.m` at its default; result files
`data/scenario/scenario_<pulsar>_array3m2_<variant>_<n>d.mat`):

| Pulsar | Variant | Days | TOAs per day | Hours of data per TOA (10–90 %) | σ median | spread of error/σ | mean red. χ² | pChi2 < 0.05 |
|---|---|---|---|---|---|---|---|---|
| B0329+54 | clean | 30 | 86.5 | 0.23 (0.17–0.30) | 705 µs | 1.10 | 0.999 | 3.9 % |
| B0329+54 | rfi | 30 | 85.0 | 0.23 (0.17–0.30) | 706 µs | 1.11 | 0.985 | 0.7 % |
| B0329+54 | rfiCal | 30 | 85.2 | 0.23 (0.17–0.30) | 707 µs | 1.10 | 0.999 | 4.3 % |
| B0329+54 | rfiCalScint | 30 | 91.4 | 0.17 (0.10–0.40) | 700 µs | 1.10 | 1.006 | 10.7 % |
| B0329+54 | rfiCalScintFixed | 30 | 98.3 | 0.17 | 784 µs | 5.54 | 1.006 | 10.3 % |
| B0329+54 | clean | 100 | 86.5 | 0.23 (0.17–0.30) | 705 µs | 1.12 | 0.999 | 4.3 % |
| B0329+54 | rfiCalScint | 100 | 87.7 | 0.20 (0.10–0.40) | 701 µs | 1.11 | 1.006 | 10.4 % |
| J0437−4715 | clean | 30 | 6.9 | 1.75 (1.25–2.33) | 12.8 µs | 1.33 | 0.999 | 3.4 % |
| J0437−4715 | rfi | 30 | 6.7 | 1.75 (1.25–2.49) | 12.8 µs | 1.18 | 0.984 | 0.0 % |
| J0437−4715 | rfiCal | 30 | 6.8 | 1.75 (1.25–2.47) | 12.8 µs | 1.19 | 0.998 | 3.4 % |
| J0437−4715 | rfiCalScint | 30 | 10.7 | 0.33 (0.08–3.42) | 12.2 µs | 1.13 | 0.998 | 1.6 % |
| J0437−4715 | rfiCalScintFixed | 30 | 11.2 | 0.33 | 15.1 µs | 2.55 | 0.998 | 2.1 % |
| J0437−4715 | clean | 100 | 7.1 | 1.67 (1.17–2.33) | 12.8 µs | 1.23 | 0.999 | 4.2 % |
| J0437−4715 | rfiCalScint | 100 | 10.1 | 0.33 (0.08–3.60) | 12.2 µs | 1.13 | 0.999 | 3.7 % |

Not looked into: J0437−4715's error bars 23–33 % small in the clean runs; pChi2 < 0.05 for 10–11 % of B0329+54's TOAs
with scintillation.

**Position of a fixed site** (B0329+54 alone, followed, the true site 100 km east; `expNavFit`, `navFit`):

| Run | Good TOAs | east ("self") | radial ("self") | calibrated: ENU errors, polar |
|---|---|---|---|---|
| 30 d, perfect clock | 2634 | 73.9 ± 9.4 km | +7.1 ± 11.9 km | ± [9.4 10.9 7.1] km; polar −4.6 ± 5.4 km |
| 30 d, SDR clock tracked daily | 2631 | 70.8 ± 9.4 km | +9.2 ± 11.9 km | – |
| 30 d, SDR clock not tracked | 1866 | 104.3 ± 9.5 km | −3.0 ± 12.0 km | – |
| 100 d, perfect clock | 8764 | 91.2 ± 5.1 km | +5.0 ± 6.5 km | ± [5.1 6.0 3.9] km; polar −4.6 ± 3.0 km |
| 30 d, truth 1000 km east | 2463 | 972.8 ± 9.4 km | −3.2 ± 11.9 km | ± [9.4 11.0 7.1] km |

The precision is ±9.4 km east and ±11.9 km radial after 30 days (was ±29.6 and ±37.2 km), ±5.1 and ±6.5 km after 100.
The low estimates of the standard runs (73.9 for 100 km, 972.8 for 1000 km: the same 26–27 km short) are one shared
noise draw, as the eight independent draws showed (section "Scratch runs" above). `expNavPlan` at the design point,
30 d: B0329+54 alone 9.5 / 11.9 / 5.4 km (east / radial / polar, calibrated), polar with a crystal clock 28.2 km; with
B0950+08, B1642−03 and B1929+10 8.8 / 9.3 / 5.3 km and 13.5 km. It assumes every pulsar followed in its own beam.

**A single pass through a fixed beam** (scratch copy; a window at full gain with the pass's SNR²):

| Beam width | Window | SNR per pass (10–90 %; lowest of 120) | TOA error per pass | Both beams on B0329+54, 30 d |
|---|---|---|---|---|
| 6.3° (follows from 3 m²) | 0.545 h | 11.5 (10.4–12.6; 9.0) | 459 µs | 108 km along one axis |
| 15° (as Jasper wrote) | 1.299 h | 17.4 (16.6–18.7; 15.6) | 303 µs | 71 km along one axis |

One beam alone gives a singular fit at either width. A 15° beam and 3 m² do not belong to one filled aperture.

**How it ran.** One background MATLAB batch, about 2.4 h of computing, paused for an hour while Jasper drove (the
process stopped with SIGSTOP and resumed by a watcher when the lid opened, so nothing ran in the sleeve). Two
launches of mine failed within seconds before that: `runScenario`'s own variables `R` and `S` overwrote the batch
script's paths in the shared base workspace (review item D2-6 in practice), and `run` with a full path changed the
folder so `pipelineParams` was not found. One step was skipped by my mistake (the motion switch not set) and run
afterwards. Files: `.claude/runs/2026-11-23_array3m2/rerun/` (copies `rr*.m`, `rerunMain.log`, `post.log`), and
`results_for_docs.md` beside it, the table of numbers given to the docs subagents.

## The docs moved to the new design point (evening of 23 Nov)

Two docs subagents, split by page so that no file had two writers, both working from one table of numbers that I
wrote from the run logs (`results_for_docs.md`) with Jasper's four decisions. One took the result pages (scenario,
home, results, navigation, simulator, limits, localisation, scripts, tests, parameters), the other the stage and
overview pages (roadmap, physics, receiver, excision, channels, fold, toa, fullband, pipeline, glossary, formulas,
lessons, history; no change needed on detection, sky, conventions, validation). I wrote the changelog entry in
`maintenance.html`.

- Scenario page: rewritten except R3. R1 is the receiver as it is now; new R1.2 "Two figures of merit"; new R7.4 "A
  single pass through a fixed beam"; the old receiver and the superseded tables are in History boxes.
- The stage pages' −54 dB became −44.2 dB (J0437−4715) and −38.0 dB (B0329+54); numbers that follow by formula were
  recomputed (for example one radar pulse against 80,000 pulsar pulses is now against 8,300).
- 16 figures regenerated from their generators (scenario R.1, R.2; navigation 10.1, 10.3, 10.4; results 1; limits 1;
  localisation L.3; parameters 2; scripts 1; roadmap 1; physics P.5; receiver 3.3; excision E.2, E.6; fold 7.4).
- My checks: all 29 pages parse, every internal link and anchor resolves, no duplicate ids; the old values that remain
  are in History boxes, dated records or runs at their own level; I read the new R1.2 and R7.4 in full and checked the
  rescaled 8,300.

Not replaced, and named as such on the pages:
- the locked-radar biases in the excision tables (runs at −54 dB; `runLockedRadar.m` has its levels typed in, so the
  new levels need a code change and a rerun);
- the fractions of TOAs beyond 3σ and 4σ (October runs; `navFit`'s 2 % rests on them). They can be read from the new
  100-day result files without a new run;
- the realistic RFI list's link budget on the receiver page (made for the earlier array, not re-estimated);
- the archived ship results; the 6-day SDR tracking table; the unit tests that carry 734 µs as a setting.

For Jasper to look at:
- The roadmap's "Physics notes from the revision" (23 Oct) now carry the 23 Nov planning-study numbers. They were a
  dated note; the subagent kept a backup of the page as it was.
- Five new rows in the history page's decision log, in the subagent's wording of his decisions of 20 and 23 Nov.
- `expNavFit` and `expNavPlan` went from `validated` to `written` on the tests page: the scripts were edited and the
  new numbers are Claude's.
- Three times in figures P.5, 3.3 and 7.4 (41 min, 2.3 min) are for the voltage chain's 10 ms test pulsar in 400 MHz
  with one polarization; the captions say so.
- The pages contain derived values that the subagents computed with the pages' own formulas (listed in their
  reports); they are not from runs.
- `api.html` is to be regenerated (`docs/tools/make_api.py`) after the commit of `runScenario.m`.

## Four inputs instead of one (Jasper: "can you confirm whether the chain still works for that")

The four analogue outputs are two beams × two polarizations. What the code does today (my reading, 23 Nov):

- `processIQ` takes one IQ stream and gives one fold. Four outputs are four independent runs of it; nothing in it
  needs the other streams (the excision, the dedispersion and the detection are per stream).
- From the fold onward the chain already handles two polarizations: the fast simulator's 'NPol' 2 gives the
  total-intensity fold, the observer passes Bnoise × 2, and every scenario run since 20 Oct uses it.
- Missing at voltage level: one step that adds the two polarization folds of a beam into one. Sums and weights of a
  fold add (`mergeSubints` uses that across time), so adding `sum`, `weight`, `weight2` and `weightX` of the two folds,
  with `prof = sum ./ weight`, is exact, and the noise model then needs the plain Bnoise, not Bnoise × 2. It needs the
  two folds on the same sub-int grid (the same predictor and time base; one board gives that) and each stream in its
  own noise units, which the chain's normalization already does.
- The two beams are two observations, not one: a different part of the sky, so a different hour angle of the same
  pulsar, or another pulsar. They are not added. Their TOAs go into `navFit` as TOAs of the same pulsar only if the
  delay between the two beam paths is known (action R12); a free constant per beam would absorb the very difference
  that carries the position.

Scratch test (Claude's run, `pol2test.m`): fast simulator, the new default design point, 200 passes of 16 sub-ints.

| Case | SNR per pass (median) | σ (median) | spread of error / σ |
|---|---|---|---|
| one stream alone, Bnoise | 8.01 | 657 µs | 1.09 ± 0.05 |
| two streams with independent noise, folds added, Bnoise | 11.34 | 465 µs | 1.10 ± 0.06 |
| the simulator's 'NPol' 2, Bnoise × 2 | 11.25 | 469 µs | 1.12 ± 0.06 |

Added over one stream: 1.414 in SNR (√2). Added over 'NPol' 2: 1.008. So the TOA stage takes an added fold as it is.
Not tested: two voltage streams through `processIQ` and then added (the generator makes one stream); two streams with
different blanking masks; unequal gains or T_sys in the two polarizations. No function for the adding exists yet; it
would be one new function with its unit test, for Jasper's OK. On the to-do list at his "put the addition of this
function on the todo list": notes §10 item 8 and the roadmap's open item "Second polarization".

## Open (for the team, through Jasper)

1. Fixed panel or a mount, and which kind. Position needs several hour angles.
2. Splitter or switch: assumed a splitter, beams live continuously (above). To confirm when the chance comes: are all
   four outputs live at the same time, and what does the SPDT switch?
3. A measured SEFD or G/T per beam (for example on the Sun): it confirms the 3 m² and the 70 K in one number.
4. Where the LNAs sit (at each bowtie, per tile, or after the last combiner).
5. What the second beam is for.
6. Which 150 MHz is clean, and whether the rest is usable at a higher T_sys.

## Files

- Notes: `logs/current_project_notes.md` §11, "Array, update 23 Nov 2026" (after the 20 Nov paragraph).
- Docs, later on 23 Nov at Jasper's "set up a subagent to change the docs and their outdated numbers": one docs
  subagent, briefed to change only what the new array makes wrong and to leave every simulated number at its design
  point. 14 pages, 163 lines added, 38 removed; a new Section R1.1 on `scenario.html` ("The hardware as of 23 November
  2026"); design-point labels on the simulated results; SEFD, sensitivity and G/T in the glossary; changelog entries in
  `maintenance.html` and `history.html`. I checked its diff on six pages and corrected one sentence (the 30 / 47 km
  provenance). Left for Jasper: the two figures that still state the old array (home page hero, roadmap Figure 1, and
  `gen_roadmap.py`, which still says 18 streams), steps E3, E3b and E4 (annotated, not rewritten), the design points,
  and the outdated numbers from the 21 Nov research (±2 ppm, the 4.15 threshold, ±7.1 µs), which were out of scope.
- Not changed: the Task 2.4 deliverable; the code; `api.html`.
- The repository is public; this log describes the team's hardware. Not committed.
