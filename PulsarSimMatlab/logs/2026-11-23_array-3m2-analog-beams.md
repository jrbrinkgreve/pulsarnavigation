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
30 / 47 km budget follows the pulsar through the day.

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
- Risk register, when Jasper asks for it: FM-43 changes character (hour-angle coverage instead of scan loss); FM-14
  shrinks; FM-15 has a size for the second beam; FM-18 gets the measured 150 MHz; FM-1 and FM-7 matter more, because a
  single pass on a low-flux day is at SNR 5; FM-46 (ADC calibration, alias band) and FM-45 (reference lock) stay.
- The fast simulator has no drift through a beam: a gain that changes within a pass, and one pass per day.

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
- Not changed: the docs (`scenario.html`, `receiver.html`, `roadmap.html`), the risk register, the Task 2.4 deliverable.
- The repository is public; this log describes the team's hardware. Not committed.
