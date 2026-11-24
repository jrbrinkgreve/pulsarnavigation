# 2 Nov 2026: risk register (draft), and two excision questions

Session log. No change to the pipeline code; one new test file (`tests/testNavFitFaults.m`). Committed at Jasper's "Commit now".

## The session

Jasper started the read-through of the paragraph-flow docs (`d5a89c0`) in his own browser. Two questions came up on the
excision page, and then the request for this register.

**Can excision blank every channel?** Continuous RFI does not get blanked: it raises the channel's own baseline
(median/ln 2 per 16 ms block), so nothing rises above the thresholds, and the optimal weighting gives the channel little
weight instead. A flag blanks a few µs, not a channel: at most 16 samples (3.8 µs) plus a 10-sample guard on each side.
Dense RFI pulls the median up, so the detector then flags less, not more. **Correction, given in the session:** my
first answer missed review item D2-22. In `testDetectRFI` test 6 (27 Oct, validated), strong bursts at independent times
in every channel made the frequency guard blank 99 % of the data, against 33 % from the own detections. D2-22 is in
`tests.html`, but not in the excision page's limits (E9). Jasper: "Not now" to adding the answer to the excision page.

**What does `Passes` do?** Each pass computes the baseline without the previous pass's flagged samples, tests the windows,
and rebuilds the coverage. The mask is the last pass's coverage only (`detectRFI.m:233`).
- `Passes = 1`: the baseline includes the RFI (5.4 % high in test 6), and the per-channel thresholds of E2.2 never
  switch on, because they start at pass 2.
- `Passes = 2` (the default): the baseline is 0.5 % high, and 58 % of the weak bursts are wholly blanked, against 51 %.
- `Passes = 3` has never been run. My estimate is a small gain on the remaining 0.5 %, at extra CPU only (each channel is
  read once).

**The request.** Jasper: "just drafting: can you think of all scenarios in which something goes wrong in this pipeline,
rate their probability x impact, and put tests + developments on the to-do list if they require action".

## Risk register (draft)

### Scope and scales

The scope is the planned use: block E (a fixed site in the Netherlands, in simulation), then block F (real recordings at
that site). Ship-only scenarios are rated for the fixed site. Each scenario is rated as the pipeline stands today, so a
known gap that is already planned still scores high; the Action column says where it is planned.

| P | Probability |
|---|---|
| 1 | rare: needs an unusual combination |
| 2 | unlikely |
| 3 | possible during a long campaign |
| 4 | likely |
| 5 | certain, or it happens in today's runs |

| I | Impact, if it happens and nothing catches it |
|---|---|
| 1 | negligible: below 1 % on the error bars, a bias far below the error |
| 2 | minor: a few % on the error bars, or fewer TOAs |
| 3 | moderate: error bars off by more than 10 %, a bias about the size of the error, or many TOAs lost |
| 4 | major: a bias of several errors, or no fix within the campaign |
| 5 | critical: no pulsar at all, or a confident wrong position, without a warning |

The score is P × I: 15–25 high, 8–12 medium, 6 and below low. "Seen" says whether the pipeline would notice by itself
(an error, a warning, χ², missing TOAs). Numbers marked "estimate" are mine. The four estimates of the first draft
were checked afterwards, and two ratings changed (see "The estimates, checked").

### What stands out

Three of the high risks concern real data entering the processing. `processIQ` has no predictor (FM-22). Sample drops
and restarts would shift the time stamps unseen (FM-9). An inverted or mis-set frequency mapping would smear the pulse
away (FM-8). Together they decide whether F detects anything at all, and two of them give no reason when it does not.

The one high risk that can be tested now, in the simulator alone, is the oscillator's daily temperature term (FM-11).
It lands on nearly the same period as the position signal, so the fit would read it as position with honest-looking
error bars. After the check it ties FM-22 at the top: with the planned ±2 ppm TCXO, a room's daily temperature swing
puts a bias of roughly 2 to 25 times the fit's own error on the position.

The timing-model terms (FM-5, FM-30) are high too, but they are already step E2.

**21 Nov, after the literature research:** a new entry joins the top, the array's effective area against direction
(FM-43, 20). A zenith-facing panel of patches gives B0329+54 about the "worst case" sensitivity as a day average, and it
costs the position more than the sensitivity, because the pulsar is seen well only near culmination. FM-11 keeps its
score, but its answer is now sized: an oven oscillator on the reference input, where fitting the daily line instead
costs a factor 6.7 at 30 days. Details in the section "21 Nov 2026: the register against the literature".

**23 Nov, the 3 m² array:** the hardware changed (about 3 m², two hardwired analogue beams, four outputs; SEFD
6.4·10⁴ Jy, about 150 MHz clean). Sensitivity is no longer the first worry: B0329+54 should show in a single pass. The
top entry beside FM-22 and FM-11 is now what FM-43 has become: with fixed beams a pulsar is seen at one hour angle a day,
and that gives no position. The mount decides it. Details in the section "23 Nov 2026: the register for the 3 m² array".

### The register, highest score first

| ID | Stage | Scenario | P | I | P×I | Seen | Action |
|---|---|---|---|---|---|---|---|
| FM-22 | 7 fold, `processIQ` | `processIQ` folds with the linear f0/F1 model only, without a predictor (D1-2). Real data drift by up to 10⁻⁴ (Doppler), about 0.36 s per hour, so any fold longer than seconds smears the pulse away. The noise calibration's on-pulse window and `rfiProfile` use the same model (D1-3, P3-3). **21 Nov:** every reference backend (DSPSR, PRESTO) folds with a polynomial or TEMPO2 predictor, so a 300 s fold of real data without one will not match theirs. Public voltages to try this on before the hardware exists: B1642−03 from Dwingeloo (438 MHz, 12 GB, with the observatory's own DSPSR fold of the same samples) and B1133+16 from Green Bank (1125–1875 MHz). | 5 | 5 | 25 | partly: no detection, without the reason | existing: D1-2, D1-3, P3-3 (F); the first thing F needs |
| FM-11 | 4a clock | The oscillator follows its temperature, which cycles daily. The clock gets a rate term at the solar day, 4 min per day from the sidereal day of the position signal: over 30 days the two drift apart by only 30° of phase (99 % of the term looks like position), over 100 days by 99° (88 %). A rate amplitude *a* gives a daily phase of *a* · 86400 s / 2π, which is 140 µs at 10⁻⁸. On B0329+54 that is 72 km, since 100 km east gives a 193 µs sine (`navigation.html`), or 2.4 times the fit's ±30 km equatorial error. A ±2 ppm TCXO changes by roughly 10⁻⁸ to 10⁻⁷ per K (estimate: the ±2 ppm band over a rated range of about 100 K, steeper in places; the SDR's datasheet decides). A room's daily swing of ±1 K then gives 70–700 km, 2–25 σ; an outdoor box more. An OCXO (about 10⁻⁹ per K or better) keeps it below 7 km. The simulated crystal has no such term (D1-17). **21 Nov (data sheets; the slopes are range totals divided by the range):** the TCXO slope of 10⁻⁸ to 10⁻⁷ per K is confirmed. An OCXO is 2·10⁻¹¹ to 1.1·10⁻¹⁰ per K, so "below 7 km" is below 0.8 km at ±1 K and below 4 km at ±5 K; a double oven, a rubidium or a chip-scale atomic clock stays below 0.1 km at ±1 K. The favoured RFSoC 4x2 has no TCXO but a bare 48 MHz crystal (±15 ppm, ±25 ppm over temperature, distributor listings; about 1–5·10⁻⁷ per K, estimate): 700–3500 km at ±1 K. Fitting a free solar-day line instead costs a factor 1/√(1 − ρ²) on the position error: 6.7 at 30 days and 2.1 at 100 days for one pulsar; 1.40 with a second, equally strong pulsar 6 h away in right ascension and 1.00 at 12 h; a faint second pulsar helps little (researcher's derivation; it reproduces the 4 Nov table within 5 %). Measured crystal clocks also show the second and third solar harmonics, which cost extra with part-day tracks. **23 Nov:** unchanged in kind with the new array. With both fixed beams on B0329+54 the position comes from the difference of two TOAs 2.2 h apart, and a daily clock line of amplitude A puts 0.57 A into that difference: about 0.5 km east per µs, the same scale as before (my arithmetic). | 5 | 5 | 25 | no: the fit reads it as position (30 days: χ² unchanged even at 10⁻⁷; 100 days: χ² 1.38 at 10⁻⁷) | existing: D1-17 (E5); **new: test, before E5**; test written: `testNavFitFaults` 1; **21 Nov:** a single-oven 10 MHz OCXO on the reference input as the minimum (R4) |
| FM-43 | E3, E4 array | *(21 Nov)* The array's effective area depends on direction and frequency, but the budgets assume 0.61 m² for every pulsar, 24 h a day. A 4×4 lattice of 7–9 dBi patches at 0.195 m pitch realizes 0.38–0.49 m² near broadside (band mean, 90 % efficiency) and 0.14–0.20 m² at 44° from its axis; at 1.6 GHz grating lobes take 30–60 % even at broadside (researcher's pattern model, ±1 to ±1.5 dB; no measurement exists for such an array). A zenith-facing panel at 52° N then has a day-averaged SEFD of about 2.7 times nominal for B0329+54, which is the "worst case" design point, and 5–8 times for the pulsars at declination 8–22°. The position suffers more than the sensitivity, because the data bunch around culmination, where the daily cosine is hard to tell from the pulsar's constant: a one-day fit weighted by that area has the east error ×3.8 and the north error ×6.1 (my scale estimate: white noise, constant T_sys), so 30 / 47 km becomes roughly 110 / 290 km. With tracks of ±3 h around culmination only, the north component is nearly lost (118 / 558 km, researcher's model). A panel that keeps the pulsar near its axis costs ×1.45. **23 Nov (the 3 m² array):** the cause above belongs to the 0.61 m² lattice at one wavelength pitch and is gone: the new array has half-wavelength pitch and no grating lobes for its beams. What stays is the hour-angle coverage, and it is sharper. The two beams are hardwired, (0°, 0°) and (20°, 45°), about 6° wide, with delays "fixed for good", so a pulsar is seen only while it drifts through a beam: 0.72 h a day for B0329+54. Seen at one fixed hour angle, the site's delay is a constant that merges with the pulsar's constant, so one beam gives no position at all. Both beams on B0329+54 (the second catches it 2.2 h from culmination, with the panel turned about its normal) give one combination, mostly east–west: about 120 km after 30 days with 150 MHz. A tracking mount gives about 10 / 16 km east / north (all three scaled from the simulated budget, not simulated). The mount is "still open to decide". Pointing a 6° beam within 0.5 dB needs ±1.3°. | 5 | 4 | 20 | partly: fewer TOAs and a larger covariance; the forecasts do not show it | **new: development** (the area against direction and frequency in the scenario runs, R1), then a decision on the panel: zenith, tilted, pointed or sub-arrays (E4); **23 Nov:** the decision is now the mount: fixed, one polar axis, or full tracking (R14); the simulator needs a drift through a fixed beam (R1, changed) |
| FM-44 | 7a prediction, E1 | *(21 Nov)* Catalogue spin parameters in mixed time-scale units. ATNF v2.7.0 gives B0329+54 and B2016+28 in TDB, and B0950+08, B1133+16, B1642−03, B1749−28, B1929+10 and B1937+21 in TCB. A TCB frequency read as TDB is wrong by 1.55·10⁻⁸: 1.34 ms per day between the two groups, which a common clock rate cannot absorb, and it breaks B1937+21's turn count within a day. It is not an error today, because the simulated truth and the observer share the same values (`runScenario.m:112`). It becomes one with real data or the PINT/TEMPO2 check. | 4 | 4 | 16 | yes: a residual trend of 1.3 ms per day | **new: development** (E2: one par file per pulsar from one release, converted once to TDB; R3) |
| FM-5 | 7a prediction | The pulsar's catalogue position is at an old epoch; proper motion and parallax are not modelled. For B0329+54 that is about 47 µs per year since the position epoch, as an annual term (notes §10 item 13, to check). A 30-day fit absorbs most of it into the clock; the rest biases a multi-pulsar solution. **21 Nov:** 47 µs per year is the upper bound (a pulsar in the ecliptic); for B0329+54's geometry it is 37.9 µs per year (researcher's calculation). ATNF v2.7.0 gives its position at MJD 56000 (VLBI), so leaving proper motion out today gives about 557 µs (167 km) of annual amplitude, with a slope of up to 9.6 µs per day; B1133+16 about 8 ms, B1929+10 about 2.5 ms. "A 30-day fit absorbs most of it into the clock" holds for one pulsar only: the slope differs per pulsar, and `navFit` has one common rate (`navFit.m:158`). | 5 | 3 | 15 | no | existing: E2 (item 13, first in its proposed order); **21 Nov:** proper motion and position epoch in E2 (R3), a rate per pulsar in `navFit` (R2) |
| FM-8 | 4 IQ, 5a | The real receiver's frequency mapping differs from the processing's: an inverted band (a mixer with its LO above the signal), or another LO, band or sample rate. The dedispersion chirp then has the wrong sign or size, and the pulse smears over up to twice the band's sweep (31.5 ms at DM 25, against B0329+54's 6.6 ms pulse). At the weak target, "inverted" looks the same as "too weak". The mixer type is not known yet. **21 Nov:** a real example: the Breakthrough Listen paper says Green Bank L-band voltages are band-inverted, but a 2016 L-band file header has a positive bandwidth and the 2017 C-band files a negative one, so the sense must be read per file and checked on the pulsar. A ZCU216 can sample 1.2–1.6 GHz only in its second Nyquist zone, which inverts the band. | 3 | 5 | 15 | no | **new: development + test** (before F); test written 3 Nov (`testFrontEndFaults` 1: 4.7 % of the amplitude left, half a turn off, all sub-ints still "good"; B0329+54: SNR 7 → 2.9) |
| FM-20 | excision | `detectRFI`, `periodicRFI` and `blankingWeights` work on whole files, and the weight file is 4.8 GB per second of data (D2-9, D2-23). Hours of real data cannot run. | 5 | 3 | 15 | yes: a memory error | existing: F (a design exists: time blocks with the filter reach as overlap) |
| FM-29 | 10 navigation | With one pulsar, the clock's wander looks like a polar offset: the polar component stays at about 64 km after 30 and after 100 days (nominal flux). The fit's clock covariance also leaves out flicker (D1-19). **21 Nov:** pulsars at other declinations give the polar component only together with known absolute phases; with a free constant per pulsar, declination spread alone gives nothing (Thompson, Moran & Swenson ch. 12; `2026-11-04_static-redesign.md` item 3 has it right). The other route is a known height. With an OCXO the unmodelled 30-day wander falls from 12.5 ms to about 0.3 ms, with a rubidium to about 13 µs (estimates on assumed noise floors). **23 Nov:** with the new array's fixed beams, pulsars at other declinations are out of reach from one site; only a mount brings them in. | 5 | 3 | 15 | yes: in the fit's covariance | existing: E3 (pulsars at other declinations), E5 (oscillator) |
| FM-30 | 7a, simulator | Truth and observer share `delayModel`, so a term missing from it is missing on both sides, and no simulation shows what it costs (notes §10 item 13). **21 Nov:** `delayModel`'s four terms agree in form and sign with TEMPO2 and PINT (a reading of the code and the papers, not a numerical comparison). Missing against the full model: proper motion and position epoch, the binary orbit, every dispersive term, the troposphere, planetary Shapiro, polar motion, a UT1 series, full nutation, the clock chain, F2. Expected difference against PINT as the code stands: 0.1–0.2 µs daily (estimate). PINT does not run on native macOS ARM Python (no extended precision), so the check needs a Linux container or Rosetta. TEMPO2 and PINT cannot fit a ground site's position, so the site fit stays outside them. | 5 | 3 | 15 | no | existing: E2 (each term with a switch; the PINT/TEMPO2 end-to-end check) |
| FM-4 | sky, 7a | The spin ephemeris ages, normal pulsars have timing noise, and B0329+54 has glitched. A drift in one pulsar's phase is not the common clock, so it leaks into the clock rate and the position. The sizes need literature numbers per pulsar (item 13). **21 Nov:** "B0329+54 has glitched" is not supported: neither the Jodrell Bank nor the ATNF glitch table lists one. Of the northern set only B1642−03 has slow glitches (77–228 µs per day). B0329+54's timing noise is 7.46 ms rms over 36.5 years (Hobbs et al. 2010), which scales to below 0.2 µs over 30 days after a ν, ν̇ fit. The real term is the age of the catalogue's spin parameters (B0329+54: epoch 1986; the other northern normal pulsars: MJD 59400): a frequency error of 10⁻¹⁰ to 10⁻⁹, 10–150 µs per day (estimate), different per pulsar. About a third of an unfitted drift leaks into the daily sine: 1–14 km per pulsar, while `navFit` has one common rate. A phase and a rate per pulsar cost the position under 0.1 % beyond three days (researcher's model). | 5 | 3 | 15 | partly: a residual trend | existing: E1, E2 (numbers, a red-noise truth); **new: test**; **21 Nov:** a phase and a rate per pulsar in `navFit` (R2) |
| FM-18 | 3 RFI | The site's real RFI differs from the 16-source model (`rfiRealistic`): more radars, DME pulses just below 1.2 GHz at the band edge, the SDR's own spurs and clock harmonics. A strong steady line costs most of its channel, one of 128 (0.78 % of the band). A comb of self-made lines costs one channel per line in the band: 40 channels (31 %) for the harmonics of a 10 MHz reference, 4 (3 %) for 100 MHz. **21 Nov:** the model is now known to differ. DME is inside the band, not below it: the ground replies of channels 113X–126X lie on 1200–1213 MHz (Pampus near Muiden on 1212 MHz; pulse pairs of 3.5 µs, 12 µs apart). Solid-state radars send chirps of 17.5–410 µs at 4–10 % duty on two or more carriers (ITU-R M.1463-3), where the two model radars have 2 µs pulses. The mobile downlink blocks are 1452–1467, 1467–1482 and 1482–1492 MHz (Staatscourant 2020, 41318), not 15/10/15 MHz, and within 1–10 km of a mast they are 33–63 dB stronger than modelled (free-space estimate). Missing from the model: Galileo E5b, BeiDou B2b and GLONASS L3 at 1197–1217 MHz, GLONASS L1 at 1598–1605 MHz, the Cabauw wind profiler (1290 MHz), radar satellites in 1215–1300 MHz (about +46 dB for seconds, a few times a month; estimate) and amateur repeaters in 1240–1300 MHz. Observatories keep 68–91 % of 1.2–1.6 GHz; Westerbork treats most of the range below 1300 MHz as unusable. **23 Nov:** measured on the new array, per Jasper: about 150 MHz of 1.2–1.6 GHz is clean at 70 K. That is SNR × 0.61, or 2.7 times the integration time, against the 400 MHz of every budget. Which 150 MHz, and whether the rest is usable at a higher T_sys, is not known. The analogue beamformer leaves no way to null a transmitter spatially. | 5 | 3 | 15 | partly | **new: development** (a site survey); **21 Nov:** update `rfiRealistic`, and survey for at least two weeks (R6) |
| FM-9 | F ingest | The SDR drops samples (overflow at 18 × 400 MS/s) or restarts with a new time origin. The processing counts samples from `info.t0`, so every later time stamp is off by the lost time. That is a clock step the fit's clock model does not expect, and many small drops add up like extra clock noise. **21 Nov:** 16 elements since 20 Nov (32 streams if dual-polarized). GUPPI RAW, VDIF and UHD mark lost samples explicitly, and published backends zero-fill and flag them, as planned here. An RFSoC design has no such layer by default: the firmware must put a sample counter latched on the PPS into every packet. Timing arrays fit instrumental steps (up to about 1 µs at Parkes) as jumps; `navFit` has no segment per restart. **23 Nov:** four streams on one board with the new array (two beams × two polarizations): 26–64 Gbit/s against the 80 Gbit/s published for the RFSoC 4x2's link, so an overflow is less likely. **24 Nov:** time-stamp faults published by observatories: 97 ms and 409 ms (Effelsberg: a clock source switched, a power cut), 0.67 or 1.34 s (uGMRT: one data buffer), over 1 ms after converter power cycles (Apertif), jumps of 4–18 µs (UTMOST-NS), and a slip of one clock period from the PPS edge on an RFSoC back end in laboratory stress tests (21CMA). All were found by timing or on the bench, none by a detection, and a pulsar shows a step only modulo its period (a whole second is 0.40 turn of B0329+54). The published defence is a pulse train locked to the PPS that arms the acquisition, measured at every start, with every re-synchronisation logged (R19). | 3 | 4 | 12 | no | **new: development + test**; processing half tested 3 Nov (`testFrontEndFaults` 3: a joined drop is an exact clock step; a zero-filled gap blanked through the excision path is harmless); **24 Nov:** R19 |
| FM-1 | sky | The pulsars are fainter than the worst-case design. Catalogue fluxes at our band are uncertain, and B0329+54 alone carries most of the information. At worst-case flux, B1642−03 and B1929+10 give no TOA in 30 days. **21 Nov:** each catalogue flux is one paper's number. B0329+54 is 200 mJy ± 25 % as a mean (203 ± 57 is Lorimer et al. 1995; no L-band series newer than the Urumqi one was found). Its refractive modulation index is 0.31–0.47 (measured at 1540 MHz), so the 5th-percentile flux is about 0.45 of the mean for a 2 h average and 0.75 for 30 days (researcher's Monte Carlo): ×5 and ×1.8 in time, on top of the SEFD worst case. B1929+10's 29 mJy is the lowest published value (others 41 ± 6 to about 90 mJy); B1642−03 is 11–55 mJy. No pulsar detection at 1.4 GHz with less than about 3.6 m² of effective area was found in the literature; well-run small systems obtain 0.5–1.0 of the ideal radiometer SNR. **23 Nov (the 3 m² array: SEFD 6.4·10⁴ Jy at the measured 70 K, about 150 MHz clean):** the margin is now set by the team's goal, a detection within a single pass through the beam. B0329+54 through the beam centre gives SNR about 11 per pass at average flux, 7 if the panel is not tilted towards it (2.6° off centre), and 5 on a low-flux day (scaled from the simulated 1.93 h per TOA, not simulated); real small systems obtain 0.5–1.0 of the ideal. SNR 7 per pass needs SEFD ≤ 1.0·10⁵ Jy. No other northern pulsar detects in one pass (1.4–9.6 h on source needed, against 0.55 h per pass), and with fixed beams only pulsars at declination 32–72° pass through a beam at all from 52° N. **24 Nov:** where surveys compare redetected pulsars with expectation, measured over expected flux has a median of 0.75–1.0 and a spread of a factor 1.6–2 either way (LOTAAS: 0.74 +0.99/−0.42), and an SEFD derived from pulsars reads high (Apertif: 85 Jy on calibrators, 130–200 Jy from pulsars, B0329+54 among them). One pass therefore tests the instrument to a factor 1.5–2 at best, ten passes to about 30 % (estimate). UTMOST-NS, a fixed array with the same rule (a TOA counts above S/N 7), gets one from 53 % of its observations. The simulated 11.5 per pass is the matched-filter figure and equals the radiometer equation with W10; the same pass reads 25 with W50, so a comparison of measured with expected S/N has to name the convention (R21). | 3 | 4 | 12 | yes: few TOAs | existing: E1 (fluxes with uncertainties), E3b (below-threshold processing), E6 (levers); **24 Nov:** R16, R21 |
| FM-10 | 4a clock, 7a | A cold start without a clock calibration. The SDR's ±2 ppm moves the pulse 7.2 ms per hour, more than B0329+54's 6.6 ms pulse, so hour-long folds smear. The runs assume a port calibration to 10⁻⁹ in rate, from an external reference such as GNSS. Without one, as in a GNSS-denied start, the first folds may detect nothing, and without TOAs there is no clock estimate. **21 Nov:** the RFSoC 4x2's reference is a bare crystal, ±15 ppm and ±25 ppm over temperature (the reference manual's bill of materials; the ppm values from distributor listings). At 25 ppm that is 90 ms per hour, 3 ms of smear in a 120 s sub-int, and about 500 rate trials instead of 40. All three candidate boards accept an external 10 MHz (FM-45). **24 Nov:** over one pass a 20 ppm offset moves the pulse 40 ms, six times its width. A constant rate offset is a period offset, so a period search recovers the detection: about 60 trials over ±25 ppm, which raises the unknown-phase threshold from 4.5 to about 5.3 (estimate), and the test of the period against the Doppler prediction (±82 ppm over the year) is lost. The one paper with this board on a telescope (Effelsberg, 2025) does not say what clocked it; the 21CMA and BURSTT back ends use an external 10 MHz and PPS (R19). | 3 | 4 | 12 | yes: no TOAs, but no way out | **new: test, then development**; **24 Nov:** R19 |
| FM-42 | 8–9 TOA, detection | *(3 Nov)* Common-mode gain noise: the receiver's gain fluctuates (1/f, temperature, supply) in every channel at once. The noise model has independent radiometer noise only, so the per-channel χ² sees 1/128 of the excess while detection and TOAs take all of it, most at the low harmonics that detection uses. B0329+54 with a 1/f knee at 10 Hz: detection variance ×2.7, unknown-phase false alarms 10⁻³ → 0.36 per sub-int, TOA error +10 %; knee 1 Hz: false alarms ×4, TOA +1 %. The knee of the real chain is not known (1–10 Hz plausible, estimate). **21 Nov:** no measurement exists for an uncooled L-band chain with an RFSoC. MeerKAT's L-band knee, scaled to 400 MHz, is about 2 Hz or more (estimate), so 1–10 Hz is the right order. **24 Nov:** the literature has no measured red-noise loss for a fold with a known ephemeris. The published factors (PALFA 1.1–2 at periods of 0.1–2 s, the CHIME stacking search an order of magnitude at 1 s and longer) are for blind Fourier searches and do not transfer. The numbers of 3 Nov above are therefore the only ones for this case, and the knee has to be measured: the off-source total-power spectrum from 0.01 to 100 Hz, per beam and polarization (R17). If the LNAs sit before the splitter, their gain fluctuations are common to both beams, and the beam without the pulsar could serve as a reference (my idea, not checked; the LNA placement is not known). | 3 | 4 | 12 | partly: `noiseRatio` rises, nothing acts on it | **new: measure** (Allan variance), then a per-harmonic noise check in `detectPulsar`; test written (`testFoldFaults` 2; log 2026-11-03); **24 Nov:** R17 |
| FM-45 | 4a clock, F | *(21 Nov)* The board's reference is not the one assumed, or is not locked. The RFSoC 4x2 runs on its ±25 ppm crystal unless an external 10 MHz is locked (FM-10). With the reference missing, the clock chip's second loop keeps every output running from the on-board VCXO, so a run can be taken on the free oscillator without notice (forum reports). The one published 10 MHz register file gives the right frequency but no deterministic phase (a researcher's decoding, to be checked). | 3 | 4 | 12 | no, unless the first loop's lock status is read | **new: development (small)** (log the lock status and the board temperature with the data, R5; a bench test of the 10 MHz lock) |
| FM-14 | E4 array | The 18 streams are not phase-coherent or not calibrated (several SDRs need one reference and a PPS), so the beam gain falls below its budget. **21 Nov:** 16 elements since 20 Nov (32 streams if dual-polarized). No restart-to-restart phase repeatability is published for any RFSoC board, and nothing at all for several RFSoC 4x2 boards (AMD guarantees ±1 sample between tiles), so a phase calibration at the start of every session is part of the design. That calibration is thermally easy on the Sun: a signal-to-noise of about 300 per baseline, per 3.125 MHz channel, per second (estimate); BURSTT does it on ZCU216 boards. 0.5 dB of beam gain allows 20° rms of phase. Cas A, Cyg A and the Moon drown in coupled receiver noise. **23 Nov:** the new array combines its 256 elements in a hardwired analogue network and delivers four outputs, so the digital side cannot calibrate, re-phase or re-point anything, and one board suffices (the multi-board question falls away). A delay error, a dead element or tile, or a beam that points off its nominal direction costs gain for good and shows only as a lower SNR. The beams have to be measured as built: drift scans of the Sun give each beam's direction, width and SEFD. **24 Nov:** every new array with a measured figure came in above its ideal SEFD: Apertif 1.2 times on calibrators, BURSTT 1.8–2.0 times (11–12 kJy against 6 kJy), the CHIME Pathfinder a few times. Two arrays with analogue sums found path errors after assembly, up to 15 cm (EMBRACE) and 20 cm (the MWA's Engineering Development Array); at 1.4 GHz an rms of 1 cm costs about 4 % of the SNR and 2 cm about 16 % (estimate, exp(−σ²/2)). Composed from the published loss terms, the panel delivers 0.72 (0.37–0.98) of its nominal sensitivity if the 3 m² and 70 K were not measured for the whole beamformed panel on the sky, and 0.90 (0.74–0.98) if they were: SNR 8.3 (4.3–11.3) or 10.4 (8.6–11.3) per pass (estimate; before red noise and the pulsar's own ±50 % per pass). SNR 7 per pass needs an SEFD of at most 1.05·10⁵ Jy (0.026 m²/K). The Sun crosses neither fixed beam at 52° N, so R12 changes. | 4 | 3 | 12 | partly: a lower SNR | existing: E4; **21 Nov:** a two-board bench test before buying four or eight (R8); **23 Nov:** measure both beams on the Sun (R12); R8 dropped; **24 Nov:** R12 changed (the panel has to be tilted to the Sun), R16 |
| FM-2 | sky, receiver | The Sun and the sky raise Tsys. The quiet Sun at 1.4 GHz is about 4·10⁵ Jy at solar minimum (a disk of about 10⁵ K), close to the nominal SEFD of 3.2·10⁵ Jy. Near the 2024–25 maximum it is higher, perhaps twice that (estimate). In the ~15° beam it cuts the SNR to 0.3–0.45 of its value. Through the sidelobes it adds 1.3–2.5 times the sidelobe level to Tsys: a few % at −20 dB, up to about 25 % at −10 dB. The array's sidelobe level is not known. Flares add bursts of minutes. The noise is measured per sub-int, so the error bars stay honest, but the forecasts, with a constant Tsys, are optimistic. B0950+08 lies 4.6° from the ecliptic, so the Sun is inside its half-beam for about 12 days around 20 August. **21 Nov (NOAA daily values at 1415 MHz, 1979–2026):** the yearly medians are 45–55 sfu at solar minimum and 110–150 sfu at maximum (monthly medians up to 185 sfu), so 4·10⁵ Jy is 10–30 % low and the maximum is 2.1–3.3 times the minimum, not "perhaps twice". With the Sun at beam centre the SNR factor is 0.36–0.41 at minimum and 0.17–0.22 at maximum. The first sidelobe of a uniform 4×4 array is −11.3 dB: about 10 K from the Sun there (estimate). Bursts above 500 sfu take at most 0.9 % of daytime in maximum years. The sky in a 15° beam is 3.2–6.3 K for every pulsar of the set except B1749−28 (about 11 K). **23 Nov:** with 3 m² the Sun is stronger against the noise: 500–1600 K in a main beam (0.5–1.5·10⁶ Jy against an SEFD of 6.4·10⁴ Jy), 25–80 K in a first sidelobe at −13 dB, a few K in the far sidelobes (my arithmetic). A zenith-facing panel at 52° N never has the Sun in a beam (the Sun stays south of declination +23.4°, the beams cover 32–72°), but the Sun crosses sidelobes every day, and a mount that follows a low-declination pulsar can meet it. | 5 | 2 | 10 | partly | **new: development** (Tsys(t) in the simulator) |
| FM-21 | 5, 5b | The last dispersion sweep of each file is unusable: 9–16 % of a 0.2 s file (D1-24, D2-14). Streaming must carry the overlap across file boundaries. | 5 | 2 | 10 | yes | existing: F |
| FM-24 | 7a | Real data processed with the analytic orbit, which is the default. Its Earth is up to 3.2·10⁵ km off, and B0329+54's Roemer delay up to 0.66 s. No pulsar is detected, and nothing says why. | 2 | 5 | 10 | no | **new: development** (a guard) |
| FM-28 | 9 TOAs | The good-TOA rule asks for detection at unknown phase (4.15) even when the predictor knows the phase (3.09), so it discards TOAs it could keep (D1-7). | 5 | 2 | 10 | yes | existing: D1-7 (review day 3) |
| FM-31 | 10 navigation | Earth orientation: UT1 − UTC is one number (up to 0.9 s, about 260 m east), polar motion is about 10 m, and the rows use the rotation angle only (below 6 km on the offsets). **21 Nov:** the 0.9 s bound ends. By Resolution 5 of the 28th CGPM (October 2026) UTC is continuous from 20 May 2027, with UT1 − UTC allowed to grow to 3600 s either way. `delayModel`'s default (UT1 = UTC) is 37 m east–west today (UT1 − UTC = −0.129 s) and grows without a ceiling. IERS Bulletin A predicts UT1 to 1.4 ms at 10 days and 4 ms at 40 days (about 1 m), so a file carried offline suffices. Polar motion is about 40 ns and the truncated nutation up to 50 ns. | 5 | 2 | 10 | no | existing: parked at block D (IERS data for F); **21 Nov:** UT1 from a Bulletin A file (R3) |
| FM-41 | 8 TOA | *(3 Nov)* The bandpass against equal channel amplitudes: `foldToTOA` leaves `ChannelGain` at 1, so the 'optimal' fit expects the same absolute pulse in every channel and weights low-gain channels most. SNR ×0.98 / 0.93 / 0.84 for a 3 / 6 / 10 dB tilt; error bars stay honest. | 5 | 2 | 10 | no | **new: development (small)**: `ChannelGain` from the channel mean powers; test written (`testFoldFaults` 1) |
| FM-7 | sky | Scintillation fades last longer than modelled: Gaussian ACFs, no refractive variation (S2-5). Low-DM pulsars such as B0950+08 get long gaps. **21 Nov:** the missing refractive variation is measured for B0329+54 itself, which carries the position: a modulation index of 0.31–0.47 at 1540 MHz, on a timescale of 8 h (one doubtful measurement) to 2–3.5 days (three scalings); FM-1 has the flux percentiles. The diffractive values (7 MHz, 13 min) sit in the middle of five studies (2.5–14 MHz, 6–17 min). B0950+08 is in weak scintillation above 100–300 MHz, so its long gaps have no confirmed mechanism, though its published fluxes scatter by 30–40 %. J0437−4715, B1929+10 and B1133+16 fade across the whole 400 MHz at once. **23 Nov:** with about 150 MHz clean the band holds about 5 scintles of B0329+54 instead of 12 (estimate), so the flux of a single 40-minute pass varies more. The single-pass goal is the first to feel it (FM-1). | 5 | 2 | 10 | yes: gaps | existing: S2-5; **21 Nov:** a slow refractive gain per pulsar in the simulator (R10) |
| FM-48 | 7 fold, 9 TOAs, 10 navigation, simulator | *(23 Nov)* A drift scan is not what the pipeline assumes. Through a fixed beam the pulsar's amplitude follows the beam gain: it changes by a factor 2 within the 0.72 h of a pass. Merging sub-ints adds them with their data weights only (`mergeSubints`), and `navFit`'s rows weigh a TOA's sub-ints equally ("equal length, constant flux", `navFit.m:15`). A TOA formed over a pass then loses 2 % of its SNR within the half-power points and 15 % over twice that span (my arithmetic), and its effective time is not the mean time when the pass is cut off-centre. The fast simulator has no gain that changes within a pass and no one-pass-a-day schedule, so no budget for the new array can be simulated yet. **24 Nov:** no transit instrument read (Molonglo, CHIME, Pushchino) weights sub-ints by the beam when forming a TOA; they fold the half-power window with a predictor. For a Gaussian beam an unweighted fold between the half-power points loses 6.6 % of the SNR against the weighted full pass, and 15 % over twice that window (estimate); a pass is worth 0.753 of the half-power crossing time at full gain, which is the 0.55 h used here. No paper defines the reference epoch of a transit TOA, so ours has to be defined with R1. | 5 | 2 | 10 | no | **new: development** (a beam gain against time in `simulateFold` and the scenario runs; sub-ints weighted by the expected gain; R1) |
| FM-49 | E4 array, mount | *(24 Nov)* The panel is mounted a few degrees off its intended tilt or rotation. The beam is 6.3° wide, so the gain falls to 0.93 at 1°, 0.76 at 2°, 0.62 at 2.6° and 0.53 at 3° off (Gaussian beam, estimate); 2.6° is simulated: SNR 11.5 → 7.2, with 57 % of the passes at SNR 7 or more. BURSTT found its array 3.0° off in azimuth after deployment, and EMBRACE's rotation was first known only to a few degrees. Two requirements meet here: a tilt of 2.6° north puts B0329+54 on the first beam's centre, and the second beam lies on the pulsar's track only if the panel's diagonal points about 69° east or west of north (with the diagonal at 45° the pulsar misses that beam by about 8°). A cross-track error does not move the time of the peak, it only lowers the SNR, so the pulsar cannot tell it from a weak array (FM-14). | 3 | 3 | 9 | partly: a lower SNR; an along-track error shifts the time of the peak (0.33° per pass at SNR 11.5, estimate) | **new: check** (survey the tilt and the rotation as mounted, to 1° or better; R20) |
| FM-50 | 8 detection, first light | *(24 Nov)* A first detection that is not one, or a real one that cannot be shown to be one. The thresholds (3.09 and 4.5, false-alarm probability 10⁻³) hold per profile and in the simulator's white noise only: the 16 sub-ints of a pass plus the total already give 1.7 % per pass, two beams, a period search (FM-10) and several ephemerides multiply it, and gain noise raises it much further (FM-42: 0.36 per sub-int with a knee at 10 Hz). For B0329+54 in 1300–1450 MHz the usual test against DM is weak: the sweep is 12.9 ms against a 6.6 ms pulse, and a fold at DM 0 keeps about half the SNR (estimate). No second pulsar and no single pulse is within reach (0.22 per pulse, estimate). The published false detections (a guiding camera's noise, an error in the Earth's orbit, a microwave oven, a mispointed band at CHIME in 2026) fell to null tests and repetition, not to a higher threshold. | 3 | 3 | 9 | no, by the pipeline; yes with folds of off-source hours (about 43 slices a day), the beam without the pulsar, and the return 3 min 56 s earlier each day | **new: test** (the confirmation set for first light, R18), then development (an empirical null, after Jasper's OK) |
| FM-33 | simulator | The fast simulator and the voltage chain drift apart unnoticed: their agreement rests on experiments without PASS/FAIL (T2-1). Every E3 forecast rests on the simulator. | 3 | 3 | 9 | no | **new: test** |
| FM-37 | F port | The GPU/FPGA port differs from `processIQ` (precision, fixed point). **21 Nov:** the GMRT pipeline needed double precision for its long unchannelized transform; per-channel 4096-point transforms are in the regime where others run single precision. | 3 | 3 | 9 | partly | existing: F (checked stage by stage against `processIQ`) |
| FM-46 | 4 IQ, receiver | *(21 Nov)* The ADC's calibration and the sampling plan. (a) The gain and time-skew background calibrations of the RFSoC ADC need an input of at least −40 dBFS (AMD PG269). The planned receiver noise, about 50 dB below full scale, is under that, so either the calibration is frozen (its interleaving spurs then grow with temperature) or the gain goes up by about 10 dB (less headroom, FM-12). (b) At 5.0 GS/s the band that folds onto 1.2–1.6 GHz is 3.4–3.8 GHz, the European 5G band (at the PYNQ default of 4.9152 GS/s: 3.32–3.72 GHz), and the 4x2's input balun passes it. A carrier there covers 5–25 % of the band and cannot be flagged away. A rate near 4.0–4.1 GS/s folds 2.4–2.9 GHz instead; a ZCU216 near 2 GS/s has five alias bands. About 70 dB of rejection is needed, against 60 dB in the front-end study. | 3 | 3 | 9 | partly: lines, or a raised noise floor in part of the band | **new: decision** (sampling rate, alias filter, gain mode; R7) |
| FM-16 | 5a–6b | Lopsided channels (narrow carriers, LTE edges, a carrier across a channel edge) break the flat-channel noise model: the error bars are 0.4–3.4 % too small, and `noiseCal` is off by default (D1-29). | 4 | 2 | 8 | partly: χ² | existing: D1-29 (decision, day 5) |
| FM-25 | 4a, 7a | Visibility gaps: for a millisecond pulsar, the clock's wander over a gap must stay below half a period (J0437: 2.9 ms), or the turn count is lost. B0329+54 has 357 ms of margin. **21 Nov:** across the daily gap the site term must also stay within a fraction of a turn: the site known to about 75 km for a 5 ms pulsar and 23 km for B1937+21, and the clock rate to 5.8·10⁻⁹ and 1.8·10⁻⁹ (0.1 turn; researcher's derivation). The ephemeris consumes the margin as well as the clock: a TCB frequency read as TDB (FM-44) breaks B1937+21's turn count within a day. **23 Nov:** with a fixed panel every pulsar has a 23-hour gap every day. | 2 | 4 | 8 | partly | none now (E1's choice of pulsars, E5) |
| FM-27 | 9, 10 | At SNR ~7 the error bars are about 5–10 % small, and 1.9 % of B0329+54's TOAs lie beyond 3σ. `navFit`'s mixture fixes the outlier fraction at 2 % and the scale at 4×; real data may hold more outliers. **21 Nov:** consistent with the literature, which puts the breakdown of the Gaussian error at S ≈ 8–10 and cuts TOAs below 8. A researcher's white-noise check with the five-component template gives 7 % small and 1.4 % beyond 3σ at S = 7. | 4 | 2 | 8 | partly: `chi2Core` | **new: development** (low) |
| FM-34 | simulator | Pulls of 1 show that the simulator agrees with the observer's model, not that the model is real (T2-3). Only real noise tests it. | 4 | 2 | 8 | no | **new: test** (the first F test) |
| FM-40 | 4 IQ | *(3 Nov)* Quantization by the ADC or an 8-bit transport: the noise model stays honest, the SNR drops by 1/(1+ε): ×0.995 at 8 bits, rms 4 LSB; ×0.93 at rms 1 LSB (room for a +40 dB radar); ×0.99 at 4 bits, rms 3 LSB. **21 Nov:** BURSTT-256 attributes a lower than expected array sensitivity to its 4-bit requantization: a caution against early requantization at weak signals. | 4 | 2 | 8 | no | **new: decision** (sc16, or sc8 at ≥ 2–4 LSB with a clip flag); test written (`testFrontEndFaults` 2) |
| FM-12 | 3–4 receiver | A strong nearby radar clips the ADC or compresses the LNA. The clipping is strong and gets blanked, but the gain dip after a pulse is a power drop, which `detectRFI` cannot see: it flags high power only (D1-11). **21 Nov:** with about 10.5 effective bits and the noise 40 dB below full scale, a signal up to about +40 dB (against the noise in 400 MHz) passes unclipped (estimate). Radar 1 of the model sits on that limit; a line-of-sight radar, an in-band DME or a mobile mast within 1–10 km exceed it. The RFSoC's Gen 3 ADC then inserts its own attenuator (15 dB) and holds it until the logic clears a flag, so the gain is lower for an unknown time after a pulse. Over-voltage and over-range flags are available to the firmware (the clip flag). | 4 | 2 | 8 | no | **new: test** (E4/F); 3 Nov: a clipped pulse also suppresses the noise in every channel, unflagged outside the frequency guard → a clip flag from the SDR (log 2026-11-03) |
| FM-38 | F | The real-time processing of 18 streams cannot keep up, and data are dropped. **21 Nov:** 16 or 32 streams since 20 Nov. Published backends run 400 MHz of dual-polarization coherent dedispersion and folding on one graphics card per beam (CHIME/Pulsar, MeerKAT); channelizing and beamforming are FPGA work in every system read. **23 Nov:** the new array delivers two beams of two polarizations, already formed in the analogue network: no beamformer to compute, and four streams instead of 16 or 32. | 2 | 3 | 6 | partly | existing: F (GPU/FPGA); FM-9's ingest makes the drops visible |
| FM-17 | excision | The frequency guard has no cap: dense, independent strong bursts blank up to 99 % (D2-22). Real radars (duty ~0.1 %) cost far less, so this is a stress case. **21 Nov:** "duty ~0.1 %" holds for the old 2 µs radars only; solid-state radars run at 4–10 % and more, with pulses of 400–1700 channel samples against the detector's longest window of 16. The 10-sample guard (2.4 µs) is shorter than one DME pulse; Green Bank blanks 12 µs per DME pulse and loses 0.5–2 % of the data. | 2 | 3 | 6 | partly: `blankedFraction` is reported, nothing warns | **new: development** (a warning); D2-22 stays |
| FM-26 | 9, 10 | A false TOA: noise or RFI passes the detection threshold and gives a TOA at a random phase. Noise alone rarely does: `adaptiveSubints` closes a group at T0 ≥ 7 at the predicted phase, and only the last group of each pass, below the target, is decided at P_FA 10⁻³. That is about one test per pulsar per pass, so about 0.1 false TOAs in 30 days, and a pulsar with no other TOAs absorbs one in its constant c_k. RFI residue is the likelier source. `navFit` puts a false TOA in its 4× wider component, but a residual of up to half a period (about 500 σ for B0329+54) still pulls with 1/16 weight. One among ~300 shifts the position by about 1–2 of its σ, depending on its phase. The detection tails below 10⁻³ are extrapolated (limits page). **21 Nov:** a search over site and clock rate adds a trials factor that `detectPulsar` does not have: about 10⁵ cells, so a threshold near 6.3σ for a total false-alarm probability of 10⁻³ (researcher's estimate). | 2 | 3 | 6 | partly: flagged as an outlier, still pulls | **new: test, then development** |
| FM-35 | drivers | `main.m` and the drivers are not run by `runAllTests` (bug #12 went unseen for two days), and scripts share the base workspace (D2-6). | 3 | 2 | 6 | partly | existing: D1-6, D2-6 (review days 4–5) |
| FM-36 | run control | A stale file is reused: `checkFold` does not check `noiseCalOnPulse` (D1-5), and RFI profile files do not record their excision settings (P3-4). | 2 | 3 | 6 | no | existing (minor) |
| FM-47 | 4 IQ, receiver | *(23 Nov)* The two beams are assumed to come from a power splitter and to be live all the time (Jasper: "i think assuming a decent power splitter is okay, beams are on continuously without switching"). He saw "SPDT", a switch, in a diagram; what it switches is not known. If the beams alternate instead, each gets half the time: SNR × 0.71 per beam, twice the integration time, a single pass of B0329+54 at SNR 8 instead of 11, and sidebands of strong RFI at the switching rate inside the clean 150 MHz. | 2 | 3 | 6 | yes, if looked for: periodic gaps in the raw voltage, lines at the switching rate | **new: check** (ask the team; look at a first recording; R13) |
| FM-3 | sky, 8 TOAs | The profile changes, but the template is fixed. B0329+54 is in its abnormal mode about 15 % of the time. The profile also evolves across 1.2–1.6 GHz, and scintillation changes which frequencies dominate. Each gives a TOA offset that lasts hours; the size needs the mode profiles. **21 Nov:** the occupancy holds (13–18 %, measured), but the abnormal spells average 26 min, not hours (a third are under 10 min, the longest seen is 2.5 h; Chen et al. 2011). From the 1400 MHz component table of Brinkman, Mitra & Rankin 2019, a fit of the normal template to the abnormal profile comes out 0.15 ± 0.04 ms early (researcher's fit; the table's signs are to be checked against the PDF): about 35 µs of scatter and 25 µs of bias per 2 h TOA, under 5 % of σ_TOA. 0.8 ms stays as the pessimistic bound (a core shift of 0.45°, which the tables neither show nor exclude; +3–5 % on σ_TOA). The abnormal mode is about 1.7 times brighter, so it carries about 23 % of the signal. Profile evolution across 1.2–1.6 GHz is 0.1–0.2 ms (estimate), so one template suffices. | 5 | 1 | 5 | no: at SNR ~7, χ² cannot see it | **new: test** (with E1's profiles) |
| FM-6 | propagation | DM error and DM variations (2.1 µs per 10⁻³), the Doppler-shifted dispersion (±7 µs annual for B0329+54), the ionosphere (3–70 ns), the troposphere and the solar wind. All are small against B0329+54's ~700 µs TOAs; they grow for millisecond pulsars. **21 Nov:** the numbers are confirmed (2.12 µs per 10⁻³ pc cm⁻³, 0.69 ns per TECU, troposphere 7.7–22 ns). The Doppler-shifted dispersion is ±9.3 µs at 1.4 GHz for B0329+54; ±7 µs is the 1.6 GHz value. Additions: Vela's DM drifts by 10–17 µs per year, the solar wind reaches about 3 µs for B0950+08 at 4.6° from the Sun, and Vela's scattering is about 23 µs at 1.4 GHz (estimate). | 5 | 1 | 5 | no | existing: E2 |
| FM-13 | 4 IQ | A real zero-IF receiver adds I/Q imbalance, a DC spike and LO leakage at the in-band LO (1.3 GHz). They act as extra carriers, which the weighting and the excision handle (D1-10, R-I3). **21 Nov:** an interleaved RFSoC ADC adds an offset spur at multiples of Fs/8: 1500 MHz at 4.0 GS/s, 1250 MHz at 5.0 GS/s, about −86 dBFS typical. It is a steady line that costs one channel. | 5 | 1 | 5 | partly | existing: F calibration; DC tested 3 Nov (`testFrontEndFaults` 2: the LO sits on a channel edge, where the dedispersion taper removes it) |
| FM-39 | numerics | A time in seconds since J2000, held in one double, resolves only about 0.12 µs in 2026. The predictor works around it (phase from the span's start); wherever absolute times are still one double, ns timing would need a two-part time. **21 Nov:** confirmed: 0.119 µs today, 0.238 µs from January 2034. | 5 | 1 | 5 | no | none (a note for millisecond pulsars); tested 3 Nov (`testFrontEndFaults` 4) |
| FM-15 | E4, sky | The pulsar's components are polarized, with different angles. With unequal X and Y gains, or one polarization, the received profile changes with the parallactic angle and with Faraday rotation. The TOA then moves on the sidereal day, which is the position's own signature, and time does not average it away. The size is unknown. **21 Nov:** the size is now estimated. For a 10 % differential power gain, a template fit on published Stokes profiles gives 11–13 µs band-averaged for B0329+54 (23 µs per channel), 23–28 µs for B0950+08, 31–50 µs for B1929+10, 1.65 µs for J0437−4715 and 0.4 µs for B1937+21 (researcher's calculation; the millisecond-pulsar values agree with van Straten 2006 within a factor 2–2.6). If all of it fell in the first sidereal harmonic, that is 6–12 km for B0329+54. Ideal crossed dipoles already differ by ±20 % in projected gain at 35° off axis. A second, cruder scale estimate gives five to ten times more; the two were not reconciled. **23 Nov:** smaller with the new array. The boresight beam has equal gains and orthogonal polarizations. In the second beam, tilted 20° along the grid's diagonal, the gains stay equal and the polarizations are 93.6° apart (ideal crossed dipoles; my geometry): about 6 % of the linear flux leaks into the summed power, a few µs for B0329+54. With a fixed panel a pulsar is always seen at the same parallactic angle, so the term is a constant per pulsar and beam. It returns with an alt-azimuth mount, where the parallactic angle turns; a polar mount keeps it constant. | 2 | 2 | 4 | no | **new: test** (E4); **21 Nov:** can run now on the EPN Stokes profiles (R9) |
| FM-32 | F time | Absolute time: UTC → TT and leap seconds. A constant error is just the clock offset; a step (a leap second, a restart) is not. **21 Nov:** no leap second can occur after 20 May 2027 (FM-31), so that step ends; restarts remain (FM-9). | 1 | 4 | 4 | no | existing: D1-15 (F); FM-9's ingest |
| FM-23 | 7 fold | A long fold needs a position prior: below 4 km for J0437, about 100 km for normal pulsars. At the fixed site the prior is known. **21 Nov:** since 2 Nov the site is not known beforehand, so this applies at acquisition. A linearized fit holds only once every block is within about 5.5 ms (1650 km of path) for B0329+54, so a global search comes first (D6c). | 1 | 4 | 4 | partly | none (the ship, parked) |
| FM-19 | excision | RFI locked to the pulsar's period leaves a 2 µs bias at −54 dB with the guard alone (0.41 µs with the periodic mask). That is small against B0329+54's TOAs. | 1 | 2 | 2 | no | none (the periodic-mask fix stays parked) |

### New actions, put on the to-do list (notes §10 item 16, roadmap "Open items")

In order of priority. None is built; each follows the usual rule (physics and plan first, then Jasper's OK).

| # | Risk | Kind | What | Fits in |
|---|---|---|---|---|
| 1 | FM-11 | test | `makeClock` with a daily rate term at the solar day (10⁻⁹ to 10⁻⁷); `expNavFit` at 30 and 100 days; the position bias against the amplitude. Simulator only. Then decide on the oscillator (OCXO) and on temperature as a regressor. For F, log the SDR's board temperature with the data either way. | before E5 |
| 2 | FM-8 | development + test | A switch in `processIQ` for an inverted band, and a test that an inverted IQ file with the switch gives the same TOAs. For F, a check of the mapping, band edges and sample rate with a known line (GPS L1 at 1575.42 MHz, or a signal generator) before the first observation. | before F |
| 3 | FM-9 | development + test | The F ingest: read the SDR's time stamps and overflow flags, fill a gap with zeros and blank it (W = 0), keep absolute time. The test can run now: clock steps, and many small ones, in `trueResiduals`, through `expNavFit`, with and without step states. | F (test in E) |
| 4 | FM-10 | test, then development | `runScenario` with an uncalibrated SDR clock (rate prior ±2 ppm). If no TOAs come: a rate search at acquisition (folds at trial rates about 10⁻⁷ apart, so about 40 trials over ±2 ppm), with the trials counted in the false-alarm budget. | E3, E5 |
| 5 | FM-33 | test | `expSimVsChain` and `expSimRFI` as PASS/FAIL tests in `runAllTests`, with tolerances from the logged numbers. | review day 5, or E3 |
| 6 | FM-3 | test | B0329+54's two modes in `simulateFold` (the abnormal one 15 % of the time, in spells of hours); the TOA offset and the position bias with one template. | E1 (profiles), then E3 |
| 7 | FM-4 | test | One pulsar with a spin-rate error, and one with a glitch, in a multi-pulsar `navFit`: does the position move, or does a residual trend show? | E2, E3 |
| 8 | FM-18 | development | An RFI survey at the site with the real hardware: spectra and time series over a day. It extends or replaces `rfiRealistic`, and settles D2-22, FM-12 and FM-17. | F, early (needs no code) |
| 9 | FM-2 | development | Tsys(t) in the simulator: the Sun, in the main beam and the sidelobes, and the sky temperature toward each pulsar. First a number per pulsar and month. | E3 |
| 10 | FM-24 | development (small) | `processIQ` warns, or refuses, when real data (a flag in `info_IQ`) come without `'Ephemeris'`. | before F |
| 11 | FM-17 | development (small) | A warning in `processIQ` when the blanked fraction is above a limit (for example 5 %). | any time |
| 12 | FM-26 | test, then development | Inject TOAs at random phases (0.1–1 %) into `testNavFit` and `expNavFit`, and measure the bias. A long noise-only run with `simulateFold` checks the detection tail at 10⁻⁴ to 10⁻⁶. If the bias matters: a uniform-in-phase outlier component in `navFit`, or a gate on the residual. | E3 |
| 13 | FM-12 | test | A receiver stage with ADC bits, clipping and LNA compression (D1-11), with a strong radar through the excision. | E4, F |
| 14 | FM-15 | test | A polarized profile with parallactic rotation and an X/Y gain imbalance in the simulator; the position bias. | E4 |
| 15 | FM-34 | test | Pulsar-free real recordings through `processIQ`: the noise model's χ² and the false-alarm rate on real noise. | F, the first test |
| 16 | FM-27 | development (low) | `navFit` estimates an error scale per pulsar (an EFAC, as pulsar timing does) and the outlier fraction from the data. | E3 |

### What the draft checked

- `navFit`'s outlier model: a fixed mixture of N(0, σ²) and N(0, (4σ)²) with 2 %, iterated by EM (`navFit.m:51–68`). There
  is no gate and no uniform component.
- Time precision: `makePredictor` takes the exact difference from `TRef` first and counts phase from the span's start
  (`makePredictor.m:47, 83`).
- The frequency mapping: `applyIQmodulation` and `channelizeIQ` assume RF = fLO + f_bb, no inversion; there is no switch.
- The clock: the runs use a port calibration to 1 ms and 10⁻⁹ (`receiver.html` 4a.3); nothing searches over clock rates.
- The SEFD design points: 3.2·10⁵ Jy nominal, 7.4·10⁵ Jy worst case, Tsys 70–80 K with ~3–5 K of sky (`scenario.html`).
  The Sun is not in Tsys.

### The estimates, checked

Jasper chose to review the ratings before the commit. Before that review, the four estimates of the first draft were
checked against the code, the docs and a short calculation. Two ratings changed and two held.

| Risk | First draft | Check | Rating |
|---|---|---|---|
| FM-11 | a rate swing of 10⁻⁸, 140 µs, "about 40 km" | 10⁻⁸ is the low end for a ±2 ppm TCXO (10⁻⁸ to 10⁻⁷ per K, estimate). On B0329+54, 140 µs is 72 km, not 40 (40 km holds for a pulsar on the equator). The fit's equatorial error is ±30 km (`navigation.html`), so the bias is 2–25 σ. `navFit`'s clock covariance (random walk at the crystal's level) does not absorb a steady daily line. | 4 × 4 = 16 → **5 × 5 = 25** |
| FM-26 | P_FA 10⁻³ per test, so false TOAs are possible | `adaptiveSubints` closes groups at T0 ≥ 7 at the predicted phase. Only the last group of a pass is decided at 10⁻³, so a 30-day run expects about 0.1 false TOAs from noise. The bias per false TOA holds: 1–2 σ, from (r/16σ)·√(p/N). | 3 × 3 = 9 → **2 × 3 = 6**; its action moves from 5 to 12 |
| FM-2 | quiet Sun 4·10⁵ Jy; main beam "about halves"; sidelobes 5–15 %; B0950+08 two weeks in August | 4·10⁵ Jy holds at solar minimum: 2kTΩ/λ² with T 10⁵ K over the disk gives 41 sfu. It is higher near maximum. Main beam: SNR × 0.3–0.45. Sidelobes: a few % to 25 %, by the unknown sidelobe level. B0950+08 lies at ecliptic latitude −4.6°, so the Sun is inside its half-beam for ~12 days around 20 August. | 5 × 2 = 10, unchanged |
| FM-18 | a self-made comb could cost tens of % | Channels are 3.125 MHz, 128 in the band. A 10 MHz comb hits 40 channels (31 %), a 100 MHz one 4 (3 %). | 4 × 3 = 12, unchanged |

Still estimates: the TCXO's temperature slope and the room's daily swing (the SDR's datasheet, or a day's measurement,
decides FM-11); the solar flux near maximum; the array's sidelobe level.

### Tests for the risks

Jasper: "you can write tests for all those cases you think could lead to failures". The project's rules still hold: the
physics and the plan of each test first, then his OK, and one new function at a time. Six of the 16 actions can get a
test now.

| Risk | Test | What it needs |
|---|---|---|
| FM-11 | the daily clock term against the fit's bias | nothing: the fault goes in through `trueResiduals`' `'Clock'` |
| FM-9 | a clock step, and many small ones | nothing (the same route) |
| FM-4 | a spin-rate error and a glitch in one of two pulsars | nothing (that pulsar's residuals) |
| FM-26 | residuals replaced by random phases | nothing |
| FM-33 | `expSimVsChain` and `expSimRFI` as PASS/FAIL tests | the saved chain sweeps; tolerances from the logged numbers |
| FM-8 | a conjugated (inverted) IQ file through `processIQ` | the test alone; the switch is a development |

The first four become cases of one new file, `tests/testNavFitFaults.m`, one case at a time. The others wait for a development
(FM-10: `simulateFold` has no clock inside a sub-int; FM-12, FM-15, FM-2, FM-24, FM-17, FM-27) or for inputs we do not have
(FM-3: E1's mode profiles; FM-18: a site survey; FM-34: real recordings).

**Case 1, FM-11 (written; Claude's run, 2.6 s; not yet validated).** Jasper OK'd the design ("Yes, write it"). The clock error
*A* sin(2π*t*/86400 s + φ) goes into `trueResiduals` with the true site equal to the assumed one; `navFit` "self", no
clock model, noise-free. The prediction *c A* sinc(Δω*T*/2) / cos δ (70.3 km at 30 days, 62.7 km at 100, for *a* = 10⁻⁸)
holds to 0.1 % at three phases and both spans (tolerance 2 %). Shown: the bias is 0.23, 2.3 and 22.6 σ at 30 days for
*a* = 10⁻⁹, 10⁻⁸ and 10⁻⁷, and 0.37, 3.7 and 36.8 σ at 100 days. One noisy run at 10⁻⁷: at 30 days a 659 km bias with
core χ² 1.120 (1.131 for the same noise alone) and no TOA flagged; at 100 days χ² 1.381 (1.002 alone). The crystal clock
model (`'Clock'`) leaves the bias unchanged at both spans. `runAllTests` runs the new file at the end, with its note for a
test not in its list. Docs: `tests.html` (table row, section, Figure 2 regenerated by `gen_testmatrix.py`, 32 tests).

**Navigation tests on hold.** Jasper: "let's hold off on the navigation tests, because we'll redesign that part next week,
I have some changes in the project final deliverable". Cases 2–4 (FM-9, FM-4, FM-26) are not written. Case 1 stays as it
is, uncommitted; it may change with the redesign. Case 2 was designed before the hold, and its physics is kept here for
after the redesign:

- A dropped buffer of *n* samples makes every later time stamp *n*/*f*s early, a clock step of −*n*/*f*s (41 µs for 16 k
  samples at 400 MS/s; the buffer size is an estimate). A restart is a step of any size; only its value modulo the period
  counts.
- One step Δ biases the position little: after the constant c_k takes its mean, its projection onto the sidereal sine is at
  most (2Δ/*N*) / sin(ωΔ*t*/2), 2.2 % of Δ over 30 days (1 ms gives ≤ 21.5 µs, 11 km, 0.36 σ).
- But it shows in χ²: +(Δ/σ)² *f*(1 − *f*), with *f* the fraction after the step (+0.46 for 1 ms at *f* = ½, against a
  scatter of 0.08). A single step large enough to matter is therefore seen; the register's "Seen: no" for FM-9 is too strict.
- Many small drops are the silent part: drops of δ at a Poisson rate λ give the clock a rate −δλ and white frequency noise
  h₀ = 2δ²λ (one 41 µs drop per hour: 1.1·10⁻⁸ and 200 µs per day, against the crystal's 76 µs). A fit that assumes the
  crystal only then states too small an error.
- A known restart is handled today by a second constant: the pulsar listed twice, `toa.k` switching at the restart.
- The planned case: (a) a 1 ms step at five random times, noise-free, against an ideal least-squares fit of [1, cos ω*t*,
  sin ω*t*] (within 1 % of the 21.5 µs maximum); (b) a 50 ms restart with the second constant, bias below 1 m; shown
  (c) steps of 0.3, 1 and 10 ms with noise: bias, χ² against 1 + (Δ/σ)²/4, flags, with and without the crystal model;
  shown (d) 41 µs drops at one per hour, 100 realizations: the crystal model with and without h₀ = 2δ²λ, pull spread.

### The new objective

Jasper, after the commit, while still working through the math of RFI excision: "note that the new objective is to
estimate the location of a static array somewhere on earth, locking onto multiple pulsars using digital beamforming, and
estimating the location of the antenna in this static setting (but the Earth moves, of course)". The test site stays
the Netherlands, probably (Jasper: "the test will probably be done in the netherlands, yes"). The dead-reckoning
navigation filter (`navFilter`, `runShip`), parked since 23 Oct, stays parked (Jasper: "but this means we can later
park the DR nav filter"). The ship is the next goal after the static field test: "let's not develop for that now,
as we first want static positioning as a field test. Do keep all of the documents for it though, it is the next goal in
the succession of this project." Recorded at the top of
notes §11, with three points for the redesign: a pulsar set that depends on the latitude, a large position prior at
acquisition (FM-23 rises), and the polar component, which needs the pulsars' absolute phases (FM-29). The docs still give
the 23 Oct endpoint (a fixed site in the Netherlands: `index.html`, `roadmap.html` block E, `scenario.html`); they change
with the redesign.

### 21 Nov 2026: the register against the literature

Jasper: "update the risk register with the findings". The findings are those of the literature research of 21 Nov
(`2026-11-21_deep-research.md`; the report and the notes behind it are in `reports/` and `research_notes/`, ignored by
git, on this Mac only). Each affected row above has a paragraph marked **21 Nov**. Every number in those paragraphs is
from a paper, a data sheet, a catalogue or a researcher's own estimate, as marked. None was simulated in the pipeline,
and none is a project result. The ratings stay a draft for Jasper's review.

**Ratings changed (six).**

| Risk | Was | Now | Why |
|---|---|---|---|
| FM-4 | 4 × 3 = 12 | 5 × 3 = 15 | No glitch of B0329+54 is catalogued and its timing noise is negligible inside a campaign, but the catalogue's spin parameters are certainly stale (epoch 1986 for B0329+54), by 10–150 µs per day, and `navFit` has no rate per pulsar. |
| FM-18 | 4 × 3 = 12 | 5 × 3 = 15 | The model is now known to differ: DME in the band, long chirped radar pulses, misplaced mobile blocks, five missing emitter families. |
| FM-7 | 3 × 2 = 6 | 5 × 2 = 10 | The refractive variation is measured for B0329+54 itself (index 0.31–0.47), not only expected for low-DM pulsars. |
| FM-12 | 3 × 2 = 6 | 4 × 2 = 8 | Headroom is about +40 dB; a mast within 10 km, an in-band DME or a line-of-sight radar exceed it, and the ADC then holds an attenuator until it is cleared. |
| FM-3 | 5 × 2 = 10 | 5 × 1 = 5 | The mode offset is 0.15 ms and the spells last 26 min on average: under 5 % of σ_TOA. At the pessimistic bound of 0.8 ms the impact stays 2. |
| FM-32 | 2 × 4 = 8 | 1 × 4 = 4 | Leap seconds end on 20 May 2027. |

**Rating unchanged, content changed:** FM-11 (the oscillator classes sized; the cost of fitting the line), FM-5 (37.9 µs
per year; about 557 µs today without proper motion; absorbed only per pulsar), FM-10 (the 4x2's crystal: 500 rate trials,
not 40), FM-1 (flux spread and the 5th percentiles), FM-2 (the Sun over the cycle; the sidelobe level), FM-15 (a size:
11–13 µs for B0329+54 at 10 % gain imbalance), FM-31 (UT1 unbounded from May 2027), FM-14 (calibration on the Sun; no
published restart repeatability), FM-29 (declination spread alone gives no polar component). Smaller notes: FM-6, FM-8,
FM-9, FM-13, FM-17, FM-22, FM-23, FM-25, FM-26, FM-27, FM-30, FM-37, FM-38, FM-39, FM-40, FM-42.

**Confirmed, no change:** the daily-sine sensitivity behind the 30 km east (Zhou et al. 2026 give cos δ / (c√2) per
equatorial component: 27.8 km with our TOAs); T_sys 60–80 K near broadside (measured ambient L-band arrays: 45–70 K at
half-wavelength pitch); zero-fill with zero weight for lost samples (FM-9); the form and sign of `delayModel`'s terms
(FM-30).

**New scenarios (four):** FM-43, the array's effective area against direction (20); FM-44, mixed TCB and TDB units in
the catalogue's spin parameters (16); FM-45, the board's reference not the one assumed, or not locked (12); FM-46, the
ADC's calibration threshold and the alias of the 3.4–3.8 GHz band (9).

**Where the 18 streams stood** (FM-9, FM-14, FM-38): the rows keep their text of 2 Nov, and their 21 Nov paragraphs give
the 4×4 array's 16 elements (32 streams if dual-polarized, not yet confirmed).

**New actions from the research.** None is built; each follows the usual rule (physics and plan first, then Jasper's OK).
They are not yet in notes §10 or in `roadmap.html`.

| # | Risk | Kind | What | Fits in |
|---|---|---|---|---|
| R1 | FM-43 | development, then a decision | The effective area against direction and frequency as a time-varying weight in the scenario runs (the fast simulator has one A_eff); rerun the position budget; then choose the panel: zenith, tilted, pointed, or a dense lattice with analogue sub-arrays. The real element's embedded pattern (a simulation or a measurement) replaces the researcher's model. | E3, before E4 |
| R2 | FM-4, FM-5, FM-9 | development | `navFit`: a phase and a rate per pulsar (a quadratic for a crystal-class clock), and a clock segment per restart. Test: the position error unchanged beyond three days; a stale frequency no longer biases the fit. | D6 / E3b |
| R3 | FM-44, FM-5, FM-31 | development | E2 inputs: one par file per pulsar from a single release, converted once to TDB; proper motion and position epoch in `delayModel`; UT1 from an IERS Bulletin A file. The PINT check in a Linux container (PINT does not run on native macOS ARM Python). | E1, E2 |
| R4 | FM-11 | decision, test | The oscillator: a single-oven 10 MHz OCXO on the reference input as the minimum, a rubidium for margin and for integrations beyond a day. Test in `testNavFitFaults`: the cost of a free solar-day line, for one pulsar and for two at different right ascensions. A day's bench measurement of the real slope. | E5 |
| R5 | FM-45, FM-10 | development (small) | Log the clock chip's lock status and the board temperature with the data; refuse or flag a run without lock. | F |
| R6 | FM-18, FM-17 | development | `rfiRealistic`: the licensed mobile blocks, DME pairs at 1200–1213 MHz, one long-pulse radar, the band-edge navigation signals. `detectRFI`: longer windows and a guard that covers a DME pulse. A desk check (antenna register, aeronautical publication), then a survey of at least two weeks with the real elements over 1.0–1.8 GHz. | F, early; the model before E3 |
| R7 | FM-46, FM-12 | decision | The sampling rate (near 4.0–4.1 GS/s), the alias filter (about 70 dB; the filter data sheets read by hand), the gain mode (noise near −38 dBFS, or train and freeze), overload flags and a PPS-latched sample counter in every packet. | Task 2.4, before F |
| R8 | FM-14 | test | Two boards on the bench: the phase between them after each restart, and its drift over hours. A phase calibration on the Sun or an injected noise source at the start of every session. | E4, before buying four or eight boards |
| R9 | FM-15 | test | Register action 14 on the published Stokes profiles (EPN): a gain imbalance with parallactic rotation, the term at once and twice per sidereal day. In E4: gain-equalize the two polarizations per direction, and keep the folded cross product. | E4 |
| R10 | FM-7, FM-1 | development | A slow refractive gain per pulsar in the simulator (a measured index, the timescale as a wide-range parameter). At first light: a few weeks of our own flux monitoring of B0329+54. | E3 |
| R11 | FM-22, FM-8, FM-9 | test | Public voltages through `processIQ` before the hardware exists: B1642−03 from Dwingeloo against its published fold, then one Green Bank coarse channel of B1133+16. The downloads (12 and 17 GB) need Jasper's go-ahead. | before F |

**Still estimates, to be replaced by a measurement or a first-hand reading:** the element's embedded pattern and the
panel's area against scan (FM-43); the phase between boards after a restart (FM-14); the refractive timescale of B0329+54
(FM-7); the signs in the mode-profile table (FM-3); the crystal's and the LNAs' data sheets (FM-10, the front-end study);
every RFI level (FM-18), which is a free-space estimate with assumed powers.

### 23 Nov 2026: the register for the 3 m² array

Jasper: "update the risk register for the new array". The hardware is the one he reported on 23 Nov
(`2026-11-23_array-3m2-analog-beams.md`; notes §11, "Array, update 23 Nov 2026"): about 3 m² of effective area, a 4×4
grid of tiles with 4×4 dual-polarized bowties each, analogue delay-and-sum hardwired to two beams, (0°, 0°) and
(20°, 45°), four analogue outputs, about 70 K measured at the output, about 150 MHz of the band clean. Each affected row
above has a paragraph marked **23 Nov**. Every number in them is a scaling of the simulated budget (1.93 h per SNR 7 TOA
of B0329+54 at 3.2·10⁵ Jy and 400 MHz) or my arithmetic, as marked. Nothing was simulated for the new array. The ratings
stay a draft for Jasper's review.

**In the team's terms:** sensitivity A_eff / T_sys = 0.043 m²/K (was 0.0087), SEFD = 6.4·10⁴ Jy (was 3.2·10⁵; the old
"worst case" was 7.4·10⁵). No worst case is defined for the new array yet.

**Ratings changed (three).**

| Risk | Was | Now | Why |
|---|---|---|---|
| FM-9 | 4 × 4 = 16 | 3 × 4 = 12 | Four streams on one board instead of 16 or 32: an overflow is less likely. What a drop does is unchanged. |
| FM-38 | 4 × 3 = 12 | 2 × 3 = 6 | Two beams, formed in the analogue network: no beamformer to compute. |
| FM-15 | 3 × 3 = 9 | 2 × 2 = 4 | Equal gains in both beams, 3.6° off orthogonal in the second; a constant per pulsar on a fixed panel. It returns with an alt-azimuth mount. |

**Rating unchanged, content changed:** FM-43 (20: the scan loss is gone; fixed beams give one hour angle a day and so no
position; the mount decides), FM-14 (nothing can be calibrated or re-pointed digitally; the beams must be measured as
built), FM-1 (the single-pass goal sets the margin: SNR 11, 7 or 5 per pass), FM-18 (150 MHz clean, measured), FM-2 (the
Sun against a lower SEFD), FM-11 (the same scale with two fixed beams), FM-7, FM-25, FM-29.

**New scenarios (two):** FM-48, a drift scan is not what the pipeline assumes (10); FM-47, the beams alternate instead
of being split (6).

**Not affected by the array:** the data entering the processing (FM-22, FM-8), the timing model and the ephemerides
(FM-4, FM-5, FM-30, FM-31, FM-44), the oscillator classes (FM-10, FM-11, FM-45), the ADC's calibration and sampling
plan (FM-46, now for one board), the excision (FM-12, FM-17, FM-19, FM-20).

**Actions, changed and new.** None is built; each follows the usual rule (physics and plan first, then Jasper's OK).
They are not yet in notes §10 or in `roadmap.html`.

| # | Risk | Kind | What | Fits in |
|---|---|---|---|---|
| R1 (changed) | FM-43, FM-48 | development | Was: the effective area against direction as a weight. Now: a drift through a fixed beam in the fast simulator (a beam gain against time, one pass a day per beam), sub-ints weighted by the expected gain, and the position fit for one beam, two beams and a mount. Only then is there a simulated budget for the new array. | E3 |
| R8 (dropped) | FM-14 | — | The two-board bench test: one board suffices for four outputs. | — |
| R12 | FM-14, FM-1 | test | Measure both beams as built with drift scans of the Sun: direction, width, SEFD per beam, and the delay between the two beam paths. It confirms the 3 m² and the 70 K in one number, in the team's own terms. | first light |
| R13 | FM-47 | check | Splitter or switch: ask the team what the SPDT does; in a first recording look for periodic gaps and for lines at a switching rate. | first light |
| R14 | FM-43 | decision | The mount: a fixed panel (no position with one beam; east–west only with both beams on one pulsar), one motor on a polar-aligned axis (B0329+54 through the day), or full tracking (other pulsars). For a fixed panel: a tilt of 2.6° north for B0329+54 and the rotation that puts the second beam on its track. | with the team |
| R15 | FM-1 | decision, then development | The reference scenario for the new array: its design points (nominal and worst case, in SEFD and usable bandwidth), then the scenario runs. Every budget in the notes and the docs is still at 3.2·10⁵ Jy and 400 MHz. | E3 |

R4 (the oscillator), R6 (the RFI model and the survey; the measured 150 MHz belongs in it) and R7 (the sampling plan, now
for one board) stand as they are.

**Later on 23 Nov.** R15 is done for one design point: `runScenario.m`'s default is `"array3m2"` (6.44·10⁴ Jy, 150 MHz), and
the scenario and navigation results were rerun for it (Claude's runs; no worst case, at Jasper's "one design point"). R1 is
not done: the repository script still follows the pulsar, and the single-pass numbers come from a scratch copy. Simulated,
with those limits: SNR 11.5 per pass for a 6.3° beam (FM-1); ±9.4 km east after 30 days when the pulsar is followed, no
position from one fixed beam, 108 km along one axis from both (FM-43). Log: `logs/2026-11-23_array-3m2-analog-beams.md`.

**24 Nov.** R14 has an answer for now (Jasper: "its fixed for now, but we can still show tracking as a what-if scenario"),
with a beam of about 6° and a splitter as working assumptions. FM-43 therefore applies in full: with the panel fixed the
array gives no position from one beam. R1 (the drift through a fixed beam in the simulator) becomes the first development.
Two additions from an analysis of the reruns: J0437−4715 has wrong-peak fits about ±57 µs off in 2.5 % of its TOAs at
SNR 7.3, with an ordinary χ² (FM-26, FM-27); and under scintillation pChi2 < 0.05 rings for 10–11 % of B0329+54's TOAs,
because the channels differ in pulse amplitude (FM-16, FM-41). Log: `logs/2026-11-24_fixed-panel-and-two-findings.md`.

**24 Nov, FM-26 corrected and half fixed.** The row says a false TOA up to 500σ off "still pulls with 1/16 weight". In fact
any TOA more than about 154σ from the fit made `navFit` return NaN for everything (both exponentials of the mixture
underflow: 0/0). Fixed the same day with a one-line guard and `testNavFit` check 12 (Jasper: "apply the guard and the unit
test"); a clock step (FM-9) had the same effect. What the row says now holds: such a TOA pulls with weight 1/k², 0.5 of
the east error at 150σ and 1.9 at 500σ for one TOA in 300. The gate or uniform outlier component of action 12 is still
to do.

### 24 Nov 2026: the register against the detection literature

Jasper, asked whether to add the findings of the day's literature research to the register: "Yes, update it". The
research is on pulsar detection and how it performs in practice (`2026-11-24_detection-literature.md`; report and notes
in the git-ignored `reports/` and `research_notes/`). Each affected row above has a paragraph marked **24 Nov**. The
sources were read through a fetch tool's summaries; the claims used here were re-opened in a verification round.
"Estimate" marks arithmetic of the researchers or mine; "simulated" is the project's own runs (Claude's runs). The
ratings stay a draft for Jasper's review.

**In the team's terms:** nominal sensitivity 0.043 m²/K, SEFD 6.4·10⁴ Jy. SNR 7 per pass on B0329+54 needs an SEFD of
at most 1.05·10⁵ Jy (0.026 m²/K); a detection in half the passes with the phase unknown needs at most 1.65·10⁵ Jy
(0.017 m²/K). Beam centre, mean flux, 150 MHz (estimates from the simulated 11.5).

**Rating changed (one).**

| Risk | Was | Now | Why |
|---|---|---|---|
| FM-14 | 3 × 3 = 9 | 4 × 3 = 12 | Every new array with a measured figure came in at 1.2 to 3 times its ideal SEFD, and two analogue arrays found path errors of 15–20 cm after assembly. Here nothing can be corrected afterwards, and whether the 3 m² was measured on the sky is not known. |

**Rating unchanged, content changed:** FM-1 (one pass tests the instrument to a factor 1.5–2 only; the S/N convention),
FM-42 (no published red-noise loss for a fold with an ephemeris: the numbers of 3 Nov are the only ones, and the knee
has to be measured), FM-10 (40 ms of drift per pass at 20 ppm; a period search as the way out and its price), FM-9
(the published time-stamp faults, all found by timing), FM-48 (nobody weights by the beam; the loss is 6.6 %).

**New scenarios (two):** FM-49, the panel is mounted a few degrees off its tilt or rotation (9); FM-50, a first
detection that is not one, or a real one that cannot be shown to be one (9).

**Correction to R12.** The action of 23 Nov reads "drift scans of the Sun" through both beams. At 52° N the Sun never
comes closer than 28.6° to the zenith, so it crosses neither fixed beam: the second beam, 20° from the panel's normal,
gets to 8.6° from it at midsummer at best. Cas A passes 6.8° north of the zenith (beam response 0.04) and Cyg A 11.3°
south. An SEFD on the sky therefore needs the panel tilted for the measurement.

**Not affected:** the timing model and the ephemerides (FM-4, FM-5, FM-30, FM-31, FM-44), the excision rows, the
navigation fit (FM-26, FM-27), FM-43 (one fixed beam gives no position; the literature has no case of two fixed beams
hours apart on one pulsar) and FM-47.

**On adding passes** (the route to weaker pulsars in FM-1's actions): the published stacking results do not show √N
over many days (measured on quiet sky up to about 49 days at CHIME; the Pushchino gains are planning values). Adding
folded profiles in phase gives √N by construction in white noise, but needs the predicted phase good to about 1 ms
from day to day, which is the clock again (FM-10, FM-11).

**Actions, changed and new.** None is built or decided; each follows the usual rule (physics and plan first, then
Jasper's OK). They are not yet in notes §10 or in `roadmap.html`.

| # | Risk | Kind | What | Fits in |
|---|---|---|---|---|
| R12 (changed) | FM-14, FM-49 | test | Was: drift scans of the Sun through both beams. Now: tilt the panel for a day so that the Sun drifts through a beam's centre (power ratio to cold sky about 11 at 64 sfu, estimate), or 6.8° north for Cas A (a step of about 2.6 %, limited by gain stability). It gives the SEFD independently of the pulsar, to compare with 6.4·10⁴ Jy and with the 1.05·10⁵ Jy that SNR 7 needs. GNSS satellites cross the zenith and give beam shape and pointing near 1575 MHz, not the absolute gain. | first light |
| R16 | FM-14, FM-1 | question | For the team: were the 3 m² and the 70 K measured for the whole beamformed panel on the sky, or derived per tile; were the path lengths of the 16 tiles measured on the bench, and to what tolerance (1 cm rms costs about 4 %, 2 cm about 16 %). The answer moves the central estimate between 0.72 and 0.90 of nominal. | with the team |
| R17 | FM-42 | measure | The off-source total-power spectrum, 0.01–100 Hz, per beam and polarization, from the first recording: the knee for FM-42. It makes that row's "measure (Allan variance)" specific. | first light |
| R18 | FM-50 | test, then development | The confirmation set for the first detection: the same chain on off-source hours and with wrong periods (an empirical null, and a check of the threshold), the beam without the pulsar, the return on following days at the sidereal time, the same recording through DSPSR or PRESTO, and a count of everything that was tried. | first light |
| R19 | FM-10, FM-9, FM-45 | decision | The clock for the first detection: an external 10 MHz and PPS (the known-phase test, the period test and adding passes in phase are then available), or the bare crystal with a period search (threshold about 5.3, no period test). In both cases a PPS-locked pulse measured at every start and a log of every re-synchronisation. The first-light part of R4. | before first light |
| R20 | FM-49 | check | Survey the panel's tilt and rotation as mounted, to 1° or better, against R14's values (2.6° north; the diagonal about 69° from north). | at mounting |
| R21 | FM-1 | rule | Every S/N, measured or expected, names its convention. The project's is the matched filter on the full profile (11.5 per pass, equal to the radiometer equation with W10); the W50 figure for the same pass is 25. | reports, docs |

### Open

- The ratings are a draft for Jasper's review. The four estimates were checked (above); FM-11 rose to 25, FM-26 fell to 6.
- `testNavFitFaults` case 1 waits for Jasper's run. Cases 2–4 (FM-9, FM-4, FM-26) are on hold until the navigation
  redesign (week of 9 Nov), together with the project's changed final deliverable. FM-33 and FM-8 are not navigation
  tests; each still needs its physics and plan, and Jasper's OK.
- The excision page's limits (E9) still lack D2-22 (Jasper: "Not now").
- Committed at Jasper's "Commit now" (2 Nov): this log, notes §10 item 16, the roadmap row, `tests/testNavFitFaults.m`, `tests.html`, `gen_testmatrix.py` (+ its README row), and the history and maintenance entries. The ratings stay a draft.
- **3 Nov:** three scenarios added from the existing stages' real-hardware review (FM-40 quantization, FM-41 the bandpass against equal channel amplitudes, FM-42 common-mode gain noise), and FM-8, FM-9, FM-12, FM-13 and FM-39 annotated with their tests (`testFrontEndFaults`, `testFoldFaults`). Details, numbers and four new actions: `logs/2026-11-03_real-hardware-faults.md`; actions in notes §10 item 17.
- **21 Nov:** the register checked against the literature research of the same day (the section above): 31 rows
  annotated, six ratings changed, four scenarios added (FM-43 to FM-46), eleven actions R1–R11. The actions are not yet in
  notes §10 or in the roadmap. Log: `logs/2026-11-21_deep-research.md`.
- **23 Nov:** the register updated for the 3 m² array with two hardwired analogue beams (the section above): 12 rows
  annotated, three ratings changed, two scenarios added (FM-47, FM-48), action R1 changed, R8 dropped, R12–R15 new.
  Nothing was simulated for the new array. Log: `logs/2026-11-23_array-3m2-analog-beams.md`.
- **24 Nov:** the register checked against the detection literature of the same day (the section above): six rows
  annotated, one rating changed (FM-14, 9 → 12), two scenarios added (FM-49, FM-50), action R12 changed, R16–R21 new.
  Nothing was simulated for it. Log: `logs/2026-11-24_detection-literature.md`.
