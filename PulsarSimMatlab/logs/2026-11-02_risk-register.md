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

### The register, highest score first

| ID | Stage | Scenario | P | I | P×I | Seen | Action |
|---|---|---|---|---|---|---|---|
| FM-22 | 7 fold, `processIQ` | `processIQ` folds with the linear f0/F1 model only, without a predictor (D1-2). Real data drift by up to 10⁻⁴ (Doppler), about 0.36 s per hour, so any fold longer than seconds smears the pulse away. The noise calibration's on-pulse window and `rfiProfile` use the same model (D1-3, P3-3). | 5 | 5 | 25 | partly: no detection, without the reason | existing: D1-2, D1-3, P3-3 (F); the first thing F needs |
| FM-11 | 4a clock | The oscillator follows its temperature, which cycles daily. The clock gets a rate term at the solar day, 4 min per day from the sidereal day of the position signal: over 30 days the two drift apart by only 30° of phase (99 % of the term looks like position), over 100 days by 99° (88 %). A rate amplitude *a* gives a daily phase of *a* · 86400 s / 2π, which is 140 µs at 10⁻⁸. On B0329+54 that is 72 km, since 100 km east gives a 193 µs sine (`navigation.html`), or 2.4 times the fit's ±30 km equatorial error. A ±2 ppm TCXO changes by roughly 10⁻⁸ to 10⁻⁷ per K (estimate: the ±2 ppm band over a rated range of about 100 K, steeper in places; the SDR's datasheet decides). A room's daily swing of ±1 K then gives 70–700 km, 2–25 σ; an outdoor box more. An OCXO (about 10⁻⁹ per K or better) keeps it below 7 km. The simulated crystal has no such term (D1-17). | 5 | 5 | 25 | no: the fit reads it as position (30 days: χ² unchanged even at 10⁻⁷; 100 days: χ² 1.38 at 10⁻⁷) | existing: D1-17 (E5); **new: test, before E5**; test written: `testNavFitFaults` 1 |
| FM-9 | F ingest | The SDR drops samples (overflow at 18 × 400 MS/s) or restarts with a new time origin. The processing counts samples from `info.t0`, so every later time stamp is off by the lost time. That is a clock step the fit's clock model does not expect, and many small drops add up like extra clock noise. | 4 | 4 | 16 | no | **new: development + test** |
| FM-5 | 7a prediction | The pulsar's catalogue position is at an old epoch; proper motion and parallax are not modelled. For B0329+54 that is about 47 µs per year since the position epoch, as an annual term (notes §10 item 13, to check). A 30-day fit absorbs most of it into the clock; the rest biases a multi-pulsar solution. | 5 | 3 | 15 | no | existing: E2 (item 13, first in its proposed order) |
| FM-8 | 4 IQ, 5a | The real receiver's frequency mapping differs from the processing's: an inverted band (a mixer with its LO above the signal), or another LO, band or sample rate. The dedispersion chirp then has the wrong sign or size, and the pulse smears over up to twice the band's sweep (31.5 ms at DM 25, against B0329+54's 6.6 ms pulse). At the weak target, "inverted" looks the same as "too weak". The mixer type is not known yet. | 3 | 5 | 15 | no | **new: development + test** (before F) |
| FM-20 | excision | `detectRFI`, `periodicRFI` and `blankingWeights` work on whole files, and the weight file is 4.8 GB per second of data (D2-9, D2-23). Hours of real data cannot run. | 5 | 3 | 15 | yes: a memory error | existing: F (a design exists: time blocks with the filter reach as overlap) |
| FM-29 | 10 navigation | With one pulsar, the clock's wander looks like a polar offset: the polar component stays at about 64 km after 30 and after 100 days (nominal flux). The fit's clock covariance also leaves out flicker (D1-19). | 5 | 3 | 15 | yes: in the fit's covariance | existing: E3 (pulsars at other declinations), E5 (oscillator) |
| FM-30 | 7a, simulator | Truth and observer share `delayModel`, so a term missing from it is missing on both sides, and no simulation shows what it costs (notes §10 item 13). | 5 | 3 | 15 | no | existing: E2 (each term with a switch; the PINT/TEMPO2 end-to-end check) |
| FM-1 | sky | The pulsars are fainter than the worst-case design. Catalogue fluxes at our band are uncertain, and B0329+54 alone carries most of the information. At worst-case flux, B1642−03 and B1929+10 give no TOA in 30 days. | 3 | 4 | 12 | yes: few TOAs | existing: E1 (fluxes with uncertainties), E3b (below-threshold processing), E6 (levers) |
| FM-4 | sky, 7a | The spin ephemeris ages, normal pulsars have timing noise, and B0329+54 has glitched. A drift in one pulsar's phase is not the common clock, so it leaks into the clock rate and the position. The sizes need literature numbers per pulsar (item 13). | 4 | 3 | 12 | partly: a residual trend | existing: E1, E2 (numbers, a red-noise truth); **new: test** |
| FM-10 | 4a clock, 7a | A cold start without a clock calibration. The SDR's ±2 ppm moves the pulse 7.2 ms per hour, more than B0329+54's 6.6 ms pulse, so hour-long folds smear. The runs assume a port calibration to 10⁻⁹ in rate, from an external reference such as GNSS. Without one, as in a GNSS-denied start, the first folds may detect nothing, and without TOAs there is no clock estimate. | 3 | 4 | 12 | yes: no TOAs, but no way out | **new: test, then development** |
| FM-18 | 3 RFI | The site's real RFI differs from the 16-source model (`rfiRealistic`): more radars, DME pulses just below 1.2 GHz at the band edge, the SDR's own spurs and clock harmonics. A strong steady line costs most of its channel, one of 128 (0.78 % of the band). A comb of self-made lines costs one channel per line in the band: 40 channels (31 %) for the harmonics of a 10 MHz reference, 4 (3 %) for 100 MHz. | 4 | 3 | 12 | partly | **new: development** (a site survey) |
| FM-38 | F | The real-time processing of 18 streams cannot keep up, and data are dropped. | 4 | 3 | 12 | partly | existing: F (GPU/FPGA); FM-9's ingest makes the drops visible |
| FM-2 | sky, receiver | The Sun and the sky raise Tsys. The quiet Sun at 1.4 GHz is about 4·10⁵ Jy at solar minimum (a disk of about 10⁵ K), close to the nominal SEFD of 3.2·10⁵ Jy. Near the 2024–25 maximum it is higher, perhaps twice that (estimate). In the ~15° beam it cuts the SNR to 0.3–0.45 of its value. Through the sidelobes it adds 1.3–2.5 times the sidelobe level to Tsys: a few % at −20 dB, up to about 25 % at −10 dB. The array's sidelobe level is not known. Flares add bursts of minutes. The noise is measured per sub-int, so the error bars stay honest, but the forecasts, with a constant Tsys, are optimistic. B0950+08 lies 4.6° from the ecliptic, so the Sun is inside its half-beam for about 12 days around 20 August. | 5 | 2 | 10 | partly | **new: development** (Tsys(t) in the simulator) |
| FM-3 | sky, 8 TOAs | The profile changes, but the template is fixed. B0329+54 is in its abnormal mode about 15 % of the time. The profile also evolves across 1.2–1.6 GHz, and scintillation changes which frequencies dominate. Each gives a TOA offset that lasts hours; the size needs the mode profiles. | 5 | 2 | 10 | no: at SNR ~7, χ² cannot see it | **new: test** (with E1's profiles) |
| FM-21 | 5, 5b | The last dispersion sweep of each file is unusable: 9–16 % of a 0.2 s file (D1-24, D2-14). Streaming must carry the overlap across file boundaries. | 5 | 2 | 10 | yes | existing: F |
| FM-24 | 7a | Real data processed with the analytic orbit, which is the default. Its Earth is up to 3.2·10⁵ km off, and B0329+54's Roemer delay up to 0.66 s. No pulsar is detected, and nothing says why. | 2 | 5 | 10 | no | **new: development** (a guard) |
| FM-28 | 9 TOAs | The good-TOA rule asks for detection at unknown phase (4.15) even when the predictor knows the phase (3.09), so it discards TOAs it could keep (D1-7). | 5 | 2 | 10 | yes | existing: D1-7 (review day 3) |
| FM-31 | 10 navigation | Earth orientation: UT1 − UTC is one number (up to 0.9 s, about 260 m east), polar motion is about 10 m, and the rows use the rotation angle only (below 6 km on the offsets). | 5 | 2 | 10 | no | existing: parked at block D (IERS data for F) |
| FM-14 | E4 array | The 18 streams are not phase-coherent or not calibrated (several SDRs need one reference and a PPS), so the beam gain falls below its budget. | 3 | 3 | 9 | partly: a lower SNR | existing: E4 |
| FM-15 | E4, sky | The pulsar's components are polarized, with different angles. With unequal X and Y gains, or one polarization, the received profile changes with the parallactic angle and with Faraday rotation. The TOA then moves on the sidereal day, which is the position's own signature, and time does not average it away. The size is unknown. | 3 | 3 | 9 | no | **new: test** (E4) |
| FM-33 | simulator | The fast simulator and the voltage chain drift apart unnoticed: their agreement rests on experiments without PASS/FAIL (T2-1). Every E3 forecast rests on the simulator. | 3 | 3 | 9 | no | **new: test** |
| FM-37 | F port | The GPU/FPGA port differs from `processIQ` (precision, fixed point). | 3 | 3 | 9 | partly | existing: F (checked stage by stage against `processIQ`) |
| FM-16 | 5a–6b | Lopsided channels (narrow carriers, LTE edges, a carrier across a channel edge) break the flat-channel noise model: the error bars are 0.4–3.4 % too small, and `noiseCal` is off by default (D1-29). | 4 | 2 | 8 | partly: χ² | existing: D1-29 (decision, day 5) |
| FM-25 | 4a, 7a | Visibility gaps: for a millisecond pulsar, the clock's wander over a gap must stay below half a period (J0437: 2.9 ms), or the turn count is lost. B0329+54 has 357 ms of margin. | 2 | 4 | 8 | partly | none now (E1's choice of pulsars, E5) |
| FM-27 | 9, 10 | At SNR ~7 the error bars are about 5–10 % small, and 1.9 % of B0329+54's TOAs lie beyond 3σ. `navFit`'s mixture fixes the outlier fraction at 2 % and the scale at 4×; real data may hold more outliers. | 4 | 2 | 8 | partly: `chi2Core` | **new: development** (low) |
| FM-32 | F time | Absolute time: UTC → TT and leap seconds. A constant error is just the clock offset; a step (a leap second, a restart) is not. | 2 | 4 | 8 | no | existing: D1-15 (F); FM-9's ingest |
| FM-34 | simulator | Pulls of 1 show that the simulator agrees with the observer's model, not that the model is real (T2-3). Only real noise tests it. | 4 | 2 | 8 | no | **new: test** (the first F test) |
| FM-7 | sky | Scintillation fades last longer than modelled: Gaussian ACFs, no refractive variation (S2-5). Low-DM pulsars such as B0950+08 get long gaps. | 3 | 2 | 6 | yes: gaps | existing: S2-5 |
| FM-12 | 3–4 receiver | A strong nearby radar clips the ADC or compresses the LNA. The clipping is strong and gets blanked, but the gain dip after a pulse is a power drop, which `detectRFI` cannot see: it flags high power only (D1-11). | 3 | 2 | 6 | no | **new: test** (E4/F) |
| FM-17 | excision | The frequency guard has no cap: dense, independent strong bursts blank up to 99 % (D2-22). Real radars (duty ~0.1 %) cost far less, so this is a stress case. | 2 | 3 | 6 | partly: `blankedFraction` is reported, nothing warns | **new: development** (a warning); D2-22 stays |
| FM-26 | 9, 10 | A false TOA: noise or RFI passes the detection threshold and gives a TOA at a random phase. Noise alone rarely does: `adaptiveSubints` closes a group at T0 ≥ 7 at the predicted phase, and only the last group of each pass, below the target, is decided at P_FA 10⁻³. That is about one test per pulsar per pass, so about 0.1 false TOAs in 30 days, and a pulsar with no other TOAs absorbs one in its constant c_k. RFI residue is the likelier source. `navFit` puts a false TOA in its 4× wider component, but a residual of up to half a period (about 500 σ for B0329+54) still pulls with 1/16 weight. One among ~300 shifts the position by about 1–2 of its σ, depending on its phase. The detection tails below 10⁻³ are extrapolated (limits page). | 2 | 3 | 6 | partly: flagged as an outlier, still pulls | **new: test, then development** |
| FM-35 | drivers | `main.m` and the drivers are not run by `runAllTests` (bug #12 went unseen for two days), and scripts share the base workspace (D2-6). | 3 | 2 | 6 | partly | existing: D1-6, D2-6 (review days 4–5) |
| FM-36 | run control | A stale file is reused: `checkFold` does not check `noiseCalOnPulse` (D1-5), and RFI profile files do not record their excision settings (P3-4). | 2 | 3 | 6 | no | existing (minor) |
| FM-6 | propagation | DM error and DM variations (2.1 µs per 10⁻³), the Doppler-shifted dispersion (±7 µs annual for B0329+54), the ionosphere (3–70 ns), the troposphere and the solar wind. All are small against B0329+54's ~700 µs TOAs; they grow for millisecond pulsars. | 5 | 1 | 5 | no | existing: E2 |
| FM-13 | 4 IQ | A real zero-IF receiver adds I/Q imbalance, a DC spike and LO leakage at the in-band LO (1.3 GHz). They act as extra carriers, which the weighting and the excision handle (D1-10, R-I3). | 5 | 1 | 5 | partly | existing: F calibration |
| FM-39 | numerics | A time in seconds since J2000, held in one double, resolves only about 0.12 µs in 2026. The predictor works around it (phase from the span's start); wherever absolute times are still one double, ns timing would need a two-part time. | 5 | 1 | 5 | no | none (a note for millisecond pulsars) |
| FM-23 | 7 fold | A long fold needs a position prior: below 4 km for J0437, about 100 km for normal pulsars. At the fixed site the prior is known. | 1 | 4 | 4 | partly | none (the ship, parked) |
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

### Open

- The ratings are a draft for Jasper's review. The four estimates were checked (above); FM-11 rose to 25, FM-26 fell to 6.
- `testNavFitFaults` case 1 waits for Jasper's run. Cases 2–4 (FM-9, FM-4, FM-26) are on hold until the navigation
  redesign (week of 9 Nov), together with the project's changed final deliverable. FM-33 and FM-8 are not navigation
  tests; each still needs its physics and plan, and Jasper's OK.
- The excision page's limits (E9) still lack D2-22 (Jasper: "Not now").
- Committed at Jasper's "Commit now" (2 Nov): this log, notes §10 item 16, the roadmap row, `tests/testNavFitFaults.m`, `tests.html`, `gen_testmatrix.py` (+ its README row), and the history and maintenance entries. The ratings stay a draft.
