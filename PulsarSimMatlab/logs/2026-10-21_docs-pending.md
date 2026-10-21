# Docs changes on hold (21 Oct 2026)

Jasper is rewriting the docs site for readability in a separate session. Claude's doc changes wait here
and go into the pages afterwards, in the tone of the rewrite. Each item: the page, where, and the
content (facts and numbers; wording to be adapted). Remove an item once it is in the page.

## 1. limits.html, "Estimation and detection", after the error-bar entry (was commit 008cfed, reverted)

**Rare outliers at SNR ~7 (20–21 Oct).** With adaptive groups about 0.5 % of B0329+54's TOAs land one core
width off (0.005–0.009 turns, 5–9σ), at SNR 7.3–7.9; χ² cannot flag them. The template's autocorrelation
has no second peak there: the narrow core sits on a broad component, and noise now and then moves the best
match onto that shoulder, which the error bar (from the curvature at the peak) does not cover. Fixed groups
reach SNR ~5.5 and add true jumps to the side peak at ±0.025 turns (19σ). None at SNR 10 (0 of 600).

**Decision (21 Oct, Jasper): keep SNR 7 for more TOAs, and use redundancy.** Block D's navigation filter
checks each TOA against the prediction from the other TOAs and the ship's dead reckoning (B0329+54 nominal:
TOA error 220 km, prediction ~70 km; outliers 1100–2000 km = 4.8–8.6σ from the prediction; a 4σ gate
rejects good TOAs about once in 4–5 years), with an error model that has a small outlier part (rate and
width from the simulator). No fixed groups for multi-component pulsars. Unhandled, one 6σ outlier among a
day's 10 TOAs would bias that day by ~130 km.

## 2. roadmap.html, "Block D: motion and barycentric phase", after its first paragraph (was 008cfed)

**Agreed for D (21 Oct): an outlier-tolerant fit.** The adaptive SNR-7 groups give the most TOAs, but
~0.5 % of B0329+54's land one core width off (link to limits#estimation). The navigation fit uses
redundancy: each TOA is checked against the prediction from the others plus the ship's dead reckoning
(innovation gating at ~4σ), with an error model that includes a small outlier part. D's first deliverable
— a 1-D range from B0329+54 recovering a 100 km site offset — is tested on simulated TOAs that contain
these outliers.

## 3. sky.html, stage 1 generatePulsarSignal: the 'Components' option (21 Oct, code a 'components' commit)

New option 'Components' [nComp x 3] rows [amplitude, centre, FWHM] (turns): the POWER profile is the sum
of the Gaussians, scaled to peak 1 (continuous maximum: 2^16-point phase grid plus every component's
centre), the same convention as simulateFold / gaussianTemplate 'Components'. Components' phase 0 at
(k+0.5)T, so pulseCenters are the reference-point times. 'power' mode only (error otherwise);
dutycycle ignored. One centred component [1 0 duty/100] = the single Gaussian bit for bit (same
operations). Each component is summed around its own nearest image, as the single Gaussian.
tests.html, testGeneratePulsarSignal test 7: three components (0.6/0.4/0.3 at 0/-0.03/+0.06 turns,
FWHM 0.02/0.05/0.015): residual/sigma mean +0.004 (+- 0.073), std 1.025; centroid 0.490974 vs model
0.490967 (+- 2.6e-5); width 1.21120e-3 vs 1.21139e-3 (+- 1.3e-6); bit-for-bit single case; amplitude
mode refused. Tests 1-6 unchanged. Own temp folder per run. Claude's run 21 Oct: all passed (4.5 s).

## 4. parameters.html (pulsar rows) and pipeline/sky pages: profileComponents (21 Oct)

New parameter profileComponents ([] default = one Gaussian of dutycycle): [amplitude, centre, FWHM] rows
(turns) of a multi-component power profile, for the generator (main, runSNRSweep) and the observer's
template (ephem.profileFWHM then holds the component matrix; foldToTOA, gaussianTemplate and processIQ take
either). checkSkyFiles compares it too (sky files made before count as []), so reused sky files with
another profile are flagged.

## 4b. tests.html (experiments) and simulator.html: expComponentsChain (21 Oct)

Multi-component pulse on the voltage chain vs simulateFold('Components'): B0329+54's five components on the
10 ms pulsar, DM 2.6, -5 dB, channel path with excision, 6 seeds x 0.2 s (Claude's run, 28 min):
voltage chain 108 TOAs, median SNR 48.1, median sigma 1.505 us, rms/pred 1.026 (+- 0.068), mean red. chi2
1.0027, pChi2 < 0.05 0.083; simulateFold 1800 TOAs: 48.4, 1.495 us, 1.010 (+- 0.017), 1.0030, 0.082 -> agree
within ~1 %: the components are checked end to end (as C2 / C4a for the single Gaussian).

## 4c. New page section (block D, D1a): makeTrajectory — where the receiver is (21 Oct)

makeTrajectory(type, ...) -> struct with pos(t), vel(t) ([3 x N], m, m/s), frame ('ECEF' Earth-fixed or 'SSB'
barycentric equatorial), tRange, params. Time t = s since J2000.0 (2000-01-01 12:00 TT). Types: 'fixed'
(WGS84 lat/lon/h, 'OffsetENU' e.g. the true site 100 km east of the assumed one), 'track' (ship waypoints ->
ECEF, 'Interp' linear (default) / pchip (turns) / spline (orbits); across 180 deg fine), 'barycentric' (a
probe's tabulated SSB positions), 'custom' (PosFcn / VelFcn handles: the entry point for a Kalman filter's
propagated state; default velocity by central difference). testMakeTrajectory (7 tests, closed forms and an
independent Bowring inverse): fixed site to 1e-9 m / 1e-10 deg, 100 km east = 100.000000 km horizontal,
linear legs exact (10.02 m/s), pchip velocity = derivative (2e-8 m/s) and continuous at turns, 180 deg
crossing (midpoint 180.00000000 deg, chord 8.6 m below), 1 AU orbit with hourly spline samples to 1 mm,
custom velocity exact, 4 bad inputs refused. Claude's run 21 Oct: all passed (0.5 s); in runAllTests (17).

## 4d. Block D page (new): delayModel — receiver -> solar-system barycentre (21 Oct, D1b)

delta = t_SSB - t_obs = roemer (r_obs.n/c) + curvature (-(r^2 - (r.n)^2)/(2cD), with psr.distance) + einstein
(TDB - TT: 1.657 ms sin g + 14 us sin 2g + v_E.r_site/c^2) + shapiro (+2 T_sun ln((|r| + r.n)/AU), r from the Sun:
minus the delay -2 T ln(|r| + r.n) of light passing the Sun; T_sun = 4.925 us). Earth model analytic (decision 1):
EMB Keplerian (Standish elements and rates), Sun around the SSB from Jupiter and Saturn (0.004-0.005 AU in 2026),
Earth around the EMB (low-precision lunar series, 4424-4932 km), ecliptic -> equator J2000; Earth rotation angle
from UT1 = TT - 69.184 s, no precession / nutation / polar motion (consistent between the simulator's truth and
prediction; DE440 and IAU precession-nutation at the end of D). Frames: 'ECEF' trajectories rotated by the ERA,
'SSB' trajectories (probes) directly (no Einstein term: the probe's own clock later). Options 'Terms', 'Moon',
'Planets'. Outputs per term plus rObs, vObs, rEarth, vEarth, rSun, nHat.
Figures worth showing: annual Roemer for B0329+54 (+-412 s), daily term at 52 N (7.6 ms), the 100 km measurement
equation (193 us daily sine), Shapiro over a year (-17 ... +6 us).
testDelayModel (9 checks; Claude's run 21 Oct: all passed, 1.5 s): orbit a(1-e) / a(1+e) to 2e-7 AU; annual
peak-to-peak towards perihelion 998.013534 s vs 998.013536 s (e and a at the two dates; extremes on a 10 s grid -
the first draft's mid-year e and hourly grid were off by 1.1e-4 / 3e-5 s); daily term vs closed form 2e-13 s,
amplitude 7.607 ms, maximum at transit; 100 km east = (R dr.n)/c to 1.5e-13 s, amplitude 193.33 us, mean ~0,
Einstein site change 31 ns; ships (still = fixed bit for bit; moving: rate = v.n/c to 3e-7); probe closed forms;
Einstein 3.3145 ms peak-to-peak, site term 1.86 us; Shapiro B0329+54 -17.12 ... +5.80 us.

## 5. roadmap.html (block D / open items) and limits.html: clock jitter (Jasper, 21 Oct: "we completely forgot about clock jitter")

Listed long ago under "Block 2 – robustness" but never planned. Three effects, all for block D:
1. Time-base error of the station clock (oscillator offset, rate, wander; Allan deviation σ_y(τ)): shifts every
   time stamp, so every TOA, directly. Offset = the clock unknown already in the fit (degenerate with the polar
   position component for one pulsar); rate and wander accumulate: ~σ_y(τ)·τ over a TOA's integration — TCXO
   (~1e-9 at 1 h) 3.6 µs per hour, OCXO (~1e-11–1e-12) 4–36 ns per hour, rubidium ~ns. B0329+54's TOAs
   (734 µs) tolerate a TCXO; J0437's (13 µs over ~17 h) need ≲ 1e-10, i.e. an OCXO or better. Without GNSS (the
   point of pulsar navigation) the clock is a state of the navigation filter: offset and rate (and aging), with
   process noise from the oscillator's Allan variance (the standard two-state clock model, h0 / h-2).
2. Within an integration a fractional frequency error y shifts the phase by y·T, like a Doppler error: smearing
   when it drifts (1e-9 over 1 h = 3.6 µs, below B0329+54's bins, comparable to J0437's 5.6 µs bins only over
   many hours).
3. ADC sampling (aperture) jitter σ_j: phase noise on the RF samples, noise fraction (2π f σ_j)² — 1 ps at
   1.6 GHz: 1e-4 (-40 dB) of the input power: negligible for the noise-like pulsar, but strong RFI (a +40 dB
   radar) then brings a noise floor near the noise level during its pulses (blanked anyway).
Plan: D3 truth gets a clock process (offset + rate + random-walk frequency) on the true time stamps; D4's
state includes clock offset and rate (Kalman two-state clock model); ADC jitter in the voltage chain
(applyIQmodulation, truth side) as a later item.

## Still to come (results of the runs of 21–22 Oct; added here as they land)

- Validation run (Claude's run of the frozen 171a716): tests.html / results.html status.
- expSimVsChain at DM 25 / L 0.2 s (after main on the full band): tests.html.
- New RFI profiles (per-channel detectRFI thresholds, 30 lags): simulator.html.
- testRFIChain on the new code (narrow Inmarsat): tests.html, simulator.html assumptions b and d.
- 100-day runs: scenario.html (final budgets), roadmap status.
