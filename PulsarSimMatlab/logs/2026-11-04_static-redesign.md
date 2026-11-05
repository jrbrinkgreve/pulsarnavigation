# 4 Nov 2026: navigation redesign for a static site (proposal, nothing built)

Jasper: "Let's continue the development, starting with a revision of phase D: the main objective changed to estimating the
position of the antenna array static on earth, by locking onto multiple pulsars. This means we can omit all the runship tests
etc, and re-do the navigation part. Let's re-design it, without dead reckoning to get a static position on earth".

**Status: proposal, waiting for Jasper's OK. No code or docs changed.** Python only, no MATLAB: covariance study and a toy search
in `.claude/runs/2026-11-04_static-study/` (scripts and outputs).

## The physics that sets the design
A site x (Earth-fixed) gives a TOA of pulsar k the delay r = -[cos(dec_k) (X cos H + Y sin H) + sin(dec_k) Z] / c + clock + J_k,
H the pulsar's hour angle at Greenwich. The Roemer part is linear in x however large x is.
1. No assumed site is needed. Fold against the geocentre: every site is within 21.3 ms, under 0.1 turn for the normal pulsars
   (B0329+54 0.03), so no ambiguity. J0437 stays ambiguous beyond 863 km.
2. Long folds smear before the site is known: up to 1.55 us/s (B0329+54 from 52 N: 0.55 us/s), and a +-2 ppm clock adds 2 us/s.
   A 120 s sub-int stays sharp (66 + 240 us against a 6.6 ms pulse), a 10 h TOA does not. So TOAs cannot start a cold fix at worst flux.
3. Z looks like a constant per pulsar, the same as J_k. Two routes to it: the site lies on the Earth's surface (north error = radial
   error / sin lat, no absolute phases, one pulsar is enough, a weak band near the equator, a hemisphere choice), or absolute pulsar
   phases with pulsars at different declinations (needs J_k to about 100 us).
4. A random clock wander hardly touches the position. A daily clock line does (FM-11). There the pulsars' different right ascensions
   matter: the position sine has a different phase per pulsar, the clock the same phase for all.
5. An absolute time error acts as a longitude error: 286 m per second at 52 N. Assumption: time known to a few seconds at the start.

## Numbers (1 sigma, 52 N 5 E, 20 deg mask; east / north km; height known; crystal clock; covariance only)
Check: B0329+54 alone, 30 d nominal, equatorial form: east 29.9, radial 37.6 km, against navFit's 29.6 / 37.2 on the simulated run.

| | TOAs (SNR 7), B0329+54 | all data, B0329+54 | all data, four pulsars |
|---|---|---|---|
| nominal 1 d | 188 / 455 | 185 / 349 | 158 / 323 |
| nominal 7 d | 62 / 100 | 61 / 98 | 57 / 96 |
| nominal 30 d | 30 / 48 | 30 / 47 | 28 / 46 |
| nominal 100 d | 16 / 26 | 16 / 26 | 15 / 26 |
| worst 1 d | no fit | 421 / 796 | 361 / 738 |
| worst 7 d | 226 / 824 | 139 / 223 | 130 / 220 |
| worst 30 d | 100 / 194 | 67 / 107 | 63 / 106 |
| worst 100 d | 60 / 94 | 37 / 59 | 34 / 58 |

Free daily clock line in the fit (all data; east / north km):

| | B0329+54 | + B0950+08 | four pulsars |
|---|---|---|---|
| nominal 30 d | 239 / 244 | 82 / 142 | 52 / 99 |
| nominal 100 d | 40 / 48 | 30 / 43 | 24 / 39 |
| worst 30 d | 545 / 555 | 193 / 331 | 120 / 230 |
| worst 100 d | 90 / 108 | 69 / 99 | 54 / 89 |

Clock noise 30 times the crystal's variance (about the TCXO estimate of 22 Oct): 30 d nominal 29.7 / 47.3 km, worst 67.4 / 107.3 km, so
random wander does not matter. Absolute phases without a height (30 d nominal, crystal clock; north / up km): J_k known to 100 us: 59 / 91 with two
pulsars, 43 / 74 with four; to 1 ms: 388 / 500 and 284 / 367.

Cold-start toy (`cold_start.py`, 8 seeds per row): B0329+54 alone, no position prior, clock rate unknown within +-2 ppm, phase unknown.
Grid search over X, Y and rate, 19k-137k trials. No RFI, no scintillation, no smear inside a sub-int, Gaussian pulse; its pass is 10 deg of
hour angle off the true transit (immaterial). The 600 km grid limits the best case to 181 km.

| flux, data | expected SNR | position error, median (max) |
|---|---|---|
| worst, 1 d | 9.9 | 599 (2367) km |
| worst, 2 d | 13.9 | 345 (1028) km |
| worst, 3 d | 17.0 | 181 (509) km |
| nominal, 0.5 d | 15.9 | 2135 (4065) km: half a day of sine looks like clock rate |
| nominal, 1 d | 22.4 | 471 (846) km |
| nominal, 2 d | 31.5 | 181 (726) km |

## The proposed chain
sub-int folds against the geocentre -> `siteSearch` (cold start, grid, no TOAs) -> no-threshold residuals per block at the current
site (first derivative of the matched filter over its expected curvature, so a fixed sigma at any SNR) -> `navFit` with priors (height, per-pulsar
constants) -> repeat from the new site. The TOA path stays as an independent check at nominal flux.

Steps, one new function each: D6a geocentre reference in makeTrajectory; D6b navFit priors and lat/lon output; D6c siteSearch; D6d
the no-threshold residuals; D6e driver runStatic. `navFilter`, `runShip` and their docs stay as they are (Jasper, 2 Nov: keep them for the ship),
`testNavFilter` leaves the default list of runAllTests but can be run by name.

## Open
- Pointing: a tracked beam needs the site. The 15 deg beam allows about +-330 km for under 0.5 dB; "anywhere" needs a grid of beams at the
  start, plus the array's heading and level. For E4; the fast simulator keeps a constant gain for now.
- Idea, not simulated: the beam itself measures direction. SNR 22 after a day with a 15 deg beam gives about 0.5 deg, about 50 km, better than
  timing at that stage, limited by how well level and heading are known. For E4.
- Older statements this changes: "the polar component needs absolute phases" (not with a known height); "extra pulsars add little" (true only
  without a daily clock term); FM-23 rises and FM-26 / FM-27 shrink on the no-threshold path; FM-11 stays the top risk.

## Later on 4 Nov: a simpler predictor for the fold (proposal, nothing built)
Jasper: "can we simplify the trajectory predictor for the folding, as the location does not need the input from the navigation filter".
Proposal: one predictor per observation, built for the geocentre and never rebuilt, so the fold takes no site or clock estimate. foldProfile
needs no change, as it only evaluates the polynomials; the filter's input came through the drivers (daily rebuilds in runShip and runScenario).
The Earth's orbit stays in the predictor (Doppler up to 1e-4, 0.36 s per hour). The site and clock enter after the fold as a shift per sub-int.
- Smear inside a 120 s sub-int, B0329+54 from 52 N: site 66 us + clock 240 us, against a 6.6 ms pulse: negligible.
- The shift must use the exact delay at the clock-corrected time: a +-2 ppm clock is up to 5.2 s off after 30 days, which is up to 518 us through the
  Doppler and 1.5 km east through the Earth's rotation.
- Limit: J0437 (0.14 ms) would smear (about 110 us from the site at 30 S, plus 240 us from the clock) and needs shorter sub-ints or a second fold. No
  millisecond pulsar is feasible from 52 N.
- Nothing is deleted: 'track', 'barycentric', 'custom', 'Clock' and 'Tol' stay for the ship and for testFoldPredictor and testTrueResiduals (probe trajectories).
