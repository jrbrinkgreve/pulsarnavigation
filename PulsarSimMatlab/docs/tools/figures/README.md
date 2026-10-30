# Computed figures of the docs

Each `gen_*.py` computes one figure from the code's formulas and constants (or from the
numbers in a results table or a logged experiment) and writes the `<figure>` snippet as `<name>.svg` in the
current folder (e.g. `python3 docs/tools/figures/gen_toa.py` from anywhere). To update a figure,
run its script and replace the matching `<figure>`
block in the page. Each caption starts with its figure label as the page numbers it (`<b>Figure 8.1.</b>`, `<b>Figure E.3.</b>`;
reference pages `<b>Figure 1.</b>`); when a figure moves or the page is renumbered, change the label in its script too. The
scripts' captions match the pages word for word (checked 30 Oct by running every script: all 73 figures identical). Needs numpy; `plotsvg.py` maps data to SVG coordinates. `pfb.py` holds the prototype of `channelizeIQ` (same grid and design formula), shared by the four 5a figures.

| Script | Figure | Page |
|---|---|---|
| `gen_navgeom.py` | navigation geometry: a pulse front and an offset station; three lines of position | physics.html (Figure P.1) |
| `gen_whysquare.py` | why square: average voltage over turns is zero, average power shows the pulse (toy) | physics.html (Figure P.2) |
| `gen_radiometer.py` | averaged power scatter vs samples per average; half-band noise; Isserlis sums | physics.html (Figure P.4) |
| `gen_foldgrow.py` | folded profiles after 1–4096 turns; SNR vs observing time for several ρ | physics.html (Figure P.5) |
| `gen_flanks.py` | timing information on the pulse flanks; the bound √2·σ_t/SNR | physics.html (Figure P.6) |
| `gen_timeprec.py` | float32 vs double time resolution; a clock built by adding 1/fs drifts | conventions.html (Figure C.3) |
| `gen_datarate.py` | bytes written per second of signal by each stage, against the 24 GB memory | pipeline.html (Figure O.4) |
| `gen_freqplan.py` | the frequency plan: RF band, LO, baseband, IQ low-pass, 128 channels | parameters.html (Figure 1) |
| `gen_timescales.py` | every time scale of the pipeline on one log ruler (J0437's TOA: the catalogue estimate and the simulated 82–95 h) | parameters.html (Figure 2) |
| `gen_visibility.py` | hours per day above 20° elevation vs station latitude, candidate pulsars | scenario.html (Figure R.2) |
| `gen_experiments.py` | what each run script varies and keeps fixed (the voltage-chain scripts, runRFIProfile, and the fast-simulator scripts runScenario and runShip) | scripts.html (Figure 1) |
| `gen_testmatrix.py` | regression / algebra / physics checks per unit test, run times (the 31 tests of `runAllTests` and experiment 3a) | tests.html (Figure 2) |
| `gen_honesty.py` | measured ÷ stated error bar for every check of the validation record, blocks 1 to D (logged numbers) | results.html (Figure 1) |
| `gen_effects.py` | size of each effect if left out, as a timing error, with what block D models and what is parked | limits.html (Figure 1) |
| `gen_neighbour.py` | lesson #8: linear assignment shares noise; variance of sums up to 1.5× | lessons.html (Figure 1) |
| `gen_timeline.py` | the blocks on the calendar, 29 Sep – 22 Oct, with session days | history.html |
| `gen_roadmap.py` | the roadmap: blocks on the calendar, block C's steps by status, what D, E and F hold | roadmap.html (Figure 1) |
| `gen_voltage.py` | the generated voltage: noise louder once per turn; power follows G (toy rate) | sky.html (Figure 1.1) |
| `gen_profile.py` | power profile G(t), neighbour sum at 20 % duty | sky.html (Figure 1.2) |
| `gen_rfspec.py` | RF spectrum: white to f_in/2, the band the dispersion stage keeps | sky.html (Figure 1.3) |
| `gen_waterfall.py` | pulse in 16 sub-bands before and after dispersion | sky.html (Figure 2.1) |
| `gen_disp.py` | dispersion delay across the band, taper W(f) | sky.html (Figure 2.2) |
| `gen_kernel.py` | the kernel as a chirp; the design grid and its empty zone | sky.html (Figure 2.3) |
| `gen_rxorder.py` | what dedispersion does to the pulsar, an impulse and a carrier (why RFI is added after dispersion) | receiver.html (Figure 3.2) |
| `gen_snrscale.py` | SNR of one pulse against snrDB, design target and test runs | receiver.html (Figure 3.3) |
| `gen_rfiscene.py` | the realistic RFI list sketched in time and frequency | receiver.html (Figure 3.4) |
| `gen_radar.py` | radar pulse spectrum with instant / 0.1 µs edges; rotating-antenna gain | receiver.html (Figure 3.5) |
| `gen_iq.py` | spectra through the IQ stage | receiver.html (Figure 4.1) |
| `gen_iqfilter.py` | the IQ low-pass response, the band and the alias zones | receiver.html (Figure 4.2) |
| `gen_smear.py` | smear left inside a channel by incoherent dedispersion, DM 25 and 100 | fullband.html (Figure 5.3) |
| `gen_support.py` | fully supported data and complete turns of the 0.1 s file, both paths | fullband.html (Figure 5.2) |
| `gen_fout.py` | the f_out sweep of 9 Oct (logged numbers): SNR and error-bar honesty vs bin width | fullband.html (Figure 6.2) |
| `gen_dmmem.py` | memory of the dispersion filters vs DM (each function's sizing rule) | channels.html (Figure 5b.1) |
| `gen_powcov.py` | powerCovariance for the 3.125 MHz channels: power correlation, V and X(L) | channels.html (Figure 6b.1) |
| `gen_proto.py` | the 5a prototype: a tapered sinc in time (zero crossings every 192 samples), a block with sloped edges in frequency | channels.html (Figure 5a.1) |
| `gen_pfbbank.py` | the 256 DFT slots and the 128 kept rows; six channel responses around the LO | channels.html (Figure 5a.2) |
| `gen_pfbfold.py` | folding keeps every 16th bin: a weighted frame with three tones; correlation of two outputs against their spacing | channels.html (Figure 5a.3) |
| `gen_pfbalias.py` | the copies of a channel response when the frame slides by D = 192 and by 256 | channels.html (Figure 5a.4) |
| `gen_blanksub.py` | energy of a radar pulse vs SNR; what subtraction would leave | excision.html (Figure E.2) |
| `gen_windows.py` | detectRFI thresholds vs window length (needs scipy) | excision.html (Figure E.4) |
| `gen_blankvar.py` | experiment 3a: blanked-bin variance vs the W and W² rules (reads `data/mc/expBlankingVariance.mat`, local only; needs scipy) | excision.html (Figure E.5) |
| `gen_locked.py` | locked-radar TOA bias vs SNR per frequency guard (table numbers) | excision.html (Figure E.6) |
| `gen_periodic.py` | periodicRFI: found pulses, gaps, predicted pulses (schematic) | excision.html (Figure E.7) |
| `gen_nearest.py` | nearest-bin rounding errors repeat every turn | fold.html (Figure 7.3) |
| `gen_subint.py` | turns per sub-int for SNR 6–7 vs snrDB | fold.html (Figure 7.4) |
| `gen_chi2.py` | reduced χ² distribution, measured cases, threshold 2 (needs scipy) | toa.html (Figure 8.2) |
| `gen_chanweights.py` | expected channel weights of the optimal fit with the illustrative RFI | toa.html (Figure 8.3) |
| `gen_rice.py` | detection thresholds vs PFA: known, unknown (Rice), naive counts (needs scipy) | detection.html (Figure D.1) |
| `gen_pd.py` | detection probability vs SNR and the measured good-TOA fractions (needs scipy) | detection.html (Figure D.4) |
| `gen_qq.py` | reading validateTOA: Q-Q plots, uncertainty of the rms ratio (illustrative; needs scipy) | validation.html (Figure 9.2) |
| `gen_excision.py` | one channel through detectRFI, blanking, dedispersion, W | excision.html (Figure E.3) |
| `gen_fold.py` | folding in three steps | fold.html (Figure 7.1) |
| `gen_toa.py` | FFTFIT: profile, cross-correlation, weights | toa.html (Figure 8.1) |
| `gen_navsignal.py` | the navigation signal: three days of simulated B0329+54 TOAs and 30 days by hour angle against the 193 µs daily sine (reads `data/scenario/pre-22oct-tag/…_30d_motion_perfect.mat`, local only; needs scipy) | navigation.html (Figure 10.1) |
| `gen_navship.py` | the ship: one ship's drift (truth against the smoothed filter, east and north); the track errors of eight ships, dead reckoning against both starts (reads `data/scenario/ship_14d_nominal*.mat`, local only; needs scipy) | navigation.html (Figures 10.5, 10.6) |
| `gen_navfit.py` | navFit converging: the east offset fitted to the first N of 100 days (reads `data/scenario/scenario_B0329_nominal_clean_100d_motion_perfect.mat`, local only; needs scipy) | navigation.html (Figure 10.3) |
| `gen_navavg.py` | a TOA is a time-average: the daily signal of a 100 km offset, a 10 h TOA and a two-pass TOA against one time per TOA | navigation.html (Figure 10.2) |
| `gen_navplan.py` | the planning study: the polar error with a crystal clock per pulsar set, 30 and 100 days (logged numbers) | navigation.html (Figure 10.4) |
| `gen_det.py` | T0 and Tmax under H0 / H1 (simulated, ~1 min) | detection.html (Figure D.2) |
| `gen_results.py` | phase D sweep and threshold sweep | results.html (Figure 2) |
| `gen_simvschain.py` | fast simulator vs voltage chain: median σ and P_D against ρ (C2) | simulator.html (Figure F.3) |
| `gen_scint.py` | scintillation as simulateFold draws it: dynamic spectra (J0437, 50 MHz scintles), J0437's band-averaged flux | physics.html (Figure P.7) |
| `gen_scen.py` | time per TOA vs precision of the candidate pulsars | scenario.html (Figure R.1) |
| `gen_clock.py` | the receiver clock's "crystal" preset: σ_y(τ) and its three noise parts from `makeClock`'s formula; the time error τ·σ_y(τ) against the TCXO and OCXO estimates and the pulse widths | receiver.html (Figure 4a.1) |
| `gen_clockdraw.py` | what `makeClock` draws: one day of each noise (white, flicker, random walk) as time error; eight crystal clocks over 30 days against τ·σ_y(τ) (a port of its three generators; numpy's random numbers) | receiver.html (Figure 4a.2) |
| `gen_delays.py` | the delay terms for B0329+54 from 52° N: orbital Roemer over 2026, rotational Roemer over two days, Einstein, Shapiro; every term's size on one log scale against the TOA noise (uses `delaymodel.py`) | localisation.html (Figure L.3) |
| `gen_timesys.py` | TAI, UTC, TDB and TCB as offsets from TT, 1972–2030 (leap-second table of IERS Bulletin C; IAU L_B) | localisation.html (Figure L.2) |
| `gen_dopsmear.py` | why the phase must be predicted: B0329+54's pulse stack in a 120 s sub-int with a Doppler factor of 10⁻⁴, and the smeared fold | fold.html (Figure 7a.2) |
| `gen_predseg.py` | `makePredictor`'s segments around a ship's turn: testMakePredictor test 5 replayed with a port of its loop (15 segments, the shortest 7 s; uses `delaymodel.py`) | fold.html (Figure 7a.3) |
| `gen_simdraw.py` | how `simulateFold` draws a sub-int: the mean with its noise band, three drawn sub-ints, the neighbour correlation ¼ (illustrative: 128 bins) | simulator.html (Figure F.1) |
| `gen_adaptive.py` | `adaptiveSubints`' rule under scintillation: the flux factor and the T0 of each group closing at 7 (illustrative, SNR 1.2 per sub-int) | simulator.html (Figure F.2) |
| `gen_pulls.py` | why long runs are judged by fixed-time pulls over many runs: eight correlated runs, and the rms each method measures over 2000 repetitions (illustrative) | validation.html (Figure 9.1) |

`delaymodel.py` is not a figure: it ports `delayModel`'s analytic branch and `makeTrajectory`'s `'fixed'` and linear `'track'` types
to Python, for `gen_delays.py` and `gen_predseg.py`. It reproduces the tests' numbers: 7.607 ms daily at 52° N, Shapiro −17.1 to
+5.8 µs, and test 5 of `testMakePredictor` (15 segments, the shortest 7 s).

The flowcharts (pipeline run control, ground-truth rule, info chain, receiver blocks,
overlap-save, overlap-add (sky.html), excision chain (excision.html), bin stamps, filterbank, good TOAs, test data, the chain of stage 7a
(fold.html#prediction), long runs and a moving receiver (pipeline.html), the two chains of localisation.html) are hand-written SVG
in the pages themselves. They use the classes in `style.css` (`.box`, `.line`, `.curve`, ...),
so they follow the light and dark themes.
