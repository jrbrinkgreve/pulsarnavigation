# Computed figures of the docs

Each `gen_*.py` computes one figure from the code's formulas and constants (or from the
numbers in a results table or a logged experiment) and writes the `<figure>` snippet as `<name>.svg` in the
current folder (e.g. `python3 docs/tools/figures/gen_toa.py` from anywhere). To update a figure,
run its script and replace the matching `<figure>`
block in the page. Needs numpy; `plotsvg.py` maps data to SVG coordinates.

| Script | Figure | Page |
|---|---|---|
| `gen_voltage.py` | the generated voltage: noise louder once per turn; power follows G (toy rate) | sky.html (stage 1) |
| `gen_profile.py` | power profile G(t), neighbour sum at 20 % duty | sky.html (stage 1) |
| `gen_rfspec.py` | RF spectrum: white to f_in/2, the band the dispersion stage keeps | sky.html (stage 1, design choices) |
| `gen_waterfall.py` | pulse in 16 sub-bands before and after dispersion | sky.html (stage 2) |
| `gen_disp.py` | dispersion delay across the band, taper W(f) | sky.html (stage 2) |
| `gen_kernel.py` | the kernel as a chirp; the design grid and its empty zone | sky.html (stage 2) |
| `gen_rxorder.py` | what dedispersion does to the pulsar, an impulse and a carrier (why RFI is added after dispersion) | receiver.html (stage 3) |
| `gen_snrscale.py` | SNR of one pulse against snrDB, design target and test runs | receiver.html (stage 3) |
| `gen_rfiscene.py` | the realistic RFI list sketched in time and frequency | receiver.html (stage 3) |
| `gen_radar.py` | radar pulse spectrum with instant / 0.1 µs edges; rotating-antenna gain | receiver.html (stage 3) |
| `gen_iq.py` | spectra through the IQ stage | receiver.html (stage 4) |
| `gen_iqfilter.py` | the IQ low-pass response, the band and the alias zones | receiver.html (stage 4) |
| `gen_smear.py` | smear left inside a channel by incoherent dedispersion, DM 5 and 100 | fullband.html (stage 5) |
| `gen_support.py` | fully supported data and complete turns of the 0.1 s file, both paths | fullband.html (stage 5) |
| `gen_fout.py` | the f_out sweep of 9 Oct (logged numbers): SNR and error-bar honesty vs bin width | fullband.html (stage 6) |
| `gen_dmmem.py` | memory of the dispersion filters vs DM (each function's sizing rule) | channels.html |
| `gen_powcov.py` | powerCovariance for the 3.125 MHz channels: power correlation, V and X(L) | channels.html |
| `gen_blanksub.py` | energy of a radar pulse vs SNR; what subtraction would leave | excision.html |
| `gen_windows.py` | detectRFI thresholds vs window length (needs scipy) | excision.html |
| `gen_blankvar.py` | experiment 3a: blanked-bin variance vs the W and W² rules (reads `data/mc/expBlankingVariance.mat`, local only; needs scipy) | excision.html |
| `gen_locked.py` | locked-radar TOA bias vs SNR per frequency guard (table numbers) | excision.html |
| `gen_periodic.py` | periodicRFI: found pulses, gaps, predicted pulses (schematic) | excision.html |
| `gen_nearest.py` | nearest-bin rounding errors repeat every turn | fold.html |
| `gen_subint.py` | turns per sub-int for SNR 6–7 vs snrDB | fold.html |
| `gen_chi2.py` | reduced χ² distribution, measured cases, threshold 2 (needs scipy) | toa.html |
| `gen_chanweights.py` | expected channel weights of the optimal fit with the illustrative RFI | toa.html |
| `gen_rice.py` | detection thresholds vs PFA: known, unknown (Rice), naive counts (needs scipy) | detection.html |
| `gen_pd.py` | detection probability vs SNR and the measured good-TOA fractions (needs scipy) | detection.html |
| `gen_qq.py` | reading validateTOA: Q-Q plots, uncertainty of the rms ratio (illustrative; needs scipy) | validation.html |
| `gen_excision.py` | one channel through detectRFI, blanking, dedispersion, W | excision.html |
| `gen_fold.py` | folding in three steps | fold.html |
| `gen_toa.py` | FFTFIT: profile, cross-correlation, weights | toa.html |
| `gen_det.py` | T0 and Tmax under H0 / H1 (simulated, ~1 min) | detection.html |
| `gen_results.py` | phase D sweep and threshold sweep | results.html |
| `gen_scen.py` | time per TOA vs precision of the candidate pulsars | scenario.html |

The flowcharts (pipeline run control, ground-truth rule, info chain, receiver blocks,
overlap-save, overlap-add (sky.html), excision chain (excision.html), bin stamps, filterbank, good TOAs, test data, roadmap) are hand-written SVG
in the pages themselves. They use the classes in `style.css` (`.box`, `.line`, `.curve`, ...),
so they follow the light and dark themes.
