# Computed figures of the docs

Each `gen_*.py` computes one figure from the code's formulas and constants (or from the
numbers in a results table) and writes the `<figure>` snippet as `<name>.svg` in the
current folder (e.g. `python3 docs/tools/figures/gen_toa.py` from anywhere). To update a figure,
run its script and replace the matching `<figure>`
block in the page. Needs numpy; `plotsvg.py` maps data to SVG coordinates.

| Script | Figure | Page |
|---|---|---|
| `gen_profile.py` | power profile G(t), neighbour sum at 20 % duty | sky.html (stage 1) |
| `gen_disp.py` | dispersion delay across the band, taper W(f) | sky.html (stage 2) |
| `gen_iq.py` | spectra through the IQ stage | receiver.html (stage 4) |
| `gen_excision.py` | one channel through detectRFI, blanking, dedispersion, W | excision.html |
| `gen_fold.py` | folding in three steps | fold.html |
| `gen_toa.py` | FFTFIT: profile, cross-correlation, weights | toa.html |
| `gen_det.py` | T0 and Tmax under H0 / H1 (simulated, ~1 min) | detection.html |
| `gen_results.py` | phase D sweep and threshold sweep | results.html |
| `gen_scen.py` | time per TOA vs precision of the candidate pulsars | scenario.html |

The flowcharts (pipeline run control, ground-truth rule, info chain, receiver blocks,
overlap-save, bin stamps, filterbank, good TOAs, test data, roadmap) are hand-written SVG
in the pages themselves. They use the classes in `style.css` (`.box`, `.line`, `.curve`, ...),
so they follow the light and dark themes.
