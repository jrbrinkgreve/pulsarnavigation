# 2026-10-13 – docs: design reasoning, readability, figures (sky page)

Jasper's request: more reasoning in the docs on *why* each algorithm is built the way it is;
less dense wording; more visualisations. Plan agreed in the chat: a "Design choices" section
per stage page (choice, alternative, why, evidence), one stage at a time in the review-week
order. Sky page (stages 1–2) done first. No code change.

## Done

- `docs/sky.html` rewritten. Same facts, new layout: "In short", "The physics", "How it
  works" in plain steps, exact formulas folded into `<details>`, then "Design choices".
  Options/outputs, validation, history and limits kept (lightly reworded; validation now
  points to D1, the missing unit test).
- Design choices, stage 1: voltage at RF (not detected power or baseband); 4 GHz sampling;
  random noise, not a pulse waveform; one Gaussian component; `'power'` mode; first pulse
  at T/2 (now mentions R3); absolute time + private random stream.
- Design choices, stage 2: exact chirp on the voltage (not sub-band shifts or a Taylor
  expansion); band-limited filter; frequency taper, not a time window; designed in
  frequency on a grid ≥ 2·Nk; overlap-add in blocks; fRef at the top, bulk delay recorded
  (now mentions D2); output as long as the input.
- New figures (computed from the code's constants; scripts in `docs/tools/figures/`, README
  updated): `gen_voltage.py`, `gen_rfspec.py`, `gen_waterfall.py`, `gen_kernel.py`; an
  overlap-add flowchart written into the page. Checked by rendering the page (Quick Look).
- FFT block-size table for the 1 s reference run (computed with `chooseBlockSize`'s rule):
  2²⁸ chosen; 2²⁹ would be 2.5 % less work but needs 23.6 GB (> 16 GB budget).
- `docs/maintenance.html` changelog row.

## Confirmed (Jasper, 13 Oct: "the reasons are correct")

These reasons were not recorded in the notes or logs; Claude reconstructed them, Jasper
confirmed them, and they are now in notes §12:
- why simulate a voltage at RF rather than detected power or baseband;
- why f_in = 4 GHz (Nyquist > 3.2 GHz, margin, whole-number decimation to 800 / 500 MHz);
- why a single Gaussian component;
- why a time-domain chirp formula was not used (stationary-phase approximation);
- why overlap-add rather than overlap-save in the forward stage;
- why the late tail is dropped (mirrors a real observation; `fullySupported` downstream).

## Next

Receiver page (stages 3–4) in the same format (started the same day, below).

## Receiver page (stages 3–4)

Jasper confirmed the sky-page reasons and asked to go on. Done the same way:
- `docs/receiver.html` rewritten: "In short", physics, plain steps, details folded, the RFI
  model with a new "in time and frequency" column, "Design choices" for each stage. All
  anchors other pages link to are kept (`#noise`, `#rfi`, `#risetime`, `#rotation`,
  `#noisetype`, `#realistic`, `#iq`, `#lo`).
- Design choices, stage 3: noise and RFI at RF after dispersion; SNR as a ratio of spectral
  densities at the pulse peak; INR relative to all in-band noise; a few simple RFI types,
  not recorded RFI; own random streams tied to the absolute sample; two RFI lists. The
  type-specific reasons (edges, rotation, `'noise'`, PRF 373 Hz, 1351.3 MHz) stay with each type.
- Design choices, stage 4: an analog I/Q receiver model (mix, low-pass, decimate); LO at
  1.3 GHz; 800 MHz output; the low-pass (360 MHz, 80 dB, 2049 taps); odd length and exact
  delay removal; the factor 2.
- New figures: `gen_rxorder.py` (dedispersion turns an impulse into an inverse sweep, leaves a
  carrier), `gen_snrscale.py` (SNR per pulse vs snrDB; −54 dB → ~62 h per TOA at SNR 7),
  `gen_rfiscene.py` (realistic list in time and frequency, continuous vs intermittent),
  `gen_radar.py` (pulse spectrum 1/f² vs 1/f⁶; rotating-antenna gain, sidelobes 98.8 % of
  the scan), `gen_iqfilter.py` (low-pass gain, flat to ±355 MHz, −80 dB beyond ±365 MHz,
  alias zones). Checked by rendering.
- Fixed: the −20 dB prediction (docs table and notes §5.3) showed 3.88 / 12.3 / 38.8 / 77 µs,
  the values without self-noise; `predictPerformance` includes self-noise and prints 3.85 /
  12.2 / 38.5 / 77.8 µs (recomputed with its formula). −5 dB (97.9, 2.85 µs) was right.

### Confirmed (Jasper, 13 Oct), added to notes §12

Reasons not recorded in the notes or logs, reconstructed by Claude:
- adding RFI after stage 4 instead would mean the I/Q filter never sees it;
- the SNR alternatives (per sample, per bin, after folding) depend on later choices;
- INR as total power: the one measure that means the same for a carrier, an LTE block and a
  burst, and what a link budget gives;
- simple RFI types instead of a recording: a recording has no ground truth; the types map
  onto blank (intermittent) vs down-weight (continuous);
- the I/Q stage as mix / low-pass / decimate instead of a whole-file FFT conversion
  (memory; hardware-like timing conventions);
- the low-pass: 80 dB because the receiver noise is not band-limited before it, and four
  alias copies fold onto each in-band frequency (≤ 10⁻⁸ each); 2049 taps for a narrow
  9.8 MHz transition, cheap with FFT convolution.

Added to notes §12 on 13 Oct.

## Remaining algorithm pages (Jasper: "continue working like this for the full catalog")

Same format on fullband, channels, excision, fold, toa, detection and validation: summaries
and plain step lists, details folded where long, design-choice boxes, and new figures
(scripts in `docs/tools/figures/`, all listed in its README):
- fullband: `gen_smear.py` (incoherent smear per channel width, DM 5 / 100), `gen_support.py`
  (fully supported data, 9 vs 8 TOAs), `gen_fout.py` (the 9 Oct f_out sweep, logged numbers).
- channels: `gen_dmmem.py` (filter memory vs DM, each function's sizing rule), `gen_powcov.py`
  (powerCovariance recomputed: matches V 0.8877, X 0.0458/0.0051/0.0020, Bnoise 3.0596 MHz).
- excision: `gen_blanksub.py`, `gen_windows.py` (Gamma thresholds within 0.5 % of the exact
  ones; η₂₅₆ ≈ 1.33 vs the −5 dB pulse peak 1.32), `gen_blankvar.py` (experiment 3a's saved
  results), `gen_locked.py` (B5b table), `gen_periodic.py` (schematic); an excision-chain flowchart.
- fold: `gen_nearest.py` (rounding errors repeat every turn), `gen_subint.py` (turns per sub-int).
- toa: `gen_chi2.py`, `gen_chanweights.py` (expected 1/level² weights with the illustrative RFI:
  carrier channel ~6·10⁻⁴; all three sources cost ~3 % SNR).
- detection: `gen_rice.py` (thresholds vs PFA; 2048 independent bins would give 4.90, 5.5
  trials 3.57, Rice 4.15), `gen_pd.py` (P_D vs SNR with the 6 Oct good-TOA fractions).
- validation: `gen_qq.py` (Q-Q reading, ratio uncertainty 1/√(2n)); boxes on predicted σ and on
  keeping the truth out of the processing (both recorded decisions).

## Site check (all 26 pages)

- Structure (script): tags balanced, no duplicate ids, every internal link and anchor resolves,
  every page in `nav.js`, every CSS class defined, every SVG has role/aria-label. Clean (the
  home page's hero figure has no caption, by design).
- Option defaults in the docs vs the code's `arguments` blocks: all agree (differences only in
  notation, e.g. `[]` that resolves to fHigh).
- Recurring numbers across pages (hours per TOA, 373·ρ / 388·ρ, Bnoise, 80,000×): consistent.
- Fixed: forward-dispersion memory at DM 100 was "~94 GB (extrapolated ×2)" in sky, limits,
  roadmap and notes §10 item 5; the sizing rule gives 2³⁰ → 47 GB for DM ~54–106 (as DM 60).
  The −20 dB prediction (3.88 → 3.85 etc., see above). Channel filter at DM 100 ~120 MB (was
  ~100 MB). Excision page: experiment 3a errors −50 % / +41 % (the "±45 %" summary stays
  elsewhere). Roadmap code-review row: "commit pending" → committed as `0499da7`. All listed
  in maintenance.html's discrepancy table.
- Every new figure was rendered and checked by eye (Quick Look); label overlaps fixed.

## Confirmed (Jasper, 13 Oct: "the reasons are correct"), added to notes §12

Receiver page: RFI added after stage 4 would bypass the I/Q filter; SNR per sample / bin / fold
depends on later choices; INR as total in-band power; simple RFI types instead of a recording
(no ground truth); I/Q as mix–filter–decimate instead of a whole-file FFT; the low-pass's 80 dB
(white receiver noise, four alias copies) and 2049 taps. Full band: partly supported samples
would bias the profile. Channels: the 128 × 3.125 MHz trade-off. Fold: cutting turns halfway
between pulses; storing the weights but not σ² (one fold serves both noise models). Validation:
"rely on Monte Carlo runs and sweeps, not one main.m run".
Added to notes §12 on 13 Oct and committed with the docs.
