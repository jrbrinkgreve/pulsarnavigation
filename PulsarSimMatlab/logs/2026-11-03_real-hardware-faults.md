# 3 Nov 2026: what a real receiver does to the existing stages, and two fault tests

Session log. Jasper, while studying the RFI excision math: "can you look for possible problems after real implementation
in the design of the existing stages (except the nav filter) and write tests for them? Take your time, go autonomously
for now, I am studying". No change to `functions/`. Two new test files: `tests/testFrontEndFaults.m` (the voltage chain)
and `tests/testFoldFaults.m` (the fold, with the fast simulator). Both pass in Claude's runs; neither is validated yet.
Navigation (`navFit`, `navFilter`) is left out, because it is on hold for the redesign (week of 9 Nov). Committed
after a second run (Jasper: "run the tests, then commit if good").

## How I looked

I read each processing stage with one question: what does a real receiver or recording do that the simulated one does
not? The stages read were `processIQ`, `channelizeIQ`, `dedisperseChannels`, `detectChannels`, `detectRFI`,
`blankChannels`, `blankingWeights`, `foldProfile`, `estimateTOA` (the 'optimal' fit, its χ² and error bars),
`detectPulsar` (thresholds, noise ratio), `foldToTOA` and `simulateFold`. I then compared each finding with the 2 Nov
register (`logs/2026-11-02_risk-register.md`). For every finding that a test could settle without a code change, I put
the fault into the data on the truth side, ran the observer's code unchanged, and checked the result against a
prediction. The rules of `testNavFitFaults` apply: the test checks the response the physics predicts, not whether the
pipeline is "good".

## What came out

Three scenarios are new (FM-40 to FM-42). Five from the register now have a test or a sharper description: FM-8, FM-9,
FM-12, FM-13 and FM-39. The ratings use the register's scales (P × I, 1–5 each).

| ID | Stage | Scenario | P | I | P×I | Seen | Action |
|---|---|---|---|---|---|---|---|
| FM-42 *(new)* | 8–9 TOA, detection | **Common-mode gain noise.** The LNA's and the receiver's gain fluctuate (1/f, temperature, supply), and every channel's power moves with it at once. The noise model is radiometer noise only, independent between channels. A common fluctuation adds its full variance to the combined profile, but only 1/128 of it to each channel. So the per-channel χ² is blind to it, while detection and TOAs take all of it. It weighs most at the low harmonics, which the detection statistic uses and the TOA hardly does. For B0329+54 with a 1/f gain spectrum: with the knee at 1 Hz, false alarms rise 4× and the TOA error rises 1 %. With the knee at 10 Hz, the detection statistic's variance is ×2.7, unknown-phase false alarms go from 10⁻³ to 0.36 per sub-int, and the TOA error rises 10 %. The knee of an uncooled L-band chain over 400 MHz is not known; 1–10 Hz is plausible (estimate). | 3 | 4 | 12 | partly: `noiseRatio` rises (+1 % at 1 Hz, +10 % at 10 Hz) and is printed, but nothing acts on it; the per-channel χ² is blind | **new: measure, then development** (below) |
| FM-41 *(new)* | 8 TOA (`estimateTOA` 'optimal') | **The bandpass against equal channel amplitudes.** The fit models the pulse as b·s_c·T with s_c = `ChannelGain`, which `foldToTOA` (and so the whole chain) leaves at 1. That means the same absolute pulse in every channel. A real band is not flat (antenna match, LNA and cable slopes, standing-wave ripple, the SDR's filter), and its gain G_c scales the pulse and the noise of channel c alike. The fit weights each channel by 1/G_c², so the low-gain channels count most, where the pulse is smallest. The SNR keeps Σ(1/G)/√(nΣ1/G²): ×0.98 for a 3 dB tilt, ×0.93 for 6 dB and ×0.84 for 10 dB, which is +4 %, +16 % and +41 % integration time. The error bars stay honest. The pulsar's own spectrum (index −1.6) costs ×0.99 the same way. | 5 | 2 | 10 | no: only the SNR is lower | **new: development (small)**: `ChannelGain` from each channel's mean power, which is observer information; the test shows it gives the flat band's SNR exactly |
| FM-40 *(new)* | 4 IQ (ADC, transport) | **Quantization.** An SDR streams 16 or 8 bits, and 18 streams at 400 MS/s favour 8. The noise model stays honest (χ² unchanged to 0.002), but the SNR drops by 1/(1+ε), with ε = the quantization noise in the band. At 8 bits with the noise at rms 4 LSB, the loss is 0.5 %. At rms 1 LSB, the level that leaves room for a +40 dB radar, it is 7 % (×0.93). At 4 bits and rms 3 LSB it is 1 %; at 2 bits, 11 %. | 4 | 2 | 8 | no | **new: decision before F**: sc16, or sc8 with the noise at ≥ 2–4 LSB plus a clip flag (FM-12) |
| FM-8 | 4 IQ, 5a | The mirrored band (rated 15). **Test written** (`testFrontEndFaults` 1). On the 10 ms test pulsar (DM 5, LO at the band centre), 4.7 % of the amplitude is left (predicted 4.3 %), half a turn off. All 8 sub-ints still pass the good-TOA rule there, with red. χ² 1.11–1.13 and none flagged: a bright pulsar is "found" in the wrong place. For B0329+54, 41 % is left, at the right phase (the smear is symmetric): an SNR-7 TOA becomes 2.9, below every threshold. An LO error of 1 MHz shifts the TOAs by +15.5 µs at DM 5 (predicted +15.8; 85 µs per MHz at B0329+54's DM). That is a constant per pulsar, which `navFit` "self" absorbs; the TCXO's ±2.8 kHz is 0.24 µs. | 3 | 5 | 15 | no | unchanged: the switch and a known-line check before F (action 2) |
| FM-9 | F ingest | Lost samples (rated 16). **Processing half tested** (`testFrontEndFaults` 3). A drop joined without a gap is an exact clock step: every later TOA moves by −n/fs (−14σ per sub-int for 41 µs on the test pulsar), and nothing in the processing sees it. A zero-filled gap without a mask moves the sub-int it hits by 1.0σ, with χ² +0.12, below the flag. The same gap blanked through the existing excision path (zeros, plus mask rows over the gap and the channel filter's reach) moves it 0.17σ, with χ² unchanged. So the planned F ingest works with today's `blankChannels` and `blankingWeights`. | 4 | 4 | 16 | no (joined); partly (zero-filled, χ² alarm) | the ingest design is confirmed; the navigation half stays on hold |
| FM-12 | 3–4 receiver | A strong radar clips the ADC (rated 6). **Found while reading, not tested.** While a clipped pulse lasts, the saturated ADC suppresses the noise in *every* channel, not only near the radar. `detectRFI` flags high power only, and the frequency guard reaches 7 channels, so the other ~110 channels keep unflagged dips at the radar's PRF. It acts like the unmasked gap of FM-9, scattered in phase. | 3 | 2 | 6 | no | **new: development** with the F ingest: the SDR's overflow/clip flag (or the fraction of full-scale samples) blanks all channels for that interval |
| FM-13 | 4 IQ | DC offset, LO leakage, I/Q imbalance (rated 5). **DC tested** (`testFrontEndFaults` 2). The LO sits on a channel edge (the grid is offset half a channel from it), and the dedispersion's sin² taper is zero there. A two's-complement offset (½ LSB on I and Q) appears in the channelizer files at +0.499 LSB² (predicted 0.5) and is gone after dedispersion (1.8·10⁻⁵). Only an LO moved to a channel centre would let it through. I/Q imbalance was not tested: its image at −25 to −40 dB is a mirrored, smeared ghost of the pulsar at 0.01–0.3 % of its amplitude (FM-8's physics), and an image of each RFI source, which is handled like RFI. | 5 | 1 | 5 | n/a | none; keep the LO on a channel edge |
| FM-39 | numerics | Absolute time (rated 5). **Tested** (`testFrontEndFaults` 4). At t0 = 8.3·10⁸ s, the TOAs agree with the t0 = 0 run to 0.115 µs: one ULP of the double that holds them, and nothing more is lost inside the chain. | 5 | 1 | 5 | n/a | none |

Found while reading, with no test needed:
- **The SDR's sample rate fixes the channel width.** `channelizeIQ` needs fs/ChanWidth to be an integer and the band to
  be a whole number of channels. Many SDRs run at 491.52 MS/s, where 3.125 MHz channels are impossible. 3.84 MHz
  channels work there (K = 128, D = 96), with a 104-channel band of 399.36 MHz. This is a configuration step for F, not
  a fault.
- **The full-band path's noise bandwidth** assumes a flat band. This is already in notes §10, "Real data: the receiver
  bandpass is not ideal". The channel path has no such problem, because each channel gets its own level and noise.
- **On a bright pulsar, noise costs less than 1/(1+ε).** Its own noise dominates on the pulse, and the radiometer model
  includes it. This mattered for the test's prediction (test 2), not for the weak target.

### New actions (notes §10 item 17)

| # | Risk | Kind | What | Fits in |
|---|---|---|---|---|
| 1 | FM-42 | measurement, then development | Measure the receiver's total-power stability (an Allan variance on a load or a cold-sky stretch, over 0.1–100 Hz) to get the knee. If it is above ~1 Hz: give `detectPulsar` a noise check per harmonic band, and leave out or whiten the lowest harmonics (a narrow pulse keeps most of its SNR without them). For E4: element LNAs with independent gain noise average down in the beam, common parts (LO, clock, supply) do not. | before F; E4 |
| 2 | FM-41 | development (small) | `ChannelGain` = each channel's mean power (normalized), × the pulsar's spectrum when known, passed by `foldToTOA` to `estimateTOA` and `detectPulsar` | any time |
| 3 | FM-40 | decision | The transport format: sc16, or sc8 with the noise at ≥ 2–4 LSB, and a clip flag | before F |
| 4 | FM-12 | development | The SDR's overflow and clip flags into the F ingest (with FM-9's gaps): blank all channels while the ADC is saturated | F |

## Tests

### `testFrontEndFaults` (voltage chain; ~7 min; ~3 GB in tempdir, deleted)

The frozen seed-43 IQ file (data/test/: the 10 ms pulsar at DM 5, −5 dB, 100 ms at 800 MS/s, LO 1.3 GHz) is changed
the way a real receiver would change it. Then `processIQ` (channel path) and `foldToTOA` run unchanged.

1. **FM-8, the frequency mapping.** The file is re-centred on 1.4 GHz first, so that the mirror maps the band onto
   itself.
   - (a) The re-centred file gives the frozen processing's TOAs within 8.8·10⁻⁷ σ.
   - (b) Mirrored: amplitude 0.047 of the true one, against a predicted 0.043 (+1.4 of its error bars). Phase −0.4983
     turns, against a predicted −0.5000. The prediction is built from each channel's delays: channel j holds the true
     channel 129 − j, so a delay of τ(f₁₂₉₋ⱼ − b) − τ(fⱼ + b). The fit is the mean over channels of ⟨Pⱼ, c_τ⟩/⟨c, c⟩.
     Each channel's own baseline takes the smeared pulse's mean, which is why 0.043 and not the raw overlap of 0.115.
   - (c) LO given 1 MHz high: +15.50 µs, against a predicted +15.77 µs (−1.7 %).
2. **Quantization**, mid-rise, step = rms/level. The SNR factor is checked against a white-error Bussgang model that
   includes the pulse's own noise; the in-band PSD is 1.12× the file's mean.

   | Case | SNR factor | Model | Result |
   |---|---|---|---|
   | 8-bit, rms 4 LSB | 0.9961 | 0.9962 | PASS |
   | 8-bit, rms 1 LSB | 0.9422 | 0.9430 | PASS |
   | 4-bit, rms 3 LSB | 0.9892 | 0.9916 | PASS |
   | 2-bit, rms 1 LSB | 0.8992 | 0.9109 | shown only: the error follows the input there; a power-compression model gave 0.70, also wrong |

   Red. χ² changes by at most 0.002 in every case. Two's complement: +0.4991 LSB² in channels 32 and 33 before
   dedispersion, 1.8·10⁻⁵ after.
3. **FM-9**: 33024 samples (41.28 µs) lost at 45.3 ms.
   - (a) Joined: 2 sub-ints before the drop are unchanged (1.4·10⁻⁷ σ); 3 after it move by −14.0, −14.1 and −14.1σ,
     against a predicted −14.0, −14.1 and −14.0.
   - (b) Zero-filled, without a mask: 1.02σ and χ² +0.124 in the sub-int it hits. Blanked (193 channel samples per
     channel): 0.165σ and χ² −0.002.
4. **FM-39, absolute time.** With t0 = 8.3·10⁸ + 0.123456789 s and TRef moved with it, the TOAs minus t0 differ from the
   t0 = 0 run by up to 0.115 µs. The check is 4 ULP of t0, 0.48 µs; the difference is 0.04 of a sub-int's TOA error. The
   phases agree to 2.1·10⁻⁶ turns. The chain loses no more than the double that holds the TOA.

All four cases passed in one run (13.4 min that time on a busy machine; 6.5 min for cases 1–3 alone before).

### `testFoldFaults` (fast simulator; ~3 min, ~2 GB of memory; no files)

The folds have 128 channels at the frozen power file's noise coefficients, f0 = 100 Hz, −20 dB, 300 sub-ints of 20
turns.

1. **FM-41.** The bandpass is a 6 dB tilt plus a ±1 dB ripple of 20 MHz period (gains 0.41–2.47).
   - Optimal fit: SNR ×0.9190, predicted 0.9166. Equal weighting: ×0.9146, predicted 0.9166.
   - `ChannelGain` = G gives exactly the flat fit (TOAs identical, SNR to 10⁻¹⁴). `ChannelGain` = the measured channel
     power gives ×1.0000.
   - With `ChannelGain` 1, the pulls have rms 0.990 (flat: 0.977) and the TOA error is ×1.087.
2. **FM-42**, white gain noise with R = 1.

   | Quantity | Measured | Predicted |
   |---|---|---|
   | Noise ratio | +0.989 | +1 |
   | T0 variance | +0.729 | +0.667 |
   | TOA variance | +0.739 σ² | +0.667 |
   | Per-channel red. χ² | +0.0075 | +0.0078 (R/128) |

   The variances are within the test's 33 % (4√(2/300)). False alarms at unknown phase: 1 → 8 in 300 sub-ints (Rice:
   0.3 → 9.5). The B0329+54 1/f numbers are in the FM-42 row above.

## Open

- Both tests are Claude's runs; Jasper's runs validate them. Second run, 3 Nov 16:16, on a busy machine, together with
  `testNavFitFaults`: all three pass with the same numbers as before (`testFrontEndFaults` 443 s, `testFoldFaults`
  156 s, `testNavFitFaults` 2.9 s).
- Not tested:
  - FM-12's clipped radar. The finding above stands as a reasoning result.
  - I/Q imbalance (low).
  - Excision on quantized data. The channel voltages are Gaussian again after the filterbank, so no change is expected.
- The register ratings of FM-40 to FM-42 are a draft for Jasper's review, like the 2 Nov ones.
- `runAllTests` runs both new files at the end, with its note for tests not in its list. `testFrontEndFaults` adds
  ~7 min. It could join the slow tests instead (only when named, like `testRFIChain`); that is Jasper's choice.
