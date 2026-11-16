# 16 Nov 2026 — design note: what a direct-sampling SDR would change

Jasper: "how would the code change if we use a direct sampling SDR?", then "Yes, note it".
A design note, not a decision: the receiver is still antenna → mixer → SDR (memory of 21 Oct; the
mixer type is not known). No code change. Read for this note: `pipelineParams.m`,
`applyIQmodulation.m`, `channelizeIQ.m`, `processIQ.m`, `addNoiseAndRFI.m`, `makeClock.m`,
`tests/testFrontEndFaults.m`; nothing was run.

## Main point

The processing (`processIQ`, `foldToTOA`) does not change. The simulator is already an ideal
direct-sampling receiver: `addNoiseAndRFI` writes a real RF stream at `f_in` = 4 GS/s, and
`applyIQmodulation` mixes it with an exact complex exponential, low-pass filters and decimates,
which is the digital downconverter (NCO, filter, decimator) of a direct-sampling SDR. A direct-sampling
receiver still delivers an IQ file with an equivalent `fLO` and an `fs`, all that `processIQ` takes.
What is missing is on the truth side: a real ADC at RF (finite bits, clipping, aliasing).

## Two ways to sample 1.2–1.6 GHz

Bandpass sampling in Nyquist zone n needs 2 f_H / n ≤ fADC ≤ 2 f_L / (n − 1).

| Option | fADC | Band after sampling | Remarks |
|---|---|---|---|
| Zone 1 | ≥ ~3.5 GS/s | in place | simplest; RFSoC-class ADC; everything up to fADC/2 reaches the ADC (headroom) |
| Zone 3 | 1.067–1.2 GS/s | at 1.1 GS/s: 0.1–0.5 GHz, not inverted | guard bands 1.1–1.2 GHz (100 MHz) and 1.6–1.65 GHz (50 MHz): a steep band-pass in front; the ADC's analog bandwidth must reach 1.6 GHz |
| Zone 4 | exactly 0.8 GS/s | 0.4–0 GHz, inverted | no guard band: not practical; an even zone mirrors the band (FM-8) |

With undersampling every other zone inside the ADC's analog bandwidth folds onto the band. Example
at 1.1 GS/s: GSM at 900 MHz lands on 200 MHz, the same place as 1.3 GHz.

## Stage by stage

| Stage | Change |
|---|---|
| `generatePulsarSignal`, `applyDispersionStream` | none |
| `addNoiseAndRFI` | undersampling: out-of-band RFI (GSM 900, LTE 800/1800, Wi-Fi 2.4) and out-of-band noise, because they fold in. Sources up to `f_in`/2 = 2 GHz are accepted today; Wi-Fi at 2.4 GHz needs a higher `f_in`. Zone 1: they count for ADC headroom only |
| new truth-side stage (e.g. `sampleADC`) | analog band-pass (its rejection sets what folds in), resampling `f_in` → fADC, quantization with a full scale and clipping at RF. This is the receiver stage FM-12 already asks for ("ADC bits, clipping and LNA compression", `pipelineParams.m:69`: "would clip a real ADC (not modelled)"). Sampling jitter optional: −20 log10(2π · 1.4 GHz · 100 fs) = −61 dBc, negligible for the pulsar, larger with strong RFI |
| `f_in` | an integer multiple of fADC: 4 GHz for fADC = 4 GS/s (no resampling at all), 4.4 GHz for 1.1 GS/s |
| `applyIQmodulation` | is the DDC; same operation. Undersampling: input at fADC, NCO at the LO's alias (1.4 − 1.1 = 0.3 GHz), an inversion switch for even zones (none today, FM-8) |
| IQ rate `fs` | Zone 1 at 4 GS/s: 800 MS/s (D = 5) and the 1.3 GHz LO stay as they are. At 1.1 GS/s: 550 MS/s (D = 2); 550/3.125 = 176 = K and 176/(4/3) = 132 = D, so the 3.125 MHz channels stay, but ±275 MHz holds the band only with the LO at 1.4 GHz (a digital NCO has no DC or LO leakage, the reason for the off-centre LO) |
| `processIQ` … `foldToTOA` | none |
| `makeClock` | none: one oscillator clocks the ADC and the NCO, so a rate y shifts the RF mapping by y · 1.4 GHz (2.8 kHz at 2 ppm), as an analog LO on the same reference does |
| risks and tests | FM-13 (I/Q imbalance, DC spike, LO leakage of a zero-IF mixer) no longer applies. `testFrontEndFaults` 2 quantizes the IQ file (the transport, FM-40); an RF ADC quantizes before the DDC and clips on the total input power, which no test covers (FM-12). New: ADC spurs (harmonics and, for interleaved ADCs, interleaving spurs fold in), leakage through the band-pass |
| `pipelineParams` / `rx` | truth side: `fADC`, bits, full scale, Nyquist zone, band-pass edges and rejection. The observer's `rx` keeps fLow, fHigh, and the IQ info its equivalent fLO and fs |

## Effect on the design point

- Bits: a weak pulsar needs few; the rest is headroom for RFI. Gain so that the thermal noise spans
  several LSB (FM-40's numbers apply per real sample), with the ADC's own noise well below k T_sys B.
- T_sys with undersampling: each folded zone the LNA still amplifies adds about another T_sys if it
  is not filtered (zones 1–4 at 1.1 GS/s: up to ×4); 30 dB of rejection leaves 0.1 % per zone.
- Hardware: 18 coherent channels at ≥ 1.1 GS/s is one or two RFSoC-type boards on one clock.

## Suggested order (if it is ever built)

1. Zone 1 at fADC = `f_in` = 4 GS/s: only the ADC stage (quantization and clipping at RF) between
   `addNoiseAndRFI` and `applyIQmodulation`. It is FM-12's stage, useful for any receiver.
2. Undersampling only if the hardware goes that way: aliasing, out-of-band RFI, the inversion switch.

## Corrections to the chat answer of this session

- It said the zone-1 variant needs the LO moved to 1.4 GHz. It does not: at 4 GS/s the present
  D = 5, 800 MS/s and 1.3 GHz LO stay. The LO at band centre is needed only for 550 MS/s.
- It said a weak band-pass "can double T_sys" with undersampling. Unfiltered it is up to one more
  T_sys per folded zone, so more than double at 1.1 GS/s.
- It presented the ADC stage as new. The clipping part is FM-12 (13th of the risk register's
  actions) and quantization at IQ level is tested (FM-40); new is only doing both at RF.
- It said the mixer faults drop out of `testFrontEndFaults`. The test has no I/Q-imbalance case;
  what drops out is risk FM-13. The DC check of case 2 still holds for a two's-complement transport.

## Docs

`docs/roadmap.html`: one bullet under "Still open, in no block yet"; `docs/maintenance.html`:
changelog row.
