# 2026-10-14 – review week, continued

Continues `2026-10-13_review-week.md` (stages 2–4). Same method: physics, code walk, tests,
refinements noted in notes §10 item 11 and the docs roadmap.

## Stage 2 again: Nk and Nf in `makeDispersionKernel`, from scratch (Jasper's question)

Explained in five steps with one figure (two bars, Nf = Nk vs Nf = 2·Nk):
1. h = the answer to a click (a 6.3 ms descending chirp); Nk = its length
   (10,000 + 25,209,771 + 10,000 + 1 = 25,229,772).
2. h is built from its spectrum H(f) (W, φ) with an inverse FFT of Nf points.
3. An FFT's time axis is a circle of Nf samples: anything outside 0…Nf−1 wraps and adds on top.
4. Nf ≥ Nk or the chirp folds onto itself; the band-edge tails before 0 and after Nk need an
   empty stretch → Nf ≥ 2·Nk (2^26 = 67.1 M), where both tails land and are cut.
5. `leak` = energy in the cut stretch (3.5e-13); trustworthy because nothing wraps back unseen.
Frequency view: phase step per bin = d/Nf turn; > ½ turn is indistinguishable from a negative
delay → again Nf > 2·Nk (0.376 turn today). Not to be confused with `Nfft` (applies h to data
blocks: Nfft ≥ n + Nk − 1).

## I2 fixed (Jasper: "fix I2 now")

Notes §5.4 numbers updated to the current set-up (fLO 1.3 GHz, fs 800 MHz: D = 5, cutoff
360 MHz, flat to ±355 MHz, stopband from ≈365 < 400 MHz, band −100…+300 MHz, mirror +1.1…+1.5 GHz,
r = 0.325 = 13/40, 800,000,001 samples per s); the 30 Sep set-up named as the one that would not
fit. Text only. Docs already correct (`receiver.html` callouts "LO at 1.3 GHz", "800 MHz output");
roadmap row and notes §10 item 11 marked done.

## Stage 5 – channelizer (`channelizeIQ`, notes §5.15)

Walk-through given (code unchanged since `f5ddaf5`): 128 channels of 3.125 MHz, K = 256-point
FFT, D = 192 (oversampling 4/3, 4.17 MHz per channel), half-slot offset β = ½ because fLO lies
on a channel edge; per output time: frame of P·K samples × β-modulated prototype, fold to K,
IFFT, phase correction e^(−i2π·k·c/K); edges outside `fullySupported`. Polyphase diagram for one
output time drawn in the chat. Tests: `testChannelizeIQ` (tones to 2.7e-7, passband 1.2e-4,
stopband −80.2 dB, noise power 1.0007, Parseval on real data 1.00007, 6 blocks) – only on
today's grid.

Noted (notes §10 item 11 "Stage 5", roadmap row "Channelizer test"): **C1** extend the test
(β = 0 grid, oversampling 2, `BlockOut` bit for bit, the refusals). Links only.

## Stage 4 completeness check (Jasper: "did you write down all things you found?")

Compared the stage-4 walk-through with what was recorded. Three gaps, fixed:
- I3 was in §10 and the roadmap but marked "for §9" without being added: now in §9 "Receiver"
  and on `docs/limits.html` (Receiver bullet).
- The legacy `'Gain','unity'` option was called "legacy" in the code table but not noted: now
  **I4** (remove; no caller passes `'Gain'`; keep `info.gainFactor` for `expectedPowerModel`).
- The 256 ns filter transient at both file ends: written down as checked and negligible.

## Stage 6a – RFI detection (`detectRFI`, notes §5.20)

Walk-through given (code unchanged since B6 `2082540`): why before dedispersion (a 1 µs pulse
would be smeared over the channel's own 32–75 µs sweep at DM 5); exponential power, thresholds
13.82 … 2.69 for windows of 1–16 samples at PFA 1e-6, exact tail via eigenvalues (ρ(1) = 0.05);
median / ln 2 baseline per 16 ms block (66,667 samples, ~0.6 % per block); time guard ±10
samples (prototype half-length 2.4 µs); frequency guard ±7 channels for events ≥ 100×; what it
deliberately does not catch. Mask picture (channel × time) drawn in the chat.

Noted 15 Oct (notes §10 item 11 "Stage 6", roadmap row "RFI-detection test"): **E1** extend `testDetectRFI` with (a) a gain ramp / step (the
reason for 16 ms blocks; the simulation has no gain drift; estimate: a baseline 1.6 % off raises
false flags ~25 % for 1-sample windows, ~60 % for 16-sample windows) and (b) Passes 1 vs 2 with
dense RFI. Remark for the long-data item: the memory check counts 48 B per sample; by count
(not measured) the 1-sample window step makes several file-length temporaries, so the peak
is likely higher.
