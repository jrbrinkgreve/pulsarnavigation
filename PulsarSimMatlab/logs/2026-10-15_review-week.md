# 2026-10-15 – review week, continued: stage 6

Continues `2026-10-14_review-week.md` (stages 5 and 6a). Same method: physics, code walk,
tests, refinements noted (not implemented) in notes §10 item 11 and the docs roadmap.

## E1 noted (Jasper: "note E1, then continue to periodicRFI")

Notes §10 item 11 "Stage 6" (6a `detectRFI`: E1 gain ramp/step and passes 1 vs 2; remark on
the 48 B/sample memory check for the long-data item; links), roadmap row "RFI-detection test".

## Stage 6b – periodic RFI (`periodicRFI`, notes §5.21; off by default since 9 Oct)

Walk-through given (code unchanged since B6 `2082540`): why (locked radar piles up in the fold;
weak sidelobe pulses are missed per pulse but predictable); steps: events → pulses (same time,
neighbouring channels) → groups by strongest channel → period from gaps ÷ 1…8, best phase,
chance n·P(Bin(n−1, q) ≥ m−1), q = 9/P, × candidates < 1e-6, m ≥ 6 → least-squares refit →
mask of every predicted pulse (data-independent, so `blankingWeights` stays exact). Comb
picture drawn in the chat. Why off: 9 Oct false emitters with dense flags (grid of pulse times,
true period never a candidate, no plausibility check); fix a/c/b designed in
`2026-10-09_handover-B-closed.md` open item 1.

Noted (Jasper: "note P1 and P2"; notes §10 item 11 "6b", §9 pointer, roadmap parked-fix entry;
additions to the parked fix, nothing to implement now):
**P1** test cases beyond the fix design's: two emitters in one channel group, a staggered PRF
(→ `unfitted`), more than `MaxEmitters`; **P2** (fix a, by reasoning, not tested) candidate
periods come from pulse pairs, so the pair's second pulse fits by construction, but the chance
formula counts it as a chance hit: factor ≈ (m − 1)/(n(n − 1)q) against the exact count at
that phase (~1.5 for 20 pulses at q ≈ 0.01; ~50 for 10 pulses at P = 2 ms, where both
chances are far below the threshold anyway).

## Stage 6c – blanking (`blankChannels`, notes §5.19)

Walk-through given (code unchanged since B1 `45d95bc`): zero, not subtract (at −54 dB one radar
pulse ≈ 80,000 pulsar pulses); zero, not delete (dedispersion needs the regular grid); per
channel: normalize the mask (clip, sort, merge), stream blocks of 2^20, zero via a difference
array, write a copy (raw files never overwritten); output = `info_chan` with new files + the
normalized mask, which `blankingWeights` uses (one source). Each blank later spreads over the
dedispersion filter's reach (130–310 samples = 31–74 µs per channel). Tests: `testBlankChannels`
(exactness, merge/clip, block size, empty mask, interface with `dedisperseChannels`) – complete
for what it does; only the refusals are unexercised (trivial).

Checked, not candidates: copies of all 128 channels are needed anyway (noise false flags alone,
5.4e-5, hit every channel: ~22 samples per channel per 0.1 s); blanking inside
`dedisperseChannels` would save one copy (~4.3 GB per second of data) but voltage simulation
cannot reach hours anyway (real hardware blanks in the stream).

Proposed, waiting for Jasper: **B-text** (minor, text only) the help of `blankChannels` says the
in-channel sweep is "~60 us per channel at DM 5"; it is 32 µs (1.6 GHz) to 75 µs (1.2 GHz),
as notes §5.19 already has (130–310 samples).
