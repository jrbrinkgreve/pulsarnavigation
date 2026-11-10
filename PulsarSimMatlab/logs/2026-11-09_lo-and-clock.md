# 9 Nov 2026: questions on the LO and the receiver clock (no code change)

Jasper: "is it true that the drift of an LO manifests itself as a clock bias later on", then "how would clock bias manifest itself
then?", then "Yes, add it" (a correction to item 5 of the 4 Nov proposal). No code, docs or tests changed.

## The LO
- An LO error shifts every frequency by df but delays nothing. Detection ignores the shift and the LO's phase.
- It reaches the TOA only through the channel frequencies: dedispersion removes slightly the wrong delay, dt = 2 t_DM(f) df / f, with
  t_DM(f) the dispersion delay at f.
- `testFrontEndFaults` 1c (3 Nov) already checks this: the LO given 1 MHz high at DM 5 moved the TOA +15.50 us against +15.77 us predicted.
  At 2 ppm, about 2.8 kHz at 1.4 GHz, that is 0.04 us at DM 5 and about 0.2 us at DM 25.
- It does not accumulate, it scales with DM, and a constant part goes into navFit's per-pulsar constants c_k. It becomes a clock bias only
  when the LO and the sample clock share one reference, and then through the sample clock. The mixer is not chosen yet. The simulator's LO
  is exact, which is fine at this size.

## The clock
- Every TOA gets the clock error b(t), the same for all pulsars at that moment. A position error differs per pulsar and follows the daily
  sine, which is how the fit separates the two.
- In navFit's "self" mode the constant part goes into the c_k, the rate is one column, and the wander is correlated noise in the GLS. A
  daily clock line (FM-11) and pulsars observed at different times (hence D4c's simultaneous beams) are where clock and position mix.
- In the fold the rate smears the pulse, 240 us in a 120 s sub-int at 2 ppm, hence makePredictor 'Clock'.
- Longitude: first answered as "right for UT1, not for the receiver clock", from a toy model folded against the true site. A two-pass toy
  model then showed that a clock offset does move the first fit against the geocentre by 286 m per second, through the TOA times, and that
  a second pass removes it. The correction under item 5 of `logs/2026-11-04_static-redesign.md` gives the numbers.

## Not committed
This note, the correction in the 4 Nov log, and the 5 Nov real-time section of `logs/2026-11-05_static-stage7.md`.
