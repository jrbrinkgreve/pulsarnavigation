# 13 Nov 2026 — general questions; T_sys 70 K, and whether it includes the sky

No code or docs change. Jasper asked general questions (why oscillator stability is given in ppm, offset
vs. fluctuation, a list of common oscillators, the Allan deviation, how system temperatures add, and
where the 70 K comes from). Only the last answer touches the project; it is recorded here.

## Where 70 K comes from

- It is not derived. It is the midpoint of Jasper's hardware spec of 6 Oct 2026: T_sys 60–80 K
  (`current_project_notes.md` §11, "Reference scenario – draft"; `2026-10-06_phase-e.md`).
- Design points: nominal 0.61 m², 70 K → SEFD 3.2·10⁵ Jy; worst case 0.3 m², 80 K → 7.4·10⁵ Jy
  (`docs/scenario.html` R1).
- Where the code uses it: the scenario's ρ = S_peak/SEFD (fast-simulator signal levels), and the RFI link budget in
  `pipelineParams.m:66`: kT_sys·B = 1.38·10⁻²³ · 70 · 4·10⁸ = 3.9·10⁻¹³ W = −124.1 dBW (checked). The
  voltage-chain default `snrDB = -5` is set directly, not from T_sys. No T_sky term exists in the code.

## Plausibility (estimate, uncooled L-band, LNA at each element)

| Contribution | Estimate |
|---|---|
| Sky (CMB 2.7 K + Galactic, off the plane, 1.4 GHz) | 3–5 K |
| Atmosphere (zenith → low elevation) | ~2–8 K |
| Ground through sidelobes / back lobes of a small array | ~5–20 K |
| Loss before the LNA (0.2–0.5 dB at 290 K) | ~15–35 K |
| LNA (NF 0.3–0.6 dB) | ~20–45 K |
| **Total** | **~45–110 K** |

60–80 K is mid-range for this budget. Caveats: it needs an LNA at each element (a bare SDR front end, NF
~5–8 dB, is ~600–1500 K); the worst case allows only +14 % on T_sys (most of its factor 2.3 in SEFD is the
area); elevation-dependent ground pickup and element-to-element noise coupling are not in it; the Sun is
not in it (FM-2).

## Inconsistency: is the sky inside the 60–80 K?

- `current_project_notes.md:4178` (decision of 6 Oct): "T_sys = T_rx (60–80 K, assumed to exclude sky) +
  T_sky(direction)", with T_sky per pulsar from a sky map.
- `docs/scenario.html` R1 and `2026-11-02_risk-register.md:148`: T_sys 70–80 K *includes* the
  beam-averaged sky (~3–5 K off the plane, ~5–20 K for plane sources).
- Effect: off the plane ≤ 5 K (< 7 %); for plane sources (Vela, B1749−28) +5–20 K on T_sys, up to ~30 %,
  so up to ~1.7× in time per TOA (t ∝ T_sys²).
- Open for Jasper: which reading holds. If the 6 Oct one: `scenario.html` R1 and the risk register's line
  need correcting, and the per-pulsar T_sky belongs with FM-2 (Tsys(t) in the simulator).

## Correction to this session's own answers

Early answers took "400 MHz" as the centre frequency; it is the bandwidth (band 1.2–1.6 GHz, LO in band
~1.3–1.4 GHz). So ±2 ppm is ±2.8 kHz at the LO (not ±800 Hz), T_sky at 1.4 GHz is 3–5 K off the plane
(not 20–50 K), and at L-band the receiver, not the sky, dominates T_sys. Corrected in the chat.
