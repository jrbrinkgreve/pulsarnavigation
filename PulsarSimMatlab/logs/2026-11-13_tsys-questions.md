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

## LNA shortlist for the 70 K design point (web search, 13 Nov)

Jasper: "for this design point, can you look online for some suitable LNA's?", then "Yes, log it".
Not a decision; hardware choice is open.

**Requirements derived.** Band 1.2–1.6 GHz, flat over 400 MHz. NF ≤ 0.4–0.5 dB at 1.4 GHz
(70 K − sky ~4 − atmosphere ~3 − ground ~10 − input loss ~15 leaves ~35–40 K for the electronics).
Gain ≥ 35–40 dB in two stages before the SDR (bare SDR NF 6–8 dB ≈ 865–1500 K: +9 K behind 20 dB,
< 0.2 K behind 40 dB). OP1dB ≳ +17 dBm, OIP3 ≳ +33 dBm (realistic-scenario radar, INR +40 dB:
−54 dBm at the LNA input, ≈ −14 dBm after 40 dB). Single supply, small, cheap (×18), stable into a
non-50 Ω array element.

**Candidates.** NFs are datasheet values at the frequency given (eval board or de-embedded); few are
specified at 1.4 GHz; expect +0.1–0.2 dB on a real board. Values from distributor pages and datasheet
summaries, not all checked against the current datasheet revision.

| Part | Type | Band | NF | Gain | OP1dB | OIP3 | Role |
|---|---|---|---|---|---|---|---|
| Guerrilla RF GRF2071 | MMIC, external match | 0.7–2.7 GHz | 0.35 dB (1.9 GHz) | 19 dB | +21 dBm | +36–38 dBm | stage 1, first choice |
| Qorvo TQP3M9037 | MMIC | 0.7–6 GHz (older rev 1.5–2.7) | 0.4 typ / 0.65 max (1.95 GHz) | ~20 dB | +20 dBm | +35 dBm | stage 1, second source |
| Qorvo QPL9057 | MMIC | 0.6–4.2 GHz | 0.54 dB (3.5 GHz) | 22.8 dB | +17 dBm | +32 dBm | stage 2 |
| Mini-Circuits PMA2-162LN+ / ZX60-P162LN+ | MMIC / coax module | 0.7–1.6 GHz | 0.5 (1 GHz), ~0.8 dB (1.6 GHz) | ~20–24 dB | ~+19 dBm | ~+30 dBm | prototype, first T_sys measurements |
| Mini-Circuits TAV-541+ | discrete E-PHEMT | 0.45–6 GHz | 0.4 (0.9 GHz), 0.5 dB (2 GHz) | design | — | high | own design, noise-matched to the element |

Rejected: ADL5523 (0.8 dB at 0.9 GHz, ~1.0 at 1.95), MACOM MAAM12000 (1.35 dB), GRF2070 (stops at
1.5 GHz), NXP GPS1502 (1164–1299 MHz only), Fairview FMAM1028 (0.7 dB coax, no gain over the ZX60),
Skyworks SKY67151 (< 0.5 dB, but sources disagree on its linearity).

**Budget check** (GRF2071 → QPL9057 → filter/cable/SDR, 0.2 dB loss before the LNA): input loss ~14 K
+ GRF2071 ~24 K + QPL9057 ~0.5 K + rest < 0.2 K → T_rx ≈ 39 K; + sky/atmosphere/ground ~17 K →
T_sys ≈ 56 K, ~15 K under 70 K. Reference: the uncooled Apertif PAF reaches T_sys ≈ 70 K
(arXiv:2109.14234; LNA + second stage 28 K).

**Design notes.** LNA at the element (input loss is the largest term after the LNA). Band-pass filter
after stage 1 (1 dB in front costs ~75 K), unless a strong nearby radar forces it forward. Same part
and layout on all 18 channels, with calibration (gain/phase match and its temperature drift matter for
beamforming). Check stability with the real element, not only into 50 Ω (mutual coupling).

Sources: Guerrilla RF GRF2071 via rfmw.com and richardsonrfpd.com; Qorvo product pages and datasheets
(TQP3M9037, QPL9057); Mini-Circuits datasheets (PMA2-162LN+, ZX60-P162LN+, TAV-541+); Skyworks
SKY67151 product page; Analog Devices ADL5523 datasheet; MACOM MAAM12000, Fairview FMAM1028 and NXP
GPS1502 datasheets/pages; Apertif arXiv:2109.14234.
