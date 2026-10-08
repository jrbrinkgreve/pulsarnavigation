# 2026-10-08 – handover: block B (RFI excision) nearly complete — start here next session

Written at the end of 8 Oct for a fresh session. The detailed record of the day is
`2026-10-08_excision.md`; the full reference is `current_project_notes.md` (§5.3, §5.19–5.21,
§7, §9, §10); the docs site (`docs/index.html`) explains the same in readable form.
Previous handover (block A): `2026-10-08_overview.md`.

## First things next session

1. **Open decision (Jasper asked to be reminded): measure first, or fix first?**
   In plain words: inside a channel that a strong wideband transmitter only partly
   covers (the edge of an LTE band), the noise is lopsided — strong in one slice of the
   channel, weak in the rest. Two consequences there: the RFI detector flags too much
   (harmless, a little data lost), and the noise is ~4 % larger than the pipeline
   thinks (reduced χ² 1.04), so TOA error bars come out slightly too small. TOAs
   themselves were fine. Options:
   - **(a) measure first** — build B5e, the realistic L-band scenario (rotating radar,
     LTE, Inmarsat, GNSS L1/L2/E6) and see how many such edge channels occur and how
     big the error-bar problem gets. Fixes nothing yet.
   - **(b) fix now** — let the pipeline measure each channel's real noise level from
     the data instead of assuming it (changes the validated TOA / detection noise
     model → design note first).
   Claude's suggestion: (a) then (b) — the scenario shows how much the fix matters and
   is its test case.
2. Jasper's run of `run('tests/testRfiNoise.m')` (B5d; passes in Claude's run).
3. Push to GitHub when convenient (local main ahead of origin by 7 commits).
4. Still unread: the 7 Oct Q&A (`2026-10-07_qa.md`).

## What was done on 8 Oct (block B)

| Step | What | Commit |
|---|---|---|
| B1 | `blankChannels`: mask → blanked copies of the channel files | 45d95bc |
| B2 | `detectRFI`: multi-scale power threshold per channel, exact thresholds (window power = Σλ·Exp, survival via expm), median/ln2 baseline, guard 10 | 45d95bc |
| B3 | excision in `main.m` (`excision = true`, default, Jasper) | f4181f5 |
| B4 | `runRFITest.m`: each RFI type with / without excision | f4181f5 |
| B5a | radar `RiseTime` (raised-cosine edges, default 0.1 µs = realistic; 0 = old) | 724501d |
| B5b | `runLockedRadar.m`: radar locked to the pulsar; leftover in baseline units, bias vs SNR | 724501d |
| B5b-2 | `detectRFI` frequency guard (strong events ≥ 100× blanked ±`FreqGuard` channels, default 7) | 724501d |
| – | ionosphere noted as an open item (notes §9, §10) | 2500223 |
| B5c | rotating antenna for any RFI type (`ScanPeriod`, `BeamTime`, `BeamWidth`, `SidelobeDB`) | 9b959f0 |
| B6 | `periodicRFI`: find periodic emitters from detections (significance-tested period fit) and blank every predicted pulse; `periodicMask = true` in main | 2082540 |
| B5d | `'noise'` RFI type (band-limited Gaussian, LTE-like); LTE cases in `runRFITest` | 9750461 |
| docs | site created (by a separate session) and updated for B6 / B5d | 9da0f49, a0f9198, 5a6ae63 |

Validated by Jasper's runs: B1–B4, B5a, B5b-2, B5c, B6. Pending: B5d (`testRfiNoise`).

## Where things stand

```
main.m, channel path (default):
channelizeIQ → detectRFI → periodicRFI → blankChannels → dedisperseChannels
→ detectChannels → blankingWeights → foldProfile('DataWeights') → estimateTOA / detectPulsar ('optimal')
```

Defaults: `frontEnd 'channels'`, `weighting 'optimal'`, `excision true`, `periodicMask true`,
detectRFI windows 1–16 samples, PFA 1e-6, guard 10, FreqGuard 7; rfiSource radar
RiseTime 0.1 µs. `rfiOn = false` (RFI studied on purpose with the run scripts).

Key results (−5 dB, seed 43; details notes §7):
- No RFI: excision costs 0.0125 % of the data, TOAs unchanged (0.07σ); periodicRFI finds
  no emitter (the pulsar is dispersed, so it is never taken for periodic RFI).
- Carrier / GNSS / LTE: harmless on the channel path even without excision — the
  'optimal' fit weights channels by 1/level² (a built-in notch). LTE costs 2–3 % SNR.
- Radar: breaks TOAs without excision (ratio ~5, χ² 34); with excision as noise only.
- Radar locked to the pulsar (worst case), bias at −54 dB: excision only 147 µs → with
  frequency guard 7: 2 µs → with the periodic mask: 0.41 µs. A locked radar seen only
  through −35 dB sidelobes: excision alone fails (locks on the radar); with the periodic
  mask leftover 1.2e-10 of the baseline, 0.001 µs.
- Insight: express RFI results in units that scale (leftover / noise baseline vs the
  pulse height ρ); −5 dB TOA tests alone understate RFI by ~80,000× relative to −54 dB.

## Next steps (after the decision above)

- B5e realistic scenario and/or per-channel noise calibration (see above).
- Later in B: staggered-PRF radars in `periodicRFI` (v1: fixed PRF only); a frequency
  guard scaled with the event's strength; `periodicRFI` for long data (it compares
  pulse pairs and writes one row per pulse and channel); spectral kurtosis /
  whole-channel flags only if needed.
- Not modelled yet, possibly as important as excision for real hardware: receiver
  dynamic range (a strong radar can clip a few-bit ADC and spread over the band);
  ionosphere (delay ~3–70 ns at 1.4 GHz, Faraday rotation, scintillation; notes §9).
- Then C (fast simulator), D (barycentric phase, Doppler, clock), E.

## How to run

Tests (from `PulsarSimMatlab`), all passing on 8 Oct:
```matlab
run('tests/testBlankChannels.m')   % ~10 s
run('tests/testDetectRFI.m')       % ~33 s, writes ~2.5 GB to tempdir
run('tests/testRadarEdges.m')      % ~5 s
run('tests/testRfiGating.m')       % ~10 s
run('tests/testPeriodicRFI.m')     % ~35 s
run('tests/testRfiNoise.m')        % ~5 s  (Jasper's run pending)
```
(plus the block A tests listed in `2026-10-08_overview.md`). The RFI regression
reference is `data/mc/rfiRef_preB5a.mat` (`tests/makeRfiReference.m`, code of f4181f5).

Scripts: `runRFITest` (cases 1–13; `runCases` picks a subset; needs `frontEnd = 'channels'`;
~45 s per case; files in `data/rfi`), `runLockedRadar` (cases 1–12, `runCases`; ~9 min
for all; files in `data/locked`).

## Working conventions (learned / agreed; keep)

- CLAUDE.md: explain the physics and the planned change first, wait for OK; one logical
  change per edit; Jasper runs MATLAB himself (Claude may run tests per unit, but
  "validated" only after Jasper's run); ground-truth rule; dated log after each
  session; **update `docs/` alongside the code** (checklist on `docs/maintenance.html`;
  regenerate `docs/api.html` with `python3 docs/tools/make_api.py` after commits that
  touch functions or scripts).
- Plain labels A–E and B1…B6; explain choices in plain words (Jasper found an
  abbreviated choice unclear on 8 Oct).
- One new function at a time; tests as functions; never run tests in parallel with
  Jasper's (shared `data/chan/`).
- When Jasper has local edits in a file Claude also changed, commit only Claude's lines.
- Commit straight to main (solo project); commit when Jasper says so.
