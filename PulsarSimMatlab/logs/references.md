# Reference register of the docs site

The source of truth for every literature reference the docs cite. `docs/references.html` is built from this file by
hand. A reference that is added, replaced or dropped is changed here first.

Started 24 Nov 2026 (Jasper: "find, and double-check those found references for the information on the pages … add the
citation as [1] etc, with all citations used in the docs in a new page. Also keep track yourself of these references").
Session log: `logs/2026-11-24_docs-references.md`.

## Rules

- **Scope** (Jasper, 24 Nov): every named method, theorem or data source, and every number taken from the literature.
  Not general textbook statements.
- **Numbers.** `[n]` was given on 24 Nov 2026, when `references.html` was first built, in order of first appearance
  along the sidebar (table "Numbers" below). Later references are appended. A number never changes once given.
- **States.**
  - **A, checked in the source:** the text was opened and the statement the docs rely on was found in it.
  - **B, record checked, content through a second source:** the bibliographic record is confirmed; the content was
    confirmed in another source that quotes or cites it.
  - **C, record checked, named method only:** the record is confirmed and the docs cite it only as the origin of a
    named method or data set; the text was not opened.
  - **D, not confirmed:** record or content still to be found. Not cited as checked.
- **Records** were confirmed in Crossref (by DOI) or the arXiv listing unless another route is named.
- **Mismatches** between a source and a page are listed in the session log; the pages are not changed for them
  (Jasper, 24 Nov).

## Numbers (fixed 24 Nov 2026)

| No. | id | No. | id |
|---|---|---|---|
| [1] | `lorimerkramer2005` | [26] | `3gpp_band32` |
| [2] | `hankinsrickett1975` | [27] | `maron2000` |
| [3] | `dicke1946` | [28] | `isgps200` |
| [4] | `isserlis1918` | [29] | `smith_sasp` |
| [5] | `reed1962` | [30] | `barnes1971` |
| [6] | `neymanpearson1933` | [31] | `allan1966` |
| [7] | `taylor1992` | [32] | `vandierendonck1984` |
| [8] | `rice1945` | [33] | `harris2003` |
| [9] | `rickett1990` | [34] | `amari1997` |
| [10] | `sosafiscella2021` | [35] | `naif_de440s` |
| [11] | `keith2013` | [36] | `fairhead1990` |
| [12] | `gwinn2006` | [37] | `hobbs2006` |
| [13] | `kaplan2005` | [38] | `shang2017` |
| [14] | `bulletinc` | [39] | `massey1951` |
| [15] | `iau2006b3` | [40] | `white1980` |
| [16] | `edwards2006` | [41] | `kalman1960` |
| [17] | `shapiro1964` | [42] | `rts1965` |
| [18] | `park2021` | [43] | `lange2001` |
| [19] | `jplapprox` | [44] | `chen2011` |
| [20] | `capitaine2003` | [45] | `manchester2005` |
| [21] | `stockham1966` | [46] | `manchester2001` |
| [22] | `matsumoto1998` | [47] | `brown2025` |
| [23] | `lecuyer1999` | [48] | `wang2005` |
| [24] | `kaiserschafer1980` | [49] | `dai2015` |
| [25] | `bellanger1976` | [50] | `lorimer1998` |

Appended on 24 Nov 2026, after the rollout, when nine of the open sources were confirmed:

| No. | id | No. | id | No. | id |
|---|---|---|---|---|---|
| [51] | `rao1945` | [54] | `wgs84` | [57] | `itu_res222` |
| [52] | `cramer1946` | [55] | `meeus1998` | [58] | `kaiser1974` |
| [53] | `kay1993` | [56] | `horizons` | [59] | `huber1967` |

Next free number: [60].

Not numbered, because no page cites them directly: `kulkarni2020` (named inside [1]), `pennucci2014` (named inside [7]),
`rainal1988`, `smirnova2006`, `gangadhara2001`, `price2016`. They stay in the tables below as helpers.

All 50 are cited since the rollout of 24 Nov 2026: 116 citations on 19 pages. Each entry on `references.html` lists the
pages that cite it ("Cited in"); that list is the record of where a citation stands, the "Cited for" column below is the
plan it came from. The first three numbers were swapped once on that day, before anything was committed, so that the
physics page cites [1], [2], [3] in order.

## A · Checked in the source (17)

| id | Reference | What was confirmed, and where | Cited for (page: claim) |
|---|---|---|---|
| `rice1945` | S. O. Rice, "Mathematical analysis of random noise", Bell Syst. Tech. J. 23(3), 282–332 (1944) and 24(1), 46–156 (1945). doi:10.1002/j.1538-7305.1944.tb00874.x, 10.1002/j.1538-7305.1945.tb00453.x | Text of the 1945 instalment (Internet Archive, `bstj24-1-46`): Part III, Section 3.3 "Expected number of zeros per second", eq. (3.3-11); eq. (3.3-14): passages through the value I with positive slope per second = exp(−I²/2ψ₀) × ½ × zeros per second. | detection D3 (Eq. D.3), physics P6, formulas, glossary: Rice's formula |
| `vandierendonck1984` | A. J. Van Dierendonck, J. B. McGraw, R. G. Brown, "Relationship between Allan variances and Kalman filter parameters", Proc. 16th Annual PTTI Meeting, 273–293 (1984). NASA N85-29238, https://ntrs.nasa.gov/citations/19850020926 | Scan read page by page. Eqs. (6)–(8), p. 276: σ_y² = (h₀/2)τ⁻¹, (2 ln 2)h₋₁, ((2π)²/6)h₋₂τ. Eq. (60), p. 284: Cov[x(Δt), ȳ(Δt)], which gives q₁ = h₀/2 and q₂ = 2π²h₋₂. Figure 3, p. 288: "a 3 parts in 10¹¹ crystal oscillator", h = 1.8e-21, 6.492e-22, 1.37e-24, about 11 ns at 200 s. Own arithmetic with eq. (60): 8.50 ns at 200 s without flicker. | receiver 4a.2 (preset, formulas, test 4), formulas (Allan variance, Kalman clock noise), glossary |
| `kaplan2005` | G. H. Kaplan, "The IAU Resolutions on Astronomical Reference Systems, Time Scales, and Earth Rotation Models", U.S. Naval Observatory Circular 179 (2005). arXiv:astro-ph/0602086 | Eq. (2.6): TDB ≈ TT + 0.001657 sin(628.3076T + 6.2401) + 0.000022 sin(575.3385T + 4.2970) + 0.000014 sin(1256.6152T + 6.1969) + …; "maximum error … is about 10 µs from 1600 to 2200". TT = TAI + 32.184 s. Quotes IAU 2000 Resolution B1.5: L_C = 1.48082686741·10⁻⁸. | fold 7a.3–7a.4 (Einstein term, test 4), localisation L3 and L6, formulas |
| `shang2017` | L. H. Shang et al., "Investigating the multifrequency pulse profiles of PSRs B0329+54 and B1642−03 in an inverse Compton scattering model", MNRAS 468, 4389–4398 (2017). arXiv:1703.03582 | Eq. (1): Gaussians with w_i the standard deviation. Table 2, 1410 MHz: (position, intensity, w) = (−11.427, 0.159, 1.479), (−1.123, 0.439, 3.281), (0, 0.585, 1.055), (7.719, 0.152, 2.751), (9.482, 0.17, 0.847), equal to `runScenario.m:118`. The 1410 MHz profile is archival, not their own observation. | toa 8.7, scenario R6–R7, simulator: B0329+54's five components |
| `sosafiscella2021` | V. Sosa Fiscella et al., "PSR J0437−4715: the Argentine Institute of Radioastronomy 2019–2020 observational campaign", ApJ 908, 158 (2021). arXiv:2010.00010 | §3: "We scale the scintillation parameters given at the frequency of 1.5 GHz by [Keith et al. 2013] to match our observations centered at 1.4 GHz and obtain the scintillation bandwidth Δν_d = 740 MHz and scintillation timescale Δt_d = 2290 s." | parameters §4 (`scint`), physics (scintillation), limits. On the pages now as "IAR 2020 (arXiv:2010.00010)" |
| `keith2013` | M. J. Keith et al., "Measurement and correction of variations in interstellar dispersion in high-precision pulsar timing", MNRAS 429, 2161–2174 (2013). arXiv:1211.5887 | Table 2 ("at a reference wavelength of 20 cm"): J0437−4715 ν₀ = 1000 MHz, τ_diff = 2486 s. Own arithmetic: scaled from 1.5 to 1.4 GHz with ν^4.4 and ν^1.2 this is 738 MHz and 2288 s. | same places as `sosafiscella2021` |
| `wang2005` | N. Wang et al., "Long-term scintillation observations of five pulsars at 1540 MHz", MNRAS 358, 270–282 (2005). arXiv:astro-ph/0501218 | Table 2: B0329+54 Δt_d = 16.9(5) min, Δν_d = 14(1) MHz at 1540 MHz; eq. (14) defines them as the 1/e and half-maximum half-widths; bandwidth index about 3.6. | scenario R6: B0329+54's scintles (see the mismatch list) |
| `gwinn2006` | C. R. Gwinn, C. Hirano, S. Boldyrev, "Interstellar scintillation of PSR J0437−4715 on two scales", A&A 453, 595–599 (2006). arXiv:astro-ph/0204367 | Conclusions (327 MHz): broad scale 16 (+8, −3) MHz and 1000 (+350, −250) s; narrow scale 0.5 ± 0.1 MHz and 90 ± 20 s. Own arithmetic: × (1400/327)^1.2 gives 5730 s and 516 s, the page's "~5700 s and ~500 s". | physics (scintillation) |
| `chen2011` | J. L. Chen et al., "Long-term monitoring of mode switching for PSR B0329+54", ApJ 741, 48 (2011). arXiv:1107.4676 | Abstract: "84.9% of the time … in the normal mode and 15.1% … in the abnormal mode". | scenario R6, limits: abnormal mode about 15 % |
| `edwards2006` | R. T. Edwards, G. B. Hobbs, R. N. Manchester, "TEMPO2, a new pulsar timing package – II. The timing model and precision estimates", MNRAS 372, 1549–1574 (2006). arXiv:astro-ph/0607664 | Eq. (8): the delay terms. §2.5.1 Einstein, §2.5.3 Roemer (eq. 13), §2.5.5 Shapiro (eq. 33). §2.5.2: ionosphere 5–100 TECU, 7–130 ns·(f/1 GHz)⁻², which is 3.6–66 ns at 1.4 GHz. | fold 7a.1 and 7a.3, localisation L4–L6, limits (ionosphere), formulas |
| `lange2001` | Ch. Lange et al., "Precision timing measurements of PSR J1012+5307", MNRAS 326, 274–282 (2001). arXiv:astro-ph/0102309 | Appendix A: "This new model has been implemented in the tempo timing software as binary model ELL1." | limits, navigation 10.5: ELL1 |
| `kulkarni2020` | S. R. Kulkarni, "Dispersion measure: confusion, constants & clarity", arXiv:2007.02886 (2020), preprint | Quotes the Handbook's a = 4.148808(3) GHz² cm³ pc⁻¹ ms; from the SI constants a = 4.1488064239(11) (eq. 8); the timing convention DM = 2.41000·10⁻¹⁶ DC. | second source for `lorimerkramer2005`; the convention note |
| `jplapprox` | JPL Solar System Dynamics, "Approximate Positions of the Planets" (Keplerian elements after E. M. Standish), https://ssd.jpl.nasa.gov/planets/approx_pos.html | Table 1 "valid for the time-interval 1800 AD – 2050 AD"; nominal errors of the Earth–Moon barycentre: 20″, 8″ and 6000 km. 20″ at 1 au is 14,500 km (own arithmetic). The page does not name Standish; `delayModel.m:50` does. | fold 7a.3, localisation L6: the analytic orbit |
| `manchester2001` | R. N. Manchester et al., "The Parkes multi-beam pulsar survey – I", MNRAS 328, 17 (2001). arXiv:astro-ph/0106522 | §2 and Table 1: system temperature 21 K, gain 0.735 K/Jy, "an equivalent system flux density of 28.6 Jy" for the central beam. | scenario R1: "Parkes has an SEFD of about 30 Jy" |
| `dai2015` | S. Dai et al., "A study of multifrequency polarization pulse profiles of millisecond pulsars", MNRAS 449, 3223–3262 (2015). arXiv:1503.01841 | Read first hand on 21 Nov (`research_notes/…/pulsar_set_catalogue_values.md:167`): J0437−4715 W50 = 8.9° (0.142 ms), W10 = 63.4° (1.01 ms) at 1400 MHz; S1400 = 150.2 ± 1.6 mJy. Today only the record and the width definitions (Table 3) were re-read. | scenario R2 and R7.3: J0437's profile and widths |
| `naif_de440s` | NAIF, generic kernels, `de440s.bsp` and `aa_summaries.txt`, https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/ | Quoted verbatim by a fetch: coverage 1849 DEC 26 to 2150 JAN 22; EARTH (399) and MOON (301) relative to EARTH BARYCENTER (3), that and SUN (10) relative to the barycentre (0); listed size "31M". | fold 7a.4: the file |
| `park2021` | R. S. Park, W. M. Folkner, J. G. Williams, D. H. Boggs, "The JPL planetary and lunar ephemerides DE440 and DE441", AJ 161, 105 (2021). doi:10.3847/1538-3881/abd414 | First 40 % of the article read by a fetch, quoted verbatim: DE440 covers 1550–2650; orbits tied to ICRF3 to about 0.2 mas; the barycentre moved about 100 km against DE430; DE440 and DE441 differ by less than 1 m. An accuracy for the Earth was **not** found (mismatch list). | fold 7a.3–7a.4, localisation, glossary, formulas: DE440 |

## B · Record checked, content through a second source (9)

| id | Reference | Second source and what it confirms | Cited for |
|---|---|---|---|
| `taylor1992` | J. H. Taylor, "Pulsar timing and relativistic gravity", Phil. Trans. R. Soc. Lond. A 341, 117–134 (1992). doi:10.1098/rsta.1992.0088 | Pennucci, Demorest & Ransom 2014 (ApJ 790, 93): the algorithm "is laid out in detail in Taylor (1992, Appendix A)"; Liu et al. 2014 (arXiv:1407.3827): "Appendix A of Taylor (1992) prescribes a Fourier frequency-domain technique for measuring the phase shift". The paper itself is behind a paywall. | toa 8.2, physics P5, glossary, formulas: FFTFIT |
| `lorimerkramer2005` | D. R. Lorimer, M. Kramer, *Handbook of Pulsar Astronomy*, Cambridge University Press. ISBN 0 521 82823 6. Copyright page 2005; often cited as 2004 | `kulkarni2020` quotes its dispersion constant, 4.148808(3)·10³ MHz² pc⁻¹ cm³ s. The equation number was not confirmed. | sky 2.1, fullband 5.1, physics P3, formulas: the dispersion constant |
| `smith_sasp` | J. O. Smith, *Spectral Audio Signal Processing*, online book, CCRMA, Stanford (2011), section "Under the hood of kaiserord", https://ccrma.stanford.edu/~jos/sasp/ | Opened: "the filter order M is estimated from stop-band attenuation A and desired transition width Δω using the empirical formula M = (A − 8)/(2.285·Δω)", given as Kaiser's formula. Kaiser's own 1974 paper is in state D. | channels 5a.2, receiver 4.2: the Kaiser rule |
| `iau2006b3` | IAU 2006 Resolution B3, "Re-definition of Barycentric Dynamical Time, TDB", https://iauarchive.eso.org/static/resolutions/IAU2006_Resol3.pdf | Search summary of the resolution text: L_B = 1.550519768·10⁻⁸, T₀ = 2443144.5003725, TDB₀ = −6.55·10⁻⁵ s. Not opened by me. | localisation L3: TCB against TDB |
| `brown2025` | G. C. Brown, S. Z. Sheikh et al., "Long-term monitoring of scintillation in the pulsar J0332+5434", ApJ 995, 120 (2025). doi:10.3847/1538-4357/ae0fff | Abstract through a search: Allen Telescope Array, 900–1956 MHz; mean scintillation bandwidth 2.6 MHz (sd 1.3) at 1140 MHz and 14.7 MHz (sd 9.1) at 1668 MHz; index 3.1 ± 1.0. The 780 s time scale of notes §11 was not found. | scenario R6. On the page now as "ATA 2025" |
| `fairhead1990` | L. Fairhead, P. Bretagnon, A&A 229, 240 (1990); title not confirmed | `kaplan2005`: eq. (2.6) "is a truncated form of a much longer and more precise series given by Fairhead & Bretagnon (1990)", with the journal reference in its list. | fold 7a.4 (History), `loadEphemeris.m:37` |
| `bulletinc` | IERS Bulletin C (leap seconds), https://hpiers.obspm.fr/iers/bul/bulc/ | Search result quoting the bulletin: "from 2017 January 1, 0h UTC, until further notice: UTC − TAI = −37 s". | localisation L3: TAI − UTC = 37 s, and so 69.184 s |
| `3gpp_band32` | 3GPP TS 37.104, operating band table: band 32 | Search result quoting the table: band 32 is downlink only, 1452–1496 MHz. | receiver 3.5, parameters: LTE band 32 |
| `isgps200` | IS-GPS-200, NAVSTAR GPS Space Segment / Navigation User Interfaces | A Keysight paper citing it: −158.5 dBW for L1 C/A at the antenna. The other GNSS levels are in state D. | receiver 3.5: "published minimum received power" |

## C · Record checked, named method or data set only (30)

Cited as the origin of a name; the text was not opened.

| id | Reference | Cited for |
|---|---|---|
| `isserlis1918` | L. Isserlis, Biometrika 12, 134–139 (1918). doi:10.1093/biomet/12.1-2.134 | Isserlis' theorem: channels 6b, excision E4, physics, formulas, glossary |
| `reed1962` | I. S. Reed, "On a moment theorem for complex Gaussian processes", IRE Trans. Inf. Theory 8, 194–195 (1962) | the complex form the pages use, cov(\|z_a\|², \|z_b\|²) = P²\|R\|² |
| `neymanpearson1933` | J. Neyman, E. S. Pearson, Phil. Trans. R. Soc. A 231, 289–337 (1933). doi:10.1098/rsta.1933.0009 | Neyman–Pearson test: detection, physics, fullband, glossary |
| `white1980` | H. White, Econometrica 48, 817 (1980). doi:10.2307/1912934 | sandwich covariance: toa 8.6, formulas, glossary |
| `massey1951` | F. J. Massey, J. Am. Stat. Assoc. 46, 68–78 (1951) | Kolmogorov–Smirnov test: toa 8.4, scripts |
| `kalman1960` | R. E. Kalman, J. Basic Eng. 82, 35–45 (1960). doi:10.1115/1.3662552 | Kalman filter: navigation, glossary |
| `rts1965` | H. E. Rauch, F. Tung, C. T. Striebel, AIAA J. 3, 1445–1450 (1965). doi:10.2514/3.3166 | the smoother: navigation 10.4, glossary |
| `dicke1946` | R. H. Dicke, Rev. Sci. Instrum. 17, 268–275 (1946). doi:10.1063/1.1770483 | radiometer equation |
| `hankinsrickett1975` | T. H. Hankins, B. J. Rickett, "Pulsar signal processing", Methods in Computational Physics 14, 55–129 (1975) | coherent dedispersion: sky, fullband, physics |
| `rickett1990` | B. J. Rickett, Annu. Rev. Astron. Astrophys. 28, 561–605 (1990) | scintillation scalings (ν^4.4, ν^1.2): physics |
| `kaiserschafer1980` | J. F. Kaiser, R. W. Schafer, IEEE Trans. ASSP 28, 105–107 (1980). doi:10.1109/TASSP.1980.1163349 | the Kaiser window |
| `bellanger1976` | M. Bellanger, G. Bonnerot, M. Coudreuse, IEEE Trans. ASSP 24, 109–114 (1976) | polyphase filterbank: channels 5a, glossary |
| `harris2003` | F. J. Harris, C. Dick, M. Rice, IEEE Trans. MTT 51, 1395–1412 (2003) | polyphase channelizer with a non-critical slide |
| `price2016` | D. C. Price, "Spectrometers and polyphase filterbanks in radio astronomy", arXiv:1607.03579 | the same, for radio astronomy |
| `stockham1966` | T. G. Stockham, Proc. AFIPS Spring Joint Computer Conf., 229 (1966). doi:10.1145/1464182.1464209 | FFT convolution in blocks (overlap-add, overlap-save) |
| `matsumoto1998` | M. Matsumoto, T. Nishimura, ACM Trans. Model. Comput. Simul. 8, 3–30 (1998) | `mt19937ar`: sky, receiver, conventions |
| `lecuyer1999` | P. L'Ecuyer, Oper. Res. 47, 159–164 (1999) | `mrg32k3a`: receiver, conventions |
| `allan1966` | D. W. Allan, Proc. IEEE 54, 221–230 (1966) | Allan variance |
| `barnes1971` | J. A. Barnes et al., IEEE Trans. Instrum. Meas. IM-20, 105–120 (1971) | power-law clock noise (h coefficients) |
| `shapiro1964` | I. I. Shapiro, Phys. Rev. Lett. 13, 789–791 (1964) | Shapiro delay |
| `capitaine2003` | N. Capitaine, P. T. Wallace, J. Chapront, A&A 412, 567–586 (2003) | precession angles (`delayModel.m:232`) |
| `hobbs2006` | G. B. Hobbs, R. T. Edwards, R. N. Manchester, MNRAS 369, 655–672 (2006) | TEMPO2 overview; predictors ("polycos") |
| `manchester2005` | R. N. Manchester, G. B. Hobbs, A. Teoh, M. Hobbs, AJ 129, 1993–2006 (2005). doi:10.1086/428488 | the ATNF pulsar catalogue: scenario R2 |
| `lorimer1998` | D. R. Lorimer et al., A&AS 128, 541–544 (1998) | the EPN profile database and its format |
| `gangadhara2001` | R. T. Gangadhara, Y. Gupta, ApJ 555, 31–39 (2001) | B0329+54's emission components (notes §11) |
| `amari1997` | S. V. Amari, R. B. Misra, IEEE Trans. Reliab. 46, 519–522 (1997) | the closed form for a sum of exponentials: excision E2.1 |
| `maron2000` | O. Maron, J. Kijak, M. Kramer, R. Wielebinski, A&AS 147, 195–203 (2000) | power-law pulsar spectra: sky 1.1, limits |
| `rainal1988` | A. J. Rainal, "Origin of Rice's formula", IEEE Trans. Inf. Theory 34, 1383–1387 (1988) | history of Rice's formula |
| `pennucci2014` | T. T. Pennucci, P. B. Demorest, S. M. Ransom, ApJ 790, 93 (2014) | second source for `taylor1992` |
| `smirnova2006` | T. V. Smirnova, C. R. Gwinn, V. I. Shishov, A&A 453, 601–607 (2006) | companion to `gwinn2006` ("Paper 2") |

## Added after the rollout (9), 24 Nov 2026

| id | Reference | State and what was confirmed | Cited for |
|---|---|---|---|
| `rao1945` | C. R. Rao, Bull. Calcutta Math. Soc. 37, 81–89 (1945) | C; record from Wikipedia's reference list for the Cramér–Rao bound, not from Crossref | physics P5, receiver 3.2, formulas, glossary: Cramér–Rao bound |
| `cramer1946` | H. Cramér, *Mathematical Methods of Statistics*, Princeton University Press (1946) | C; same route | same |
| `kay1993` | S. M. Kay, *Fundamentals of Statistical Signal Processing, Vol. I: Estimation Theory*, Prentice Hall (1993), ISBN 0-13-345711-7 | B; publisher's listing and table of contents: the Cramér–Rao lower bound is chapter 3 | physics P5 |
| `wgs84` | NGA.STND.0036_1.0.0_WGS84, *Department of Defense World Geodetic System 1984* (8 July 2014) | B; a = 6378137.0 m, 1/f = 298.257223563 in a copy of the specification (UN OOSA), equal to `makeTrajectory.m:86` | localisation L2, glossary |
| `meeus1998` | J. Meeus, *Astronomical Algorithms*, 2nd ed., Willmann-Bell (1998), ISBN 0-943396-61-1 | C; record from a secondary list; chapter 22 is nutation by a software library's documentation | localisation L6, fold 7a.3: nutation (`delayModel.m:233`) |
| `horizons` | J. D. Giorgini et al., "JPL's on-line solar system data service", BAAS 28, 1158 (1996) | C; record from NASA's report server and a later abstract; ten authors | fold 7a.3 (test 11), localisation L6 |
| `itu_res222` | ITU, Resolution 222 (Rev. WRC-23) | B; search summary of the resolution and of NTIA's compendium: 1525–1559 MHz mobile-satellite, space to Earth, used by Inmarsat. The full title was not confirmed | receiver 3.5: Inmarsat downlink |
| `kaiser1974` | J. F. Kaiser, Proc. 1974 IEEE Int. Symp. Circuits and Systems, 20–23 (April 1974) | C; record from MathWorks' reference list | receiver 4.2, channels 5a.2, next to [29] |
| `huber1967` | P. J. Huber, Proc. Fifth Berkeley Symp. Math. Statist. Probab. 1, 221–233 (1967) | C; record from a bibliography service | toa 8.6, formulas, glossary, next to [40] |

Settled without a number:
- `psrcat_version`: version 2.6.5 exists (26 Aug 2025, 4303 pulsars) and 2.7.0 followed on 5 Oct 2025, by the catalogue's
  own history page. Noted inside [45].
- `irwin1999`: A. Irwin, T. Fukushima, A&A 348, 642–652 (1999), as quoted inside `kaplan2005`. No page names it.
- `iers2010`: G. Petit, B. Luzum (eds.), *IERS Conventions (2010)*, IERS Technical Note 36, Verlag des BKG, Frankfurt am
  Main, ISBN 3-89888-989-6. Record confirmed; no page names it, and its content was not opened.
- `oppenheimschafer`, `kay1998`: dropped, not needed ([29] and [21] cover the rule and the block convolution; [6] the test).

## D · Still not confirmed (4)

| What the docs need | State on 24 Nov 2026 |
|---|---|
| GNSS received levels and carrier frequencies (receiver 3.5) | GPS L1 C/A is [28]. Found in secondary sources only: GLONASS G2 = 1246 + 0.4375·K MHz (ITU-R M.1904-1), Galileo E6 band 1258.29–1299.21 MHz (centre 1278.75 MHz by arithmetic), BeiDou B1I 1561.098 MHz at −163 dBW (its ICD, version 2.0), B3I 1268.52 MHz. Levels for E6, B3I and G2 on the ground: not found |
| Radar characteristics (receiver 3.4: "still to be checked against ITU emission masks") | not looked up |
| Troposphere 8 ns and 44 ns (localisation L6); sky temperature 3–5 K (scenario R1); the constant 40.3 (limits, formulas) | not found; the ionospheric delay per TEC unit is confirmed through `edwards2006` |
| A textbook source for the matched filter and for first-order error propagation | Kay's volume II was not confirmed for this |
