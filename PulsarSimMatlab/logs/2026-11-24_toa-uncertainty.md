# 24 Nov 2026: the TOA error bar derived step by step, and the reference mismatches noted (docs only; no code change)

Jasper: "note the mismatches for later. Can you go to the TOA uncertainty estimation part, and make the derivations there
clearer? I have a hard time following 100%". Asked which steps lose him, he picked all four: why C′ is a weighted sum,
why δτ = −δC′/C″, what dᵀ·Cov·d means, and where √2·σ_t/SNR comes from. Scope: Section 8.3 only, with one computed
figure. The plan was approved before any edit. MATLAB was not run; nothing is committed.

## The mismatches
- `docs/roadmap.html`, "Open items": a `todo` row "Reference mismatches (24 Nov)" with the six that matter most and a
  pointer to `logs/2026-11-24_docs-references.md`, which keeps all of them with page and line.
- Not added to `logs/current_project_notes.md` §10: that file holds another session's uncommitted changes. The item
  can be copied there once it is free.

## Section 8.3 of `docs/toa.html`
The four opening lines became five steps, each one a paragraph that opens from the one before. The derivation follows
`fitOne` in `functions/estimateTOA.m`.
1. The correlation as a sum over phase bins, C(τ) = c₀·Σ p_j·s̃(φ_j − τ) with c₀ = N/2, Eq. (8.7). Its slope is linear in
   the data, C′ = Σ d_j·p_j with d_j = −c₀·s̃′(φ_j − τ), Eq. (8.8) (`dvec`).
2. At the true shift the model's part of the slope is zero, so δC′ = Σ d_j·n_j, with variance
   Σ d_j²·Cov_jj + 2·Σ_L Σ_j d_j·d_{j+L}·Cov_{j,j+L} = dᵀ·Cov·d, Eq. (8.9) (`varC1`). The second sum is the shared noise
   of neighbouring phase bins; without it the error bar is about 18 % too small.
3. The peak is where the slope is zero, so the lift δC′ moves it by −δC′/C″, Eq. (8.10).
4. var(τ̂) = var(δC′)/C″² = dᵀ·Cov·d/C″², Eq. (8.11) (the old Eq. 8.7), and σ_TOA = σ_τ/f_ref. C″ is taken as a constant (`c2`).
5. White noise: var(τ̂) = σ²/(b²·Σ s̃′²). For a Gaussian pulse Σ s̃′² = N·√π/(2σ_φ) and SNR² = b²·N·σ_φ·√π/σ², hence
   σ_TOA = √2·σ_t/SNR. The fitted SNR is 0.962 of that one, which is the page's "≈".
The amplitude error is written out with weights g_j = c₀·s̃(φ_j − τ̂) (`varC0`, `bErr`), in place of "and similarly".

Kept: the id `err`, the paragraphs on the two noise models, the callout, every code name and link of the old text.
The equation numbers in this note are the final ones (see "Step four and the numbering" below). The section's prose
(without equations and symbol lists) grew from 362 to 1105 words. The plan had estimated about 350 more; it became
about 740 more.

Tightened the same day (Jasper: "try to tighten it a bit"): the section's prose is now about 925 words, of which about
645 are the derivation and the rest the unchanged paragraphs on the noise models. Every equation, the figure and the
five steps stayed. Shortened: each step's lead-in and closing sentence, the hill picture (two short sentences), the
delta-method paragraph, the white-noise check and the amplitude paragraph. The symbol lists lost the entries for n_j
(defined in the opening sentence) and for "d, C″ at the peak" (now part of the C′, C″ entry); s̃ and s̃′ share one entry.
Rechecked afterwards: every code name, link and number of the section is still there (`dvec`, `varC1`, `weightX`, `c2`,
`varC0`, `bErr`; the links to fold, physics and detection; 0.962; 18 %), 1492 links with none broken, equation and
figure labels as before, no KaTeX or console errors in the browser, no sideways scroll.

## New Figure 8.2
`docs/tools/figures/gen_toaerr.py`: the profile of Figure 8.1 with more noise (SNR 14.6). (a) C(τ) near its peak
without noise and for one draw. (b) C′(τ) of both, with the lift δC′ and the two zero crossings. (c) the fitted shift
over 4000 draws against the Gaussian of Eq. (8.11). The script prints: Eq. (8.11) 20.52 µs, scatter 20.99 µs (ratio 1.023);
for the drawn case a shift of 31.35 µs against −δC′/C″ = 34.33 µs.
Renumbered: the χ² figure 8.2 → 8.3 and the notch figure 8.3 → 8.4, in `toa.html` (captions and two references),
`gen_chi2.py`, `gen_chanweights.py` and the figures' `README.md` (plus a row for the new script).

## Numerical check of the derivation (Python with numpy, a port of `fitOne`; not MATLAB)
Profile: 2048 bins, 5 % Gaussian, b = 0.5, noise 0.12 per bin (SNR 36.6), 4000 draws per case.
- Step 1: `dvec` = −(N/2)·s̃′ to 1.1·10⁻¹⁴ (relative); C = (N/2)·Σ p·s̃ to rounding; the amplitude weights likewise
  (2·10⁻¹⁶).
- Step 2: Σ d_j = 0 exactly; the model's part of C′(τ₀) is 1.3·10⁻¹⁰ on a scale of 5·10⁶.
- Step 3: τ̂ − τ₀ against −δC′/C″, draw by draw: correlation 0.9989, difference 4.7 % of the scatter (white noise).
- Step 4: scatter of τ̂ divided by Eq. (8.11): 0.988 ± 0.011 for white noise, 0.985 ± 0.011 with a neighbour
  correlation of ¼. With the lag terms left out the second ratio is 1.206 (expected √1.5 = 1.22).
- Step 5: C″ of the model equals −(N/2)·b·Σ s̃′² to all digits; σ/(b·√Σ s̃′²) = √2·σ_φ/SNR = 8.2088·10⁻⁴ turns, both to
  all digits, against 8.2212·10⁻⁴ from Eq. (8.11); the fitted SNR is 0.9616 of the known-baseline one.
- Amplitude: scatter of b̂ divided by its formula 0.983 (white), 0.990 (correlated).

## Checks
1492 links, none broken; equation labels on the page in order (then 8.1 to 8.12 with 8.7a–d; now 8.1 to 8.16); figure labels 8.1 to 8.4;
the captions of all four figure scripts equal the page word for word; every code name and link of the old section is
still there. In the browser (local docs server, dark theme): the section renders without KaTeX or console errors;
at phone width the page does not scroll sideways, and the wide equations scroll inside their own box as (8.6) and
(8.10) already did. The light theme was not looked at.

## Step four and the numbering (same day)
Jasper: "The fourth step takes the variance of Eq. (8.7d), with C″ as a fixed number (c2 in the code)." is not clear
yet. Also, the equation numbering seems off".
- Step four now says why C″ may be taken as a constant: from one noise draw to the next δC′ changes completely (pure
  noise, mean zero, variance of Eq. 8.9), while C″ is set by the pulse and changes by a fraction of order 1/SNR. Then
  "a random number divided by a constant has its variance divided by the constant squared", and the equation shows the
  step between: var(τ̂) = var(δC′)/C″² = dᵀ·Cov·d/C″². Checked for the Gaussian of the numerical check: the relative
  noise of C″ is σ·√(Σ s̃″²)/(b·Σ s̃′²) = √3/SNR, both 0.03086 at SNR 56.1 (2048 bins, 5 %, b = 0.5, noise 0.12).
- Numbering: the first version labelled the steps (8.7a) to (8.7d) and kept the result as (8.7), so that a plain (8.7)
  came after (8.7d). Now one sequence: (8.7) to (8.10) for the steps, (8.11) for the result, and the old (8.8) to
  (8.12) are (8.12) to (8.16). Followed in `toa.html` (28 mentions, among them the caption of Figure 8.1 and the page's
  code table), `gen_toaerr.py`, `gen_toa.py`, the figures' `README.md` and today's changelog row. No other page, no code
  and no note cites these numbers (searched); older changelog rows keep the numbers of their day. The two white-noise
  equations stay without a number, as side results do on other pages.
- Checks: labels on the page (8.1) to (8.16) in order; both figure scripts rerun in the scratchpad and their output
  equals the page's figure blocks; 1492 links, none broken; no KaTeX or console errors.
