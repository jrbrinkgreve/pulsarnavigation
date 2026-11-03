# 3 Nov 2026: questions on detectRFI's exact thresholds (E2.1); derivation added to the docs

Jasper asked about the excision page, section E2.1 (Eqs. E.1 and E.2), during a read of the docs. No code change.

## Questions and answers
- **Where does E2.1 come from?** Three standard results put together, not one paper:
  1. Each channel sample is a linear combination of Gaussian receiver noise, so a window is jointly Gaussian, zero mean and
     circular, with covariance m·toeplitz(rho). For white input, rho(d) = sum h(r) h(r+dD) / sum h^2.
  2. Whitening with the eigenvectors of R turns the power sum into sum lambda_k E_k with E_k unit exponentials (Eq. E.1).
  3. The tail of a sum of exponentials is a phase-type (hypoexponential) survival function, [1 0 .. 0] expm(Q x) 1 (Eq. E.2).
     Check: n = 1 gives e^-x = 1e-6 at x = 13.82 = eta_1.
- **Why the mixing phase drops out:** with mix-then-filter and the LO phase tied to the absolute input sample (channelizeIQ),
  both outputs see each shared input sample at the same LO phase, so the phase cancels inside the sum. In the chat, Claude
  first gave a different reason (a unitary diagonal similarity, which would apply to a filter-then-mix channelizer) and corrected it.
- **Window vs block:** a block (16 ms, about 66 700 channel samples, fixed tiles) sets the baseline m from the median; a window
  (1-16 samples, sliding by half) is the test statistic. Each window is compared with eta_n times the baseline of the block it
  starts in. The block is long for a precise, RFI-robust m; the window is short to match a burst.
- **Cost of the exact form vs the matched Gamma:** only the threshold numbers change; the per-window test is identical.
  `windowThresholds` runs once for the white thresholds and once per channel in pass 2 (129 calls for 128 channels).
  Python stand-in (not MATLAB): about 0.8 ms per call exact vs 0.06 ms Gamma, so about 0.1 s per run, independent of data length,
  against tens of seconds per second of data for the rest of detectRFI. Jasper: no MATLAB timing now (machine busy).
- **The "chain of exponential stages":** Jasper read it as the chance of passing the excision stage. It is a calculation device
  only: each lambda_k E_k is read as a random waiting time, with the power level x in the role of time, and "still waiting at x"
  is "the sum exceeds x". The thresholds describe noise alone, without RFI.

## Docs change (Jasper: "update the wording and add a short but clear derivation part")
- excision.html E2.1: one sentence saying the thresholds describe noise alone; Q's "(a chain of exponential stages)" became
  "(the rates of a chain of random waiting times, see the derivation below)", because "stages" means pipeline stages elsewhere.
- New collapsed box "Derivation of Eqs. (E.1) and (E.2)" after the symbol list: covariance of a window, whitening to E.1,
  the waiting-time chain to E.2, the n = 1 check, and why the partial-fraction closed form cancels. Equations unnumbered, so
  E.3 and E.4 keep their numbers.
- Checked in the browser: no KaTeX or console errors; the four new equations fit the column (the chain equation split over two lines).
  details tags balanced (14/14).
- maintenance.html changelog row; history.html session entry and log row.

## Follow-up: why the closed form fails (Jasper: "investigate to be sure, report back and correct then")
The first version of the derivation said the closed form cancels because "the eigenvalues all lie close to 1". At n = 16 they span
0.61-1.05, so that was imprecise. Checked in Python (not MATLAB) with the real prototype from data/chan/test_rx_IQ_chan_info.mat
(3853 taps, D = 192); scripts closed_form_check.py and eig_vs_spectrum.py in the session scratchpad.
- Same rho as detectRFI: the expm thresholds come out 13.816, 8.391, 5.387, 3.681, 2.686 (table: 13.82 8.39 5.39 3.68 2.69).
- P at the table thresholds, three independent ways: the closed form in 80-digit decimal, uniformization (a sum of positive terms
  only) and scipy expm. Uniformization matches the 80-digit value to 1e-14; expm to 3e-9 (n = 8) and 7e-10 (n = 16).
- The closed form in double precision: right at n = 2 and 4; 7.57 at n = 8 and -3.1e25 at n = 16 (true values 1.005e-6 and
  9.56e-7). Largest terms 2.5e15 and 8.7e39.
- The cause is near-equal eigenvalues, not eigenvalues near 1. At n = 16, with the same range 0.614-1.051 but evenly spaced, the
  double closed form is right to 2e-12; evenly spaced over 0.95-1.05 it loses about 2 digits (5e-3); the real eigenvalues
  (closest pair 1.4e-8 apart) give nonsense.
- Why near-equal: the channel's power spectrum is flat at 1.0513-1.0515 across the whole channel (|f| <= 1.5625 MHz, spread
  4.2e-4) and falls off only in the oversampling margin (1.037 at 1.67 MHz, 0.53 at 2.08 MHz). The eigenvalues roughly follow
  the spectrum: at n = 16, ten of sixteen lie within 1e-4 of the flat level; at n = 8, four of eight (closest pair 1.4e-9).
- excision.html, last paragraph of the derivation box, rewritten with these facts. maintenance.html changelog row.

## Left open
- detectRFI.m comments still say "chain of exponential stages" (lines ~31 and ~440); not changed (not asked).
- Eq. (E.1) already overflowed the column at a narrow window before this change (it scrolls sideways); not changed.
