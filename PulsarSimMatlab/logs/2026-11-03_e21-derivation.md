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

## More questions on E2 (answered in chat, no change)
- **Why is the channel correlation (E2.2) measured only from pass 2?** It must leave out pulsed RFI (a +20 dB radar pulse is
  coherent over ~8 samples and would make the channel look correlated, raising its own thresholds) but keep a carrier's own peaks
  (pass 1 flags ~1400/s of them with the white thresholds). Pass 1's flagged intervals, sorted by peak (>= 30 x baseline = pulsed
  RFI), make that split. The second pass itself exists for the baseline (E2.3: 5.4 % high under dense bursts, 0.5 % after pass 2);
  the channel thresholds reuse it (20 Oct, logs/2026-10-20_blockC.md).
- **Does pass 2 use pass 1's mask?** Only to leave samples out of its estimates: the correlation (strong intervals only), the
  baseline (all pass-1 flags; strong intervals only for own-threshold channels; blocks with < 1000 unflagged samples keep pass 1's
  value). The window test runs over every sample again, and the final mask is pass 2's alone (detectRFI.m:239-282).

## To read on 4 Nov (Jasper: "remind me to read this last part tomorrow"): E4, why blanked parts need correcting after dedispersion
Before dedispersion a blank is easy: a radar pulse blanked with its guard is ~50 channel samples of exact zeros (52 in Figure E.3: 32 flagged, because the windows reach past the pulse, plus a guard of 10 on each side; this line first said ~28, corrected in the 4 Nov audit), and a fold could
give exactly those samples weight 0.

Dedispersion spreads every blank. It is a convolution with a chirp filter h whose sweep is ~660-1560 channel samples (159-374 us)
at the default DM 25, and ~130-310 (32-75 us) at the DM 5 of the docs' figure, y(a) = sum h(l) x(a-l), so every output sample within reach of the blank misses part of its inputs. In a chirp each lag belongs to
one frequency, so a blank removes a thin slice of the channel's band from each output it touches, a different slice per output.
No output is empty or complete: each holds a fraction W of valid data (Figure E.3d). Its detected power averages
W x (noise + pulsar), so every blank becomes a dip in the dedispersed power, as wide as that sweep.

A naive fold counts every sample as complete, so the dips go into the profile:
- blanks at random phases (RFI unrelated to the pulsar): the dips spread over all phases and lower the baseline almost uniformly;
  harmless (the baseline is fitted out), apart from a little extra noise the error model does not know;
- blanks that follow pulse phase (a radar locked to the period, E6; periodicRFI blanking predicted pulses, E8; blanks triggered by
  the pulse): the dips stack at the same phase, distort the profile and bias the TOA. At -54 dB the pulse peak is ~4e-6 of the
  baseline, so a 1e-5 dip at a fixed phase already beats the pulse. Fake blanking at phase 0.30-0.36: naive fold -0.574
  (= 1 - 0.425, the mean W there), weighted fold -0.0019 +- 0.0009 (fold.html 7.5).

What blankingWeights corrects:
- the mean: the fold divides by sum a W instead of sum a, so a W = 0.5 bin counts as half a sample with the right mean. A blank
  removes the same frequency slice from noise and pulsar (both fill the band), so dividing by W restores both;
- the noise: V and X(L) per bin (variance, shared noise with bin k+L) set the error bars, chi2 and optimal weights; no simple rule
  (W and W^2 both off by -50 % to +40 %, expBlankingVariance), so they are exact from Eq. (E.3);
- empty bins (E4.3): W = 0 bins add nothing; filter-tail power there biased nearly empty phase bins by up to +285 %.

Not recorded as the reason, Claude's reasoning: throwing away every output sample a blank touches would cost the whole sweep per
blank, most samples nearly complete; fractional weights keep all valid data and still remove the bias.

## 4 Nov: more questions on E4 (answered in chat); DM correction
Checked in Python (not MATLAB), noise-free, one profile, equal weights; FWHM 5 % of a 10 ms period (sigma 212 us), NBin 2048; scripts
dip_bias.py, relative_peak.py, w_energy.py and w_energy_ends.py in the session scratchpad.
- **How does missing data bias the TOA?** It does not by itself: left out of the fit, a gap on the flank gives 0.000 us (error bar
  x 1.16). Counting it as zero power gives -467 us at -20 dB. The dip is almost all missing noise: d = (1 - e)(1 + rho t), and to first
  order tau = (1/rho) sum e t' / sum t'^2, confirmed by the fit. A dip of depth 1e-4, width 0.005 turns, at -30 dB: 0 at the peak,
  -+3.33 us on the flanks, 0 far off-pulse or uniform. On the trailing flank: -0.335, -3.33, -32.7, -300 us at -20, -30, -40, -54 dB
  (first order fails at -54 dB, where the dip is 25 x the peak); only the pulsar lost: -0.002 us at every strength.
- **Does the data left act as a relative peak?** Only when it is about pulse-sized: all blanked but a window of 0.2 turns or less, or
  dips every 2.5 FWHM or less, and the fit centres on the gap. Wider stretches are a plateau with no position (a dip 25 x the peak far
  off-pulse: 0 us; the same dip with no pulsar: 84 % of shifts within 0.1 % of the best). Then a dip only pushes the pulse away.
- **Is W energy-weighted or a count of keep(i)?** Energy-weighted: w(a) = sum |h(l)|^2 keep(a-l) / sum |h|^2, W its bin mean; the default
  'channelizer' form adds the R_x cross terms (blankingWeights.m:202-221). For a 28-sample blank at DM 25, 1.402 GHz: deepest 1 - W is
  0.0330, against 0.0294 for a count over the sweep and 0.0005 for a count over all 54 348 taps the filter carries (bulk delay). The 52 samples of Figure E.3 give a deepest 1 - W of 0.058 at DM 25, 1.4 GHz (4 Nov audit).
  |h|^2 is flat to ~0.3 % mid-sweep, up to +12.7 % near the ends. The area sum(1 - W) is exactly the blank's length.
- **DM correction:** the "130-310 samples" quoted earlier in this log is the sweep at DM 5, which E1 gave without naming its DM. E1 now
  gives both: 660-1560 at DM 25 (same definition, 2 K DM dF / f^3; 99 % of the energy lies in 644-1504), 130-310 at DM 5. Figure E.3
  stays at DM 5 (Jasper's rule of 21 Oct); its caption does not name the DM, not changed.
