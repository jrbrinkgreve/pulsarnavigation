function testFrontEndFaults()
%TESTFRONTENDFAULTS  The channel path under the faults of a real receiver: each case puts one fault into the IQ file.
%{
Run from the PulsarSimMatlab folder: run('tests/testFrontEndFaults.m'), or runAllTests("testFrontEndFaults").
~7 min (13 min on a busy machine); writes ~3 GB of temporary files to tempdir (deleted at the end). The simulated
receiver writes an ideal IQ file, a real one does not. Each case changes the frozen seed-43 IQ
file (data/test/, testParams: the 10 ms pulsar at DM 5, -5 dB, 100 ms at 800 MS/s, LO 1.3 GHz)
the way a real receiver would (truth side, in this test), runs the observer's chain unchanged
(processIQ, channel path without excision, then foldToTOA) and checks its response against a
prediction. Risk register: logs/2026-11-02_risk-register.md; this test: logs/2026-11-03_real-hardware-faults.md.

  1. FM-8, the frequency mapping. processIQ assumes RF = fLO + f_bb. A mixer with its LO above
     the band, or I and Q swapped, gives RF = fLO - f_bb: the spectrum is mirrored about the LO.
     The file is first re-centred on an LO of 1.4 GHz (x(n) exp(-i pi n/4): the same RF, the band
     at -200..+200 MHz), so that the mirror maps the band onto itself. Observer channel j then
     holds the true channel 129 - j, reversed: its content at in-channel offset b is true RF
     f_(129-j) - b, which arrives
       tau(f_(129-j) - b) - tau(f_j + b),   tau(f) = 4.148808e3 s MHz^2 * DM / f^2,
     after the time the observer's dedispersion assigns to it: -6.3 to +6.3 ms across the band
     at DM 5 (0.63 turns of the 10 ms pulsar), ~0.1 ms within a channel. Prediction: channel j's
     profile is the template moved by those delays (averaged over b), and the fit's amplitude is
     the mean over channels of <P_j, c_tau> / <c, c> at the best tau (c = template - its mean:
     each channel's own baseline takes the smeared pulse's mean), as a fraction of the true one.
     (a) The re-centred file, with the LO given as 1.4 GHz, gives the frozen processing's TOAs
         (data/test/chan/test_dedisp_chan_power.dat, folded the same way) within 1e-3 sigma:
         the re-centring changes nothing.
     (b) The mirrored file, processed as if it were not: the total fold's amplitude within 3
         of its error bars of the predicted fraction of the re-centred one (0.043), and its
         phase within 0.01 turns of the predicted one (half a turn: the channels moved by
         +-0.5 turns pile up there).
     (c) The LO given 1 MHz above the true one (1.401 GHz): channel j holds true RF f_j + b - 1 MHz,
         which arrives tau(f_j + b - 1 MHz) - tau(f_j + b) late (10-24 us across the band); the fit
         takes the mean over the channels. The total fold's time offset against (a) within 5 %
         of that mean.
     Shown: the sub-ints of (b) - how many pass the good-TOA rule, at which phase, with which
     red. chi^2; for B0329+54 (P 714.5 ms, DM 26.76, five components) the same prediction.
  2. Quantization. A real receiver delivers b bits of I and of Q, from its ADC or after the
     transport (an SDR streams 16 or 8 bits; 18 streams at 400 MS/s favour 8). The file is
     quantized here (mid-rise, step = rms / level) and processed. Bussgang: Q(x) = alpha x + e,
     e uncorrelated with x and, for 4 bits or more, white over the sample rate, so in the band the
     quantization adds the noise epsQ = var(e) / (alpha^2 kappa sig^2) (kappa = in-band PSD / mean
     PSD of the file, measured: 1.12). A weak pulsar's SNR keeps 1 / (1 + epsQ). This one is bright
     (pulse peak ~0.3 x the baseline), its own noise dominates on the pulse and it loses a little
     less: the fit's radiometer model (variance (1 + epsQ + beta T)^2 per bin) gives the prediction.
     Checked, for 8 bits at rms 4 LSB and at 1 LSB (the level that leaves room for a +40 dB radar)
     and 4 bits at rms 3 LSB: the total fold's SNR within 0.005 of the prediction, and red. chi^2
     within 0.01 of the unquantized fold (the noise model stays honest). Shown: 2 bits at rms 1 LSB,
     beyond the white-error model (there the error follows the input).
     Two's complement (floor instead of mid-rise: -1/2 LSB on I and on Q) puts a tone of 1/2 LSB^2
     at the LO, as LO leakage and a DC offset do. The LO sits on a channel edge (the grid is offset
     by half a channel from it). Checked: the two channels that meet there get +0.5 LSB^2 in their
     channelizeIQ files (within 5 %), and their dedispersed power changes by < 1e-3: the
     dedispersion's sin^2 taper is zero at the channel edge.
  3. FM-9, a lost buffer: n = 33024 samples (41.28 us, 43 detected bins) at 45.3 ms.
     (a) The stream is joined without the gap: every later sample is labelled n/fs early, a clock
         step. Checked: the sub-ints whose data all come before the drop (ending a band sweep,
         6.3 ms, before it) keep their TOAs within 1e-3 sigma; those after it move by -n/fs within
         0.1 sigma.
     (b) The gap is filled with zeros (time stamps kept, as the F ingest plans). Shown: without a
         mask, the TOA and red. chi^2 of the sub-ints whose channels hold the hole. Checked: with
         the gap blanked through the excision path (detectRFI's mask plus rows for the gap, widened
         by the channel filter's half-length; blankChannels, blankingWeights), every TOA within 0.3
         sigma of the same chain without the gap (the gap's sub-int loses ~1 % of its information
         in the channels where the hole meets the pulse: an expected change of ~0.1 sigma), and
         red. chi^2 within 0.005.
  4. FM-39, absolute time. The simulations start at t0 = 0; real data carry TT seconds since J2000
     (8.3e8 s in 2026), where one double resolves 0.12 us. The same file with t0 = 8.3e8 +
     0.123456789 s and the ephemeris' TRef moved with it. Checked: every TOA minus t0 within 4 ULP
     of t0 (0.48 us) of the t0 = 0 run: the chain loses no more than the double that holds the TOA.
Errors at the end if any check fails. (3 Nov 2026, risk register FM-8, FM-9, FM-13, FM-39, FM-40.)
%}

root = fullfile(fileparts(mfilename('fullpath')), '..');
addpath(fullfile(root, 'functions'));
oldDir = cd(root); restoreDir = onCleanup(@() cd(oldDir));
run(fullfile(fileparts(mfilename('fullpath')), 'testParams.m'));   % frozen test settings, data/test/ (not pipelineParams)
tmp = fullfile(tempdir, 'testFrontEndFaults');
if ~isfolder(tmp), mkdir(tmp); end
cleanTmp = onCleanup(@() cleanDir(tmp));
fails = {};
template = gaussianTemplate(nBin, ephem.profileFWHM);
info_IQ = loadInfo(fileIQ);
nChan = round((fHigh - fLow) / chanWidth);
fC = (fLow + ((1:nChan).' - 0.5) * chanWidth) / 1e6;      % [MHz] observer's channel centres
bC = ((1:41) - 21) / 40 * chanWidth / 1e6;                % [MHz] in-channel offsets
tauD = @(f) 4.148808e3 * ephem.DM ./ f.^2;                % [s] dispersion delay, f in MHz
rx   = struct('fLow', fLow, 'fHigh', fHigh);
proc = struct('frontEnd', 'channels', 'refFreq', refFreq, 'f_out', f_out, 'nBin', nBin, ...
    'subintPeriods', subintPeriods, 'chanWidth', chanWidth, 'weighting', weighting, ...
    'excision', false, 'excisionArgs', {{}}, 'periodicMask', false, 'periodicArgs', {{}});
pF = struct('chanIQ', fullfile(tmp, 'c_chan'), 'chanBlank', fullfile(tmp, 'c_chan_blanked'), ...
    'chanDedisp', fullfile(tmp, 'c_dedisp_chan'), 'chanPower', fullfile(tmp, 'c_dedisp_chan_power.dat'), ...
    'chanDedispB', fullfile(tmp, 'c_dedispB_chan'), 'chanPowerB', fullfile(tmp, 'c_dedispB_chan_power.dat'), ...
    'chanWeight', fullfile(tmp, 'c_dedispB_chan_weights.dat'), 'foldChan', fullfile(tmp, 'c_fold_chan.mat'));
chain = @(info) runChain(info, rx, ephem, proc, pF);

% ---------------------------------------------------------------------------------
% 1. FM-8: the frequency mapping (mirrored band, LO offset)
% ---------------------------------------------------------------------------------
detR = loadInfo(fullfile(dataDir, 'chan', 'test_dedisp_chan_power.dat'));     % the frozen processing
[ifR, foldR] = foldProfile(detR, 'x', ephem.f0, 'F1', ephem.F1, 'TRef', ephem.TRef, 'NBin', nBin, ...
    'SubintPeriods', subintPeriods, 'NoiseCoeffs', [detR.noise.V, detR.noise.X], 'SaveFile', false, ...
    'Verbose', false);
[toaR, detR0] = foldToTOA(foldR, ifR, ephem.profileFWHM, detR.noise.Bnoise, 'Weighting', weighting, ...
    'Verbose', false);

fLO14 = 1.4e9;                                            % band centre: the mirror maps the band onto itself
recFile = fullfile(tmp, 'iq_lo14.dat'); mirFile = fullfile(tmp, 'iq_lo14_mirrored.dat');
shiftIQ(info_IQ.file, info_IQ.N, (fLO14 - info_IQ.fLO) / info_IQ.actualFsOut, recFile, mirFile);
infoA = info_IQ; infoA.file = recFile; infoA.fLO = fLO14;
infoB = infoA;   infoB.file = mirFile;                    % mirrored, the observer is not told
infoC = infoA;   infoC.fLO = fLO14 + 1e6;                 % the LO given 1 MHz too high
[toaA, ~]    = chain(infoA);
[toaB, detB] = chain(infoB);
[toaC, ~]    = chain(infoC);

% (a) the re-centred file
v = toaR.valid & toaA.valid;
dA = max(abs(toaA.toa(v) - toaR.toa(v)) ./ toaR.toaErr(v));
okA = isequal(toaA.valid, toaR.valid) && dA < 1e-3;
fprintf(['1a. FM-8 re-centred on 1.4 GHz: %d TOAs, max |TOA - frozen processing| %.1e sigma ' ...
         '(< 1e-3): %s\n'], nnz(v), dA, passStr(okA));
if ~okA, fails{end+1} = 'FM-8 re-centring'; end

% (b) the mirrored band
dMir = tauD(flipud(fC) - bC) - tauD(fC + bC);              % [nChan x nb] delays [s]
[fracB, tauB] = movedFit(template, ephem.f0 * dMir);
fracMeas = toaB.total.amp / toaA.total.amp;
zB = (toaB.total.amp - fracB * toaA.total.amp) / toaB.total.ampErr;
dPhase = mod(toaB.total.phase - toaA.total.phase - tauB + 0.5, 1) - 0.5;
okB = abs(zB) < 3 && abs(dPhase) < 0.01;
fprintf(['1b. FM-8 mirrored band: channel delays %.2f to %.2f ms (%.2f to %.2f turns); total fold: amplitude ' ...
         '%.3f of the re-centred (predicted %.3f; %+.1f of its error bars, < 3), phase %+.4f turns (predicted ' ...
         '%+.4f; difference %+.4f, < 0.01), SNR %.0f (re-centred %.0f): %s\n'], min(dMir(:))*1e3, ...
    max(dMir(:))*1e3, min(dMir(:))*ephem.f0, max(dMir(:))*ephem.f0, fracMeas, fracB, zB, ...
    mod(toaB.total.phase - toaA.total.phase + 0.5, 1) - 0.5, tauB, dPhase, toaB.total.snr, ...
    toaA.total.snr, passStr(okB));
if ~okB, fails{end+1} = 'FM-8 mirrored band'; end

% (c) the LO given 1 MHz too high
dLO = tauD(fC + bC - 1) - tauD(fC + bC);
predC = mean(dLO(:));
[~, tauC] = movedFit(template, ephem.f0 * dLO);
measC = toaC.total.timeOffset - toaA.total.timeOffset;
okC = abs(measC / predC - 1) < 0.05;
dSub = toaC.toa(toaC.valid) - toaA.toa(toaC.valid);
fprintf(['1c. FM-8 LO given 1 MHz high: channel delays %.1f to %.1f us; total fold %+.2f us later ' ...
         '(predicted %+.2f us, mean over channels; %+.2f us at the predicted fit''s peak), sub-ints ' ...
         '%+.2f +- %.2f us (TOA error %.2f us): %s\n'], min(dLO(:))*1e6, max(dLO(:))*1e6, measC*1e6, ...
    predC*1e6, tauC / ephem.f0 * 1e6, mean(dSub)*1e6, std(dSub)*1e6, median(toaA.toaErr(toaA.valid))*1e6, ...
    passStr(okC));
if ~okC, fails{end+1} = 'FM-8 LO offset'; end

% shown: the mirrored band's sub-ints, and B0329+54
g = toaB.good;
fprintf(['    shown, mirrored band: %d of %d sub-ints good (re-centred %d), their phases %s turns, ' ...
         'red. chi2 %.3f-%.3f (re-centred %.3f-%.3f), %d chi2-flagged; known-phase T0 at phase 0 ' ...
         'median %.1f (re-centred SNR %.0f per sub-int)\n'], nnz(g), nnz(toaB.valid), nnz(toaA.good), ...
    mat2str(round(toaB.phase(g).', 3)), min(toaB.redChi2(toaB.valid)), max(toaB.redChi2(toaB.valid)), ...
    min(toaA.redChi2(toaA.valid)), max(toaA.redChi2(toaA.valid)), nnz(toaB.flagChi2 & toaB.valid), ...
    median(detB.T0(toaB.valid)), median(toaA.snr(toaA.valid)));
comps = [0.159 -11.427 1.479; 0.439 -1.123 3.281; 0.585 0 1.055; 0.152 7.719 2.751; 0.170 9.482 0.847];
comps = [comps(:, 1), comps(:, 2) / 360, comps(:, 3) * 2*sqrt(2*log(2)) / 360];   % deg -> turns, FWHM
tauP = @(f) 4.148808e3 * 26.7641 ./ f.^2;                 % B0329+54's DM
f0P = 1.399541538720;
dP = tauP(flipud(fC) - bC) - tauP(fC + bC);
[fracP, tauPh] = movedFit(gaussianTemplate(nBin, comps), f0P * dP);
fprintf(['    shown, B0329+54 (DM 26.76, P 714.5 ms, 128 channels, LO at the band centre): channel ' ...
         'delays %.1f to %.1f ms; the fit keeps %.3f of the amplitude, at %+.4f turns (%+.1f ms): a ' ...
         'TOA of SNR 7 becomes SNR %.1f\n'], min(dP(:))*1e3, max(dP(:))*1e3, fracP, tauPh, tauPh/f0P*1e3, 7*fracP);

% ---------------------------------------------------------------------------------
% 2. Quantization: the ADC, or the transport, keeps b bits of I and of Q
% ---------------------------------------------------------------------------------
[sigIQ, kappa] = iqLevel(info_IQ, fLow, fHigh);           % rms of I (= of Q); in-band PSD / mean PSD
edgeCh = round((info_IQ.fLO - fLow) / chanWidth) + [0 1];  % the two channels that meet at the LO
qc = {'8-bit, rms 4 LSB', 8, 4, 'mid', true; '8-bit, rms 4 LSB, two''s complement', 8, 4, 'tc', true; ...
      '8-bit, rms 1 LSB', 8, 1, 'mid', true; '4-bit, rms 3 LSB', 4, 3, 'mid', true; ...
      '2-bit, rms 1 LSB', 2, 1, 'mid', false};         % last column: checked (false: shown only)
qFile = fullfile(tmp, 'iq_quantized.dat');
vR = toaR.valid;
betaR = median(toaR.amp(vR) ./ toaR.baseline(vR));       % pulse peak / baseline of the unquantized fold
[rawEdge, detMean] = deal(zeros(size(qc, 1), 2));
for i = 1:size(qc, 1)
    [lab, bits, lev, kind, chk] = qc{i, :};
    quantizeIQ(info_IQ.file, info_IQ.N, sigIQ / lev, bits, kind, qFile);
    infoQ = info_IQ; infoQ.file = qFile;
    [toaQ, detQ, outQ] = chain(infoQ);
    rawEdge(i, :) = rawChanPower(pF.chanIQ, edgeCh, outQ.info_chan.fullySupported);
    detMean(i, :) = outQ.info_det.chanMeanPower(edgeCh);
    [epsQ, alphaQ] = quantNoise(lev, bits, kappa);
    etaP = snrFactor(template, betaR, epsQ);
    etaM = toaQ.total.snr / toaR.total.snr;
    v = toaQ.valid & vR;
    dChi = mean(toaQ.redChi2(v)) - mean(toaR.redChi2(v));
    ok = abs(etaM - etaP) < 0.005 && abs(dChi) < 0.01;
    verdict = passStr(ok);
    if ~chk, verdict = 'shown'; end
    fprintf(['2. quantization, %s: SNR x %.4f (white-error model %.4f, < 0.005 apart; Bussgang gain %.4f, ' ...
             'error %.4f of the in-band noise, so x %.4f for a weak pulsar); red. chi2 %+.4f against the ' ...
             'unquantized fold (< 0.01), noise ratio %.3f (unquantized %.3f); TOAs move by up to %.2f ' ...
             'sigma, error bars x %.3f: %s\n'], lab, etaM, etaP, alphaQ, epsQ, 1 / (1 + epsQ), dChi, ...
        median(detQ.noiseRatio(v)), median(detR0.noiseRatio(v)), ...
        max(abs(toaQ.toa(v) - toaR.toa(v)) ./ toaR.toaErr(v)), median(toaQ.toaErr(v) ./ toaR.toaErr(v)), ...
        verdict);
    if chk && ~ok, fails{end+1} = ['quantization ' lab]; end %#ok<AGROW>
end
% the two's-complement offset: -1/2 LSB on I and on Q, a tone at the LO, which sits on the
% edge between channels edgeCh: in their channel files at full power (|DC|^2 = 1/2 LSB^2), and
% gone after dedispersion (its sin^2 taper is 0 at the channel edge)
dcRaw = rawEdge(2, :) - rawEdge(1, :);
dcDet = abs(detMean(2, :) ./ detMean(1, :) - 1);
ok = all(abs(dcRaw / 0.5 - 1) < 0.05) && all(dcDet < 1e-3);
fprintf(['2. two''s-complement DC at the LO (%.4g GHz, between channels %d and %d): channel files +%s LSB^2 ' ...
         '(predicted 0.5), dedispersed channels changed by %s (< 1e-3): %s\n'], info_IQ.fLO/1e9, edgeCh, ...
    mat2str(dcRaw, 4), mat2str(dcDet, 2), passStr(ok));
if ~ok, fails{end+1} = 'quantization DC offset'; end

% ---------------------------------------------------------------------------------
% 3. FM-9: a lost buffer of samples
% ---------------------------------------------------------------------------------
infoCh = loadInfo(fullfile(dataDir, 'chan', 'test_rx_IQ_chan'));
binIn = infoCh.decimation * detR.binLen;                  % input samples per detected bin (768)
nDrop = 43 * binIn;                                       % 33024 samples, 41.28 us
fsIQ  = info_IQ.actualFsOut;
gStart = binIn * round(45.3e-3 * fsIQ / binIn);           % first lost sample (0-based)
tDrop = info_IQ.t0 + gStart / fsIQ;
sweep = tauD(fLow / 1e6) - tauD(fHigh / 1e6);             % [s] across the band
dropFile = fullfile(tmp, 'iq_dropped.dat'); gapFile = fullfile(tmp, 'iq_gap.dat');
cutIQ(info_IQ.file, info_IQ.N, gStart, nDrop, dropFile, gapFile);
infoD = info_IQ; infoD.file = dropFile; infoD.N = info_IQ.N - nDrop;
infoG = info_IQ; infoG.file = gapFile;
toaD = chain(infoD);
toaG = chain(infoG);
span = @(t) [t.tRef - 1/ephem.f0, t.tRef + 1/ephem.f0];  % a sub-int's data lie within (1 turn each)

% (a) joined without the gap: a clock step
[iR, iD] = commonSub(toaR, toaD);
sR = span(toaR); sR = sR(iR, :);
pre  = sR(:, 2) < tDrop - sweep - 2e-4;                   % all channels' data before the drop
post = sR(:, 1) > tDrop + 2e-4;                           % all after it
dT = (toaD.toa(iD) - toaR.toa(iR)) ./ toaR.toaErr(iR);
step = -nDrop / fsIQ ./ toaR.toaErr(iR);
ok = any(pre) && any(post) && all(abs(dT(pre)) < 1e-3) && all(abs(dT(post) - step(post)) < 0.1);
fprintf(['3a. FM-9 %d samples (%.2f us) dropped at %.2f ms, the stream joined: %d sub-ints before the ' ...
         'drop move by up to %.1e sigma (< 1e-3), %d after it by %s sigma (predicted -n/fs = %s, within ' ...
         '0.1), %d straddling: %s %s\n'], nDrop, nDrop / fsIQ * 1e6, tDrop * 1e3, nnz(pre), ...
    max(abs(dT(pre))), nnz(post), mat2str(dT(post).', 3), mat2str(step(post).', 3), ...
    nnz(~pre & ~post), mat2str(dT(~pre & ~post).', 3), passStr(ok));
if ~ok, fails{end+1} = 'FM-9 dropped samples'; end

% (b) the gap filled with zeros: shown without a mask; checked with the gap blanked
[iR, iG] = commonSub(toaR, toaG);
sR = span(toaR); sR = sR(iR, :);
hole = sR(:, 2) > tDrop - sweep - 2e-4 & sR(:, 1) < tDrop + nDrop / fsIQ + 2e-4;   % sub-ints with the hole
dTg = (toaG.toa(iG) - toaR.toa(iR)) ./ toaR.toaErr(iR);
dCg = toaG.redChi2(iG) - toaR.redChi2(iR);
fprintf(['3b. FM-9 the same gap filled with zeros, no mask (shown): sub-ints with the hole move by %s ' ...
         'sigma, red. chi2 %s; the others by up to %.1e sigma\n'], mat2str(dTg(hole).', 3), ...
    mat2str(dCg(hole).', 3), max(abs(dTg(~hole))));
gapRows = @(ic) [(1:ic.nChan).', repmat([ceil((gStart - ic.groupDelayRemoved) / ic.decimation) + 1, ...
    floor((gStart + nDrop - 1 + ic.groupDelayRemoved) / ic.decimation) + 1], ic.nChan, 1)];
[toaE, fracE] = maskChain(info_IQ, [], rx, ephem, proc, tmp);
[toaM, fracM, nGapS] = maskChain(infoG, gapRows, rx, ephem, proc, tmp);
[iE, iM] = commonSub(toaE, toaM);
sE = span(toaE); sE = sE(iE, :);
holeE = sE(:, 2) > tDrop - sweep - 2e-4 & sE(:, 1) < tDrop + nDrop / fsIQ + 2e-4;
dTm = (toaM.toa(iM) - toaE.toa(iE)) ./ toaE.toaErr(iE);
dCm = toaM.redChi2(iM) - toaE.redChi2(iE);
ok = any(holeE) && all(abs(dTm) < 0.3) && all(abs(dCm(holeE)) < 0.005);
fprintf(['3b. FM-9 the gap blanked (%d channel samples per channel, the filter''s reach included): ' ...
         'blanked fraction %.4f %% -> %.4f %%; sub-ints with the hole move by %s sigma, red. chi2 %s ' ...
         '(< 0.005); all sub-ints within 0.3 sigma (max %.3f): %s\n'], nGapS, 100 * fracE, 100 * fracM, ...
    mat2str(dTm(holeE).', 3), mat2str(dCm(holeE).', 3), max(abs(dTm)), passStr(ok));
if ~ok, fails{end+1} = 'FM-9 blanked gap'; end

% ---------------------------------------------------------------------------------
% 4. FM-39: absolute time stamps
% ---------------------------------------------------------------------------------
tAbs = 8.3e8 + 0.123456789;                               % [s] TT since J2000, in 2026
infoT = info_IQ; infoT.t0 = tAbs;
ephT = ephem; ephT.TRef = ephem.TRef + tAbs;              % the same pulse phases
outT = processIQ(infoT, rx, ephT, proc, pF, 'Verbose', false);
toaT = foldToTOA(outT.fold, outT.info_fold, ephT.profileFWHM, outT.Bnoise, 'Weighting', weighting, ...
    'Verbose', false);
[iR, iT] = commonSub(toaR, toaT);
dAbs = (toaT.toa(iT) - tAbs) - toaR.toa(iR);
ok = numel(iT) == nnz(toaR.valid) && max(abs(dAbs)) <= 4 * eps(tAbs);
fprintf(['4. FM-39 absolute time (t0 = %.9f s): TOAs - t0 against t0 = 0 differ by up to %.3f us (< 4 ' ...
         'ULP of t0 = %.3f us; %.3f of a sub-int''s TOA error, %.2f us); phases %.1e turns apart: %s\n'], ...
    tAbs, max(abs(dAbs)) * 1e6, 4 * eps(tAbs) * 1e6, max(abs(dAbs) ./ toaR.toaErr(iR)), ...
    median(toaR.toaErr(iR)) * 1e6, max(abs(toaT.phase(iT) - toaR.phase(iR))), passStr(ok));
if ~ok, fails{end+1} = 'FM-39 absolute time'; end

if isempty(fails)
    fprintf('testFrontEndFaults: ALL PASSED\n');
else
    error('testFrontEndFaults:failed', 'FAILED: %s', strjoin(fails, ', '));
end
end


% =================================================================================
function [toa, det, out] = runChain(info_IQ, rx, ephem, proc, pF)
%RUNCHAIN  The observer's chain on one IQ file: processIQ (channel path), then foldToTOA.
out = processIQ(info_IQ, rx, ephem, proc, pF, 'Verbose', false);
[toa, det] = foldToTOA(out.fold, out.info_fold, ephem.profileFWHM, out.Bnoise, ...
    'Weighting', out.weighting, 'Verbose', false);
end

function shiftIQ(inFile, N, cyc, outFile, mirFile)
%SHIFTIQ  x(n) exp(-2 pi i cyc n) into outFile (the band moves down by cyc*fs: a new LO cyc*fs
% higher, the same RF), and its complex conjugate (the mirrored spectrum) into mirFile. cyc must
% be 1/8 here, so that the rotation repeats every 8 samples and is exact.
rot = single(exp(-2i*pi*cyc*(0:7).'));
if abs(cyc - 1/8) > 1e-12, error('testFrontEndFaults:shift', 'Only a shift of fs/8.'); end
fi = fopen(inFile, 'r', 'ieee-le'); ci = onCleanup(@() fclose(fi));
fo = fopen(outFile, 'w', 'ieee-le'); co = onCleanup(@() fclose(fo));
fm = fopen(mirFile, 'w', 'ieee-le'); cm = onCleanup(@() fclose(fm));
blk = 2^23; n0 = 0;
while n0 < N
    nb = min(blk, N - n0);
    raw = fread(fi, [2, nb], 'single=>single');
    z = complex(raw(1, :), raw(2, :)) .* rot(mod(n0 + (0:nb-1), 8) + 1).';
    fwrite(fo, [real(z); imag(z)], 'single');
    fwrite(fm, [real(z); -imag(z)], 'single');
    n0 = n0 + nb;
end
end

function [frac, tau] = movedFit(T, s)
%MOVEDFIT  The fit's response to channels whose pulse is moved: channel j's profile is the
% template T moved by s(j, :) turns (averaged over the columns, the in-channel offsets), and the
% fit, one amplitude and one phase for all channels with equal weights, correlates the mean
% channel profile with c = T - mean(T). frac = max over tau of <P, c_tau> / <c, c> (the amplitude
% as a fraction of the unmoved one), tau [turns] where it peaks (64 x finer than the bins).
N = numel(T); U = 64;
k = [0:ceil(N/2)-1, -floor(N/2):-1].';
c = T(:) - mean(T);
S = mean(exp(-2i*pi*k.' .* reshape(s, [], 1, size(s, 2))), 3);   % [nChan x N]: mean over b
Y = fft(T(:)) .* mean(S, 1).' .* conj(fft(c)) / N;                 % <P, c_tau> in Fourier space
Ypad = zeros(N*U, 1);
Ypad(mod(k, N*U) + 1) = Y;
x = real(ifft(Ypad)) * N * U;
[xm, m] = max(x);
frac = xm / sum(c.^2);
tau = mod((m - 1) / (N*U) + 0.5, 1) - 0.5;
end

function [toa, frac, nGap] = maskChain(info_IQ, extraRows, rx, ephem, proc, tmp)
%MASKCHAIN  processIQ's excision branch stage by stage, with extra mask rows: extraRows(info_chan)
% -> [M x 3] rows added to detectRFI's mask ([] = none). frac: mean blanked fraction over the
% channels; nGap: channel samples per channel in the extra rows.
info_chan = channelizeIQ(info_IQ.file, fullfile(tmp, 'm_chan'), info_IQ.actualFsOut, info_IQ.fLO, ...
    rx.fLow, rx.fHigh, 'ChanWidth', proc.chanWidth, 'T0', info_IQ.t0, 'Verbose', false);
mask = detectRFI(info_chan, 'Verbose', false);
nGap = 0;
if ~isempty(extraRows)
    rows = extraRows(info_chan);
    mask = [mask; rows];
    nGap = rows(1, 3) - rows(1, 2) + 1;
end
info_chanD = blankChannels(info_chan, mask, fullfile(tmp, 'm_chan_blanked'), 'Verbose', false);
info_dc = dedisperseChannels(info_chanD, fullfile(tmp, 'm_dedisp'), ephem.DM, 'RefFreq', proc.refFreq, ...
    'Verbose', false);
info_det = detectChannels(info_dc, fullfile(tmp, 'm_power.dat'), info_dc.fs / round(info_dc.fs / proc.f_out), ...
    'Verbose', false);
info_w = blankingWeights(info_chan, info_dc, info_det, info_chanD.mask, fullfile(tmp, 'm_weights.dat'), ...
    'Verbose', false);
[info_fold, fold] = foldProfile(info_det, fullfile(tmp, 'm_fold.mat'), ephem.f0, 'F1', ephem.F1, ...
    'TRef', ephem.TRef, 'NBin', proc.nBin, 'SubintPeriods', proc.subintPeriods, 'DataWeights', info_w, ...
    'Verbose', false);
toa = foldToTOA(fold, info_fold, ephem.profileFWHM, info_det.noise.Bnoise, 'Weighting', proc.weighting, ...
    'Verbose', false);
frac = mean(info_chanD.blankedFraction);
end

function cutIQ(inFile, N, g, n, dropFile, gapFile)
%CUTIQ  Samples g .. g+n-1 (0-based) lost: dropFile without them (the stream joined), gapFile
% with zeros in their place (time stamps kept).
fi = fopen(inFile, 'r', 'ieee-le');  ci = onCleanup(@() fclose(fi));
fd = fopen(dropFile, 'w', 'ieee-le'); cd_ = onCleanup(@() fclose(fd));
fg = fopen(gapFile, 'w', 'ieee-le');  cg = onCleanup(@() fclose(fg));
blk = 2^23; n0 = 0;
while n0 < N
    nb = min(blk, N - n0);
    raw = fread(fi, [2, nb], 'single=>single');
    lost = (n0 + (0:nb-1)) >= g & (n0 + (0:nb-1)) < g + n;
    fwrite(fd, raw(:, ~lost), 'single');
    raw(:, lost) = 0;
    fwrite(fg, raw, 'single');
    n0 = n0 + nb;
end
end

function [ia, ib] = commonSub(a, b)
%COMMONSUB  Indices of the sub-ints valid in both TOA structs, matched by their reference turn.
ta = a.turnRef; ta(~a.valid) = NaN;
tb = b.turnRef; tb(~b.valid) = NaN;
[~, ia, ib] = intersect(ta, tb);
end

function [sig, kappa] = iqLevel(info_IQ, fLow, fHigh)
%IQLEVEL  rms of I (and of Q) from the first 2^23 samples, and kappa = the mean PSD in the RF band
% over the mean PSD of the whole sample rate (averaged 1024-point periodograms).
n = min(2^23, floor(info_IQ.N / 1024) * 1024);
fid = fopen(info_IQ.file, 'r', 'ieee-le'); c = onCleanup(@() fclose(fid));
raw = fread(fid, [2, n], 'single=>double');
sig = sqrt(mean(raw(:).^2));
P = mean(abs(fft(reshape(complex(raw(1, :), raw(2, :)), 1024, []))).^2, 2);
fb = ((0:1023).' - 1024 * ((0:1023).' >= 512)) / 1024 * info_IQ.actualFsOut;   % signed bin frequency
in = fb >= fLow - info_IQ.fLO & fb <= fHigh - info_IQ.fLO;
kappa = mean(P(in)) / mean(P);
end

function quantizeIQ(inFile, N, lsb, bits, kind, outFile)
%QUANTIZEIQ  I and Q each to 2^bits levels of step lsb (input units), written in LSB units:
% 'mid' mid-rise (floor + 1/2, symmetric, clipped at +-(2^(bits-1) - 1/2)); 'tc' two's
% complement (floor, clipped to -2^(bits-1) .. 2^(bits-1) - 1: the mid-rise levels - 1/2).
L = 2^(bits - 1);
fi = fopen(inFile, 'r', 'ieee-le');  ci = onCleanup(@() fclose(fi));
fo = fopen(outFile, 'w', 'ieee-le'); co = onCleanup(@() fclose(fo));
blk = 2^23; n0 = 0;
while n0 < N
    nb = min(blk, N - n0);
    y = floor(fread(fi, [2, nb], 'single=>double') / lsb);
    if strcmp(kind, 'mid'), y = min(max(y + 0.5, -L + 0.5), L - 0.5);
    else,                   y = min(max(y, -L), L - 1);
    end
    fwrite(fo, y, 'single');
    n0 = n0 + nb;
end
end

function [epsQ, alpha] = quantNoise(sig, bits, kappa)
%QUANTNOISE  In-band quantization noise of a mid-rise quantizer (2^bits levels, step 1) on a
% Gaussian of rms sig per component, as a fraction epsQ of the receiver noise. Bussgang:
% Q(x) = alpha x + e, alpha = E[x Q] / sig^2, e uncorrelated with x, of variance
% E[Q^2] - alpha^2 sig^2 and white over the sample rate; the receiver noise has kappa times the
% mean PSD in the band. So epsQ = var(e) / (alpha^2 kappa sig^2).
L = 2^(bits - 1);
t = [-Inf, -L+1:L-1, Inf];                                % thresholds
q = (-L:L-1) + 0.5;                                       % levels
pdf = @(x) exp(-x.^2 / (2*sig^2)) / (sig*sqrt(2*pi));
cdf = @(x) 0.5 * erfc(-x / (sig*sqrt(2)));
ExQ = sum(q .* sig^2 .* (pdf(t(1:end-1)) - pdf(t(2:end))));
EQ2 = sum(q.^2 .* (cdf(t(2:end)) - cdf(t(1:end-1))));
alpha = ExQ / sig^2;
epsQ = (EQ2 - alpha^2 * sig^2) / (alpha^2 * kappa * sig^2);
end

function r = snrFactor(T, beta, epsQ)
%SNRFACTOR  The amplitude's SNR with extra white noise epsQ (in units of the baseline) over the
% SNR without: weighted least squares of y = a + b T, per-bin variance (1 + epsQ + beta T)^2 (the
% radiometer model, the pulse's own noise included; beta = pulse peak / baseline). For a weak
% pulse (beta -> 0) r = 1 / (1 + epsQ).
T = T(:);
sb = @(e) sqrt(sum(1 ./ (1 + e + beta*T).^2) / (sum(1 ./ (1 + e + beta*T).^2) * ...
    sum(T.^2 ./ (1 + e + beta*T).^2) - sum(T ./ (1 + e + beta*T).^2)^2));
r = sb(0) / sb(epsQ);
end

function p = rawChanPower(base, ch, sup)
%RAWCHANPOWER  Mean |x|^2 of channelizeIQ's channel files ch (before dedispersion), over the
% fully supported samples sup = [first last].
p = zeros(1, numel(ch));
for i = 1:numel(ch)
    fid = fopen(sprintf('%s_ch%03d.dat', base, ch(i)), 'r', 'ieee-le');
    fseek(fid, (sup(1) - 1) * 8, 'bof');
    raw = fread(fid, [2, sup(2) - sup(1) + 1], 'single=>double');
    fclose(fid);
    p(i) = mean(sum(raw.^2, 1));
end
end

function cleanDir(d)
delete(fullfile(d, '*.dat')); delete(fullfile(d, '*.mat'));
end

function s = passStr(p)
if p, s = 'PASS'; else, s = 'FAIL'; end
end
