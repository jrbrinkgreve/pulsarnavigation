function testDetectRFI()
%TESTDETECTRFI  Unit tests for detectRFI (B2): false flags, the bright pulsar, real RFI.
%{
Run from the PulsarSimMatlab folder: run('tests/testDetectRFI.m').
Makes its own noise-only and noise + RFI channel files (0.05 s through the
real receiver chain: addNoiseAndRFI on a zero sky signal -> applyIQmodulation
-> channelizeIQ; the same noise seed, so noise + RFI minus noise only is the
RFI alone, channel by channel; ~2.5 GB in tempdir, deleted on the way), and
uses the seed-43 channel files (data/test/chan/test_rx_IQ_chan_*: -5 dB pulsar,
no RFI) for test 2.

  1. Noise only:
     (a) baseline (median / ln 2 per block) = the block's mean power;
     (b) flagged windows per scale at PFA 1e-3 and 1e-4 (one pass) =
         nWindows * PFA, i.e. the thresholds (exact distribution of the
         window power, not the Gamma approximation) are right, within 4 sigma
         (sigma from the scatter of the counts between channels:
         overlapping windows cluster);
     (c) default settings (PFA 1e-6, two passes): flagged windows as
         expected, and the blanked fraction of clean data (incl. guard) at
         most sum nWindows * PFA * (n + 2*guard) / N (overlaps ignored).
  2. The -5 dB pulsar: flagged windows ON the pulse (|t - t_pulse| < FWHM)
     and OFF it (> 3 FWHM), per scale, vs the prediction from the same Gamma
     statistics (exact distribution) with the pulse power 1 + rho*p(t) (per channel: dispersion
     delay, intra-channel smearing) and the detector's own baseline; only
     channels outside the dispersion stage's 8 MHz band-edge taper (there
     the simulated pulse is weaker);
     (a) the default windows 1-16; (b) windows up to 256 samples (the
     reason for the default: the on-pulse blanked fraction is shown).
  3. RFI: radar 1300 MHz +20 dB (scenario), broadband impulses 1000/s
     +20 dB, carrier 1351.3 MHz -5 dB, GNSS L1 -10 dB; scored with the
     RFI-only signal:
     (a) radar: every pulse's samples in channels 32-33 blanked; fraction of
         the radar energy removed (channels 31-34, +-20 us) > 0.9999. Shown:
         the same over all channels: the pulse's spectrum has sidelobes
         (~1/f^2 near the carrier; with the default 0.1 us edges ~1/f^6
         beyond ~5 MHz, with instant edges 1/f^2 over tens of MHz); the
         detector finds them where they stand out per channel;
     (b) impulses: fraction of their energy removed (other channels, +-15
         samples, away from the radar pulses) > 0.95. Expected ~0.98: per channel a burst carries ~80
         times a sample's noise energy, but ~exponentially distributed, so
         the weakest ~15 % stay below the thresholds;
     (c) constant-envelope RFI is not flagged: in the carrier and GNSS
         channels (49, 120, 121) no more flags away from the impulses than
         noise would give; their baseline vs noise only shown.
     (d) frequency guard: the mask equals the own detections plus every
         strong event (peak >= FreqGuardMin x baseline) copied to +-FreqGuard
         channels at the same samples (recomputed independently from
         coverage); no strong events on noise or from the -5 dB pulsar.
     Shown: channels flagged at >= half of the radar pulses, and the radar
     energy left over all channels with and without the frequency guard
     (each channel's steady RFI level, mainly the far-reaching GNSS
     sidelobes, subtracted first).
  4. Continuous narrowband RFI, per-channel thresholds (21 Oct 2026): three
     200 kHz noise-like carriers (channels 60, 70, 80; INR -21 dB as the
     narrow-Inmarsat case: about the channel's noise power) and a +20 dB
     radar on top of the channel-70 carrier:
     (a) at PFA 1e-3 (two passes) the flagged windows in the two radar-free
         carrier channels match nWindows * PFA per scale (4 sigma) with their
         own thresholds; the white thresholds flag > 3x as many over all
         window lengths (at n = 1 a sample's power is exponential whatever
         the correlation, so there only the pass-2 baseline differs; shown:
         the counts at the default PFA, and the thresholds);
     (b) every radar pulse's samples in channel 70 are still blanked;
     (c) channels away from the RFI (outside 55-85) keep the white thresholds
         and give the same mask bit for bit as 'ChannelThresholds' false, and
         no channel of the noise-only data (test 1c) takes its own.
  5. Gain changes (27 Oct 2026, review R-E1): the noise-only channels of test 1 with their
     power x G(t), one pass, PFA 1e-3: a ramp of +10 % over the 0.05 s file (3.3 % per 16.7 ms
     baseline block, far faster than a real receiver drifts) and a step of +5 % in the middle
     of block 2.
     (a) the baselines follow the gain: gained / clean baseline (the same noise) = each
         block's effective gain, the median of the gained exponential power (solved here),
         within 4 sigma of the scatter between channels;
     (b) flagged windows per scale = sum over windows of S_n(eta * n * G_eff / G(t)), S_n the
         exact survival function (the run's eigenvalues), within 4 sigma: false flags stay
         near nWindows * PFA (the small excess is printed);
     (c) ramp: in the last tenth of each block (gain above the baseline) and the first tenth
         (below) the counts match the same prediction; their ratio is printed.
  6. One pass against two under dense bursts (R-E1), on the same noise: strong bursts (+20 dB,
     4 samples, 3 % of the samples) and weak ones (+3 dB, 16 samples, 1 %) at random times,
     'ChannelThresholds' false (the pass logic alone):
     (a) the pass-1 baseline / clean = the median of the mixture (solved here per channel and
         block from the burst samples), within 4 sigma;
     (b) the pass-2 baseline / clean = the median of what pass 1 left unflagged (the weak
         bursts it missed included), within 4 sigma;
     (c) strong bursts wholly blanked by the channel's own detections: at least 99.9 % in
         both passes (a Gaussian burst fades: its total power misses every threshold with
         probability ~1e-4, Gamma estimate printed);
     (d) weak bursts wholly blanked by the own detections: more in pass 2 than in pass 1
         (shown: the blanked fraction with and without the frequency guard, which copies
         each channel's strong bursts into its neighbours: these bursts are independent
         per channel, unlike broadband RFI);
     (e) false windows in clean stretches (no burst within the guard), each pass = the
         prediction from that pass's baseline, within 4 sigma;
     (f) the defaults ('ChannelThresholds' true): strong bursts as in (c), weak ones at
         least as well as in pass 1; shown: channels with their own thresholds, the baseline.
Uses its own temporary folder per run (21 Oct 2026).
Errors at the end if any check fails.
%}

root = fullfile(fileparts(mfilename('fullpath')), '..');
addpath(fullfile(root, 'functions'));
oldDir = cd(root); restoreDir = onCleanup(@() cd(oldDir));
run(fullfile(fileparts(mfilename('fullpath')), 'testParams.m'));   % frozen test settings, data/test/ (not pipelineParams)
tmp = tempname(tempdir);                         % own folder per run (21 Oct): two runs at once
mkdir(tmp);                                      %   (e.g. Jasper's and Claude's) cannot collide
cleanTmp = onCleanup(@() rmdir(tmp, 's'));
fails = {};

% ---------------------------------------------------------------------------------
% Data: noise only and noise + RFI through the receiver chain (0.05 s)
% ---------------------------------------------------------------------------------
tD = tic;
zeroFile = fullfile(tmp, 'zero.dat');
writeZeros(zeroFile, round(0.05 * f_in));
rfiT = [ ...
    rfiSource('pulsed',  'Freq', 1300e6, 'PulseWidth', 2e-6, 'PRF', 373, 'ChirpBW', 1e6, ...
              'INRdB', 20, 'Label', 'L-band radar'), ...
    rfiSource('impulse', 'Rate', 1000, 'Duration', 200e-9, 'INRdB', 20, 'Label', 'impulses'), ...
    rfiSource('cw',      'Freq', 1351.3e6, 'INRdB', -5, 'Label', 'carrier'), ...
    rfiSource('bpsk',    'Freq', 1575.42e6, 'ChipRate', 1.023e6, 'INRdB', -10, 'Label', 'GNSS L1')];
rx = struct('f_in', f_in, 'fs', fs, 'fLO', fLO, 'band', [fLow fHigh], 'order', filterOrder);
chan0 = makeChannels(zeroFile, fullfile(tmp, 'n'), struct([]), rx);
[chan1, info_rx1] = makeChannels(zeroFile, fullfile(tmp, 'r'), rfiT, rx);
c4  = [60 70 80];                                   % test 4: narrow carriers, radar on the middle one
fc4 = fLow + (c4 - 0.5) * chanWidth + [-1.1 0.4 1.2] * 1e6;
rfi4 = arrayfun(@(f) rfiSource('noise', 'Freq', f, 'Bandwidth', 200e3, 'INRdB', -21, ...
    'Label', 'narrow carrier'), fc4, 'UniformOutput', false);
rfi4 = [rfi4{:}, rfiSource('pulsed', 'Freq', fc4(2), 'PulseWidth', 2e-6, 'PRF', 373, ...
    'ChirpBW', 1e6, 'INRdB', 20, 'Label', 'radar on a carrier')];
[chan4, info_rx4] = makeChannels(zeroFile, fullfile(tmp, 'c'), rfi4, rx);
delete(zeroFile);
nCh = chan0.nChan; Nc = chan0.N; fsC = chan0.fs; sup = chan0.fullySupported;
fprintf('data: noise-only and noise + RFI channel files, %d x %d samples (%.0f s)\n', ...
    nCh, Nc, toc(tD));

% ---------------------------------------------------------------------------------
% 1. Noise only
% ---------------------------------------------------------------------------------
% (a) baseline vs the block mean power
[~, d1] = detectRFI(chan0, 'PFA', 1e-4, 'Passes', 1, 'Verbose', false);
e = d1.blockEdges; nBlk = numel(e) - 1;
bm = zeros(nCh, nBlk);
for j = 1:nCh
    p = readPower(chan0.chanFiles(j));
    for b = 1:nBlk
        bm(j, b) = mean(p(max(e(b) + 1, sup(1)) : min(e(b + 1), sup(2))));
    end
end
ratio = d1.baseline(:) ./ bm(:);
pass = abs(mean(ratio) - 1) < 4 * std(ratio) / sqrt(numel(ratio));
fprintf('1a. baseline / block mean power %.5f +- %.5f (scatter per block %.4f, %d blocks): %s\n', ...
    mean(ratio), std(ratio) / sqrt(numel(ratio)), std(ratio), numel(ratio), passStr(pass));
if ~pass, fails{end+1} = 'baseline'; end

% (b) the thresholds: false windows at PFA 1e-3, 1e-4 (one pass: unbiased baseline)
for pfa = [1e-3 1e-4]
    [~, d] = detectRFI(chan0, 'PFA', pfa, 'Passes', 1, 'Verbose', false);
    [ok, r, sr] = countCheck(d.nFlaggedWindows, d.expectedFalse);
    fprintf('1b. PFA %.0e, one pass: flagged / expected windows per scale %s (+- %s): %s\n', ...
        pfa, mat2str(round(r, 3)), mat2str(round(sr, 3)), passStr(ok));
    if ~ok, fails{end+1} = sprintf('false windows %.0e', pfa); end %#ok<AGROW>
end

% (c) defaults
[~, dD] = detectRFI(chan0, 'Verbose', false);
[okW, r, sr] = countCheck(dD.nFlaggedWindows, dD.expectedFalse);
bound = sum(dD.nWindows .* dD.pFA .* (dD.scales + 2*dD.guard)) / Nc;
frac  = mean(dD.flaggedFraction);
okFr  = frac <= bound * (1 + 4 / sqrt(sum(dD.expectedFalse)));
fprintf(['1c. defaults (PFA %.0e, %d passes, guard %d): flagged / expected windows %s (+- %s); ' ...
         'blanked fraction of clean data %.2e (at most %.2e expected): %s\n'], dD.pFA, ...
    dD.passes, dD.guard, mat2str(round(r, 2)), mat2str(round(sr, 2)), frac, bound, passStr(okW && okFr));
if ~(okW && okFr), fails{end+1} = 'defaults on noise'; end

% ---------------------------------------------------------------------------------
% 2. The -5 dB pulsar (seed-43 channel files)
% ---------------------------------------------------------------------------------
t2 = tic;
chanP = loadInfo(fullfile(dataDir, 'chan', 'test_rx_IQ_chan'));
info_rx = loadInfo(fileRx);
if ~(isempty(info_rx.rfi) && isfinite(info_rx.snrDB))
    error('testDetectRFI:data', 'Test 2 needs the seed-43 files without RFI (data/test/, tests/makeTestData.m).');
end
info_gen = loadInfo(fileRaw);
info_disp = loadInfo(fileDispersed);
rhoP = 10^(info_rx.snrDB / 10);
fwhm = info_gen.FWHM_power;                          % power profile (ground truth)
sig  = fwhm / (2 * sqrt(2 * log(2)));
Kdm  = 4.148808e3 * 1e12 * ephem.DM;                 % [s Hz^2]
NcP = chanP.N; fsP = chanP.fs; t0P = chanP.t0; dF = chanP.chanWidth;
% channels with the full pulse: outside the dispersion stage's band-edge taper
useCh = find(chanP.chanFreqs - dF/2 >= info_disp.fLow + info_disp.edgeWidth & ...
             chanP.chanFreqs + dF/2 <= info_disp.fHigh - info_disp.edgeWidth).';
inS = false(1, NcP); inS(chanP.fullySupported(1):chanP.fullySupported(2)) = true;
tS  = t0P + (0:NcP - 1) / fsP;
runs = cell(1, 2); masks = cell(1, 2);
[masks{1}, runs{1}] = detectRFI(chanP, 'KeepWindows', true, 'Verbose', false);
[masks{2}, runs{2}] = detectRFI(chanP, 'Scales', 2.^(0:8), 'KeepWindows', true, 'Verbose', false);
nR = numel(runs); nSmax = 9;
lsv = cell(nR, nSmax);                               % log survival of the exact distribution
for r = 1:nR
    d = runs{r};
    for s = 1:numel(d.scales)
        xs = d.threshold(s) * d.scales(s) * linspace(0.6, 1.1, 200);
        ls = arrayfun(@(x) log(sumExpSurvival(d.lambda{s}, x)), xs);
        lsv{r, s} = griddedInterpolant(xs, ls, 'pchip');
    end
end
mOn = zeros(nR, nSmax); mOff = mOn; pOn = mOn; pOff = mOn;
bOn = zeros(1, nR); bOff = bOn; nOn = 0; nOff = 0;
for j = useCh
    fc  = chanP.chanFreqs(j);
    tP0 = ephem.TRef + Kdm * (1/fc^2 - 1/refFreq^2);  % pulse 0 in this channel
    dist = @(t) abs(t - tP0 - round((t - tP0) * ephem.f0) / ephem.f0);   % to the nearest pulse
    sweep = Kdm * abs(1/(fc - dF/2)^2 - 1/(fc + dF/2)^2);
    sE  = sqrt(sig^2 + sweep^2 / 12);                % smeared within the channel
    amp = rhoP * sig / sE;
    p  = readPower(chanP.chanFiles(j));
    dS = dist(tS);
    onS = inS & dS < fwhm; offS = inS & dS > 3*fwhm;
    m  = mean(p(offS));                              % noise level of this channel
    nOn = nOn + nnz(onS); nOff = nOff + nnz(offS);
    for r = 1:nR
        d = runs{r};
        blk = blockIndex(d.blockEdges);
        for s = 1:numel(d.scales)
            n  = d.scales(s);
            i0 = 1:d.stride(s):NcP - n + 1;
            dc = dist(t0P + (i0 - 1 + (n - 1)/2) / fsP);
            rr = 1 + amp * exp(-dc.^2 / (2 * sE^2));
            P  = exp(lsv{r, s}(d.threshold(s) * n * d.baseline(j, blk(i0)) ./ (m * rr)));
            pOn(r, s)  = pOn(r, s)  + sum(P(dc < fwhm));
            pOff(r, s) = pOff(r, s) + sum(P(dc > 3*fwhm));
            ds = dist(t0P + (d.flaggedWindows{j, s} - 1 + (n - 1)/2) / fsP);
            mOn(r, s)  = mOn(r, s)  + nnz(ds < fwhm);
            mOff(r, s) = mOff(r, s) + nnz(ds > 3*fwhm);
        end
        cv = coveredOf(masks{r}, j, NcP);
        bOn(r)  = bOn(r)  + nnz(cv & onS);
        bOff(r) = bOff(r) + nnz(cv & offS);
    end
end
tol = @(pr) 4 * sqrt(2 * pr + 1);                    % overlapping windows cluster: var <= 2x
lab = {'2a. windows 1-16  ', '2b. windows 1-256 '};
for r = 1:nR
    nS = numel(runs{r}.scales);
    ok = all(abs(mOn(r, 1:nS) - pOn(r, 1:nS)) < tol(pOn(r, 1:nS))) && ...
         all(abs(mOff(r, 1:nS) - pOff(r, 1:nS)) < tol(pOff(r, 1:nS)));
    fprintf(['%s(%.1f dB, channels %d-%d): flagged windows ON the pulse %s, predicted %s; OFF %s, ' ...
             'predicted %s; blanked fraction on-pulse %.2e, off-pulse %.2e: %s\n'], lab{r}, ...
        info_rx.snrDB, useCh(1), useCh(end), mat2str(mOn(r, 1:nS)), mat2str(round(pOn(r, 1:nS))), mat2str(mOff(r, 1:nS)), ...
        mat2str(round(pOff(r, 1:nS))), bOn(r) / nOn, bOff(r) / nOff, passStr(ok));
    if ~ok, fails{end+1} = sprintf('pulsar run %d', r); end %#ok<AGROW>
end
fprintf('    (test 2: %.0f s)\n', toc(t2));

% ---------------------------------------------------------------------------------
% 3. RFI, scored with the RFI-only signal
% ---------------------------------------------------------------------------------
[mask3, d3] = detectRFI(chan1, 'KeepWindows', true, 'Verbose', false);
[mask30, d30] = detectRFI(chan1, 'FreqGuard', 0, 'Verbose', false);   % for 3d
blk3 = blockIndex(d30.blockEdges);
cov0 = false(nCh, Nc); cov5 = cov0; strongCov = cov0;
eA0k = 0;
radar = info_rx1.rfi(1).params;
sOf = @(t) round((t - chan1.t0) * fsC) + 1;          % time -> channel sample (1-based)
tR = radar.startTime + (0 : floor(Nc / fsC * radar.prf)).' / radar.prf;   % pulse starts
tR = tR(sOf(tR + radar.pulseWidth) <= Nc);
sI = sOf(info_rx1.rfi(2).evStart(:) / f_in);         % impulses
hw = round(20e-6 * fsC);
radarCh = 31:34; ceCh = [49 120 121]; skipCE = [48:50, 119:122]; skipI = [radarCh, skipCE];
eR = 0; eRk = 0; eA = 0; eAk = 0; eS = 0; eI = 0; eIk = 0; allPulses = true;
hitR = zeros(1, nCh);
inI = false(1, Nc);
for i = 1:numel(sI), inI(max(1, sI(i) - 15) : min(Nc, sI(i) + 15)) = true; end
for j = 1:nCh
    x1 = readIQ(chan1.chanFiles(j));
    Pr = abs(x1 - readIQ(chan0.chanFiles(j))).^2;
    cv = coveredOf(mask3, j, Nc);
    cov5(j, :) = cv;
    cov0(j, :) = coveredOf(mask30, j, Nc);
    % strong events of the unguarded mask: peak power / baseline >= FreqGuardMin
    p1 = abs(x1).^2;
    r0 = mask30(mask30(:, 1) == j, 2:3);
    for i = 1:size(r0, 1)
        if max(p1(r0(i, 1):r0(i, 2))) / d30.baseline(j, blk3(r0(i, 1))) >= d3.freqGuardMin
            strongCov(j, r0(i, 1):r0(i, 2)) = true;
        end
    end
    inR = false(1, Nc);
    for k = 1:numel(tR)
        core = max(1, sOf(tR(k))) : sOf(tR(k) + radar.pulseWidth);
        inR(max(1, core(1) - hw) : min(Nc, core(end) + hw)) = true;
        hitR(j) = hitR(j) + any(cv(core)) / numel(tR);
        if any(j == [32 33]), allPulses = allPulses && all(cv(core)); end
    end
    rR = inR & ~inI; rI = inI & ~inR;               % radar / impulse cells, kept apart
    if ~any(j == skipCE)                             % (carrier, GNSS: always on, not scored)
        % the steady RFI level of this channel (GNSS sinc^2 sidelobes reach far, 1/f^2)
        % is subtracted, so the all-channel numbers measure the radar alone
        Pq = Pr - mean(Pr(~inR & ~inI));
        eA = eA + sum(Pq(rR)); eAk = eAk + sum(Pq(rR & cv)); eA0k = eA0k + sum(Pq(rR & cov0(j, :)));
        if ~any(j == [32 33]), eS = eS + sum(Pq(rR)); end
    end
    if any(j == radarCh)
        eR = eR + sum(Pr(rR)); eRk = eRk + sum(Pr(rR & cv));
    elseif ~any(j == skipI)
        eI = eI + sum(Pr(rI)); eIk = eIk + sum(Pr(rI & cv));
    end
end
fR = eRk / eR; fA = eAk / eA; fA0 = eA0k / eA; fI = eIk / eI;
pass = allPulses && fR > 1 - 1e-4;
fprintf(['3a. radar: %d pulses, all blanked in channels 32-33 %d; energy removed in channels ' ...
         '31-34 %.6f (1 - %.1e): %s\n'], numel(tR), allPulses, fR, 1 - fR, passStr(pass));
fprintf(['    all channels: %.2e of the radar energy outside channels 32-33 (spectral sidelobes ' ...
         'of the pulse); left after excision %+.1e (without the frequency guard %+.1e; precision ' ...
         '~1e-6 from the level subtraction); channels flagged at >= half the pulses %s\n'], ...
    eS / eA, 1 - fA, 1 - fA0, mat2str(find(hitR >= 0.5)));
if ~pass, fails{end+1} = 'radar'; end
pass = fI > 0.95;
fprintf('3b. impulses: %d bursts, energy removed %.4f (expected ~0.98): %s\n', numel(sI), fI, passStr(pass));
if ~pass, fails{end+1} = 'impulses'; end
flagsCE = 0;
for j = ceCh
    for s = 1:numel(d3.scales)
        st = d3.flaggedWindows{j, s};
        flagsCE = flagsCE + nnz(min(abs(st - sI), [], 1) > 30);
    end
end
expCE = numel(ceCh) * sum(d3.nWindows) * d3.pFA;
pass  = flagsCE <= expCE + tol(expCE);
lvl = mean(d3.baseline(ceCh, :) ./ dD.baseline(ceCh, :), 2).';
fprintf(['3c. constant-envelope channels %s: %d flagged windows away from impulses (noise: %.1f); ' ...
         'baseline / noise-only %s; blanked %s %%, other channels %.3f %%: %s\n'], mat2str(ceCh), ...
    flagsCE, expCE, mat2str(round(lvl, 2)), mat2str(round(100 * d3.flaggedFraction(ceCh).', 3)), ...
    100 * median(d3.flaggedFraction(setdiff(1:nCh, [radarCh, ceCh]))), passStr(pass));
if ~pass, fails{end+1} = 'constant envelope'; end

% 3d. frequency guard: the mask = own detections + the strong ones (peak >= FreqGuardMin x
% baseline) in +-FreqGuard channels, same samples; none on noise or from the pulsar
G = d3.freqGuard;
covE = cov0;
for j = 1:nCh
    nb = setdiff(max(1, j - G) : min(nCh, j + G), j);
    covE(j, :) = covE(j, :) | any(strongCov(nb, :), 1);
end
okG = isequal(covE, cov5);
nsN = sum(dD.nStrong); nsP = sum(runs{1}.nStrong);
pass = okG && nsN == 0 && nsP == 0 && sum(d3.nStrong) > 0;
fprintf(['3d. frequency guard +-%d channels: mask = own + strong events copied (independent ' ...
         'coverage check) %d; strong events: RFI data %d, noise-only %d, -5 dB pulsar %d; blanked ' ...
         '%.4f %% (without the guard %.4f %%): %s\n'], G, okG, sum(d3.nStrong), nsN, nsP, ...
    100 * mean(d3.flaggedFraction), 100 * mean(d30.flaggedFraction), passStr(pass));
if ~pass, fails{end+1} = 'frequency guard'; end

% ---------------------------------------------------------------------------------
% 4. Continuous narrowband RFI: per-channel thresholds (21 Oct 2026)
% ---------------------------------------------------------------------------------
c4q = c4([1 3]);                                     % carrier channels without the radar
[m4o, d4o] = detectRFI(chan4, 'ChannelThresholds', false, 'Verbose', false);
[m4n, d4n] = detectRFI(chan4, 'Verbose', false);
[~, d4p] = detectRFI(chan4, 'PFA', 1e-3, 'Verbose', false);
[~, d4q] = detectRFI(chan4, 'PFA', 1e-3, 'ChannelThresholds', false, 'Verbose', false);
% (a) false flags on the radar-free carriers: own thresholds vs white ones
[okP, rP, sP] = countCheck(d4p.nFlaggedWindows(c4q, :), numel(c4q) * d4p.nWindows * d4p.pFA);
rQ  = sum(d4q.nFlaggedWindows(c4q, :), 1) ./ (numel(c4q) * d4q.nWindows * d4q.pFA);
rQt = sum(d4q.nFlaggedWindows(c4q, :), 'all') / (numel(c4q) * sum(d4q.nWindows) * d4q.pFA);
okA = okP && all(d4p.ownThreshold(c4q)) && rQt > 3;
fprintf(['4a. narrow carriers (channels %s), PFA %.0e: flagged / expected windows per scale %s (+- %s) ' ...
         'with own thresholds, %s with white ones (all lengths %.1f); default PFA, flagged windows per carrier channel ' ...
         '%s (white: %s); thresholds n = %d: own %s, white %.2f: %s\n'], mat2str(c4q), d4p.pFA, ...
    mat2str(round(rP, 2)), mat2str(round(sP, 2)), mat2str(round(rQ, 1)), rQt, ...
    mat2str(sum(d4n.nFlaggedWindows(c4, :), 2).'), mat2str(sum(d4o.nFlaggedWindows(c4, :), 2).'), ...
    d4n.scales(end), mat2str(round(d4n.thresholdUsed(c4, end).', 2)), d4o.threshold(end), passStr(okA));
if ~okA, fails{end+1} = 'narrow carriers: false flags'; end
% (b) the radar on the channel-70 carrier is still blanked, pulse by pulse
radar4 = info_rx4.rfi(end).params;
sOf4 = @(t) round((t - chan4.t0) * fsC) + 1;
tR4 = radar4.startTime + (0 : floor(Nc / fsC * radar4.prf)).' / radar4.prf;
tR4 = tR4(sOf4(tR4 + radar4.pulseWidth) <= Nc);
cv70 = coveredOf(m4n, c4(2), Nc);
hit4 = arrayfun(@(t) all(cv70(max(1, sOf4(t)) : sOf4(t + radar4.pulseWidth))), tR4);
okB = all(hit4) && d4n.ownThreshold(c4(2));
fprintf('4b. radar (+20 dB) on the channel-%d carrier: %d of %d pulses blanked with own thresholds: %s\n', ...
    c4(2), nnz(hit4), numel(hit4), passStr(okB));
if ~okB, fails{end+1} = 'narrow carriers: radar'; end
% (c) channels away from the RFI: as before, bit for bit; noise only: no own thresholds
white = setdiff(1:nCh, 55:85);
rowsW = @(m) m(ismember(m(:, 1), white), :);
okC = ~any(d4n.ownThreshold(white)) && isequal(rowsW(m4n), rowsW(m4o)) && ~any(dD.ownThreshold);
fprintf(['4c. channels outside 55-85: own thresholds in %d, mask = ''ChannelThresholds'' false bit for ' ...
         'bit %d; noise-only data (1c): own thresholds in %d channels: %s\n'], nnz(d4n.ownThreshold(white)), ...
    isequal(rowsW(m4n), rowsW(m4o)), nnz(dD.ownThreshold), passStr(okC));
if ~okC, fails{end+1} = 'narrow carriers: white channels'; end

% ---------------------------------------------------------------------------------
% 5. Gain changes (review R-E1)
% ---------------------------------------------------------------------------------
t5  = tic;
eB  = d1.blockEdges; nB = numel(eB) - 1;           % the clean pass-1 run of 1a: its baseline
blk = blockIndex(eB);
inSupS = false(1, Nc); inSupS(sup(1):sup(2)) = true;
ns  = round((eB(2) + eB(3)) / 2);                   % the step: middle of block 2
Gs  = {1 + 0.10 * (0:Nc-1) / Nc, 1 + 0.05 * ((1:Nc) > ns)};
gNames = {'ramp +10 % over the file', 'step +5 % in block 2'};
for c = 1:2
    G = Gs{c};
    chG = modChannels(chan0, fullfile(tmp, sprintf('g%d', c)), @(x, j) x .* sqrt(G));
    [~, dG] = detectRFI(chG, 'PFA', 1e-3, 'Passes', 1, 'KeepWindows', true, 'Verbose', false);
    gEff = zeros(1, nB);
    for b = 1:nB
        gB = G(eB(b) + 1 : eB(b + 1)); gB = gB(inSupS(eB(b) + 1 : eB(b + 1)));
        gEff(b) = mixMedian(gB, ones(size(gB)));    % the median of the gained exponential power
    end
    rB = dG.baseline ./ d1.baseline;                % the same noise: little scatter
    okA = all(abs(mean(rB, 1) - gEff) < max(4 * std(rB, 0, 1) / sqrt(nCh), 1e-6));   % constant-gain blocks: exact
    % (b) the flagged windows against the exact survival function at G_eff / G(t)
    nSc = numel(dG.scales);
    expd = zeros(1, nSc); eFirst = expd; eLast = expd; cFirst = expd; cLast = expd;
    for s = 1:nSc
        n = dG.scales(s);
        i0 = 1:dG.stride(s):Nc - n + 1;
        Gw = G(i0 + floor((n - 1) / 2));
        x  = dG.threshold(s) * n * gEff(blk(i0)) ./ Gw;
        F  = survTable(dG.lambda{s}, dG.threshold(s) * n);
        S  = exp(F(x));
        expd(s) = nCh * sum(S);
        pos = (i0 - eB(blk(i0)) - 1) ./ (eB(blk(i0) + 1) - eB(blk(i0)));   % position in the block, 0..1
        f1 = pos < 0.1; f2 = pos >= 0.9;
        eFirst(s) = nCh * sum(S(f1)); eLast(s) = nCh * sum(S(f2));
        for j = 1:nCh
            w = dG.flaggedWindows{j, s};
            pw = (w - eB(blk(w)) - 1) ./ (eB(blk(w) + 1) - eB(blk(w)));
            cFirst(s) = cFirst(s) + nnz(pw < 0.1); cLast(s) = cLast(s) + nnz(pw >= 0.9);
        end
    end
    [okB, r, sr] = countCheck(dG.nFlaggedWindows, expd);
    k = sr .* expd ./ sqrt(expd);                   % overdispersion of the counts (clustered windows)
    okC = all(abs(cFirst - eFirst) < 4 * k .* sqrt(eFirst)) && all(abs(cLast - eLast) < 4 * k .* sqrt(eLast));
    pass = okA && okB && (c == 2 || okC);
    fprintf(['5. %s: baseline / clean per block %s (effective gain %s, %d); flagged / predicted %s ' ...
             '(+- %s), predicted / nominal %s; '], gNames{c}, mat2str(round(mean(rB, 1), 5)), ...
        mat2str(round(gEff, 5)), okA, mat2str(round(r, 3)), mat2str(round(sr, 3)), ...
        mat2str(round(expd ./ (nCh * dG.nWindows * 1e-3), 4)));
    if c == 1
        fprintf('last / first tenth of the blocks measured %s, predicted %s (%d): %s\n', ...
            mat2str(round(cLast ./ cFirst, 2)), mat2str(round(eLast ./ eFirst, 2)), okC, passStr(pass));
    else
        fprintf('%s\n', passStr(pass));
    end
    if ~pass, fails{end+1} = ['gain change, ' gNames{c}]; end %#ok<AGROW>
    cellfun(@delete, cellstr(chG.chanFiles));
end
fprintf('   (test 5: %.0f s)\n', toc(t5));

% ---------------------------------------------------------------------------------
% 6. One pass against two under dense bursts (review R-E1)
% ---------------------------------------------------------------------------------
t6  = tic;
rs6 = RandStream('mt19937ar', 'Seed', 23);
slot = 32;                                          % one burst per slot of 32 samples at most
nSlot = floor((sup(2) - sup(1) + 1) / slot) - 2;
nStr = round(0.03 * Nc / 4); nWk = round(0.01 * Nc / 16);
mJ = mean(d1.baseline, 2);                          % each channel's noise power
sStart = cell(nCh, 1); wStart = cell(nCh, 1);
for j = 1:nCh
    sl = sup(1) + slot * randperm(rs6, nSlot, nStr + nWk);   % slot starts
    sStart{j} = sl(1:nStr) + randi(rs6, slot - 4, 1, nStr) - 1;
    wStart{j} = sl(nStr+1:end) + randi(rs6, slot - 16, 1, nWk) - 1;
end
cgn = @(m, n) sqrt(m / 2) * complex(randn(rs6, 1, n), randn(rs6, 1, n));   % complex Gaussian, power m
chB = modChannels(chan0, fullfile(tmp, 'b'), @(x, j) addBursts(x, sStart{j}, 4, cgn(100 * mJ(j), 4 * nStr), ...
    wStart{j}, 16, cgn(2 * mJ(j), 16 * nWk)));
isS = @(j) spanMask(sStart{j}, 4, Nc); isW = @(j) spanMask(wStart{j}, 16, Nc);
[m1, q1] = detectRFI(chB, 'Passes', 1, 'ChannelThresholds', false, 'KeepWindows', true, 'Verbose', false);
[m2, q2] = detectRFI(chB, 'Passes', 2, 'ChannelThresholds', false, 'KeepWindows', true, 'Verbose', false);
[m3, q3] = detectRFI(chB, 'KeepWindows', true, 'Verbose', false);            % defaults
g6 = q1.guard;
% (a), (b) the baselines against the median of the mixture (all samples / what pass 1 left)
pred1 = zeros(nCh, nB); pred2 = pred1;
for j = 1:nCh
    lev = ones(1, Nc); lev(isW(j)) = 3; lev(isS(j)) = 101;   % power level / noise of each sample
    cov1 = windowCover(q1.flaggedWindows(j, :), q1.scales, g6, Nc);
    for b = 1:nB
        idx = eB(b) + 1 : eB(b + 1);
        u = inSupS(idx);
        pred1(j, b) = mixMedian(ones(1, nnz(u)), lev(idx(u)));
        u2 = u & ~cov1(idx);
        pred2(j, b) = mixMedian(ones(1, nnz(u2)), lev(idx(u2)));
    end
end
rb1 = q1.baseline ./ d1.baseline; rb2 = q2.baseline ./ d1.baseline; rb3 = q3.baseline ./ d1.baseline;
dev = @(rb, pr) abs(mean(rb(:) - pr(:))) < 4 * std(rb(:) - pr(:)) / sqrt(numel(rb));
okA = dev(rb1, pred1); okB = dev(rb2, pred2);
% (c), (d) strong bursts blanked, weak bursts wholly blanked, by the channel's own detections
% (the mask also holds the frequency guard: these bursts are independent per channel, so each
% channel's strong bursts are copied into 14 neighbours and the mask covers most of the data)
[allS, fracW, fOwn, fMask] = deal(zeros(1, 3));
qs = {q1, q2, q3}; ms = {m1, m2, m3};
for r6 = 1:3
    nAll = 0; nW = 0;
    for j = 1:nCh
        cv = windowCover(qs{r6}.flaggedWindows(j, :), qs{r6}.scales, g6, Nc);
        cs = cumsum([0, cv]);
        nAll = nAll + nnz(cs(sStart{j} + 4) - cs(sStart{j}) == 4);
        nW = nW + nnz(cs(wStart{j} + 16) - cs(wStart{j}) == 16);
        fOwn(r6) = fOwn(r6) + mean(cv) / nCh;
        fMask(r6) = fMask(r6) + mean(coveredOf(ms{r6}, j, Nc)) / nCh;
    end
    allS(r6) = nAll / (nCh * nStr); fracW(r6) = nW / (nCh * nWk);
end
% a Gaussian burst fades: its 4 samples' power total is Gamma(4) x 101 noise powers, below the
% 4-sample threshold eta_4 * 4 * baseline with probability ~1e-4 (the shorter windows rarely help)
pMiss = gammainc(4 * q1.threshold(q1.scales == 4) * mean(rb1(:)) / 101, 4);
okC = all(allS(1:2) >= 0.999); okD = fracW(2) > fracW(1);
% (e) false windows in clean stretches against each pass's own baseline
qq = {q1, q2}; [okE, rE, sE] = deal(false(1, 2), cell(1, 2), cell(1, 2));
for r6 = 1:2
    q = qq{r6};
    cnt = zeros(nCh, numel(q.scales)); expd = zeros(1, numel(q.scales));
    F = arrayfun(@(s) survTable(q.lambda{s}, q.threshold(s) * q.scales(s)), 1:numel(q.scales), 'UniformOutput', false);
    for j = 1:nCh
        near = dilate(isS(j) | isW(j), g6);
        cn = cumsum([0, near]);
        for s = 1:numel(q.scales)
            n = q.scales(s);
            i0 = 1:q.stride(s):Nc - n + 1;
            cl = cn(i0 + n) - cn(i0) == 0;          % windows with no burst within the guard
            x = q.threshold(s) * n * q.baseline(j, blk(i0(cl))) / mJ(j);
            expd(s) = expd(s) + sum(exp(F{s}(x)));
            w = q.flaggedWindows{j, s};
            cnt(j, s) = nnz(cn(w + n) - cn(w) == 0);
        end
    end
    [okE(r6), rE{r6}, sE{r6}] = countCheck(cnt, expd);
end
% (f) the defaults
okF = allS(3) >= 0.999 && fracW(3) >= fracW(1);
pass = okA && okB && okC && okD && all(okE) && okF;
fprintf(['6. dense bursts (%d strong, %d weak per channel): baseline / clean, pass 1 %.4f (mixture median %.4f, %d), ' ...
         'pass 2 %.4f (median of the unflagged %.4f, %d); strong bursts wholly blanked %s (missed %s of %d; Gamma ' ...
         'estimate %.0f); weak bursts wholly blanked ' ...
         '%s (pass 1, pass 2, defaults); clean-stretch false windows / predicted, pass 1 %s (+- %s), pass 2 %s ' ...
         '(+- %s) (%d %d); defaults: own thresholds in %d channels, baseline / clean %.4f; blanked fraction, own ' ...
         'detections %s, with the frequency guard %s: %s\n'], nStr, nWk, ...
    mean(rb1(:)), mean(pred1(:)), okA, mean(rb2(:)), mean(pred2(:)), okB, mat2str(round(allS, 5)), ...
    mat2str(round((1 - allS) * nCh * nStr)), nCh * nStr, pMiss * nCh * nStr, mat2str(round(fracW, 3)), ...
    mat2str(round(rE{1}, 2)), mat2str(round(sE{1}, 2)), mat2str(round(rE{2}, 2)), mat2str(round(sE{2}, 2)), okE, ...
    nnz(q3.ownThreshold), mean(rb3(:)), mat2str(round(fOwn, 3)), mat2str(round(fMask, 3)), passStr(pass));
if ~pass, fails{end+1} = 'dense bursts, one pass against two'; end
cellfun(@delete, cellstr(chB.chanFiles));
fprintf('   (test 6: %.0f s)\n', toc(t6));

% ---------------------------------------------------------------------------------
if isempty(fails)
    fprintf('testDetectRFI: ALL PASSED\n');
else
    error('testDetectRFI:failed', 'FAILED: %s', strjoin(fails, ', '));
end
end


% =================================================================================
function [info_chan, info_rx] = makeChannels(zeroFile, base, rfi, rx)
% Receiver chain on a zero sky signal: noise (seed 11, std 1) + rfi -> IQ -> channels.
info_rx = addNoiseAndRFI(zeroFile, [base '_rx.dat'], rx.f_in, 'Band', rx.band, 'NoiseStd', 1, ...
    'RFI', rfi, 'Seed', 11, 'Verbose', false);
info_IQ = applyIQmodulation(info_rx.file, [base '_iq.dat'], rx.f_in, rx.fs, rx.fLO, ...
    'FilterOrder', rx.order, 'Band', rx.band, 'Verbose', false);
delete(info_rx.file);
info_chan = channelizeIQ(info_IQ.file, [base '_chan'], info_IQ.actualFsOut, info_IQ.fLO, ...
    rx.band(1), rx.band(2), 'T0', info_IQ.t0, 'Verbose', false);
delete(info_IQ.file);
end

function [ok, r, sr] = countCheck(counts, expected)
% Flagged windows summed over channels vs expected, per scale; sigma from the
% scatter between channels (at least Poisson).
T  = sum(counts, 1);
sg = sqrt(max(size(counts, 1) * var(counts, 0, 1), expected));
ok = all(abs(T - expected) < 4 * sg);
r  = T ./ expected; sr = sg ./ expected;
end

function P = sumExpSurvival(lam, x)
% P(sum_k lam_k E_k > x), E_k unit exponentials (as in detectRFI): phase-type, expm.
mu = 1 ./ lam(:);
Q  = diag(-mu);
if numel(mu) > 1, Q = Q + diag(mu(1:end-1), 1); end
E  = expm(Q * x);
P  = max(sum(E(1, :)), realmin);
end

function F = survTable(lam, x0)
% log P(sum_k lam_k E_k > x) as a pchip interpolant over 0.8 .. 1.25 x0 (the exact survival
% function of sumExpSurvival, tabulated once per scale).
xs = x0 * linspace(0.8, 1.25, 300);
F = griddedInterpolant(xs, arrayfun(@(v) log(sumExpSurvival(lam, v)), xs), 'pchip', 'linear');
end

function g = mixMedian(G, lev)
% The median of independent exponential powers with means G .* lev (one per sample), / ln 2:
% what the median / ln 2 baseline converges to, in units of the noise power (1 for pure noise).
s = G(:) .* lev(:);
x = fzero(@(x) mean(exp(-x ./ s)) - 0.5, log(2) * median(s));
g = x / log(2);
end

function info = modChannels(info0, base, fun)
% A copy of the channel set info0 with channel j's samples x (complex row) replaced by fun(x, j).
info = info0;
info.chanFiles = strings(info0.nChan, 1);
for j = 1:info0.nChan
    info.chanFiles(j) = sprintf('%s_ch%03d.dat', base, j);
    z = fun(readIQ(info0.chanFiles(j)), j);
    fid = fopen(info.chanFiles(j), 'w', 'ieee-le');
    fwrite(fid, [real(z); imag(z)], 'single');
    fclose(fid);
end
end

function z = addBursts(x, sSt, sLen, sVals, wSt, wLen, wVals)
% x plus bursts: sVals in the sLen samples from each start sSt, wVals likewise (column order).
z = x;
iS = sSt(:).' + (0:sLen-1).';
z(iS(:)) = z(iS(:)) + sVals;
iW = wSt(:).' + (0:wLen-1).';
z(iW(:)) = z(iW(:)) + wVals;
end

function m = spanMask(st, len, N)
% Samples st .. st+len-1 of 1..N (logical row), for every start st.
d = accumarray([st(:); st(:) + len], [ones(numel(st), 1); -ones(numel(st), 1)], [N + 1, 1]);
m = cumsum(d(1:N)).' > 0;
end

function m = dilate(m, g)
% A logical row widened by g samples on each side.
m = conv(double(m), ones(1, 2*g + 1), 'same') > 0;
end

function cv = windowCover(fw, scales, g, N)
% Samples covered by the flagged windows fw{s} (starts, scale s) widened by the guard g.
a = []; e = [];
for s = 1:numel(scales)
    a = [a, fw{s} - g]; e = [e, fw{s} + scales(s) - 1 + g]; %#ok<AGROW>
end
a = max(a, 1); e = min(e, N);
d = accumarray([a, e + 1].', [ones(1, numel(a)), -ones(1, numel(e))].', [N + 1, 1]);
cv = cumsum(d(1:N)).' > 0;
end

function blk = blockIndex(edges)
blk = zeros(1, edges(end));
for b = 1:numel(edges) - 1, blk(edges(b) + 1 : edges(b + 1)) = b; end
end

function cv = coveredOf(mask, j, N)
% Samples of channel j inside the mask rows.
r = mask(mask(:, 1) == j, 2:3);
d = accumarray([r(:, 1); r(:, 2) + 1], [ones(size(r, 1), 1); -ones(size(r, 1), 1)], [N + 1, 1]);
cv = cumsum(d(1:N)).' > 0;
end

function writeZeros(file, N)
fid = fopen(file, 'w', 'ieee-le');
c = onCleanup(@() fclose(fid));
blk = 4e6;
for n0 = 0:blk:N - 1
    fwrite(fid, zeros(1, min(blk, N - n0), 'single'), 'single');
end
end

function p = readPower(file)
x = readRaw(file);
p = x(1, :).^2 + x(2, :).^2;
end

function z = readIQ(file)
x = readRaw(file);
z = complex(x(1, :), x(2, :));
end

function x = readRaw(file)
fid = fopen(file, 'r', 'ieee-le');
c = onCleanup(@() fclose(fid));
x = fread(fid, [2 Inf], 'single=>double');
end

function s = passStr(p)
if p, s = 'PASS'; else, s = 'FAIL'; end
end
