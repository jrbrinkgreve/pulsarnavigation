function testDetectRFI()
%TESTDETECTRFI  Unit tests for detectRFI (B2): false flags, the bright pulsar, real RFI.
%{
Run from the PulsarSimMatlab folder: run('tests/testDetectRFI.m').
Makes its own noise-only and noise + RFI channel files (0.05 s through the
real receiver chain: addNoiseAndRFI on a zero sky signal -> applyIQmodulation
-> channelizeIQ; the same noise seed, so noise + RFI minus noise only is the
RFI alone, channel by channel; ~2.5 GB in tempdir, deleted on the way), and
uses the seed-43 channel files (data/chan/test_rx_IQ_chan_*: -5 dB pulsar,
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
         the same over all channels: the pulse switches on and off
         instantly, so its spectrum has sinc sidelobes (~1/f^2) over tens of
         MHz; the detector finds them where they stand out per channel;
     (b) impulses: fraction of their energy removed (other channels, +-15
         samples, away from the radar pulses) > 0.95. Expected ~0.98: per channel a burst carries ~80
         times a sample's noise energy, but ~exponentially distributed, so
         the weakest ~15 % stay below the thresholds;
     (c) constant-envelope RFI is not flagged: in the carrier and GNSS
         channels (49, 120, 121) no more flags away from the impulses than
         noise would give; their baseline vs noise only shown.
     Shown: channels flagged at >= half of the radar pulses.
Errors at the end if any check fails.
%}

root = fullfile(fileparts(mfilename('fullpath')), '..');
addpath(fullfile(root, 'functions'));
oldDir = cd(root); restoreDir = onCleanup(@() cd(oldDir));
pipelineParams;
tmp = fullfile(tempdir, 'testDetectRFI');
if ~isfolder(tmp), mkdir(tmp); end
cleanTmp = onCleanup(@() cleanFolder(tmp));
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
    error('testDetectRFI:data', 'Test 2 needs the seed-43 files without RFI (main.m, rfiOn = false).');
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
    Pr = abs(readIQ(chan1.chanFiles(j)) - readIQ(chan0.chanFiles(j))).^2;
    cv = coveredOf(mask3, j, Nc);
    inR = false(1, Nc);
    for k = 1:numel(tR)
        core = max(1, sOf(tR(k))) : sOf(tR(k) + radar.pulseWidth);
        inR(max(1, core(1) - hw) : min(Nc, core(end) + hw)) = true;
        hitR(j) = hitR(j) + any(cv(core)) / numel(tR);
        if any(j == [32 33]), allPulses = allPulses && all(cv(core)); end
    end
    rR = inR & ~inI; rI = inI & ~inR;               % radar / impulse cells, kept apart
    if ~any(j == skipCE)                             % (carrier, GNSS: always on, not scored)
        eA = eA + sum(Pr(rR)); eAk = eAk + sum(Pr(rR & cv));
        if ~any(j == [32 33]), eS = eS + sum(Pr(rR)); end
    end
    if any(j == radarCh)
        eR = eR + sum(Pr(rR)); eRk = eRk + sum(Pr(rR & cv));
    elseif ~any(j == skipI)
        eI = eI + sum(Pr(rI)); eIk = eIk + sum(Pr(rI & cv));
    end
end
fR = eRk / eR; fA = eAk / eA; fI = eIk / eI;
pass = allPulses && fR > 1 - 1e-4;
fprintf(['3a. radar: %d pulses, all blanked in channels 32-33 %d; energy removed in channels ' ...
         '31-34 %.6f (1 - %.1e): %s\n'], numel(tR), allPulses, fR, 1 - fR, passStr(pass));
fprintf(['    all channels: %.2e of the radar energy outside channels 32-33 (sidelobes of the ' ...
         'sharp-edged pulse), removed %.5f (1 - %.1e); channels flagged at >= half the pulses %s\n'], ...
    eS / eA, fA, 1 - fA, mat2str(find(hitR >= 0.5)));
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

function cleanFolder(d)
f = dir(fullfile(d, '*'));
for i = 1:numel(f)
    if ~f(i).isdir, delete(fullfile(d, f(i).name)); end
end
end
