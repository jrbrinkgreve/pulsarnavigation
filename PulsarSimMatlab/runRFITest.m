%{
RUNRFITEST - effect of each RFI type on the channel path, without and with excision (B4)

For each case (noise only; one source of the pipelineParams scenario at a
time; all five) the receiver stage runs with that RFI (main's noise seed, so
the receiver noise is the same in every case), then the channel path of
main.m twice: without excision, and with it (detectRFI -> periodicRFI if
periodicMask (B6; in runRFITest since 9 Oct 2026, before that detectRFI
only) -> blankChannels -> dedisperseChannels -> detectChannels ->
blankingWeights -> foldProfile with DataWeights). Both use 'optimal'
weighting. A table at the end: blanked
fraction, good TOAs, median SNR, median TOA error, rms / predicted error
(vs ground truth), median reduced chi^2, total-fold offset.

Reference for the full band without excision: notes §7, RFI test of 5 Oct
(carrier then at 1350 MHz; in the full band its frequency does not matter).

Rotating radar (B5c-2): the scenario radar with a rotating antenna (10 s
scan, 1.4 deg beam -> BeamWidth 39 ms): one beam passage centred in the
file, and sidelobes only (beam 5 s away) at -25 / -30 / -35 dB. For cases
with a pulsed source the table also gives, from the source's parameters
(ground truth, validation only): the fraction of its pulses blanked in its
own channels, and the fraction of its energy in the pulses that were missed.

LTE (B5d-2): band-limited Gaussian noise (rfiSource 'noise'), 20 MHz at
1472 MHz (inside LTE band 32), INR -10 dB and 0 dB of the whole-band noise
(~+3 / +13 dB in its ~6 channels): cases 12-13.

Realistic scenario (B5e, 9 Oct 2026): case 14 = rfiRealistic of
pipelineParams (16 sources: GNSS, two rotating radars, LTE band 32,
Inmarsat, impulses, a local carrier), ~2.5 min; case 15 = the same with
Inmarsat as ten narrow (200 kHz) carriers inside channels 105-114 instead
of three wide blocks (the worst case for the error bars). For every case a per-channel
noise check (evaluation only): r = relative variance of the off-pulse fold
bins (> 3 FWHM from the pulse) per channel over the median channel's
(model: all channels alike, so r = 1); a channel only partly covered by
noise-like RFI has r > 1. Valid only when little is blanked: r compares
with the median channel, not with the blanking-aware model (heavy
blanking changes a channel's noise legitimately). The table gives the number of channels with
r > 1.1 and how much the TOA error bars are too small because of it:
sqrt(sum F r / sum F) with F = 1/a_c^2 (the fit's channel weights;
independent channels). After the table, per case with excision, the
channels with r > 1.1, or with flagged windows > 10x the expected false
flags and > 3x the median channel's, with the sources covering them
(ground truth, for the labels), and the periodic emitters found.

Needs: the sky files (main.m with runStage.sky once). Files go to data/rfi
and are overwritten per case; main's files and data/chan are not touched.
~45 s per case; runCases picks a subset.
%}

% Paths and parameters
scriptDir = fileparts(mfilename('fullpath'));
addpath(fullfile(scriptDir, 'functions'));
pipelineParams;
if ~strcmp(frontEnd, 'channels')
    error('runRFITest:frontEnd', 'runRFITest is for the channel path (frontEnd = ''channels'').');
end

% B5e: Inmarsat as narrow carriers (~200 kHz, BGAN-like), one inside each of channels
% 105-114 at -21 dB: each about as strong as its channel's noise, the worst case for
% the error bars (weight still ~1/4, noise ~4-5x the model's)
inmF = fLow + ((105:114) - 0.5) * chanWidth + [-1.1 0.4 -0.6 1.2 -0.2 0.9 -1.3 0.1 0.7 -0.9] * 1e6;
inmNarrow = arrayfun(@(f) rfiSource('noise', 'Freq', f, 'Bandwidth', 200e3, 'INRdB', -21, ...
    'Label', sprintf('Inmarsat carrier %.2f MHz', f / 1e6)), inmF, 'UniformOutput', false);
inmNarrow = [inmNarrow{:}];

% cases: indices into rfiScenario ([] = noise only), rfiSource structs, or {structs, label}
rad = rfiScenario(3);                % the scenario radar, with a rotating antenna
rotRadar = @(beamTime, sl, label) rfiSource('pulsed', 'Freq', rad.freq, 'PulseWidth', rad.pulseWidth, ...
    'PRF', rad.prf, 'ChirpBW', rad.chirpBW, 'StartTime', rad.startTime, 'INRdB', rad.INRdB, ...
    'RiseTime', rad.riseTime, 'ScanPeriod', 10, 'BeamTime', beamTime, 'BeamWidth', 39e-3, ...
    'SidelobeDB', sl, 'Label', label);
cases = {[], 1, 2, 3, 4, 5, 1:5, ...
    rotRadar(L/2, -30, 'rotating radar, beam passage'), ...
    rotRadar(5, -25, 'radar sidelobes -25 dB'), rotRadar(5, -30, 'radar sidelobes -30 dB'), ...
    rotRadar(5, -35, 'radar sidelobes -35 dB'), ...
    rfiSource('noise', 'Freq', 1472e6, 'Bandwidth', 20e6, 'INRdB', -10, 'Label', 'LTE 20 MHz, -10 dB'), ...
    rfiSource('noise', 'Freq', 1472e6, 'Bandwidth', 20e6, 'INRdB', 0, 'Label', 'LTE 20 MHz, 0 dB'), ...
    {rfiRealistic, 'realistic scenario (B5e)'}, ...
    {[rfiRealistic(~startsWith({rfiRealistic.label}, 'Inmarsat')), inmNarrow], 'realistic, Inmarsat narrow'}};
runCases = 1:numel(cases);           % e.g. 8:11 rotating radar, 12:13 LTE, 14:15 realistic (B5e)
rThr     = 1.1;                      % per-channel check: list channels with r above this

rfiDir  = fullfile(dataDir, "rfi");
rRx     = fullfile(rfiDir, "rfi_rx.dat");
rIQ     = fullfile(rfiDir, "rfi_rx_IQ.dat");
rChan   = fullfile(rfiDir, "rfi_chan");
rBlank  = fullfile(rfiDir, "rfi_chan_blanked");
rDedisp = fullfile(rfiDir, "rfi_dedisp_chan");
rPower  = fullfile(rfiDir, "rfi_dedisp_chan_power.dat");
rWeight = fullfile(rfiDir, "rfi_dedisp_chan_weights.dat");

info_gen  = loadInfo(fileRaw);
info_disp = loadInfo(fileDispersed);
template  = gaussianTemplate(nBin, ephem.profileFWHM);

nC = numel(runCases);
res = struct('label', {}, 'excision', {}, 'blanked', {}, 'nGood', {}, 'nValid', {}, ...
    'snr', {}, 'err', {}, 'ratio', {}, 'redChi2', {}, 'offset', {}, 'offsetErr', {}, ...
    'pulses', {}, 'missedE', {}, 'chan', {});
fprintf('runRFITest: %d cases, SNR %.1f dB, L = %.3g s, noise seed %d\n', nC, snrDB, L, noiseSeed);
tAll = tic;
for c = runCases
    tCase = tic;
    idx = cases{c};
    if isempty(idx)
        rfiCase = struct([]); label = 'noise only';
    elseif iscell(idx)                               % {rfiSource structs, label}
        rfiCase = idx{1}; label = idx{2};
    elseif isstruct(idx)
        rfiCase = idx; label = strjoin({rfiCase.label}, ' + ');
    else
        rfiCase = rfiScenario(idx);
        label = strjoin({rfiCase.label}, ' + ');
        if numel(idx) == numel(rfiScenario), label = 'all five'; end
    end

    % Receiver (synthetic): noise + this RFI at RF, IQ conversion, filterbank
    info_rx = addNoiseAndRFI(info_disp.file, rRx, info_disp.actualFsOut, ...
        'Band', [fLow fHigh], 'SNRdB', snrDB, 'SignalInfo', info_gen, ...
        'RFI', rfiCase, 'Seed', noiseSeed, 'Verbose', false);
    info_IQ = applyIQmodulation(info_rx.file, rIQ, info_rx.actualFsOut, fs, fLO, ...
        'FilterOrder', filterOrder, 'Band', [fLow fHigh], 'Verbose', false);
    delete(info_rx.file);
    info_chan = channelizeIQ(info_IQ.file, rChan, info_IQ.actualFsOut, info_IQ.fLO, ...
        fLow, fHigh, 'ChanWidth', chanWidth, 'T0', info_IQ.t0, 'Verbose', false);
    delete(info_IQ.file);
    srcLabels = sourcesIn(rfiCase, info_chan);       % ground truth, labels of the listing

    % Processing (ephemeris + receiver settings only), as main.m's channel path
    for exc = [false true]
        prfs = [];                                   % periodic emitters found [Hz]
        if exc
            [rfiMask, info_rfi] = detectRFI(info_chan, excisionArgs{:}, 'Verbose', false);
            if periodicMask             % + every predicted pulse of periodic RFI (as main.m)
                [perRows, info_per] = periodicRFI(info_chan, info_rfi, periodicArgs{:}, 'Verbose', false);
                rfiMask = [rfiMask; perRows]; %#ok<AGROW>
                prfs = [info_per.emitters.prf];
            end
            info_chanD = blankChannels(info_chan, rfiMask, rBlank, 'Verbose', false);
        else
            info_chanD = info_chan;
        end
        info_dc  = dedisperseChannels(info_chanD, rDedisp, ephem.DM, 'RefFreq', refFreq, 'Verbose', false);
        info_det = detectChannels(info_dc, rPower, info_dc.fs / round(info_dc.fs / f_out), 'Verbose', false);
        if exc
            info_w = blankingWeights(info_chan, info_dc, info_det, info_chanD.mask, rWeight, 'Verbose', false);
            noiseArgs = {'DataWeights', info_w};
            blanked = mean(info_chanD.blankedFraction);
            [pulses, missedE] = pulseStats(rfiCase, rfiMask, info_chan);
        else
            pulses = NaN; missedE = NaN;
            noiseArgs = {'NoiseCoeffs', [info_det.noise.V, info_det.noise.X]};
            blanked = 0;
        end
        [info_fold, fold] = foldProfile(info_det, fullfile(rfiDir, "rfi_fold.mat"), ephem.f0, 'F1', ephem.F1, 'TRef', ephem.TRef, ...
            'NBin', nBin, 'SubintPeriods', subintPeriods, noiseArgs{:}, 'SaveFile', false, 'Verbose', false);
        toaArgs = {'Bnoise', info_det.noise.Bnoise, 'Weighting', weighting};
        toa = estimateTOA(fold, info_fold, template, toaArgs{:}, 'Verbose', false);
        detection = detectPulsar(fold, info_fold, template, toaArgs{:}, 'Verbose', false);
        good = toa.valid & detection.detectedUnknown & ~toa.flagChi2;
        v = toa.valid;
        if any(v)
            val = validateTOA(toa, info_gen, 'Plot', false, 'Label', sprintf('%s, excision %d', label, exc));
        else
            val = struct('ratio', NaN, 'totalOffset', NaN, 'totalOffsetErr', NaN);
        end
        nChan = numel(info_chan.chanFreqs);
        if exc
            flags = sum(info_rfi.nFlaggedWindows, 2).';
            expFlags = sum(info_rfi.nWindows) * info_rfi.pFA;
        else
            flags = nan(1, nChan); expFlags = NaN;
        end
        ch = channelNoise(fold, toa, ephem.profileFWHM);
        ch.freq = info_chan.chanFreqs(:).'; ch.flags = flags; ch.expFlags = expFlags;
        ch.src = srcLabels; ch.prf = prfs;
        res(end+1) = struct('label', label, 'excision', exc, 'blanked', blanked, ...
            'nGood', nnz(good), 'nValid', nnz(v), 'snr', median(toa.snr(v)), ...
            'err', median(toa.toaErr(v)) * 1e6, 'ratio', val.ratio, ...
            'redChi2', median(toa.redChi2(v)), 'offset', val.totalOffset * 1e6, ...
            'offsetErr', val.totalOffsetErr * 1e6, 'pulses', pulses, 'missedE', missedE, ...
            'chan', ch); %#ok<SAGROW>
    end
    fprintf('  case %d of %d (%s): %.0f s\n', c, numel(cases), label, toc(tCase));
end
fprintf('runRFITest: done in %.0f s\n\n', toc(tAll));

% Table
fprintf('%-34s %4s %9s %6s %7s %9s %8s %8s %16s %8s %9s %5s %9s\n', 'case', 'exc', 'blanked', 'good', ...
    'SNR', 'err [us]', 'rms/pred', 'redChi2', 'offset [us]', 'pulses', 'E missed', ...
    sprintf('r>%g', rThr), 'err bars');
for k = 1:numel(res)
    r = res(k);
    pc = '       -'; me = '        -';
    if ~isnan(r.pulses), pc = sprintf('%7.1f%%', 100*r.pulses); me = sprintf('%9.1e', r.missedE); end
    fprintf('%-34s %4d %8.4f%% %3d/%-2d %7.1f %9.3f %8.3f %8.2f %+8.2f +- %5.2f %s %s %5d %+8.2f%%\n', ...
        r.label, r.excision, 100*r.blanked, r.nGood, r.nValid, r.snr, r.err, r.ratio, r.redChi2, ...
        r.offset, r.offsetErr, pc, me, nnz(r.chan.r > rThr), 100 * (r.chan.errBarFactor - 1));
end
fprintf(['pulses: fraction of the pulsed source''s pulses blanked in its own channels; E missed: ' ...
         'fraction of its energy in the missed pulses (ground truth)\n' ...
         'r>%g: channels whose off-pulse fold noise is > %g x the median channel''s; err bars: ' ...
         'how much the TOA error bars are too small because of them (sqrt(sum F r / sum F) - 1)\n'], ...
         rThr, rThr);

% Per-channel listing, runs with excision
fprintf(['\nPer channel (with excision): channels with r > %g, or with flagged windows > 10x the ' ...
         'expected false flags and > 3x the median channel''s (broadband impulses flag every ' ...
         'channel); level = baseline a_c / the median channel''s (fit weight 1/level^2)\n'], rThr);
for k = find([res.excision])
    c = res(k).chan;
    prfText = 'none';
    if ~isempty(c.prf), prfText = strjoin(compose('%.3f Hz', c.prf), ', '); end
    fprintf(['%s: %d channel(s) with r > %g, mean r %.3f, error bars %.2f %% too small; ' ...
             'median flags %g; periodic emitters: %s\n'], res(k).label, nnz(c.r > rThr), rThr, ...
        c.meanR, 100 * (c.errBarFactor - 1), median(c.flags), prfText);
    for j = find(c.r > rThr | c.flags > max(10 * c.expFlags, 3 * median(c.flags)))
        fprintf('  ch %3d %9.3f MHz  r %6.2f  level %8.2f  flags %6d (exp. %.1f)  %s\n', j, ...
            c.freq(j) / 1e6, c.r(j), c.level(j), c.flags(j), c.expFlags, c.src{j});
    end
end
save(fullfile(rfiDir, "runRFITest_results.mat"), 'res', 'cases', 'runCases');


% ======================================================================
%  Local functions
% ======================================================================
function [frac, missedE] = pulseStats(rfiCase, mask, info_chan)
% For the first pulsed source of the case (ground truth): fraction of its pulses whose
% samples are blanked in all its own channels (those within ChanWidth/2 + ChirpBW/2 of
% its frequency), and the fraction of its energy (per-pulse antenna gain) in the others.
frac = NaN; missedE = NaN;
if isempty(rfiCase), return; end
k = find(strcmp({rfiCase.type}, 'pulsed'), 1);
if isempty(k), return; end
d = rfiCase(k);
fsC = info_chan.fs; Nc = info_chan.N;
ch = find(abs(info_chan.chanFreqs - d.freq) < info_chan.chanWidth/2 + d.chirpBW/2).';
sOf = @(t) round((t - info_chan.t0) * fsC) + 1;
tk = d.startTime + (0 : floor(Nc / fsC * d.prf)).' / d.prf;      % pulse starts
tk = tk(sOf(tk) >= 1 & sOf(tk + d.pulseWidth) <= Nc);
g = ones(size(tk));                                              % antenna gain per pulse
if isfinite(d.scanPeriod)
    T = d.scanPeriod;
    dt = mod(tk + d.pulseWidth/2 - d.beamTime + T/2, T) - T/2;
    g = max(exp(-4*log(2) * (dt / d.beamWidth).^2), 10^(d.sidelobeDB/10));
end
hit = true(size(tk));
for j = ch
    r = mask(mask(:, 1) == j, 2:3);
    cv = false(1, Nc);
    for i = 1:size(r, 1), cv(r(i, 1):r(i, 2)) = true; end
    for i = 1:numel(tk)
        hit(i) = hit(i) && all(cv(sOf(tk(i)) : sOf(tk(i) + d.pulseWidth)));
    end
end
frac = mean(hit);
missedE = sum(g(~hit)) / sum(g);
end


function ch = channelNoise(fold, toa, fwhm)
% Per-channel noise check (evaluation only, B5e). r: relative variance (var / mean^2) of
% the off-pulse fold bins (phase > 3 FWHM from the pulse at phase 0), averaged over the
% valid sub-ints, over the median channel's. The model is the same in every channel, so
% r = 1; a channel whose spectrum is not flat (partly covered by noise-like RFI) has
% fewer independent samples per bin -> r > 1. level: the fit's baseline a_c over the
% median channel's. errBarFactor = sqrt(sum F r / sum F), F = 1/a_c^2 (the fit's channel
% weights): true / reported TOA error when channel c's noise is r times the model's
% (independent channels, weights as used). meanR ~ red. chi^2 relative to clean data.
ph  = fold.phase(:);
off = min(ph, 1 - ph) > 3 * fwhm;
v   = find(toa.valid);
P   = fold.prof(off, v, :);                                  % [nOff x nValid x nChan]
rv  = var(P, 0, 1, 'omitnan') ./ mean(P, 1, 'omitnan').^2;  % [1 x nValid x nChan]
rv  = reshape(mean(rv, 2, 'omitnan'), 1, []);
r   = rv / median(rv, 'omitnan');
a   = median(toa.chanBaseline(v, :), 1, 'omitnan');
F   = 1 ./ a.^2;
use = isfinite(r) & isfinite(F) & F > 0;
ch  = struct('r', r, 'level', a / median(a, 'omitnan'), 'meanR', mean(r(use)), ...
    'errBarFactor', sqrt(sum(F(use) .* r(use)) / sum(F(use))));
end


function lab = sourcesIn(rfiCase, info_chan)
% Ground truth, for the per-channel listing only: labels of the sources whose band
% overlaps each channel (cw: the carrier; bpsk: main lobe +-ChipRate; pulsed: the chirp;
% noise: +-Bandwidth/2). Impulses are broadband (every channel) and not listed.
f   = info_chan.chanFreqs(:).';
lab = repmat({''}, 1, numel(f));
for d = rfiCase(:).'
    switch d.type
        case 'cw',     h = 0;
        case 'bpsk',   h = d.chipRate;
        case 'pulsed', h = d.chirpBW / 2;
        case 'noise',  h = d.bandwidth / 2;
        otherwise,     continue
    end
    hit = abs(f - d.freq) < info_chan.chanWidth / 2 + h;
    if ~any(hit), continue; end
    lab(hit) = strcat(lab(hit), {[d.label '; ']});
end
lab = regexprep(lab, '; $', '');
end
