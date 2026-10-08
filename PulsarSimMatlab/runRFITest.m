%{
RUNRFITEST - effect of each RFI type on the channel path, without and with excision (B4)

For each case (noise only; one source of the pipelineParams scenario at a
time; all five) the receiver stage runs with that RFI (main's noise seed, so
the receiver noise is the same in every case), then the channel path of
main.m twice: without excision, and with it (detectRFI -> blankChannels ->
dedisperseChannels -> detectChannels -> blankingWeights -> foldProfile with
DataWeights). Both use 'optimal' weighting. A table at the end: blanked
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

% cases: indices into rfiScenario ([] = noise only) or rfiSource structs
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
    rfiSource('noise', 'Freq', 1472e6, 'Bandwidth', 20e6, 'INRdB', 0, 'Label', 'LTE 20 MHz, 0 dB')};
runCases = 1:numel(cases);           % e.g. 8:11 rotating radar, 12:13 LTE

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
    'pulses', {}, 'missedE', {});
fprintf('runRFITest: %d cases, SNR %.1f dB, L = %.3g s, noise seed %d\n', nC, snrDB, L, noiseSeed);
tAll = tic;
for c = runCases
    tCase = tic;
    idx = cases{c};
    if isempty(idx)
        rfiCase = struct([]); label = 'noise only';
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

    % Processing (ephemeris + receiver settings only), as main.m's channel path
    for exc = [false true]
        if exc
            rfiMask = detectRFI(info_chan, excisionArgs{:}, 'Verbose', false);
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
        res(end+1) = struct('label', label, 'excision', exc, 'blanked', blanked, ...
            'nGood', nnz(good), 'nValid', nnz(v), 'snr', median(toa.snr(v)), ...
            'err', median(toa.toaErr(v)) * 1e6, 'ratio', val.ratio, ...
            'redChi2', median(toa.redChi2(v)), 'offset', val.totalOffset * 1e6, ...
            'offsetErr', val.totalOffsetErr * 1e6, 'pulses', pulses, 'missedE', missedE); %#ok<SAGROW>
    end
    fprintf('  case %d of %d (%s): %.0f s\n', c, numel(cases), label, toc(tCase));
end
fprintf('runRFITest: done in %.0f s\n\n', toc(tAll));

% Table
fprintf('%-34s %4s %9s %6s %7s %9s %8s %8s %16s %8s %9s\n', 'case', 'exc', 'blanked', 'good', ...
    'SNR', 'err [us]', 'rms/pred', 'redChi2', 'offset [us]', 'pulses', 'E missed');
for k = 1:numel(res)
    r = res(k);
    pc = '       -'; me = '        -';
    if ~isnan(r.pulses), pc = sprintf('%7.1f%%', 100*r.pulses); me = sprintf('%9.1e', r.missedE); end
    fprintf('%-34s %4d %8.4f%% %3d/%-2d %7.1f %9.3f %8.3f %8.2f %+8.2f +- %5.2f %s %s\n', r.label, ...
        r.excision, 100*r.blanked, r.nGood, r.nValid, r.snr, r.err, r.ratio, r.redChi2, ...
        r.offset, r.offsetErr, pc, me);
end
fprintf(['pulses: fraction of the pulsed source''s pulses blanked in its own channels; E missed: ' ...
         'fraction of its energy in the missed pulses (ground truth)\n']);
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
