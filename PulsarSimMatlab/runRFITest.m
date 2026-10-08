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

Needs: the sky files (main.m with runStage.sky once). Files go to data/rfi
and are overwritten per case; main's files and data/chan are not touched.
~1 min per case.
%}

% Paths and parameters
scriptDir = fileparts(mfilename('fullpath'));
addpath(fullfile(scriptDir, 'functions'));
pipelineParams;
if ~strcmp(frontEnd, 'channels')
    error('runRFITest:frontEnd', 'runRFITest is for the channel path (frontEnd = ''channels'').');
end

cases = {[], 1, 2, 3, 4, 5, 1:5};     % indices into rfiScenario; [] = noise only

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

nC = numel(cases);
res = struct('label', {}, 'excision', {}, 'blanked', {}, 'nGood', {}, 'nValid', {}, ...
    'snr', {}, 'err', {}, 'ratio', {}, 'redChi2', {}, 'offset', {}, 'offsetErr', {});
fprintf('runRFITest: %d cases, SNR %.1f dB, L = %.3g s, noise seed %d\n', nC, snrDB, L, noiseSeed);
tAll = tic;
for c = 1:nC
    tCase = tic;
    idx = cases{c};
    if isempty(idx)
        rfiCase = struct([]); label = 'noise only';
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
        else
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
            'offsetErr', val.totalOffsetErr * 1e6); %#ok<SAGROW>
    end
    fprintf('  case %d/%d (%s): %.0f s\n', c, nC, label, toc(tCase));
end
fprintf('runRFITest: done in %.0f s\n\n', toc(tAll));

% Table
fprintf('%-34s %4s %9s %6s %7s %9s %8s %8s %16s\n', 'case', 'exc', 'blanked', 'good', ...
    'SNR', 'err [us]', 'rms/pred', 'redChi2', 'offset [us]');
for k = 1:numel(res)
    r = res(k);
    fprintf('%-34s %4d %8.4f%% %3d/%-2d %7.1f %9.3f %8.3f %8.2f %+8.2f +- %5.2f\n', r.label, ...
        r.excision, 100*r.blanked, r.nGood, r.nValid, r.snr, r.err, r.ratio, r.redChi2, ...
        r.offset, r.offsetErr);
end
save(fullfile(rfiDir, "runRFITest_results.mat"), 'res', 'cases');
