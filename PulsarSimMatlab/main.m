%{
PULSAR NAVIGATION - synthetic data pipeline (main script)

Main thought: pulsar signals arrive at a known location at very accurate
times. A pulse arriving earlier means we are closer to the pulsar in that
direction; combining several pulsars gives our position.

Pipeline (status):
  Sky (synthetic)
    1. generatePulsarSignal    noise-like pulses, ground truth  [done]
    2. applyDispersionStream   interstellar dispersion          [done]
  Receiver (synthetic)
    3. addNoiseAndRFI          receiver noise + RFI at RF       [done]
    4. applyIQmodulation       downconversion to baseband       [done]
    -  clock jitter, pulsar-Earth motion, polarisation?,
       3x3 array element signals, multiple pulsars               [todo]
  Processing (frontEnd in pipelineParams: 'fullband' or 'channels')
    (  array beamforming                                         [todo] )
    full band:
    5. applyInverseDispersion  coherent dedispersion            [done]
    6. detectPower             square-law detection             [done]
    channelized (128 x 3.125 MHz):
    5a. channelizeIQ           polyphase filterbank              [done]
        (RFI excision per channel: blank samples, blankingWeights  [todo, B])
    5b. dedisperseChannels     coherent dedispersion per channel [done]
    6b. detectChannels         power per channel + exact noise   [done]
    both:
    7. foldProfile             folding with a phase model       [done]
                               (channels: exact narrow-channel noise, NoiseCoeffs)
    8. estimateTOA             template matching (FFTFIT; channels: 'Weighting') [done]
       detectPulsar            noise normalization, NP detector,
                               TOA quality (good = detected, chi2 ok)  [done]
    -  barycentric / timing corrections, residuals               [todo]
    -  navigation solution (multi-pulsar)                        [todo]
  Validation
    - plot/check functions per stage against ground truth        [done]
    9. validateTOA             TOA vs ground truth              [done]
    - Monte Carlo: runMonteCarlo (noise seeds), runSNRSweep (SNR,
      detection, H0), runH0 (long noise-only run)                [done]

Rule: the PROCESSING stages only use what an observer would know
(ephemeris + receiver settings). Ground truth (info_gen, info_disp,
info_rx) is only used by the synthetic stages and the check/plot functions.
%}

% Run control
runStage.sky      = false;   % pulsar signal + dispersion (the slow part; independent of noise)
runStage.receiver = true;    % receiver noise/RFI + IQ (rerun this alone to change SNR or RFI)
runStage.process  = true;   % dedispersion + detection
runStage.fold     = true;
runStage.toa      = true;   % TOA estimation + validation (needs the fold)

plots.dispersion = false;   % raw vs dispersed (reads the big RF files)
plots.iq         = false;   % dispersed IQ vs dedispersed IQ
plots.detected   = true;
plots.fold       = true;
plots.toa        = true;
closeFigures     = true;

if closeFigures, close all; end

% Paths
scriptDir = fileparts(mfilename('fullpath'));
addpath(fullfile(scriptDir, 'functions'));

% Parameters (pulsar, simulation, receiver, noise/RFI, processing, ephemeris, files)
pipelineParams;

% Disk estimate for the streamed files
% raw + dispersed + receiver (real float32 at f_in), IQ + dedispersed (complex at fs), power
diskGB = L * (3*f_in*4 + 2*fs*8 + f_out*4) / 1e9;
fullBand = strcmp(frontEnd, 'fullband');
if ~fullBand && ~strcmp(frontEnd, 'channels')
    error('main:frontEnd', 'frontEnd must be ''fullband'' or ''channels''.');
end
if ~fullBand                    % channel IQ + dedispersed (4/3 oversampled) + power per channel
    diskGB = diskGB + L * (2*(fHigh - fLow)*4/3*8 + (fHigh - fLow)/chanWidth*f_out*4) / 1e9;
end
fprintf('main: L = %.3g s -> about %.1f GB of data files in "%s"\n', L, diskGB, dataDir);

% Sky: pulsar signal and interstellar dispersion
if runStage.sky
    info_gen = generatePulsarSignal(fileRaw, T, f_in, A, L, dutycycle, ...
        'Seed', seed, 'EnvelopeMode', genEnvMode, 'Verbose', true);
    info_disp = applyDispersionStream(info_gen.file, fileDispersed, DM, ...
        info_gen.actualFsOut, fLow, fHigh, 'RefFreq', refFreq);
else
    info_gen  = loadInfo(fileRaw);
    info_disp = loadInfo(fileDispersed);
end

if plots.dispersion
    plotDispersionCheck(info_gen, info_disp, 0, 20e-3);
end

% Receiver: noise + interference at RF (not dispersed), then IQ conversion
if runStage.receiver
    info_rx = addNoiseAndRFI(info_disp.file, fileRx, info_disp.actualFsOut, ...
        'Band', [fLow fHigh], 'SNRdB', snrDB, 'SignalInfo', info_gen, ...
        'RFI', rfi, 'Seed', noiseSeed);
    info_IQ = applyIQmodulation(info_rx.file, fileIQ, info_rx.actualFsOut, fs, fLO, ...
        'FilterOrder', filterOrder, 'Band', [fLow fHigh]);
else
    info_rx = loadInfo(fileRx);
    info_IQ = loadInfo(fileIQ);
end
checkConsistency(info_gen, info_disp, info_rx, info_IQ, T, f_in, L, DM, fLow, fHigh, fLO, fs, snrDB, rfi);

% Processing (uses only the ephemeris 'ephem' and receiver settings)
if fullBand
    if runStage.process
        % Coherent dedispersion
        info_dedisp = applyInverseDispersion(info_IQ.file, fileDedisp, info_IQ.actualFsOut, ...
            info_IQ.fLO, ephem.DM, fLow, fHigh, 'RefFreq', refFreq);
        % Square-law detection
        info_det = detectPower(info_dedisp.file, fileEnvelope, info_dedisp.fs, info_dedisp.fLO, ...
            f_out, 'FullySupported', info_dedisp.fullySupported);
    else
        info_dedisp = loadInfo(fileDedisp);
        info_det    = loadInfo(fileEnvelope);
    end
else
    if runStage.process
        % Filterbank: 128 channels of 3.125 MHz (oversampled 4/3)
        info_chan = channelizeIQ(info_IQ.file, fileChanIQ, info_IQ.actualFsOut, info_IQ.fLO, ...
            fLow, fHigh, 'ChanWidth', chanWidth, 'T0', info_IQ.t0);
        % (RFI excision per channel goes here: blank samples; mask -> blankingWeights; block B)
        % Coherent dedispersion per channel, all channels aligned at refFreq
        info_dc = dedisperseChannels(info_chan, fileChanDedisp, ephem.DM, 'RefFreq', refFreq);
        % Square-law detection per channel; bins of a whole number of channel samples
        info_det = detectChannels(info_dc, fileChanPower, info_dc.fs / round(info_dc.fs / f_out));
    else
        info_dc  = loadInfo(fileChanDedisp);
        info_det = loadInfo(fileChanPower);
    end
end

if fullBand
    truthArgs = {'InfoDisp', info_disp, 'InfoIQ', info_IQ, 'InfoDedisp', info_dedisp, 'InfoRx', info_rx};
    if plots.iq
        [~, check_iq] = plotIQCheck(info_gen, info_disp, info_IQ, info_dedisp, 0, 20e-3, 'InfoRx', info_rx);
    end
    if plots.detected
        [~, check_det] = plotDetectedPower(info_det, info_gen, 0, 20e-3, truthArgs{:});
    end
elseif plots.iq || plots.detected || plots.fold
    fprintf('main: the iq / detected / fold check plots are for the full-band path only\n');
end

% Folding
if runStage.fold
    if fullBand
        [info_fold, fold] = foldProfile(info_det, fileFold, ephem.f0, ...
            'F1', ephem.F1, 'TRef', ephem.TRef, 'NBin', nBin, 'SubintPeriods', subintPeriods);
        if plots.fold
            [~, check_fold] = plotFoldCheck(fold, info_fold, info_gen, truthArgs{:});
        end
    else
        % narrow channels: neighbouring time bins correlate -> exact phase-bin noise
        [info_fold, fold] = foldProfile(info_det, fileFoldChan, ephem.f0, ...
            'F1', ephem.F1, 'TRef', ephem.TRef, 'NBin', nBin, 'SubintPeriods', subintPeriods, ...
            'NoiseCoeffs', [info_det.noise.V, info_det.noise.X]);
    end
elseif runStage.toa
    if fullBand, S = load(fileFold, 'fold', 'info'); else, S = load(fileFoldChan, 'fold', 'info'); end
    fold = S.fold; info_fold = S.info;                  % reuse the saved fold
end

% TOA estimation and validation
if runStage.toa
    template = gaussianTemplate(nBin, ephem.profileFWHM);
    % Noise-equivalent bandwidth for the radiometer noise model, from the
    % observer's own dedispersion taper (receiver noise only passes this one).
    % When the pulsar dominates (high SNR, in-pulse) the true value is ~0.5 %
    % lower, since the simulated signal also passed the dispersion-stage taper.
    % Channel path: Bnoise of ONE channel (its own dedispersion taper), and the
    % channels combined as set by 'weighting'.
    if fullBand
        Bnoise = noiseBandwidth(fLow, fHigh, info_dedisp.edgeWidth);
        toaArgs = {'Bnoise', Bnoise};
    else
        Bnoise = info_det.noise.Bnoise;
        toaArgs = {'Bnoise', Bnoise, 'Weighting', weighting};
    end
    [toa, info_toa] = estimateTOA(fold, info_fold, template, toaArgs{:});

    % Detection (Neyman-Pearson, P_FA 1e-3 per sub-int) and TOA quality.
    % Good TOA: fitted, pulsar detected without using the ephemeris phase
    % (unknown phase), and the fit residuals consistent with the noise (chi^2).
    [detection, info_detect] = detectPulsar(fold, info_fold, template, toaArgs{:});
    toa.good = toa.valid & detection.detectedUnknown & ~toa.flagChi2;
    fprintf('main: %d of %d fitted TOAs good (%d not detected, %d chi2-flagged)\n', ...
        nnz(toa.good), nnz(toa.valid), nnz(toa.valid & ~detection.detectedUnknown), ...
        nnz(toa.valid & toa.flagChi2));

    val = validateTOA(toa, info_gen, 'Plot', plots.toa);
    if any(toa.good) && nnz(toa.good) < nnz(toa.valid)
        toaGood = toa; toaGood.valid = toa.good;          % validate the good TOAs only
        valGood = validateTOA(toaGood, info_gen, 'Plot', false, 'Label', 'good TOAs');
    end
    % Best achievable (ground truth, exact band tapers; see expectedPowerModel;
    % full-band path; the channel path uses the same 390 MHz of band)
    if fullBand
        M = expectedPowerModel(info_gen, info_disp, info_IQ, info_dedisp, info_rx);
    else
        M = struct('hasNoise', false);
        fprintf('main: channel path, %d channels, weighting ''%s'', channels used %s\n', ...
            info_fold.nChan, weighting, mat2str(unique(toa.nChanUsed(toa.valid)).'));
    end
    if M.hasNoise
        fprintf(['main: predicted best TOA error per sub-int (%d turns) %.4g us, ' ...
                 'SNR per sub-int %.1f; whole file %.4g us\n'], subintPeriods, ...
            M.toaErrPulse/sqrt(subintPeriods)*1e6, ...
            M.snrPulse*sqrt(subintPeriods), ...
            M.toaErrPulse/sqrt(info_gen.nPulses)*1e6);
    end
end


% ======================================================================
%  Local functions
%  ======================================================================
function checkConsistency(info_gen, info_disp, info_rx, info_IQ, T, f_in, L, DM, fLow, fHigh, fLO, fs, snrDB, rfi)
%CHECKCONSISTENCY  Warn when reused files were made with other parameters.
d = {};
if abs(info_gen.T - T) > 1e-12,       d{end+1} = 'T';    end
if abs(info_gen.fs - f_in) > 1e-3,    d{end+1} = 'f_in'; end
if abs(info_gen.L - L) > 1e-12,       d{end+1} = 'L';    end
if abs(info_disp.DM - DM) > 1e-12,    d{end+1} = 'DM';   end
if abs(info_disp.fLow - fLow) > 1e-3 || abs(info_disp.fHigh - fHigh) > 1e-3
    d{end+1} = 'band';
end
if ~(isequal(info_rx.snrDB, snrDB) || (isinf(snrDB) && info_rx.noiseStd == 0))
    d{end+1} = 'snrDB';
end
if ~strcmp(info_rx.inFile, info_disp.file), d{end+1} = 'receiver input'; end
rxRFI = struct([]);                                   % RFI sources the receiver file was made with
if isfield(info_rx, 'rfi') && ~isempty(info_rx.rfi), rxRFI = [info_rx.rfi.params]; end
if ~(isempty(rxRFI) && isempty(rfi)) && ~isequaln(rxRFI, rfi)   % isequaln: NaN parameters match
    d{end+1} = 'RFI';
end
if abs(info_IQ.fLO - fLO) > 1e-3,     d{end+1} = 'fLO';  end
if abs(info_IQ.actualFsOut - f_in/round(f_in/fs)) > 1e-3
    d{end+1} = 'fs';
end
if ~isempty(d)
    warning('main:stale', ['Files on disk were generated with different ' ...
        'parameters (%s). Rerun the corresponding stage.'], strjoin(d, ', '));
end
end