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
  Processing: processIQ (IQ file in, fold out; frontEnd in pipelineParams:
  'fullband' or 'channels'), then foldToTOA (step 8)
    (  array beamforming                                         [todo] )
    full band:
    5. applyInverseDispersion  coherent dedispersion            [done]
    6. detectPower             square-law detection             [done]
    channelized (128 x 3.125 MHz):
    5a. channelizeIQ           polyphase filterbank              [done]
        RFI excision (pipelineParams 'excision'):
        detectRFI              power threshold per channel -> mask [done, B2]
        periodicRFI            periodic RFI: all predicted pulses  [done, B6]
        blankChannels          zero the masked samples           [done, B1]
    5b. dedisperseChannels     coherent dedispersion per channel [done]
    6b. detectChannels         power per channel + exact noise   [done]
        blankingWeights        mask -> exact data weights W, V, X [done, A3a]
    both:
    7. foldProfile             folding with a phase model       [done]
                               (channels: exact narrow-channel noise, NoiseCoeffs;
                                with excision the data weights, DataWeights)
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
Since 19 Oct 2026 the processing is two functions, processIQ(info_IQ, rx,
ephem, proc, files) up to the fold and foldToTOA(fold, ...) after it, which
get nothing else: the rule holds by their inputs.
%}

% Run control
runStage.sky      = false;   % pulsar signal + dispersion (the slow part; independent of noise)
runStage.receiver = true;    % receiver noise/RFI + IQ (rerun this alone to change SNR or RFI)
runStage.process  = true;   % processIQ: dedispersion + detection + fold (channel path: also
                             % filterbank, excision and blankingWeights); false = reuse the
                             % saved intermediate files and fold (checked)
runStage.toa      = true;   % foldToTOA (TOAs, detection, good TOAs) + validation

plots.dispersion = false;   % raw vs dispersed (reads the big RF files)
plots.iq         = false;   % dispersed IQ vs dedispersed IQ (full band only)
plots.detected   = true;    % detected power vs expected (full band only)
plots.fold       = true;    % folded profile vs expected (full band only)
plots.toa        = true;    % TOA errors vs ground truth (both paths)
closeFigures     = true;

if closeFigures, close all; end

% Paths
scriptDir = fileparts(mfilename('fullpath'));
addpath(fullfile(scriptDir, 'functions'));

% Parameters (pulsar, simulation, receiver, noise/RFI, processing, ephemeris, files)
pipelineParams;

% Disk estimate for the streamed files
% both paths: raw + dispersed + receiver (real float32 at f_in), IQ (complex at fs)
diskGB = L * (3*f_in*4 + fs*8) / 1e9;
fullBand = strcmp(frontEnd, 'fullband');
if ~fullBand && ~strcmp(frontEnd, 'channels')
    error('main:frontEnd', 'frontEnd must be ''fullband'' or ''channels''.');
end
if fullBand                     % dedispersed (complex at fs) + power
    diskGB = diskGB + L * (fs*8 + f_out*4) / 1e9;
else                            % channel IQ + dedispersed (4/3 oversampled) + power per channel
    diskGB = diskGB + L * (2*(fHigh - fLow)*4/3*8 + (fHigh - fLow)/chanWidth*f_out*4) / 1e9;
    if excision                 % blanked channel copies + weights (W, V, X(1..7) per bin)
        diskGB = diskGB + L * ((fHigh - fLow)*4/3*8 + (fHigh - fLow)/chanWidth*f_out*9*4) / 1e9;
    end
end
fprintf('main: L = %.3g s -> about %.1f GB of data files in "%s"\n', L, diskGB, dataDir);

% Sky: pulsar signal and interstellar dispersion
if runStage.sky
    info_gen = generatePulsarSignal(fileRaw, T, f_in, A, L, dutycycle, ...
        'Seed', seed, 'EnvelopeMode', genEnvMode, 'Components', profileComponents, 'Verbose', true);
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
checkConsistency(info_gen, info_disp, info_rx, info_IQ, ...
    T, f_in, L, A, dutycycle, genEnvMode, seed, DM, fLow, fHigh, refFreq, ...   % sky
    snrDB, noiseSeed, rfi, fLO, fs, filterOrder);                               % receiver
if ~fullBand && excision && snrDB > -5          % also noise-free (Inf); see detectRFI's help
    warning('main:brightExcision', ['Excision is on and the pulsar is bright (snrDB %g > -5): ' ...
        'detectRFI flags the pulse itself (10 %% of the peak windows at 0 dB, almost all from ' ...
        '+5 dB, every one noise-free), so amplitude, SNR and red. chi2 come out wrong. ' ...
        'Set excision = false for bright or noise-free checks.'], snrDB);
end

% Processing up to the fold: processIQ gets only the IQ file, the receiver settings rx, the
% ephemeris ephem and the processing settings proc (pipelineParams); the truth stays here
out = processIQ(info_IQ, rx, ephem, proc, files, 'Reuse', ~runStage.process);
info_det    = out.info_det;
info_dedisp = out.info_dedisp;      % full band ([] on the channel path)
info_fold   = out.info_fold;
fold        = out.fold;

% Check plots against the ground truth (full band)
if fullBand
    truthArgs = {'InfoDisp', info_disp, 'InfoIQ', info_IQ, 'InfoDedisp', info_dedisp, 'InfoRx', info_rx};
    if plots.iq
        [~, check_iq] = plotIQCheck(info_gen, info_disp, info_IQ, info_dedisp, 0, 20e-3, 'InfoRx', info_rx);
    end
    if plots.detected
        [~, check_det] = plotDetectedPower(info_det, info_gen, 0, 20e-3, truthArgs{:});
    end
    if plots.fold
        [~, check_fold] = plotFoldCheck(fold, info_fold, info_gen, truthArgs{:});
    end
elseif plots.iq || plots.detected || plots.fold
    fprintf('main: the iq / detected / fold check plots are for the full-band path only\n');
end

% TOAs, detection, good TOAs (observer), then validation against the ground truth
if runStage.toa
    [toa, detection, info_toaFit] = foldToTOA(fold, info_fold, ephem.profileFWHM, out.Bnoise, ...
        'Weighting', out.weighting);
    if any(toa.valid)
        val = validateTOA(toa, info_gen, 'Plot', plots.toa);
    else                        % e.g. empty phase bins in every sub-int (foldProfile:coverage)
        fprintf(['main: no sub-int gave a TOA (0 of %d fitted; see the warnings above), ' ...
                 'nothing to validate\n'], numel(toa.valid));
    end
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
        chUsed = 'none';
        if any(toa.valid), chUsed = mat2str(unique(toa.nChanUsed(toa.valid)).'); end
        fprintf('main: channel path, %d channels, weighting ''%s'', excision %d, channels used %s\n', ...
            info_fold.nChan, weighting, excision, chUsed);
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
function checkConsistency(info_gen, info_disp, info_rx, info_IQ, T, f_in, L, A, dutycycle, ...
    genEnvMode, seed, DM, fLow, fHigh, refFreq, snrDB, noiseSeed, rfi, fLO, fs, filterOrder)
%CHECKCONSISTENCY  Warn when reused files were made with other parameters.
% Sky files: checkSkyFiles (also used by the run scripts); then the receiver files.
d = checkSkyFiles(info_gen, info_disp, T, f_in, L, A, dutycycle, genEnvMode, seed, ...
    DM, fLow, fHigh, refFreq, profileComponents);
if ~(isequal(info_rx.snrDB, snrDB) || (isinf(snrDB) && info_rx.noiseStd == 0))
    d{end+1} = 'snrDB';
end
if ~strcmp(info_rx.inFile, info_disp.file), d{end+1} = 'receiver input'; end
rxRFI = struct([]);                                   % RFI sources the receiver file was made with
if isfield(info_rx, 'rfi') && ~isempty(info_rx.rfi), rxRFI = [info_rx.rfi.params]; end
if ~(isempty(rxRFI) && isempty(rfi)) && ~isequaln(rxRFI, rfi)   % isequaln: NaN parameters match
    d{end+1} = 'RFI';
end
% the seed only matters with receiver noise (RFI needs noise too: addNoiseAndRFI:ref)
if info_rx.noiseStd > 0 && info_rx.seed ~= noiseSeed
    d{end+1} = 'noiseSeed';
end
if abs(info_IQ.fLO - fLO) > 1e-3,     d{end+1} = 'fLO';  end
if abs(info_IQ.actualFsOut - f_in/round(f_in/fs)) > 1e-3
    d{end+1} = 'fs';
end
M = filterOrder + 1; M = M + (mod(M, 2) == 0);       % filter length as applyIQmodulation makes it
if info_IQ.filterLen ~= M,            d{end+1} = 'filterOrder'; end
if ~isempty(d)
    warning('main:stale', ['Files on disk were generated with different ' ...
        'parameters (%s). Rerun the corresponding stage.'], strjoin(d, ', '));
end
end
