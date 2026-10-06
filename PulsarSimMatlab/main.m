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
  Processing
    (  array beamforming, RFI excision - linear, before squaring  [todo] )
    5. applyInverseDispersion  coherent dedispersion            [done]
    6. detectPower             square-law detection             [done]
    7. foldProfile             folding with a phase model       [done]
    8. estimateTOA             FFT template matching (FFTFIT)   [done]
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
runStage.receiver = true;   % receiver noise/RFI + IQ (rerun this alone to change SNR or RFI)
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

truthArgs = {'InfoDisp', info_disp, 'InfoIQ', info_IQ, 'InfoDedisp', info_dedisp, 'InfoRx', info_rx};
if plots.iq
    [~, check_iq] = plotIQCheck(info_gen, info_disp, info_IQ, info_dedisp, 0, 20e-3, 'InfoRx', info_rx);
end
if plots.detected
    [~, check_det] = plotDetectedPower(info_det, info_gen, 0, 20e-3, truthArgs{:});
end

% Folding
if runStage.fold
    [info_fold, fold] = foldProfile(info_det, fileFold, ephem.f0, ...
        'F1', ephem.F1, 'TRef', ephem.TRef, 'NBin', nBin, 'SubintPeriods', subintPeriods);
    if plots.fold
        [~, check_fold] = plotFoldCheck(fold, info_fold, info_gen, truthArgs{:});
    end
elseif runStage.toa
    S = load(fileFold, 'fold', 'info');                 % reuse the saved fold
    fold = S.fold; info_fold = S.info;
end

% TOA estimation and validation
if runStage.toa
    template = gaussianTemplate(nBin, ephem.profileFWHM);
    % Noise-equivalent bandwidth for the radiometer noise model, from the
    % observer's own dedispersion taper (receiver noise only passes this one).
    % When the pulsar dominates (high SNR, in-pulse) the true value is ~0.5 %
    % lower, since the simulated signal also passed the dispersion-stage taper.
    Bnoise = noiseBandwidth(fLow, fHigh, info_dedisp.edgeWidth);
    [toa, info_toa] = estimateTOA(fold, info_fold, template, 'Bnoise', Bnoise);

    % Detection (Neyman-Pearson, P_FA 1e-3 per sub-int) and TOA quality.
    % Good TOA: fitted, pulsar detected without using the ephemeris phase
    % (unknown phase), and the fit residuals consistent with the noise (chi^2).
    [detection, info_detect] = detectPulsar(fold, info_fold, template, 'Bnoise', Bnoise);
    toa.good = toa.valid & detection.detectedUnknown & ~toa.flagChi2;
    fprintf('main: %d of %d fitted TOAs good (%d not detected, %d chi2-flagged)\n', ...
        nnz(toa.good), nnz(toa.valid), nnz(toa.valid & ~detection.detectedUnknown), ...
        nnz(toa.valid & toa.flagChi2));

    val = validateTOA(toa, info_gen, 'Plot', plots.toa);
    if any(toa.good) && nnz(toa.good) < nnz(toa.valid)
        toaGood = toa; toaGood.valid = toa.good;          % validate the good TOAs only
        valGood = validateTOA(toaGood, info_gen, 'Plot', false, 'Label', 'good TOAs');
    end
    % Best achievable (ground truth, exact band tapers; see expectedPowerModel)
    M = expectedPowerModel(info_gen, info_disp, info_IQ, info_dedisp, info_rx);
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