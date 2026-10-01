%{
RUNMONTECARLO - TOA error bars over many receiver-noise realizations

The pulsar signal is fixed (sky files from a main.m run, same 'seed');
only the receiver noise changes per run. Each realization goes through
the same receiver and processing stages as main.m, and all TOAs are
pooled in validateTOA. For correct error bars the pooled normalized
errors (TOA - truth)/toaErr are N(0,1): rms ratio ~1, |z| < 1 for ~68 %.

Check: seed i = 1 equals main's noiseSeed (seed + 1), so its TOAs must be
identical to those of main.m with the same parameters.

Needs: the sky files (runStage.sky in main.m at least once with the
current parameters). Parameters come from pipelineParams.m. Intermediate
files go to data/mc and are overwritten for every seed.
%}

nSeeds = 10;            % number of receiver-noise realizations

% Paths and parameters
scriptDir = fileparts(mfilename('fullpath'));
addpath(fullfile(scriptDir, 'functions'));
pipelineParams;

mcDir        = fullfile(dataDir, "mc");
mcRx         = fullfile(mcDir, "mc_rx.dat");
mcIQ         = fullfile(mcDir, "mc_rx_IQ.dat");
mcDedisp     = fullfile(mcDir, "mc_dedisp.dat");
mcEnvelope   = fullfile(mcDir, "mc_envelope.dat");
mcFold       = fullfile(mcDir, "mc_fold.mat");    % not saved (SaveFile false)

% Sky files from main.m; refuse if they were made with other parameters
info_gen  = loadInfo(fileRaw);
info_disp = loadInfo(fileDispersed);
if abs(info_gen.T - T) > 1e-12 || abs(info_gen.L - L) > 1e-12 || ...
        abs(info_gen.fs - f_in) > 1e-3 || info_gen.seed ~= seed || ...
        abs(info_disp.DM - DM) > 1e-12 || ...
        abs(info_disp.fLow - fLow) > 1e-3 || abs(info_disp.fHigh - fHigh) > 1e-3
    error('runMonteCarlo:stale', ['Sky files do not match pipelineParams ' ...
        '(T, L, f_in, seed, DM or band). Run main.m with runStage.sky = true first.']);
end

template = gaussianTemplate(nBin, ephem.profileFWHM);

fprintf('runMonteCarlo: %d noise realizations, SNR %.1f dB, L = %.3g s, %d turn(s) per sub-int\n', ...
    nSeeds, snrDB, L, subintPeriods);
tAll = tic;
for i = 1:nSeeds
    tSeed = tic;
    ns = seed + i;                              % i = 1 is main's noiseSeed

    % Receiver (synthetic): noise + RFI at RF, IQ conversion
    info_rx = addNoiseAndRFI(info_disp.file, mcRx, info_disp.actualFsOut, ...
        'Band', [fLow fHigh], 'SNRdB', snrDB, 'SignalInfo', info_gen, ...
        'RFI', rfi, 'Seed', ns, 'Verbose', false);
    info_IQ = applyIQmodulation(info_rx.file, mcIQ, info_rx.actualFsOut, fs, fLO, ...
        'FilterOrder', filterOrder, 'Band', [fLow fHigh], 'Verbose', false);

    % Processing (ephemeris + receiver settings only), same calls as main.m
    info_dedisp = applyInverseDispersion(info_IQ.file, mcDedisp, info_IQ.actualFsOut, ...
        info_IQ.fLO, ephem.DM, fLow, fHigh, 'RefFreq', refFreq, 'Verbose', false);
    info_det = detectPower(info_dedisp.file, mcEnvelope, info_dedisp.fs, info_dedisp.fLO, ...
        f_out, 'FullySupported', info_dedisp.fullySupported, 'Verbose', false);
    [info_fold, fold] = foldProfile(info_det, mcFold, ephem.f0, ...
        'F1', ephem.F1, 'TRef', ephem.TRef, 'NBin', nBin, 'SubintPeriods', subintPeriods, ...
        'SaveFile', false, 'Verbose', false);
    Bnoise = noiseBandwidth(fLow, fHigh, info_dedisp.edgeWidth);
    toa = estimateTOA(fold, info_fold, template, 'Bnoise', Bnoise, 'Verbose', false);

    toaAll(i) = toa; %#ok<SAGROW>
    v = toa.valid;
    fprintf('  seed %d (%d/%d): %d TOAs, median SNR %.1f, median TOA error %.3f us, %.0f s\n', ...
        ns, i, nSeeds, nnz(v), median(toa.snr(v)), median(toa.toaErr(v))*1e6, toc(tSeed));
end
fprintf('runMonteCarlo: done in %.0f s\n', toc(tAll));

% Pooled validation against the ground truth
val = validateTOA(toaAll, repmat(info_gen, 1, nSeeds), ...
    'Label', sprintf('MC %d seeds, %.1f dB,', nSeeds, snrDB));
if isfield(info_rx.prediction, 'toaErrPulse')
    fprintf(['runMonteCarlo: predicted best TOA error per sub-int (%d turns) %.3g us, ' ...
             'SNR per sub-int %.1f\n'], subintPeriods, ...
        info_rx.prediction.toaErrPulse/sqrt(subintPeriods)*1e6, ...
        info_rx.prediction.snrPulse*sqrt(subintPeriods));
end
