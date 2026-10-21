%{
RUNMONTECARLO - TOA error bars over many receiver-noise realizations

The pulsar signal is fixed (sky files from a main.m run, same 'seed');
only the receiver noise changes per run. Each realization is simulated
(addNoiseAndRFI, applyIQmodulation) and then processed exactly as in
main.m: processIQ (IQ -> fold) and foldToTOA (fold -> TOAs, detection,
good TOAs), on the front end pipelineParams chooses (default: channels with
RFI excision and optimal weighting; 'fullband' the reference path). All
TOAs are pooled in validateTOA. For correct error bars the pooled
normalized errors (TOA - truth)/toaErr are N(0,1): rms ratio ~1, |z| < 1
for ~68 %.

Also pooled: pChi2 over all sub-ints (uniform on [0, 1] when the noise
model is right: mean 0.5, 5 % below 0.05, Kolmogorov-Smirnov distance), the
detection rate (unknown phase), and on the channel path with excision the
blanked fraction (on pure noise: the false-flag rate of detectRFI).

Fixed pulsar: the pulsar's own noise is the same in every run, so its
pure self-noise part of the variance (a few % at -5 dB; at the pulse peak
s^2/(1+s)^2 ~ 6 %, s = peak signal/noise ~0.32) is not sampled: expect the
ratio slightly below 1. runSNRSweep varies the pulsar seed as well.

Check: seed i = 1 equals main's noiseSeed (seed + 1), so its TOAs must be
identical to those of main.m with the same parameters (either front end).

Needs: the sky files (runStage.sky in main.m at least once with the
current parameters). Parameters come from pipelineParams.m. Intermediate
files go to data/mc and are overwritten for every seed.
(19 Oct 2026: on processIQ / foldToTOA, both front ends - A3b, the
multi-seed Monte Carlo of the default channel path. Before: full band only.)
%}

nSeeds = 10;            % number of receiver-noise realizations

% Paths and parameters
scriptDir = fileparts(mfilename('fullpath'));
addpath(fullfile(scriptDir, 'functions'));
pipelineParams;

mcDir  = fullfile(dataDir, "mc");
mcChan = fullfile(mcDir, "chan");
mcRx   = fullfile(mcDir, "mc_rx.dat");
mcIQ   = fullfile(mcDir, "mc_rx_IQ.dat");
% processIQ's output files (own names, so main's files in data/ stay)
mcFiles = struct('dedisp', fullfile(mcDir, "mc_dedisp.dat"), ...
    'envelope', fullfile(mcDir, "mc_envelope.dat"), 'fold', fullfile(mcDir, "mc_fold.mat"), ...
    'chanIQ', fullfile(mcChan, "mc_rx_IQ_chan"), 'chanBlank', fullfile(mcChan, "mc_rx_IQ_chan_blanked"), ...
    'chanDedisp', fullfile(mcChan, "mc_dedisp_chan"), 'chanPower', fullfile(mcChan, "mc_dedisp_chan_power.dat"), ...
    'chanDedispB', fullfile(mcChan, "mc_dedispB_chan"), 'chanPowerB', fullfile(mcChan, "mc_dedispB_chan_power.dat"), ...
    'chanWeight', fullfile(mcChan, "mc_dedispB_chan_weights.dat"), 'foldChan', fullfile(mcDir, "mc_fold_chan.mat"));
if ~isfolder(mcChan), mkdir(mcChan); end

% Sky files from main.m; refuse if they were made with other parameters
info_gen  = loadInfo(fileRaw);
info_disp = loadInfo(fileDispersed);
d = checkSkyFiles(info_gen, info_disp, T, f_in, L, A, dutycycle, genEnvMode, seed, ...
    DM, fLow, fHigh, refFreq, profileComponents);
if ~isempty(d)
    error('runMonteCarlo:stale', ['Sky files do not match pipelineParams (%s). ' ...
        'Run main.m with runStage.sky = true first.'], strjoin(d, ', '));
end
fullBand = strcmp(frontEnd, 'fullband');
excisionOn = ~fullBand && excision;
if fullBand && ~isempty(rfi)
    warning('runMonteCarlo:noExcision', ['RFI is on (rfiOn) but the full-band path has no ' ...
        'excision: the RFI goes into the TOAs unremoved.']);
end

pathStr = 'full band';
if ~fullBand, pathStr = sprintf('channels, excision %d, weighting ''%s''', excisionOn, weighting); end
fprintf('runMonteCarlo: %d noise realizations, SNR %.1f dB, L = %.3g s, %d turn(s) per sub-int, %s\n', ...
    nSeeds, snrDB, L, subintPeriods, pathStr);
blanked = nan(nSeeds, 1);
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

    % Processing (ephemeris + receiver settings only), as in main.m
    out = processIQ(info_IQ, rx, ephem, proc, mcFiles, 'Verbose', false);
    [toa, detection] = foldToTOA(out.fold, out.info_fold, ephem.profileFWHM, out.Bnoise, ...
        'Weighting', out.weighting, 'Verbose', false);

    toaAll(i) = toa; %#ok<SAGROW>
    detAll(i) = detection; %#ok<SAGROW>
    if excisionOn, blanked(i) = mean(out.info_w.blankedFraction); end
    v = toa.valid;
    blankStr = '';
    if excisionOn, blankStr = sprintf(', blanked %.4f %%', 100*blanked(i)); end
    fprintf(['  seed %d (%d/%d): %d TOAs (%d good), median SNR %.1f, median TOA error %.3f us, ' ...
             'median red. chi2 %.4f, min pChi2 %.3g%s, %.0f s\n'], ns, i, nSeeds, ...
        nnz(v), nnz(toa.good), median(toa.snr(v)), median(toa.toaErr(v))*1e6, ...
        median(toa.redChi2(v)), min(toa.pChi2(v)), blankStr, toc(tSeed));
end
fprintf('runMonteCarlo: done in %.0f s\n', toc(tAll));

% Pooled validation against the ground truth: all fitted TOAs, then the good ones
label = sprintf('MC %d seeds, %.1f dB, %s,', nSeeds, snrDB, frontEnd);
val = validateTOA(toaAll, repmat(info_gen, 1, nSeeds), 'Label', label);
toaGood = toaAll;
for i = 1:nSeeds, toaGood(i).valid = toaAll(i).good; end
if nnz([toaGood.valid]) < nnz([toaAll.valid])
    valGood = validateTOA(toaGood, repmat(info_gen, 1, nSeeds), 'Plot', false, ...
        'Label', [label ' good TOAs']);
end

% Noise-model check: pChi2 uniform on [0, 1] (Kolmogorov-Smirnov, 1e-3 level D < 1.95/sqrt(n))
p  = sort(vertcat(toaAll.pChi2));
p  = p(isfinite(p)); n = numel(p);
u  = (1:n).' / n;
Dks = max(max(u - p), max(p - (u - 1/n)));
fprintf(['runMonteCarlo: pChi2 over %d sub-ints: mean %.3f (0.5 +- %.3f), fraction < 0.05 %.3f ' ...
         '(0.05 +- %.3f), KS D %.3f (1e-3 level %.3f)\n'], n, mean(p), 1/sqrt(12*n), ...
    mean(p < 0.05), sqrt(0.05*0.95/n), Dks, 1.95/sqrt(n));
vAll = vertcat(toaAll.valid);
dUnk = vertcat(detAll.detectedUnknown);
fprintf('runMonteCarlo: detected (unknown phase) %d of %d fitted sub-ints; good TOAs %d\n', ...
    nnz(dUnk & vAll), nnz(vAll), nnz(vertcat(toaAll.good)));
if excisionOn
    fprintf('runMonteCarlo: blanked fraction mean %.4f %% (per seed %.4f-%.4f %%)\n', ...
        100*mean(blanked), 100*min(blanked), 100*max(blanked));
end

% Best achievable (ground truth, exact band tapers; see expectedPowerModel; full band only:
% the channel path has no such prediction yet, notes §10 item 8)
if fullBand
    M = expectedPowerModel(info_gen, info_disp, info_IQ, out.info_dedisp, info_rx);
    if M.hasNoise
        fprintf(['runMonteCarlo: predicted best TOA error per sub-int (%d turns) %.4g us, ' ...
                 'SNR per sub-int %.1f\n'], subintPeriods, ...
            M.toaErrPulse/sqrt(subintPeriods)*1e6, ...
            M.snrPulse*sqrt(subintPeriods));
    end
end
