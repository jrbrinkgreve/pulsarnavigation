%{
RUNRFIPROFILE - per-channel RFI statistics for the fast simulator (block C4b)

Runs the voltage chain without pulsar (zero sky + receiver noise + RFI, main's noise seed)
through processIQ on the channel path with excision, then measures with rfiProfile what the
RFI leaves in each channel's detected power: level, variance and lag covariances of the time
bins (to lag 30 since 21 Oct 2026: narrow carriers correlate beyond the clean model's 7), valid
fraction. simulateFold('RFI', prof) reproduces these in folds of any length.

Cases:
  'realistic'      rfiRealistic of pipelineParams (16 sources, B5e)
  'inmarsatNarrow' the same with Inmarsat as ten narrow (200 kHz) carriers in channels
                   105-114 (runRFITest case 15, the worst case for the error bars)
  'noiseOnly'      no RFI (sanity check: level 1, noise coefficients = the clean model)
Lengths:
  1 s    steady part (production profiles): the en-route radar's beam passage is moved
         outside the run (BeamTime = Lcal + 5 s: sidelobes only, as most of its 10 s scan);
         beam passages are time structure (block C3)
  0.1 s  exactly runRFITest's settings (beam passage at L/2 inside the run), to validate the
         simulator against runRFITest cases 14 and 15
Saves data/rfi/profile/rfiProfile_<case>_<L>s.mat (prof); the large intermediate files of
each run are deleted afterwards. Needs no sky files. Run time: 0.1 s ~1.5 min (noise only)
to ~8 min (RFI: the 16 sources' generation and the blanking weights dominate); 1 s ~10x.
(20 Oct 2026, block C4b.)
%}

lengths = [0.1 1];                   % [s]
caseNames = {'realistic', 'inmarsatNarrow', 'noiseOnly'};
overwrite = false;                   % false: skip cases whose profile file exists

scriptDir = fileparts(mfilename('fullpath'));
addpath(fullfile(scriptDir, 'functions'));
pipelineParams;
pDir = fullfile(dataDir, "rfi", "profile");
if ~isfolder(pDir), mkdir(pDir); end
pFiles = struct('chanIQ', fullfile(pDir, "p_chan"), 'chanBlank', fullfile(pDir, "p_chan_blanked"), ...
    'chanDedisp', fullfile(pDir, "p_dedisp_chan"), 'chanPower', fullfile(pDir, "p_dedisp_chan_power.dat"), ...
    'chanDedispB', fullfile(pDir, "p_dedispB_chan"), 'chanPowerB', fullfile(pDir, "p_dedispB_chan_power.dat"), ...
    'chanWeight', fullfile(pDir, "p_dedispB_chan_weights.dat"), 'foldChan', fullfile(pDir, "p_fold_chan.mat"));
procP = proc; procP.frontEnd = 'channels'; procP.excision = true;
procP.subintPeriods = 10;            % the fold is not used; fewer sub-ints, less memory
procP.weightsArgs = {'MaxMemoryGB', 8};   % blankingWeights holds ~4.8 GB at 1 s (26 GB machine)

% Inmarsat as narrow carriers (as runRFITest case 15)
inmF = fLow + ((105:114) - 0.5) * chanWidth + [-1.1 0.4 -0.6 1.2 -0.2 0.9 -1.3 0.1 0.7 -0.9] * 1e6;
inmNarrow = arrayfun(@(f) rfiSource('noise', 'Freq', f, 'Bandwidth', 200e3, 'INRdB', -21, ...
    'Label', sprintf('Inmarsat carrier %.2f MHz', f / 1e6)), inmF, 'UniformOutput', false);
inmNarrow = [inmNarrow{:}];
isInm = startsWith({rfiRealistic.label}, 'Inmarsat');

tAll = tic;
for Lc = lengths
    zeroFile = fullfile(pDir, "p_zero.dat");
    info_zero = generatePulsarSignal(zeroFile, T, f_in, 0, Lc, dutycycle, ...
        'Seed', seed, 'EnvelopeMode', genEnvMode, 'Verbose', false);
    for ci = 1:numel(caseNames)
        tCase = tic;
        fOut = fullfile(pDir, sprintf('rfiProfile_%s_%gs.mat', caseNames{ci}, Lc));
        if ~overwrite && isfile(fOut)
            fprintf('runRFIProfile: %s exists, skipped\n', fOut);
            continue
        end
        switch caseNames{ci}
            case 'realistic',      src = rfiRealistic;
            case 'inmarsatNarrow', src = [rfiRealistic(~isInm), inmNarrow];
            case 'noiseOnly',      src = struct([]);
        end
        if Lc ~= L && ~isempty(src)              % steady part: beam passages outside the run
            for k = 1:numel(src)
                if isfinite(src(k).scanPeriod) && src(k).scanPeriod > 0 && ...
                        abs(mod(src(k).beamTime - Lc/2 + src(k).scanPeriod/2, src(k).scanPeriod) ...
                            - src(k).scanPeriod/2) < Lc
                    src(k).beamTime = Lc + 5;
                end
            end
        end
        info_rx = addNoiseAndRFI(info_zero.file, fullfile(pDir, "p_rx.dat"), info_zero.actualFsOut, ...
            'Band', [fLow fHigh], 'NoiseStd', 1, 'RFI', src, 'Seed', noiseSeed, 'Verbose', false);
        info_IQ = applyIQmodulation(info_rx.file, fullfile(pDir, "p_rx_IQ.dat"), info_rx.actualFsOut, ...
            fs, fLO, 'FilterOrder', filterOrder, 'Band', [fLow fHigh], 'Verbose', false);
        delete(info_rx.file);
        out = processIQ(info_IQ, rx, ephem, procP, pFiles, 'Verbose', false);
        delete(info_IQ.file);
        label = sprintf('%s, %g s, zero sky, excision (runRFIProfile)', caseNames{ci}, Lc);
        prof = rfiProfile(out.info_det, out.info_w, 'Label', label, 'Lmax', 30);
        prof.sources = src; prof.noiseSeed = noiseSeed;
        save(fOut, 'prof');
        fprintf('runRFIProfile: %s -> %s (%.0f s)\n', label, fOut, toc(tCase));
    end
    delete(info_zero.file);
end
% free the disk: the large intermediate files of the last run
delete(fullfile(pDir, "p_*.dat"));
fprintf('runRFIProfile: done in %.1f min\n', toc(tAll)/60);
