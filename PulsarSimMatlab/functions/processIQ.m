function out = processIQ(info_IQ, rx, ephem, proc, files, opts)
%PROCESSIQ  The processing chain: IQ samples in, fold out (observer side only).
%{
What an observer does with the receiver's IQ data up to the fold: front end
(full band or channelized, with RFI excision), detection, fold. The TOAs
follow in a separate step, foldToTOA (out.Bnoise and out.weighting are its
inputs), so a fold from elsewhere (block C's fast simulator) goes through the
same TOA code, and main.m can redo the TOAs alone. The input is the IQ
file the receiver writes (in simulation applyIQmodulation; later the hardware's
ADC/DDC output, a real recording), so the simulation stays outside and the
ground-truth rule holds by the signature: this function only gets the IQ file,
receiver settings, the ephemeris and processing settings. It is also the MATLAB
reference ("golden model") for a later FPGA/GPU implementation: the same stage
boundaries and intermediate files.

  out = processIQ(info_IQ, rx, ephem, proc, files)
  out = processIQ(..., 'Reuse', true)          % saved files and fold, checked
  [toa, detection] = foldToTOA(out.fold, out.info_fold, ephem.profileFWHM, ...
                               out.Bnoise, 'Weighting', out.weighting)

Inputs (pipelineParams builds rx, proc, files from its variables):
  info_IQ  the IQ file's info (applyIQmodulation or loadInfo): file,
           actualFsOut, fLO, t0 - receiver settings only
  rx       receiver: fLow, fHigh [Hz] (observed RF band)
  ephem    ephemeris: f0 [Hz], F1 [Hz/s], TRef [s], DM [pc cm^-3]
           (profileFWHM is for foldToTOA)
  proc     processing settings:
             frontEnd       'fullband' or 'channels'
             refFreq        [Hz] dedispersion reference = TOA reference frequency
             f_out          [Hz] detected-power bin rate
             nBin, subintPeriods   fold: phase bins, turns per sub-int
             chanWidth      [Hz] channel width (channels)
             weighting      'equal' or 'optimal' (channels; full band: 'equal')
             excision       RFI excision before dedispersion (channels)
             excisionArgs   name-value options for detectRFI
             periodicMask   with excision: periodicRFI's windows too
             periodicArgs   name-value options for periodicRFI
  files    output files: full band dedisp, envelope, fold; channels chanIQ,
           chanBlank, chanDedisp, chanPower (no excision), chanDedispB,
           chanPowerB (excision: own files, so the unblanked ones stay),
           chanWeight, foldChan

Name-value options:
  'Reuse'    false (default): run every step. true: reuse the saved
             intermediate files and the saved fold (main.m: runStage.process
             = false), after checks that they belong to this IQ file and
             these settings (warnings processIQ:staleProcessing,
             processIQ:staleFold).
  'Verbose'  passed to every stage (default true)

Steps (main.m before 19 Oct, same calls):
  1. Front end, IQ -> detected power per time bin
     full band: applyInverseDispersion (coherent dedispersion of the band),
                detectPower (|y|^2 in bins of f_out)
     channels:  channelizeIQ (filterbank); with excision detectRFI (mask)
                [+ periodicRFI] -> blankChannels (zeroed copies);
                dedisperseChannels (per channel, aligned at refFreq);
                detectChannels (|y|^2 per channel, exact noise V, X); with
                excision blankingWeights (data weights W, V, X per bin)
  2. foldProfile with the ephemeris phase model (channels: NoiseCoeffs, or
     DataWeights with excision)
  3. the inputs of foldToTOA: Bnoise for the radiometer noise model (full
     band noiseBandwidth of the dedispersion taper; channels BnoiseChan of
     one channel) and the weighting (full band 'equal')

Output out (struct): frontEnd; fold, info_fold; Bnoise, weighting (for
  foldToTOA); the intermediate infos for checks before and after folding:
  info_det, info_dedisp (full band) or info_chan, info_dc, info_w (channels;
  [] where not made or reused); with excision also rfiMask (detectRFI's
  rows [+ periodicRFI's]), info_rfi, info_per (periodicRFI; [] when off) and
  info_chanD (blankChannels: mask, blankedFraction). Every intermediate file
  stays on disk.

(19 Oct 2026: moved out of main.m - the §10 item 7 redesign: processing only,
cut at IQ, both front ends.)
%}

arguments
    info_IQ struct
    rx      struct
    ephem   struct
    proc    struct
    files   struct
    opts.Reuse   (1,1) logical = false
    opts.Verbose (1,1) logical = true
end

needFields(rx, {'fLow', 'fHigh'}, 'rx');
needFields(ephem, {'f0', 'F1', 'TRef', 'DM'}, 'ephem');
needFields(proc, {'frontEnd', 'refFreq', 'f_out', 'nBin', 'subintPeriods', 'chanWidth', ...
    'weighting', 'excision', 'excisionArgs', 'periodicMask', 'periodicArgs'}, 'proc');
fullBand = strcmp(proc.frontEnd, 'fullband');
if ~fullBand && ~strcmp(proc.frontEnd, 'channels')
    error('processIQ:frontEnd', 'proc.frontEnd must be ''fullband'' or ''channels''.');
end
if fullBand
    needFields(files, {'dedisp', 'envelope', 'fold'}, 'files');
else
    needFields(files, {'chanIQ', 'chanBlank', 'chanDedisp', 'chanPower', 'chanDedispB', ...
        'chanPowerB', 'chanWeight', 'foldChan'}, 'files');
end
reuse = opts.Reuse;
vb    = {'Verbose', opts.Verbose};
excision = ~fullBand && proc.excision;
out = struct('frontEnd', proc.frontEnd, 'info_det', [], 'info_dedisp', [], 'info_chan', [], ...
    'info_dc', [], 'info_w', [], 'info_fold', [], 'fold', [], 'Bnoise', [], 'weighting', '', ...
    'rfiMask', [], 'info_rfi', [], 'info_per', [], 'info_chanD', []);

% ---- 1. Front end: IQ -> detected power ---------------------------------------------------
if fullBand
    if ~reuse
        % Coherent dedispersion
        info_dedisp = applyInverseDispersion(info_IQ.file, files.dedisp, info_IQ.actualFsOut, ...
            info_IQ.fLO, ephem.DM, rx.fLow, rx.fHigh, 'RefFreq', proc.refFreq, vb{:});
        % Square-law detection
        info_det = detectPower(info_dedisp.file, files.envelope, info_dedisp.fs, info_dedisp.fLO, ...
            proc.f_out, 'FullySupported', info_dedisp.fullySupported, vb{:});
    else
        info_dedisp = loadInfo(files.dedisp);
        info_det    = loadInfo(files.envelope);
    end
    out.info_dedisp = info_dedisp;
else
    if excision                 % own files: the unblanked ones stay for comparison
        dcFile = files.chanDedispB; pwFile = files.chanPowerB;
    else
        dcFile = files.chanDedisp;  pwFile = files.chanPower;
    end
    info_w = [];
    if ~reuse
        % Filterbank: 128 channels of 3.125 MHz (oversampled 4/3)
        info_chan = channelizeIQ(info_IQ.file, files.chanIQ, info_IQ.actualFsOut, info_IQ.fLO, ...
            rx.fLow, rx.fHigh, 'ChanWidth', proc.chanWidth, 'T0', info_IQ.t0, vb{:});
        % RFI excision BEFORE dedispersion (a radar pulse is a few samples here):
        % find the RFI, zero those samples in copies of the channel files
        if excision
            [rfiMask, info_rfi] = detectRFI(info_chan, proc.excisionArgs{:}, vb{:});
            if proc.periodicMask        % + every predicted pulse of periodic RFI (radars)
                [perRows, out.info_per] = periodicRFI(info_chan, info_rfi, proc.periodicArgs{:}, vb{:});
                rfiMask = [rfiMask; perRows];
            end
            info_chanD = blankChannels(info_chan, rfiMask, files.chanBlank, vb{:});
            out.rfiMask = rfiMask; out.info_rfi = info_rfi; out.info_chanD = info_chanD;
        else
            info_chanD = info_chan;
        end
        % Coherent dedispersion per channel, all channels aligned at refFreq
        info_dc = dedisperseChannels(info_chanD, dcFile, ephem.DM, 'RefFreq', proc.refFreq, vb{:});
        % Square-law detection per channel; bins of a whole number of channel samples
        info_det = detectChannels(info_dc, pwFile, info_dc.fs / round(info_dc.fs / proc.f_out), vb{:});
        % Dedispersion spreads each blank over the filter's reach: the exact valid
        % fraction W and noise V, X of every detected bin, from the mask and each filter
        if excision
            info_w = blankingWeights(info_chan, info_dc, info_det, info_chanD.mask, ...
                files.chanWeight, vb{:});
        end
        out.info_chan = info_chan;
    else
        info_dc  = loadInfo(dcFile);
        info_det = loadInfo(pwFile);
        if excision, info_w = loadInfo(files.chanWeight); end
    end
    out.info_dc = info_dc;
    out.info_w  = info_w;
end
out.info_det = info_det;
if reuse                        % reused files: made from this IQ file, with these settings?
    if fullBand
        chain = {info_IQ.file, info_IQ, ''; files.dedisp, info_dedisp, 'inFile'; ...
                 files.envelope, info_det, 'inFile'};
        checkProcessing(chain, info_dedisp, info_det, proc.f_out, ephem.DM, proc.refFreq, ...
            rx.fLow, rx.fHigh, []);
    else
        chain = {info_IQ.file, info_IQ, ''; files.chanIQ, [], 'inFile'};
        if excision, chain(end+1, :) = {files.chanBlank, [], 'rawFile'}; end
        chain = [chain; {dcFile, info_dc, 'inFile'; pwFile, info_det, 'inFile'}];
        if excision, chain(end+1, :) = {files.chanWeight, info_w, 'powerFile'}; end
        checkProcessing(chain, info_dc, info_det, proc.f_out, ephem.DM, proc.refFreq, ...
            rx.fLow, rx.fHigh, proc.chanWidth);
    end
end

% ---- 2. Fold -------------------------------------------------------------------------------
if fullBand, foldFile = files.fold; else, foldFile = files.foldChan; end
if ~reuse
    foldArgs = {'F1', ephem.F1, 'TRef', ephem.TRef, 'NBin', proc.nBin, ...
                'SubintPeriods', proc.subintPeriods};
    if ~fullBand
        % narrow channels: neighbouring time bins correlate -> exact phase-bin noise;
        % with excision the data weights carry it per bin (and the valid fraction W)
        if excision
            foldArgs = [foldArgs, {'DataWeights', info_w}];
        else
            foldArgs = [foldArgs, {'NoiseCoeffs', [info_det.noise.V, info_det.noise.X]}];
        end
    end
    [info_fold, fold] = foldProfile(info_det, foldFile, ephem.f0, foldArgs{:}, vb{:});
else
    S = load(foldFile, 'fold', 'info');
    fold = S.fold; info_fold = S.info;                  % reuse the saved fold
    wFile = '';                                         % data-weight file it must have used
    if excision, wFile = files.chanWeight; end
    checkFold(info_fold, foldFile, info_det, ephem, proc.nBin, proc.subintPeriods, wFile);
end
out.info_fold = info_fold;
out.fold      = fold;

% ---- 3. Inputs of foldToTOA -------------------------------------------------------------------
% Noise-equivalent bandwidth for the radiometer noise model, from the observer's own
% dedispersion taper (receiver noise only passes this one). When the pulsar dominates
% (high SNR, in-pulse) the true value is ~0.5 % lower, since the simulated signal also
% passed the dispersion-stage taper. Channel path: Bnoise of ONE channel (its own
% dedispersion taper), and the channels combined as set by proc.weighting.
if fullBand
    out.Bnoise    = noiseBandwidth(rx.fLow, rx.fHigh, info_dedisp.edgeWidth);
    out.weighting = 'equal';
else
    out.Bnoise    = info_det.noise.Bnoise;
    out.weighting = proc.weighting;
end
end


% ===========================================================================================
%  Local functions (reuse checks: observer settings vs the saved info files; from main.m)
% ===========================================================================================
function needFields(s, names, label)
%NEEDFIELDS  Error when struct s lacks one of the fields it must have.
miss = names(~isfield(s, names));
if ~isempty(miss)
    error('processIQ:fields', '%s is missing field(s): %s.', label, strjoin(miss, ', '));
end
end

function checkProcessing(chain, info_dd, info_det, f_out, DM, refFreq, fLow, fHigh, chanWidth)
%CHECKPROCESSING  Warn when reused processing files do not belong to this run.
% chain: rows {file, info (or [] = load it), link field}, from the IQ file to the
% power (and weight) file. Each file must have been made from the one before it
% (its link field) and after it (info file at least as new: a rerun receiver makes
% the older products stale). info_dd: the dedispersion info (DM, refFreq, band,
% chanWidth); chanWidth [] on the full-band path.
d = {};
tPrev = infoTime(chain{1, 1});
for k = 2:size(chain, 1)
    [f, info, link] = chain{k, :};
    [~, n] = fileparts(f);
    if isempty(info)
        if isnan(infoTime(f)), d{end+1} = sprintf('%s missing', n); continue; end %#ok<AGROW>
        info = loadInfo(f);
    end
    if ~strcmp(string(info.(link)), string(chain{k-1, 1}))
        d{end+1} = sprintf('%s not made from %s', n, chain{k-1, 1}); %#ok<AGROW>
    end
    t = infoTime(f);
    if t < tPrev, d{end+1} = sprintf('%s older than its input', n); end %#ok<AGROW>
    tPrev = t;
    if strcmp(link, 'powerFile') && (info.N ~= info_det.N || info.Lmax ~= info_det.noise.Lmax)
        d{end+1} = sprintf('%s: other time bins or lags than the power file', n); %#ok<AGROW>
    end
end
if info_det.binLen ~= round(info_det.fsIn / f_out)
    d{end+1} = sprintf('f_out (power file: %d samples per bin, this run: %d)', ...
        info_det.binLen, round(info_det.fsIn / f_out));
end
if abs(info_dd.DM - DM) > 1e-12,           d{end+1} = 'DM';      end
if abs(info_dd.refFreq - refFreq) > 1e-3,  d{end+1} = 'refFreq'; end
if abs(info_dd.fLow - fLow) > 1e-3 || abs(info_dd.fHigh - fHigh) > 1e-3, d{end+1} = 'band'; end
if ~isempty(chanWidth) && abs(info_dd.chanWidth - chanWidth) > 1e-3, d{end+1} = 'chanWidth'; end
if ~isempty(d)
    warning('processIQ:staleProcessing', ['Reused processing files do not match this run ' ...
        '(%s). Set runStage.process = true.'], strjoin(d, '; '));
end
end

function checkFold(info_fold, foldFile, info_det, ephem, nBin, subintPeriods, wFile)
%CHECKFOLD  Warn when a reused fold was not made from this power file and these
% settings (phase model, bins, turns per sub-int, data weights; wFile '' = none).
d = {};
if ~strcmp(string(info_fold.sourceFile), string(info_det.file)), d{end+1} = 'power file'; end
f = dir(foldFile);
if f.datenum < infoTime(info_det.file), d{end+1} = 'older than the power file'; end
if info_fold.NBin ~= nBin,                   d{end+1} = 'nBin';          end
if info_fold.subintPeriods ~= subintPeriods, d{end+1} = 'subintPeriods'; end
if any(abs([info_fold.f0, info_fold.F1, info_fold.TRef] - [ephem.f0, ephem.F1, ephem.TRef]) > ...
        1e-12 * max(1, abs([ephem.f0, ephem.F1, ephem.TRef])))
    d{end+1} = 'ephemeris';
end
if ~strcmp(string(info_fold.dataWeights), string(wFile)), d{end+1} = 'data weights'; end
if ~isempty(d)
    warning('processIQ:staleFold', ['The reused fold does not match this run (%s). ' ...
        'Set runStage.process = true.'], strjoin(d, ', '));
end
end

function t = infoTime(dataFile)
%INFOTIME  Modification time (datenum) of <dataFile>_info.mat, NaN when missing: the
% info file is saved after the data, so it marks when the stage finished.
[d, name] = fileparts(dataFile);
f = dir(fullfile(d, name + "_info.mat"));
t = NaN;
if ~isempty(f), t = f.datenum; end
end
