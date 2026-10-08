function info = blankChannels(info_chan, mask, outBase, opts)
%BLANKCHANNELS  Zero the masked samples of the channel IQ files (RFI excision, B1).
%{
RFI excision blanks (sets to zero) the channel voltage samples hit by RFI,
BEFORE per-channel coherent dedispersion: there a radar pulse is a few
microseconds long in one or two channels; after dedispersion it would be
smeared over the in-channel sweep (~60 us per channel at DM 5). Blanking
removes signal and noise alike (weight 0, nothing subtracted): at -54 dB
one radar pulse carries the energy of ~80,000 pulsar pulses, so even
99.9 % suppression by subtraction would leave more RFI than pulsar.

This function only applies a given mask (from the RFI detector, or a
test). It writes blanked copies of the channel files (the raw files stay,
for comparison) and returns an info struct that dedisperseChannels accepts
unchanged. The same mask then goes to blankingWeights, which turns it into
the exact data weights of the detected time bins (A3a).

  info = blankChannels(info_chan, mask, outBase)

Inputs:
  info_chan  info of channelizeIQ (chanFiles, nChan, N, file).
  mask       [M x 3] blanked intervals [channel, firstSample, lastSample]
             (1-based, inclusive, on the channel sample grid; the format of
             blankingWeights). Clipped to the file; overlapping and
             unsorted rows allowed; rows with first > last (after clipping)
             are ignored. [] = nothing blanked (plain copies).
  outBase    path without extension; channel j -> <outBase>_ch001.dat ...
             (cf32), info -> <outBase>_info.mat (loadInfo(outBase) works).
             Must differ from info_chan.file: the raw files are never
             overwritten.

Name-value options:
  'BlockSize'  samples per channel per read/write block (default 2^20);
               the result does not depend on it.
  'SaveInfo'   save info to <outBase>_info.mat (default true).
  'Verbose'    print a summary (default true).

Output info: info_chan with file = outBase and chanFiles = the blanked
  files; every other field unchanged (same sample grid, prototype filter,
  channel frequencies), so dedisperseChannels and blankingWeights see the
  same channelizer. Added: rawFile (info_chan.file), mask (normalized:
  per channel sorted, merged and clipped intervals, [channel, first, last];
  blanks the same samples as the input mask), nMaskRows (rows given),
  blankedSamples and blankedFraction (per channel), blockSize, elapsed
  (of this stage).
%}

arguments
    info_chan         struct
    mask              double
    outBase                 {mustBeTextScalar}
    opts.BlockSize    (1,1) double {mustBeInteger, mustBePositive} = 2^20
    opts.SaveInfo     (1,1) logical = true
    opts.Verbose      (1,1) logical = true
end

tStart  = tic;
outBase = char(outBase);
nChan   = info_chan.nChan;
Nc      = info_chan.N;                              % complex samples per channel
if strcmp(outBase, char(info_chan.file))
    error('blankChannels:sameFile', 'outBase must differ from info_chan.file (raw files).');
end
if isempty(mask), mask = zeros(0, 3); end
if size(mask, 2) ~= 3 || any(mask(:, 1) < 1 | mask(:, 1) > nChan | mask(:, 1) ~= round(mask(:, 1)))
    error('blankChannels:mask', 'mask must be [channel, first, last] rows with valid channels.');
end
if any(mask(:, 2:3) ~= round(mask(:, 2:3)), 'all')
    error('blankChannels:mask', 'mask first/last samples must be whole numbers.');
end

outDir = fileparts(outBase);
if ~isempty(outDir) && ~isfolder(outDir), mkdir(outDir); end

chanFiles = strings(nChan, 1);
ivAll     = cell(nChan, 1);
blankedSamples = zeros(nChan, 1);
blk = opts.BlockSize;
for j = 1:nChan
    chanFiles(j) = sprintf('%s_ch%03d.dat', outBase, j);
    if strcmp(chanFiles(j), info_chan.chanFiles(j))
        error('blankChannels:sameFile', 'Channel %d: output would overwrite the raw file.', j);
    end
    iv = mergeIntervals(mask(mask(:, 1) == j, 2:3), Nc);
    ivAll{j} = [j * ones(size(iv, 1), 1), iv];
    blankedSamples(j) = sum(iv(:, 2) - iv(:, 1) + 1);

    [fidIn, msg] = fopen(info_chan.chanFiles(j), 'r', 'ieee-le');
    if fidIn == -1
        error('blankChannels:openIn', 'Could not open "%s": %s', info_chan.chanFiles(j), msg);
    end
    cleanIn = onCleanup(@() fclose(fidIn));
    [fidOut, msg] = fopen(chanFiles(j), 'w', 'ieee-le');
    if fidOut == -1
        error('blankChannels:openOut', 'Could not open "%s": %s', chanFiles(j), msg);
    end
    cleanOut = onCleanup(@() fclose(fidOut));

    for s0 = 1:blk:Nc                               % block: samples s0 .. s1 (1-based)
        s1 = min(s0 + blk - 1, Nc);
        n  = s1 - s0 + 1;
        x  = fread(fidIn, [2 n], 'single=>single'); % rows I, Q
        if size(x, 2) < n
            error('blankChannels:read', 'Unexpected end of "%s".', info_chan.chanFiles(j));
        end
        sel = iv(:, 2) >= s0 & iv(:, 1) <= s1;      % intervals reaching into this block
        if any(sel)
            a = max(iv(sel, 1), s0) - s0 + 1;       % block-relative first, last
            b = min(iv(sel, 2), s1) - s0 + 1;
            d = accumarray([a; b + 1], [ones(numel(a), 1); -ones(numel(b), 1)], [n + 1, 1]);
            x(:, cumsum(d(1:n)) > 0) = 0;
        end
        c = fwrite(fidOut, x, 'single');
        if c ~= 2*n
            error('blankChannels:write', 'Channel %d: wrote %d of %d values (disk full?).', j, c, 2*n);
        end
    end
    clear cleanIn cleanOut                          % close both files
    if opts.Verbose && (j == nChan || mod(j, 32) == 0)
        fprintf('  blankChannels: %d/%d channels\n', j, nChan);
    end
end

info = info_chan;
if isfield(info, 'infoFile'), info = rmfield(info, 'infoFile'); end
info.file            = string(outBase);
info.chanFiles       = chanFiles;
info.rawFile         = info_chan.file;
info.mask            = vertcat(ivAll{:});           % normalized [channel, first, last]
info.nMaskRows       = size(mask, 1);
info.blankedSamples  = blankedSamples;
info.blankedFraction = blankedSamples / Nc;
info.blockSize       = blk;
info.elapsed         = toc(tStart);

if opts.SaveInfo
    infoFile = [outBase '_info.mat'];
    save(infoFile, 'info');
    info.infoFile = infoFile;
end
if opts.Verbose
    fprintf(['blankChannels: %d mask rows -> %d intervals in %d of %d channels, %.4f %% of ' ...
             'samples blanked; wrote %s_ch###.dat in %.1f s\n'], info.nMaskRows, ...
        size(info.mask, 1), nnz(blankedSamples), nChan, 100*mean(info.blankedFraction), ...
        outBase, info.elapsed);
end
end


% =====================================================================================
function iv = mergeIntervals(r, Nc)
%MERGEINTERVALS  [first, last] rows -> clipped to 1..Nc, sorted, merged where they
% overlap or touch; rows empty after clipping (first > last) dropped.
r  = [max(r(:, 1), 1), min(r(:, 2), Nc)];
r  = r(r(:, 1) <= r(:, 2), :);
if isempty(r), iv = zeros(0, 2); return; end
r  = sortrows(r);
endMax = cummax(r(:, 2));                           % furthest end so far
newRun = [true; r(2:end, 1) > endMax(1:end-1) + 1]; % gap before this row -> new interval
lastRow = [find(newRun(2:end)); size(r, 1)];        % last row of each interval
iv = [r(newRun, 1), endMax(lastRow)];
end
