function info = detectPower(inFile, outFile, fs, fLO, f_out, opts)
%DETECTPOWER  Square-law detection of streamed complex baseband data.
%{
Converts complex baseband (I/Q) samples into power, |z|^2, averaged over
time bins of length binLen = round(fs/f_out) samples, optionally split
into nChan frequency sub-bands. This is where the data rate collapses
(e.g. 1e9 complex samples/s -> 1e4..1e6 power values/s), so everything
downstream (normalization, folding, matched filtering) is cheap.

  info = detectPower(inFile, outFile, fs, fLO, f_out)                 % total power
  info = detectPower(..., 'NChan', 16, 'Band', [1.2e9 1.6e9])        % sub-bands

Inputs:
  inFile  interleaved float32 I,Q (cf32), little-endian (e.g. dedispersed).
  outFile float32, little-endian power file (layout below).
  fs      [Hz] complex sample rate of inFile.
  fLO     [Hz] LO frequency (RF = fLO + baseband); used for channel edges.
  f_out   [Hz] requested output rate (time bins per second). The achieved
          rate fs/binLen is in info.actualFout.

Name-value options:
  'NChan'           number of sub-bands (default 1).
  'Band'            [fLow fHigh] RF range to split into NChan equal channels.
                    Required if NChan > 1. With NChan = 1 and no Band, the
                    total baseband power is computed directly (fastest path;
                    after applyInverseDispersion the data is already
                    band-limited to [fLow fHigh]).
  'FineBinsPerChan' FFT bins per channel for the channelizer (default 8);
                    sets the channelizer frame length.
  'T0'              [s] time of input sample 1 (default 0).
  'FullySupported'  [first last] 1-based input sample range that is fully
                    valid (info_dedisp.fullySupported); converted to a bin
                    range in info.fullySupportedBins.
  'BlockSamples'    target input samples per block (default 2^22).
  'SaveInfo'        save info to <outFile>_info.mat (default true).
  'Verbose'         print a summary (default true).

Output file layout:
  float32, [nChan x nBins] column-major: all channels of bin 1, then all
  channels of bin 2, ... With nChan = 1 this is a plain power time series.
  Value = mean power per input sample in that bin (and channel), so with
  nChan > 1 the channels add up to the in-band part of mean |z|^2.

Timing:
  Bin k (1-based) averages input samples (k-1)*binLen+1 .. k*binLen; its
  time (centroid of those samples) is
     t_k = T0 + ((k-1)*binLen + (binLen-1)/2) / fs
       = info.binTime0 + (k-1) * info.binDt.
  The boxcar average is symmetric, so it does not bias arrival times.
  A trailing partial bin is dropped (info.samplesDropped).

Method:
  NChan = 1, no Band: p = I.^2 + Q.^2 (no sqrt, no complex arrays), summed
    per bin in double precision.
  Channelized: non-overlapping rectangular frames of frameLen samples
    (frameLen divides binLen and gives >= FineBinsPerChan FFT bins per
    channel), FFT, |Z|^2 per bin, summed into channels by one small
    matrix multiply, then summed over the frames in each time bin.
    Normalization follows Parseval: sum|z|^2 = (1/L) sum|Z|^2 per frame,
    so the channels of a band that covers the whole baseband add up
    exactly to the total power. Leakage between neighbouring channels is
    limited to their edge bins (~1/FineBinsPerChan of each channel).
%}

arguments
    inFile                     {mustBeTextScalar}
    outFile                    {mustBeTextScalar}
    fs                   (1,1) double {mustBePositive, mustBeFinite}
    fLO                  (1,1) double {mustBePositive, mustBeFinite}
    f_out                (1,1) double {mustBePositive, mustBeFinite}
    opts.NChan           (1,1) double {mustBeInteger, mustBePositive} = 1
    opts.Band                  double = []
    opts.FineBinsPerChan (1,1) double {mustBeInteger, mustBePositive} = 8
    opts.T0              (1,1) double {mustBeFinite} = 0
    opts.FullySupported        double = []
    opts.BlockSamples    (1,1) double {mustBePositive} = 2^22
    opts.SaveInfo        (1,1) logical = true
    opts.Verbose         (1,1) logical = true
end

tStart  = tic;
inFile  = char(inFile);
outFile = char(outFile);
if strcmp(inFile, outFile)
    error('detectPower:sameFile', 'inFile and outFile must differ.');
end
nChan = opts.NChan;

% ---- Time binning ----------------------------------------------------------------
binLen = max(1, round(fs / f_out));
actualFout = fs / binLen;
if abs(actualFout - f_out) > 1e-6 * f_out
    warning('detectPower:rate', ...
        'f_out = %.6g Hz is not fs/integer; using %.6g Hz (binLen = %d).', ...
        f_out, actualFout, binLen);
end

% ---- Channelizer setup -------------------------------------------------------------
fastPath = (nChan == 1) && isempty(opts.Band);
if fastPath
    frameLen  = binLen;
    A         = [];
    chanEdges = [fLO - fs/2, fLO + fs/2];
    chanFreqs = fLO;
    chanBins  = binLen;
else
    if isempty(opts.Band) || numel(opts.Band) ~= 2
        error('detectPower:band', 'Band = [fLow fHigh] is required when NChan > 1.');
    end
    fLow = opts.Band(1); fHigh = opts.Band(2);
    if fLow >= fHigh || fLow < fLO - fs/2 || fHigh > fLO + fs/2
        error('detectPower:band', ...
            'Band must satisfy fLO - fs/2 <= fLow < fHigh <= fLO + fs/2.');
    end
    chanBW  = (fHigh - fLow) / nChan;
    targetL = ceil(opts.FineBinsPerChan * fs / chanBW);
    frameLen = pickFrameLen(binLen, targetL);
    if frameLen < targetL
        warning('detectPower:resolution', ...
            ['Time bin (%d samples) is shorter than the frame needed for %d ' ...
             'FFT bins per channel (%d); using %d bins per channel. Lower f_out ' ...
             'or NChan for cleaner channels.'], ...
            binLen, opts.FineBinsPerChan, targetL, floor(frameLen*chanBW/fs));
    end
    [A, chanBins] = channelMatrix(frameLen, nChan, fs, fLO, fLow, fHigh, binLen);
    if any(chanBins == 0)
        error('detectPower:emptyChannel', ...
            'Some channels contain no FFT bins; lower NChan or f_out.');
    end
    chanEdges = fLow + (0:nChan) * chanBW;
    chanFreqs = fLow + ((1:nChan) - 0.5) * chanBW;
end
framesPerBin = binLen / frameLen;

% ---- Input size -------------------------------------------------------------------
fInfo = dir(inFile);
if isempty(fInfo)
    error('detectPower:notFound', 'inFile "%s" not found.', inFile);
end
Nx = floor(fInfo.bytes / 8);
nBinsTotal = floor(Nx / binLen);
if nBinsTotal == 0
    error('detectPower:short', 'File shorter than one time bin (%d samples).', binLen);
end
binsPerBlock = max(1, round(opts.BlockSamples / binLen));
nBlocks = ceil(nBinsTotal / binsPerBlock);

if opts.Verbose
    if fastPath
        fprintf(['detectPower: %d samples @ %.4g Hz -> %d bins @ %.6g Hz ' ...
                 '(binLen %d), total power\n'], Nx, fs, nBinsTotal, actualFout, binLen);
    else
        fprintf(['detectPower: %d samples @ %.4g Hz -> %d bins @ %.6g Hz ' ...
                 '(binLen %d), %d channels of %.4g MHz over %.4f-%.4f GHz, ' ...
                 'frame %d (%.3g MHz FFT bins)\n'], Nx, fs, nBinsTotal, actualFout, ...
            binLen, nChan, (chanEdges(2)-chanEdges(1))/1e6, chanEdges(1)/1e9, ...
            chanEdges(end)/1e9, frameLen, fs/frameLen/1e6);
    end
end

% ---- Files (explicit little-endian) -------------------------------------------------
[fidIn, msg] = fopen(inFile, 'r', 'ieee-le');
if fidIn == -1
    error('detectPower:openIn', 'Could not open "%s": %s', inFile, msg);
end
cleanupIn = onCleanup(@() fclose(fidIn)); %#ok<NASGU>
outDir = fileparts(outFile);
if ~isempty(outDir) && ~isfolder(outDir)
    mkdir(outDir);
end
[fidOut, msg] = fopen(outFile, 'w', 'ieee-le');
if fidOut == -1
    error('detectPower:openOut', 'Could not open "%s": %s', outFile, msg);
end
cleanupOut = onCleanup(@() fclose(fidOut)); %#ok<NASGU>

% ---- Stream -------------------------------------------------------------------------------
chanSum = zeros(nChan, 1);
binsDone = 0;
for b = 1:nBlocks
    nb = min(binsPerBlock, nBinsTotal - binsDone);
    nS = nb * binLen;
    raw = fread(fidIn, [2 nS], 'single=>single');
    if size(raw, 2) < nS
        error('detectPower:read', 'Unexpected end of file in block %d.', b);
    end

    if fastPath
        p   = raw(1, :).^2 + raw(2, :).^2;
        out = sum(reshape(p, binLen, nb), 1, 'double') / binLen;      % 1 x nb
    else
        Z   = fft(reshape(complex(raw(1, :), raw(2, :)), frameLen, []));
        P   = real(Z).^2 + imag(Z).^2;                                  % L x nFrames
        C   = A * P;                                                    % nChan x nFrames
        out = reshape(sum(reshape(C, nChan, framesPerBin, nb), 2, 'double'), nChan, nb);
    end
    clear raw p Z P C

    c = fwrite(fidOut, single(out), 'single');
    if c ~= numel(out)
        error('detectPower:write', 'Wrote %d of %d values at block %d.', c, numel(out), b);
    end
    chanSum  = chanSum + sum(out, 2);
    binsDone = binsDone + nb;
end

% ---- Info --------------------------------------------------------------------------------------
info = struct();
info.file           = outFile;
info.inFile         = inFile;
info.format         = 'float32 [nChan x nBins], column-major (all channels per bin)';
info.precision      = 'single';
info.byteOrder      = 'ieee-le';
info.isComplex      = false;
info.detector       = 'square-law |z|^2, mean per input sample per bin';
info.nChan          = nChan;
info.N              = binsDone;                  % time bins
info.fs             = actualFout;
info.actualFout     = actualFout;
info.fsIn           = fs;
info.binLen         = binLen;
info.binDt          = binLen / fs;
info.binTime0       = opts.T0 + (binLen - 1) / (2*fs);   % centroid time of bin 1
info.T0             = opts.T0;
info.fLO            = fLO;
info.chanFreqs      = chanFreqs;
info.chanEdges      = chanEdges;
info.chanFftBins    = chanBins;
info.frameLen       = frameLen;
info.chanMeanPower  = chanSum.' / binsDone;
info.samplesUsed    = binsDone * binLen;
info.samplesDropped = Nx - binsDone * binLen;
if isempty(opts.FullySupported)
    info.fullySupportedBins = [1, binsDone];
else
    s = opts.FullySupported;
    info.fullySupportedBins = [ceil((s(1) - 1) / binLen) + 1, min(binsDone, floor(s(2) / binLen))];
end
info.elapsed        = toc(tStart);

if opts.SaveInfo
    [d, name] = fileparts(outFile);
    infoFile = fullfile(d, [name '_info.mat']);
    save(infoFile, 'info');
    info.infoFile = infoFile;
end
if opts.Verbose
    fprintf('detectPower: wrote %d bins x %d channel(s) to %s in %.1f s\n', ...
        binsDone, nChan, outFile, info.elapsed);
end
end


% ======================================================================================
function L = pickFrameLen(binLen, targetL)
%PICKFRAMELEN  Smallest divisor of binLen >= targetL, preferring FFT-friendly
% sizes (largest prime factor <= 7). Falls back to binLen itself.
d1 = 1:floor(sqrt(binLen));
d1 = d1(mod(binLen, d1) == 0);
divs = unique([d1, binLen ./ d1]);
cand = divs(divs >= targetL);
if isempty(cand)
    L = binLen;
    return
end
smooth = arrayfun(@(x) max(factor(x)) <= 7, cand);
if any(smooth)
    L = cand(find(smooth, 1));
else
    L = cand(1);
end
end


function [A, chanBins] = channelMatrix(L, nChan, fs, fLO, fLow, fHigh, binLen)
%CHANNELMATRIX  nChan x L matrix summing FFT-bin powers into channels,
% including the Parseval (1/L) and per-bin-mean (1/binLen) scaling.
k  = 0:L-1;
ks = k;
neg = k >= ceil(L/2);
ks(neg) = k(neg) - L;                          % signed bin index (FFT order)
f  = fLO + ks * (fs / L);                      % RF frequency of each bin
bw = (fHigh - fLow) / nChan;
ch = floor((f - fLow) / bw) + 1;
ok = f >= fLow & f < fHigh & ch >= 1 & ch <= nChan;
A  = zeros(nChan, L, 'single');
A(sub2ind([nChan L], ch(ok), k(ok) + 1)) = single(1 / (L * binLen));
chanBins = accumarray(ch(ok).', 1, [nChan 1]).';
end