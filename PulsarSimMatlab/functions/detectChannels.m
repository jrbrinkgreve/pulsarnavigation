function info = detectChannels(info_dc, outFile, f_out, opts)
%DETECTCHANNELS  Square-law detection of every dedispersed channel into one power file.
%{
Third stage of the channelized front end. Each dedispersed channel file is
detected with the validated detectPower (|y|^2 averaged over time bins), and
the channels are interleaved into ONE float32 file with the existing
detected-power layout [nChan x nBins] (all channels of bin 1, then bin 2,
...), so foldProfile can read it unchanged. All channels share one time
grid (dedisperseChannels aligned them at the reference frequency).

The info struct also carries what the noise model needs per channel:
  - BnoiseChan: noise bandwidth of one channel (its band-edge taper);
  - noise: exact statistics of detected time bins in a narrow channel
    (powerCovariance): variance V and covariances X(L) with the bins L =
    1..Lmax later, relative to the radiometer value m^2/(Bnoise*binDt). In a
    3.125 MHz channel with 0.96 us bins V = 0.888, X(1) = 0.046 (the power
    has a correlation time ~1/B ~ 0.33 us); V + 2*sum(X) -> 1.
    Without blanking these are the same constants in every channel and bin.
    (Blanking, step 2, will make them per-bin streams computed from the mask.)

  info = detectChannels(info_dc, outFile, f_out)

Inputs:
  info_dc  info struct of dedisperseChannels (chanFiles, chanFreqs, chanWidth,
           fs, t0, fullySupported, edgeWidth, BnoiseChan).
  outFile  power file (float32 [nChan x nBins], little-endian); info is saved
           to <outFile>_info.mat.
  f_out    [Hz] requested bin rate. The bin length is a whole number of
           channel samples, so the rate used is fs/round(fs/f_out) (e.g.
           1.0417 MHz for 1 MHz at 4.1667 MHz); one warning if it differs.

Name-value options:
  'Tolerance'  passed to powerCovariance (default 0.002).
  'SaveInfo'   save info to <outFile>_info.mat (default true).
  'Verbose'    print a summary (default true).

Output info: the fields of detectPower's info that foldProfile uses (file,
  nChan, N, binTime0, binDt, fullySupportedBins, byteOrder, chanFreqs, ...)
  plus chanWidth, BnoiseChan, BnoiseTotal, noise (powerCovariance struct),
  chanMeanPower, samplesPerBin.
%}

arguments
    info_dc         struct
    outFile               {mustBeTextScalar}
    f_out           (1,1) double {mustBePositive, mustBeFinite}
    opts.Tolerance  (1,1) double {mustBePositive} = 0.002
    opts.SaveInfo   (1,1) logical = true
    opts.Verbose    (1,1) logical = true
end

tStart  = tic;
outFile = char(outFile);
nChan   = info_dc.nChan;
fs      = info_dc.fs;
binLen  = max(1, round(fs / f_out));
fOut    = fs / binLen;
if abs(fOut - f_out) > 1e-6 * f_out
    warning('detectChannels:rate', ...
        'f_out = %.6g Hz is not fs/integer; using %.6g Hz (%d samples per bin).', ...
        f_out, fOut, binLen);
end

outDir = fileparts(outFile);
if ~isempty(outDir) && ~isfolder(outDir)
    mkdir(outDir);
end
[d0, name] = fileparts(outFile);
tmpFiles = strings(nChan, 1);

% ---- Detect each channel (validated detectPower) ----------------------------------------
chanMean = zeros(nChan, 1);
for j = 1:nChan
    tmpFiles(j) = fullfile(d0, sprintf('%s_tmpch%03d.dat', name, j));
    dj = detectPower(info_dc.chanFiles(j), tmpFiles(j), fs, info_dc.chanFreqs(j), fOut, ...
        'FullySupported', info_dc.fullySupported, 'T0', info_dc.t0, ...
        'SaveInfo', false, 'Verbose', false);
    if j == 1
        det1 = dj;
    elseif dj.N ~= det1.N || ~isequal(dj.fullySupportedBins, det1.fullySupportedBins)
        error('detectChannels:grid', 'Channel %d has a different time grid.', j);
    end
    chanMean(j) = dj.chanMeanPower;
end
nBins = det1.N;

% ---- Interleave into [nChan x nBins], block by block ------------------------------------
cleanTmp = onCleanup(@() deleteFiles(tmpFiles));
fidT = zeros(nChan, 1);
for j = 1:nChan
    fidT(j) = fopen(tmpFiles(j), 'r', 'ieee-le');
    if fidT(j) == -1
        error('detectChannels:openTmp', 'Could not open "%s".', tmpFiles(j));
    end
end
cleanFid = onCleanup(@() arrayfun(@fclose, fidT));
[fidOut, msg] = fopen(outFile, 'w', 'ieee-le');
if fidOut == -1
    error('detectChannels:openOut', 'Could not open "%s": %s', outFile, msg);
end
cleanOut = onCleanup(@() fclose(fidOut));
blk = 2^18;
for b0 = 1:blk:nBins
    nb = min(blk, nBins - b0 + 1);
    X = zeros(nChan, nb, 'single');
    for j = 1:nChan
        X(j, :) = fread(fidT(j), [1 nb], 'single=>single');
    end
    c = fwrite(fidOut, X, 'single');
    if c ~= numel(X)
        error('detectChannels:write', 'Wrote %d of %d values (disk full?).', c, numel(X));
    end
end
clear cleanOut                                       % close the output,
clear cleanFid                                       % then the temporaries,
clear cleanTmp                                       % then delete them

% ---- Noise statistics of a channel's detected bins --------------------------------------
nc = powerCovariance(info_dc.chanWidth, info_dc.edgeWidth, fs, binLen, ...
    'Tolerance', opts.Tolerance);

info = struct();
info.file           = outFile;
info.inFile         = info_dc.file;
info.format         = 'float32 [nChan x nBins], column-major (all channels per bin)';
info.precision      = 'single';
info.byteOrder      = 'ieee-le';
info.isComplex      = false;
info.detector       = 'square-law |z|^2 per channel (detectPower), mean per channel sample per bin';
info.nChan          = nChan;
info.N              = nBins;
info.fs             = fOut;
info.actualFout     = fOut;
info.fsIn           = fs;
info.binLen         = binLen;
info.samplesPerBin  = binLen;
info.binDt          = det1.binDt;
info.binTime0       = det1.binTime0;
info.T0             = det1.T0;
info.fullySupportedBins = det1.fullySupportedBins;
info.chanFreqs      = info_dc.chanFreqs;
info.chanWidth      = info_dc.chanWidth;
info.chanMeanPower  = chanMean;
info.BnoiseChan     = info_dc.BnoiseChan;
info.BnoiseTotal    = nChan * info_dc.BnoiseChan;
info.noise          = nc;
info.elapsed        = toc(tStart);

if opts.SaveInfo
    infoFile = fullfile(d0, [name '_info.mat']);
    save(infoFile, 'info');
    info.infoFile = infoFile;
end
if opts.Verbose
    fprintf(['detectChannels: %d channels x %d bins @ %.6g Hz (%d samples, %.4g us); ' ...
             'noise per bin V %.4f, X(1) %.4f, Lmax %d (captured %.4f); %.1f s\n'], ...
        nChan, nBins, fOut, binLen, info.binDt*1e6, nc.V, nc.X(1), nc.Lmax, nc.captured, ...
        info.elapsed);
end
end


% =========================================================================================
function deleteFiles(files)
for i = 1:numel(files)
    if isfile(files(i)), delete(files(i)); end
end
end
