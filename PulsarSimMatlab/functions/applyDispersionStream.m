function info = applyDispersionStream(inFile, outFile, DM, fs, fLow, fHigh, opts)
%APPLYDISPERSIONSTREAM  Apply interstellar dispersion to a streamed RF file.
%{
Streams a real float32 RF-sampled file through a dispersion FIR filter
using overlap-add FFT convolution. Memory use is bounded by the FFT size,
never by the file size. The block size is chosen automatically from the
DM and band so that the overlap-add carry always fits inside the next
block.

  info = applyDispersionStream(inFile, outFile, DM, fs)
  info = applyDispersionStream(inFile, outFile, DM, fs, fLow, fHigh)
  info = applyDispersionStream(..., 'RefFreq', 1.4e9, 'MaxMemoryGB', 8, ...)

Inputs:
  inFile   real float32, little-endian ('ieee-le') sample file, e.g. from
           generatePulsarSignal.
  outFile  output file, same format, same sample count as inFile.
  DM       [pc cm^-3] dispersion measure (0 gives a pure band-pass).
  fs       [Hz] sample rate of inFile.
  fLow     [Hz] lower band edge (default 1.2e9).
  fHigh    [Hz] upper band edge (default 1.6e9). Must be < fs/2.

Name-value options:
  'RefFreq'     [Hz] reference frequency, fLow <= RefFreq <= fHigh
                (default fHigh). Content at RefFreq keeps its original
                time; every other frequency f is shifted by
                tau(f) - tau(RefFreq), tau(f) = 4.148808e3 s MHz^2 * DM/f^2.
                With the default, the top of the band stays put and
                lower frequencies arrive later, as observed in reality.
  'EdgeFrac'    raised-cosine band-edge taper width, as a fraction of
                (fHigh - fLow), placed inside the band (default 0.02).
  'GuardTime'   [s] kernel guard before/after the delay sweep, holding the
                ringing of the band-edge taper (default 20/edge width).
  'MaxLeakage'  warn if the kernel energy lost by truncation exceeds this
                fraction (default 1e-8).
  'BlockLen'    samples per block; [] (default) = automatic. A manual value
                must be >= kernel length - 1 unless the file fits in one
                block.
  'MaxMemoryGB' memory budget for the streaming FFT buffers (default 16).
  'SaveInfo'    save info to <outFile>_info.mat (default true).
  'Verbose'     print progress (default true).

Output:
  info  struct: file format, N, fs, dispersion parameters and convention,
        kernel and block-size details, samples written, timing.

--- Physics / sign convention ---------------------------------------------
Cold plasma gives a PHASE advance and a GROUP delay. The transfer function
used here is the standard coherent-dispersion chirp

  H(f) = W(f) * exp( +i*2*pi*K*DM*(f - fRef)^2 / (f*fRef^2) )
              * exp( -i*2*pi*f*s0/fs )

whose group delay, -1/(2*pi) * dphi/df, equals

  tau_g(f) = K*DM/f^2 - K*DM/fRef^2 + s0/fs.

So low frequencies arrive later, fRef arrives exactly s0 samples into the
kernel, and the stream function removes those s0 samples again, so fRef
lands on its original time. W(f) is 1 in the band, has raised-cosine edges
of width EdgeFrac*(fHigh-fLow) just inside fLow and fHigh, and is 0
outside. The output is therefore band-limited to [fLow, fHigh]. That is
necessary (delays diverge as f -> 0) and is also what a receiver with that
band sees. Output power is about 2*(B - 1.25*edgeW)/fs of the input power
for white input; see info.kernelEnergy.

The constant bulk delay of fRef relative to infinite frequency,
tau(fRef) = K*DM/fRef^2, is NOT applied; it is reported in
info.bulkDelayRef so arrival times can be converted to infinite frequency.

To undo this stage, the inverse must use the conjugate chirp with the SAME
DM, fRef and constant K (after downconversion: baseband frequency
f - fLO, with fRef mapping to fRef - fLO).

--- Overlap-add -------------------------------------------------------------
Each block of n samples is convolved (FFT size Nfft >= n + M - 1, M =
kernel length). The first n samples are final once the previous block's
carry (M-1 samples) has been added; the last M-1 samples become the new
carry. Automatic sizing guarantees blockLen >= M-1, so a carry only ever
spills into the immediately following block. Among power-of-two FFT sizes
that fit the memory budget, the one with the lowest total FFT cost for
this file is chosen.

Edges: the carry starts at zero (no data before sample 1). After the last
block the carry is flushed; the part that maps into [1, N] after removing
s0 is written, the rest (delayed past the end of the file) is dropped.
%}

arguments
    inFile                 {mustBeTextScalar}
    outFile                {mustBeTextScalar}
    DM               (1,1) double {mustBeNonnegative, mustBeFinite}
    fs               (1,1) double {mustBePositive, mustBeFinite}
    fLow             (1,1) double {mustBePositive, mustBeFinite} = 1.2e9
    fHigh            (1,1) double {mustBePositive, mustBeFinite} = 1.6e9
    opts.RefFreq           double = []
    opts.EdgeFrac    (1,1) double {mustBePositive} = 0.02
    opts.GuardTime         double = []
    opts.MaxLeakage  (1,1) double {mustBePositive} = 1e-8
    opts.BlockLen          double = []
    opts.MaxMemoryGB (1,1) double {mustBePositive} = 16
    opts.SaveInfo    (1,1) logical = true
    opts.Verbose     (1,1) logical = true
end

tStart  = tic;
inFile  = char(inFile);
outFile = char(outFile);
fRef    = opts.RefFreq;
if isempty(fRef), fRef = fHigh; end

% ---- Parameter checks ---------------------------------------------------
if fLow >= fHigh
    error('applyDispersionStream:band', 'fLow must be below fHigh.');
end
if fHigh >= fs/2
    error('applyDispersionStream:nyquist', ...
        'fHigh (%.4g Hz) must be below fs/2 (%.4g Hz).', fHigh, fs/2);
end
if ~isscalar(fRef) || fRef < fLow || fRef > fHigh
    error('applyDispersionStream:ref', 'RefFreq must lie in [fLow, fHigh].');
end
if opts.EdgeFrac > 0.5
    error('applyDispersionStream:edge', 'EdgeFrac must be <= 0.5.');
end
if strcmp(inFile, outFile)
    error('applyDispersionStream:sameFile', 'inFile and outFile must differ.');
end

% ---- Input size -----------------------------------------------------------
fInfo = dir(inFile);
if isempty(fInfo)
    error('applyDispersionStream:notFound', 'inFile "%s" not found.', inFile);
end
if mod(fInfo.bytes, 4) ~= 0
    warning('applyDispersionStream:partial', ...
        'inFile size is not a multiple of 4 bytes; ignoring trailing bytes.');
end
Nx = floor(fInfo.bytes / 4);
if Nx == 0
    error('applyDispersionStream:empty', 'inFile "%s" is empty.', inFile);
end

% ---- Kernel -----------------------------------------------------------------
[h, kinfo] = makeDispersionKernel(DM, fLow, fHigh, fRef, fs, ...
    opts.EdgeFrac, opts.GuardTime, opts.Verbose);
M  = numel(h);
s0 = kinfo.zeroLag;                  % samples to remove so fRef stays put
if kinfo.leakage > opts.MaxLeakage
    warning('applyDispersionStream:leakage', ...
        'Kernel truncation loses %.2e of its energy; increase GuardTime.', ...
        kinfo.leakage);
end

% ---- Block size -------------------------------------------------------------
bytesPerFftSample = 44;              % see estimate in chooseBlockSize
maxBytes = opts.MaxMemoryGB * 1e9;
if isempty(opts.BlockLen)
    [Nfft, blockLen] = chooseBlockSize(M, Nx, maxBytes, bytesPerFftSample);
    autoBlock = true;
else
    blockLen = min(round(opts.BlockLen), Nx);
    if blockLen < 1
        error('applyDispersionStream:blockLen', 'BlockLen must be positive.');
    end
    if blockLen < Nx && blockLen < M - 1
        error('applyDispersionStream:blockLen', ...
            ['BlockLen (%d) is shorter than the carry (M-1 = %d), so a ' ...
             'carry would spill over more than one block. Use at least %d, ' ...
             'or leave BlockLen empty for automatic sizing.'], ...
            blockLen, M - 1, M - 1);
    end
    Nfft = 2^nextpow2(blockLen + M - 1);
    autoBlock = false;
    if Nfft * bytesPerFftSample > maxBytes
        warning('applyDispersionStream:memory', ...
            'BlockLen %d needs about %.1f GB, above MaxMemoryGB = %.1f.', ...
            blockLen, Nfft*bytesPerFftSample/1e9, opts.MaxMemoryGB);
    end
end
nBlocks = ceil(Nx / blockLen);

if opts.Verbose
    fprintf(['applyDispersionStream: N = %d, kernel M = %d, Nfft = 2^%d, ' ...
             'blockLen = %d (%.0f%% useful), %d block(s), ~%.1f GB\n'], ...
        Nx, M, log2(Nfft), blockLen, 100*blockLen/Nfft, nBlocks, ...
        Nfft*bytesPerFftSample/1e9);
end

% ---- Files (explicit little-endian) -----------------------------------------
[fidIn, msg] = fopen(inFile, 'r', 'ieee-le');
if fidIn == -1
    error('applyDispersionStream:openIn', 'Could not open "%s": %s', inFile, msg);
end
cleanupIn = onCleanup(@() fclose(fidIn)); %#ok<NASGU>

outDir = fileparts(outFile);
if ~isempty(outDir) && ~isfolder(outDir)
    mkdir(outDir);
end
[fidOut, msg] = fopen(outFile, 'w', 'ieee-le');
if fidOut == -1
    error('applyDispersionStream:openOut', 'Could not open "%s": %s', outFile, msg);
end
cleanupOut = onCleanup(@() fclose(fidOut)); %#ok<NASGU>

% ---- Overlap-add stream -------------------------------------------------------
H = fft(h, Nfft);                    % single complex, computed once
clear h
carry    = zeros(1, M - 1, 'single');
fullPos  = 0;                        % full-convolution samples finalized so far
nWritten = 0;

for b = 1:nBlocks
    xBlock = fread(fidIn, [1 blockLen], 'single=>single');
    n = numel(xBlock);
    if n == 0
        break
    end

    Y = fft(xBlock, Nfft);
    clear xBlock
    Y = Y .* H;
    y = ifft(Y, 'symmetric');        % real single
    clear Y
    y = y(1:n + M - 1);

    y(1:M-1) = y(1:M-1) + carry;     % add previous block's tail
    carry = y(n+1:end);              % tail for the next block
    % y(1:n) is final. Its full-convolution indices are fullPos+1..fullPos+n;
    % output index = full index - s0; keep 1..Nx.
    i1 = max(1, s0 + 1 - fullPos);
    i2 = min(n, Nx + s0 - fullPos);
    if i2 >= i1
        nWritten = nWritten + writeChecked(fidOut, y(i1:i2), b);
    end
    fullPos = fullPos + n;

    if opts.Verbose
        fprintf('  block %d/%d done (%.1f%%)\n', b, nBlocks, 100*fullPos/Nx);
    end
end

% Flush: carry holds full indices Nx+1 .. Nx+M-1.
i1 = max(1, s0 + 1 - Nx);
i2 = min(M - 1, s0);
if i2 >= i1
    nWritten = nWritten + writeChecked(fidOut, carry(i1:i2), nBlocks + 1);
end

if nWritten ~= Nx
    warning('applyDispersionStream:count', ...
        'Wrote %d samples but read %d.', nWritten, Nx);
end

% ---- Info ---------------------------------------------------------------------
info = struct();
info.file              = outFile;
info.inFile            = inFile;
info.precision         = 'single';
info.byteOrder         = 'ieee-le';
info.isComplex         = false;
info.N                 = Nx;
info.fs                = fs;
info.actualFsOut       = fs;
info.DM                = DM;
info.fLow              = fLow;
info.fHigh             = fHigh;
info.refFreq           = fRef;
info.dispersionConst   = kinfo.K;            % s MHz^2 pc^-1 cm^3
info.delayConvention   = ['group delay tau(f)-tau(refFreq), ' ...
                          'tau(f) = dispersionConst*DM/(f/1e6)^2 s'];
info.smearTime         = kinfo.smearTime;    % tau(fLow) - tau(fHigh)
info.bulkDelayRef      = kinfo.bulkDelayRef; % tau(refFreq), not applied
info.edgeWidth         = kinfo.edgeWidth;
info.kernelLen         = M;
info.kernelZeroLag     = s0;
info.kernelLeakage     = kinfo.leakage;
info.kernelEnergy      = kinfo.energy;       % output/input power, white input
info.Nfft              = Nfft;
info.blockLen          = blockLen;
info.nBlocks           = nBlocks;
info.autoBlockSize     = autoBlock;
info.memEstimateGB     = Nfft * bytesPerFftSample / 1e9;
info.kernelDesignGB    = kinfo.designMemGB;
info.samplesWritten    = nWritten;
info.elapsed           = toc(tStart);

if opts.SaveInfo
    [d, name] = fileparts(outFile);
    infoFile = fullfile(d, [name '_info.mat']);
    save(infoFile, 'info');
    info.infoFile = infoFile;
end

if opts.Verbose
    fprintf('applyDispersionStream: wrote %d samples to %s in %.1f s\n', ...
        nWritten, outFile, info.elapsed);
end
end


% =============================================================================
function cnt = writeChecked(fid, data, blockNum)
cnt = fwrite(fid, data, 'single');
if cnt ~= numel(data)
    error('applyDispersionStream:write', ...
        'Wrote %d of %d samples at block %d (disk full?).', ...
        cnt, numel(data), blockNum);
end
end


% =============================================================================
function [Nfft, blockLen] = chooseBlockSize(M, Nx, maxBytes, bytesPerFftSample)
%CHOOSEBLOCKSIZE  Pick the FFT size for overlap-add.
% Memory per FFT sample (bytes), all single precision:
%   H (complex, 8) + spectrum (8) + product temporary (8)
%   + ifft workspace (8) + real output (4) + input block (4) + slices (4)
%   = ~44
% Constraint: blockLen = Nfft - M + 1 >= M - 1, i.e. Nfft >= 2M - 2, so the
% carry of one block only ever lands in the next one.
% Among power-of-two sizes from that minimum up to "whole file in one
% block", pick the one with the lowest total cost nBlocks * Nfft*log2(Nfft)
% that fits the memory budget.

nOne = 2^nextpow2(Nx + M - 1);           % whole file in a single block
nMin = 2^nextpow2(max(2*M - 2, 1));

if nOne <= nMin
    % File shorter than one minimum block: single FFT, carry is only flushed
    Nfft = nOne;
    blockLen = Nx;
    if Nfft * bytesPerFftSample > maxBytes
        memError(Nfft, bytesPerFftSample, maxBytes);
    end
    return
end

p      = log2(nMin):log2(nOne);
cands  = 2.^p;
blk    = min(cands - M + 1, Nx);
nBlk   = ceil(Nx ./ blk);
cost   = nBlk .* cands .* p;
fits   = cands * bytesPerFftSample <= maxBytes;
if ~any(fits)
    memError(nMin, bytesPerFftSample, maxBytes);
end
cost(~fits) = Inf;
[~, i]   = min(cost);
Nfft     = cands(i);
blockLen = blk(i);
end


function memError(Nfft, bytesPerFftSample, maxBytes)
error('applyDispersionStream:memory', ...
    ['The smallest usable FFT (2^%d) needs about %.1f GB, above the ' ...
     'budget of %.1f GB. Raise MaxMemoryGB, lower DM/bandwidth/fs, or ' ...
     'split the band into channels.'], ...
    log2(Nfft), Nfft*bytesPerFftSample/1e9, maxBytes/1e9);
end


% =============================================================================
function [h, kinfo] = makeDispersionKernel(DM, fLow, fHigh, fRef, fs, ...
    edgeFrac, guardTime, verbose)
%MAKEDISPERSIONKERNEL  One-sided, band-limited dispersion FIR kernel.
%
% Designed in the frequency domain on a grid at least twice the kernel
% length (so nothing wraps), transformed to time, and truncated to
%   [pre-guard | delay sweep | post-guard].
% Zero lag of the reference frequency sits at sample s0 (0-based), which
% the stream function removes again. All phases are computed in double.

K   = 4.148808e3;                        % s MHz^2 pc^-1 cm^3
KDM = K * 1e12 * DM;                     % s Hz^2: tau(f) = KDM / f^2
tau = @(f) KDM ./ f.^2;

edgeW = edgeFrac * (fHigh - fLow);
if isempty(guardTime)
    guardTime = 20 / edgeW;              % leakage ~1e-12 in tests
end
g = ceil(guardTime * fs);

preSpan  = tau(fRef) - tau(fHigh);       % how far fHigh leads fRef (>= 0)
postSpan = tau(fLow) - tau(fRef);        % how far fLow lags fRef   (>= 0)
s0 = g + ceil(preSpan * fs);             % zero lag of fRef, samples
Nk = s0 + ceil(postSpan * fs) + g + 1;   % kernel length
Nf = 2^nextpow2(2 * Nk);                 % design grid, no wrap-around

% In-band bins only (0-based bin index kb, frequency kb*fs/Nf)
kb = ceil(fLow * Nf / fs) : floor(fHigh * Nf / fs);
fb = kb * (fs / Nf);

% Raised-cosine band edges inside the band
W = ones(size(fb));
lo = fb < fLow + edgeW;
W(lo) = sin(pi/2 * (fb(lo) - fLow) / edgeW).^2;
hi = fb > fHigh - edgeW;
W(hi) = sin(pi/2 * (fHigh - fb(hi)) / edgeW).^2;

% Chirp phase (group delay tau(f) - tau(fRef)) plus exact integer shift of
% s0 samples; mod() keeps the linear term exact for large bin indices.
phiChirp = 2*pi*KDM * (fb - fRef).^2 ./ (fb * fRef^2);
phiShift = 2*pi * mod(kb * s0, Nf) / Nf;
Hb = W .* exp(1i * (phiChirp - phiShift));
clear fb W lo hi phiChirp phiShift

Hfull = zeros(1, Nf);
Hfull(kb + 1)      = Hb;
Hfull(Nf - kb + 1) = conj(Hb);           % conjugate symmetry -> real h
clear Hb
hFull = ifft(Hfull, 'symmetric');
clear Hfull

E    = sum(hFull.^2);
h    = hFull(1:Nk);
leak = 1 - sum(h.^2) / E;
h    = single(h);

kinfo = struct( ...
    'K',            K, ...
    'smearTime',    tau(fLow) - tau(fHigh), ...
    'bulkDelayRef', tau(fRef), ...
    'edgeWidth',    edgeW, ...
    'guardSamples', g, ...
    'zeroLag',      s0, ...
    'kernelLen',    Nk, ...
    'designGrid',   Nf, ...
    'leakage',      leak, ...
    'energy',       E, ...
    'designMemGB',  Nf * 40 / 1e9);

if verbose
    fprintf(['makeDispersionKernel: DM = %.3f, band %.3f-%.3f GHz, ref %.3f GHz, ' ...
             'smear = %.3f ms, kernel = %d samples (%.1f MB), leakage = %.1e\n'], ...
        DM, fLow/1e9, fHigh/1e9, fRef/1e9, kinfo.smearTime*1e3, Nk, ...
        Nk*4/1e6, leak);
end
end