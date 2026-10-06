function info = applyInverseDispersion(inFile, outFile, fs, fLO, DM, fLow, fHigh, opts)
%APPLYINVERSEDISPERSION  Coherent dedispersion of a streamed complex baseband file.
%{
Removes interstellar dispersion from complex baseband (I/Q) data with a
frequency-domain chirp filter, streamed with overlap-save. This is the
exact inverse of applyDispersionStream (same constant, DM, reference
frequency and sign convention), applied after applyIQmodulation. On real
data, DM comes from the pulsar catalogue and [fLow fHigh] is the analysed
receiver band.

  info = applyInverseDispersion(inFile, outFile, fs, fLO, DM)
  info = applyInverseDispersion(inFile, outFile, fs, fLO, DM, fLow, fHigh)
  info = applyInverseDispersion(..., 'RefFreq', fHigh, 'MaxMemoryGB', 4)

Inputs:
  inFile   interleaved float32 I,Q (cf32), little-endian.
  outFile  dedispersed output, same format, same length, same time grid.
  fs       [Hz] complex sample rate of inFile.
  fLO      [Hz] LO used for downconversion (f_RF = fLO + f_baseband).
  DM       [pc cm^-3] dispersion measure to remove.
  fLow     [Hz] lower edge of the band to dedisperse (default 1.2e9).
  fHigh    [Hz] upper edge (default 1.6e9). The band must lie strictly
           inside fLO +- fs/2.

Name-value options:
  'RefFreq'      [Hz] frequency whose timing is left unchanged (default
                 fHigh). Must match the forward stage for an exact round
                 trip; any other value leaves a constant shift of
                 tau(RefFreq_fwd) - tau(RefFreq).
  'AllowRefOutsideBand'  false (default): RefFreq must lie in [fLow, fHigh].
                 true: any RefFreq > 0. Used per channel of the channelized
                 front end (dedisperseChannels): every narrow channel is
                 referenced to the same frequency (e.g. the top of the whole
                 band), so the delays between channels are removed by the
                 chirp itself. The filter formula below is exact for any
                 RefFreq; a reference above the band only adds a pure advance
                 tau(f) - tau(RefFreq) > 0 for every f, so the overlap is then
                 all on the future side (nPast = guard only).
  'EdgeFrac'     raised-cosine band-edge taper, fraction of the band,
                 inside the band (default 0.02). Keeps the filter's
                 impulse response compact; content outside the band is
                 removed (it is only noise after dispersion anyway).
  'GuardTime'    [s] extra overlap beyond the delay sweep for the taper's
                 ringing (default 20/edge width).
  'Nfft'         FFT size; [] (default) = automatic (see below).
  'MaxMemoryGB'  budget for the working buffers (default 4).
  'CheckLeakage' true (default): one extra IFFT at setup measures how much
                 filter energy falls outside the overlap region.
  'MaxLeakage'   warn above this leakage fraction (default 1e-6).
  'SaveInfo'     save info to <outFile>_info.mat (default true).
  'Verbose'      print a summary (default true).

Output:
  info  struct: format, N, fs, conventions, overlap sizes, FFT size,
        leakage, and the range of fully supported output samples.

--- Filter ------------------------------------------------------------------
Forward (applyDispersionStream):  exp(+i*2*pi*K*DM*(f-fRef)^2/(f*fRef^2))
Inverse, per baseband bin, f = fLO + f_bb inside [fLow, fHigh]:

  H(f_bb) = W(f) * exp(-i*2*pi*K*DM*(f - fRef)^2 / (f*fRef^2))
  K = 4.148808e3 s MHz^2 pc^-1 cm^3 (applied with f in Hz: K*1e12)

Group delay -(tau(f) - tau(fRef)): frequencies below fRef are advanced,
those above are delayed, so fRef keeps its time and the sweep collapses
onto it. H is 0 outside the band.

--- Overlap-save --------------------------------------------------------------
Output sample n needs input from n-nPast to n+nFuture:
  nFuture = ceil((tau(fLow)  - tau(fRef)) * fs) + guard   (largest advance)
  nPast   = ceil((tau(fRef)  - tau(fHigh)) * fs) + guard  (largest delay)
Each FFT of length Nfft yields step = Nfft - nPast - nFuture outputs. A
single preallocated buffer is shifted by step each block (no per-block
concatenation). The start is zero-padded by nPast, the end by zeros once
the file runs out. Outputs nPast+1 .. N-nFuture (1-based) see real data
over the whole filter; outside that range part of the sweep lies beyond
the file and is missing (info.fullySupported).

--- FFT size and cost -----------------------------------------------------------
Candidates are powers of two from 2^nextpow2(2*(nPast+nFuture)) (so
step >= overlap) up to "whole file in one FFT". The one with the lowest
total cost nBlocks*Nfft*log2(Nfft) within MaxMemoryGB is used. Cost scales
with fs: use the lowest IQ rate that holds the band (e.g. 500 MHz for a
400 MHz band) to halve the work of this and every later stage.
%}

arguments
    inFile                  {mustBeTextScalar}
    outFile                 {mustBeTextScalar}
    fs                (1,1) double {mustBePositive, mustBeFinite}
    fLO               (1,1) double {mustBePositive, mustBeFinite}
    DM                (1,1) double {mustBeNonnegative, mustBeFinite}
    fLow              (1,1) double {mustBePositive, mustBeFinite} = 1.2e9
    fHigh             (1,1) double {mustBePositive, mustBeFinite} = 1.6e9
    opts.RefFreq            double = []
    opts.AllowRefOutsideBand (1,1) logical = false
    opts.EdgeFrac     (1,1) double {mustBePositive} = 0.02
    opts.GuardTime          double = []
    opts.Nfft               double = []
    opts.MaxMemoryGB  (1,1) double {mustBePositive} = 4
    opts.CheckLeakage (1,1) logical = true
    opts.MaxLeakage   (1,1) double {mustBePositive} = 1e-6
    opts.SaveInfo     (1,1) logical = true
    opts.Verbose      (1,1) logical = true
end

tStart  = tic;
inFile  = char(inFile);
outFile = char(outFile);
fRef    = opts.RefFreq;
if isempty(fRef), fRef = fHigh; end

% ---- Checks -------------------------------------------------------------------
if strcmp(inFile, outFile)
    error('applyInverseDispersion:sameFile', 'inFile and outFile must differ.');
end
if fLow >= fHigh
    error('applyInverseDispersion:band', 'fLow must be below fHigh.');
end
if fLow <= fLO - fs/2 || fHigh >= fLO + fs/2
    error('applyInverseDispersion:band', ...
        ['Band %.4g-%.4g Hz does not fit inside the baseband coverage ' ...
         '%.4g-%.4g Hz (fLO +- fs/2).'], fLow, fHigh, fLO - fs/2, fLO + fs/2);
end
if ~isscalar(fRef) || ~isfinite(fRef) || fRef <= 0
    error('applyInverseDispersion:ref', 'RefFreq must be a positive scalar.');
end
if ~opts.AllowRefOutsideBand && (fRef < fLow || fRef > fHigh)
    error('applyInverseDispersion:ref', ...
        'RefFreq must lie in [fLow, fHigh] (or set AllowRefOutsideBand).');
end
if opts.EdgeFrac > 0.5
    error('applyInverseDispersion:edge', 'EdgeFrac must be <= 0.5.');
end

% ---- Overlap sizes ----------------------------------------------------------------
Kc  = 4.148808e3;                        % s MHz^2 pc^-1 cm^3
KDM = Kc * 1e12 * DM;                    % s Hz^2
tau = @(f) KDM ./ f.^2;
edgeW = opts.EdgeFrac * (fHigh - fLow);
guardTime = opts.GuardTime;
if isempty(guardTime), guardTime = 20 / edgeW; end
g = ceil(guardTime * fs);
% max(., 0): only matters for RefFreq outside the band (then one side is
% guard only); inside the band both terms are >= 0 and nothing changes.
nFuture  = max(ceil((tau(fLow) - tau(fRef)) * fs), 0) + g;
nPast    = max(ceil((tau(fRef) - tau(fHigh)) * fs), 0) + g;
nOverlap = nPast + nFuture;

% ---- Input size --------------------------------------------------------------------
fInfo = dir(inFile);
if isempty(fInfo)
    error('applyInverseDispersion:notFound', 'inFile "%s" not found.', inFile);
end
if mod(fInfo.bytes, 8) ~= 0
    warning('applyInverseDispersion:partial', ...
        'inFile size is not a multiple of 8 bytes; ignoring trailing bytes.');
end
Nx = floor(fInfo.bytes / 8);             % complex samples
if Nx == 0
    error('applyInverseDispersion:empty', 'inFile "%s" is empty.', inFile);
end

% ---- FFT size ---------------------------------------------------------------------------
% Bytes per FFT sample (complex single = 8): buffer + spectrum + H
% + product temporary + ifft output + read/write staging = ~56.
bytesPerFftSample = 56;
maxBytes = opts.MaxMemoryGB * 1e9;
if isempty(opts.Nfft)
    Nfft = chooseNfft(nOverlap, Nx, maxBytes, bytesPerFftSample);
    autoNfft = true;
else
    Nfft = opts.Nfft;
    if Nfft <= nOverlap
        error('applyInverseDispersion:nfft', ...
            'Nfft (%d) must exceed the overlap nPast+nFuture (%d).', Nfft, nOverlap);
    end
    autoNfft = false;
end
step    = Nfft - nOverlap;
nBlocks = ceil(Nx / step);

% ---- Filter on the FFT grid (in-band bins only, chunked, double phase) -------------
H = complex(zeros(1, Nfft, 'single'));
fbin  = fs / Nfft;
kLoS  = ceil((fLow  - fLO) / fbin);      % signed bin range of the band
kHiS  = floor((fHigh - fLO) / fbin);
chunk = 2^22;
for c0 = kLoS:chunk:kHiS
    ks = c0 : min(c0 + chunk - 1, kHiS);         % signed bin indices
    f  = fLO + ks * fbin;                        % RF frequency, double
    W  = ones(size(f));
    lo = f < fLow + edgeW;   W(lo) = sin(pi/2 * (f(lo) - fLow) / edgeW).^2;
    hi = f > fHigh - edgeW;  W(hi) = sin(pi/2 * (fHigh - f(hi)) / edgeW).^2;
    cyc = KDM * (f - fRef).^2 ./ (f * fRef^2);   % chirp phase in cycles
    idx = mod(ks, Nfft) + 1;                     % FFT order, 1-based
    H(idx) = single(W .* exp(-2i*pi*cyc));
end
clear f W lo hi cyc idx ks

% ---- Optional self-check: energy outside the allowed lags ----------------------------
leakage = NaN;
if opts.CheckLeakage
    h = ifft(H);
    E = sum(abs(h).^2);
    inside = sum(abs(h(1:nPast+1)).^2) + sum(abs(h(Nfft-nFuture+1:Nfft)).^2);
    leakage = max(0, 1 - inside / E);
    clear h
    if leakage > opts.MaxLeakage
        warning('applyInverseDispersion:leakage', ...
            'Filter leakage %.2e outside the overlap; increase GuardTime.', leakage);
    end
end

if opts.Verbose
    fprintf(['applyInverseDispersion: N = %d @ %.4g Hz, band %.4f-%.4f GHz, ' ...
             'ref %.4f GHz, sweep %.3f ms\n  nPast = %d, nFuture = %d, ' ...
             'Nfft = 2^%d, step = %d (%.0f%% useful), %d block(s), ~%.2f GB, ' ...
             'leakage %.1e\n'], ...
        Nx, fs, fLow/1e9, fHigh/1e9, fRef/1e9, (tau(fLow) - tau(fHigh))*1e3, ...
        nPast, nFuture, log2(Nfft), step, 100*step/Nfft, nBlocks, ...
        Nfft*bytesPerFftSample/1e9, leakage);
end

% ---- Files (explicit little-endian) -------------------------------------------------------
[fidIn, msg] = fopen(inFile, 'r', 'ieee-le');
if fidIn == -1
    error('applyInverseDispersion:openIn', 'Could not open "%s": %s', inFile, msg);
end
cleanupIn = onCleanup(@() fclose(fidIn)); %#ok<NASGU>

outDir = fileparts(outFile);
if ~isempty(outDir) && ~isfolder(outDir)
    mkdir(outDir);
end
[fidOut, msg] = fopen(outFile, 'w', 'ieee-le');
if fidOut == -1
    error('applyInverseDispersion:openOut', 'Could not open "%s": %s', outFile, msg);
end
cleanupOut = onCleanup(@() fclose(fidOut)); %#ok<NASGU>

% ---- Overlap-save stream -----------------------------------------------------------------
% In every block, buf(nPast+1) is the first input sample whose output is
% produced by that block; buf(1:nPast) is its past (zeros in block 1, the
% true start edge) and buf(nPast+step+1:end) its future.
buf = complex(zeros(1, Nfft, 'single'));
[buf, nRead] = fillBuffer(fidIn, buf, nPast, Nfft - nPast, Nx, 0);
written = 0;

for b = 1:nBlocks
    Y = fft(buf);
    Y = Y .* H;
    Y = ifft(Y);

    nV = min(step, Nx - written);
    v  = Y(nPast + 1 : nPast + nV);
    clear Y
    iq = [real(v); imag(v)];                     % 2 x nV, column-major = I,Q,...
    c  = fwrite(fidOut, iq, 'single');
    if c ~= 2*nV
        error('applyInverseDispersion:write', ...
            'Wrote %d of %d values at block %d (disk full?).', c, 2*nV, b);
    end
    written = written + nV;
    clear v iq

    if written < Nx
        buf(1:nOverlap) = buf(step+1:Nfft);      % keep the overlap
        [buf, nRead] = fillBuffer(fidIn, buf, nOverlap, step, Nx, nRead);
    end
end

% ---- Info -----------------------------------------------------------------------------------
info = struct();
info.file            = outFile;
info.inFile          = inFile;
info.format          = 'cf32 interleaved (I0,Q0,I1,Q1,...)';
info.precision       = 'single';
info.byteOrder       = 'ieee-le';
info.isComplex       = true;
info.N               = written;
info.fs              = fs;
info.actualFsOut     = fs;
info.fLO             = fLO;
info.DM              = DM;
info.fLow            = fLow;
info.fHigh           = fHigh;
info.refFreq         = fRef;
info.refOutsideBand  = fRef < fLow || fRef > fHigh;
info.dispersionConst = Kc;
info.timingConvention = 'same grid as input; content at refFreq unshifted';
info.bulkDelayRef    = tau(fRef);                % tau(refFreq) vs infinite freq, not removed
info.edgeWidth       = edgeW;
info.nPast           = nPast;
info.nFuture         = nFuture;
info.fullySupported  = [nPast + 1, max(nPast + 1, Nx - nFuture)];   % 1-based sample range
info.Nfft            = Nfft;
info.step            = step;
info.nBlocks         = nBlocks;
info.autoNfft        = autoNfft;
info.leakage         = leakage;
info.memEstimateGB   = Nfft * bytesPerFftSample / 1e9;
info.samplesRead     = nRead;
info.samplesWritten  = written;
info.elapsed         = toc(tStart);

if opts.SaveInfo
    [d, name] = fileparts(outFile);
    infoFile = fullfile(d, [name '_info.mat']);
    save(infoFile, 'info');
    info.infoFile = infoFile;
end

if opts.Verbose
    fprintf('applyInverseDispersion: wrote %d complex samples to %s in %.1f s\n', ...
        written, outFile, info.elapsed);
end
end


% =========================================================================================
function [buf, nRead] = fillBuffer(fid, buf, offset, count, Nx, nRead)
%FILLBUFFER  Read up to count complex samples into buf(offset+1 : offset+count),
% zero-filling whatever the file cannot supply (true end edge).
nWant = min(count, Nx - nRead);
nGot  = 0;
if nWant > 0
    raw  = fread(fid, [2 nWant], 'single=>single');
    nGot = size(raw, 2);
    if nGot > 0
        buf(offset + 1 : offset + nGot) = complex(raw(1, :), raw(2, :));
    end
end
if nGot < count
    buf(offset + nGot + 1 : offset + count) = 0;
end
nRead = nRead + nGot;
end


% =========================================================================================
function Nfft = chooseNfft(nOverlap, Nx, maxBytes, bytesPerFftSample)
%CHOOSENFFT  Lowest total-cost power-of-two FFT size within the memory budget.
nMin = 2^nextpow2(2 * nOverlap);          % step >= overlap
nOne = 2^nextpow2(Nx + nOverlap);         % whole file in one FFT
if nOne <= nMin
    Nfft = nOne;
else
    p     = log2(nMin):log2(nOne);
    cands = 2.^p;
    nBlk  = ceil(Nx ./ (cands - nOverlap));
    cost  = nBlk .* cands .* p;
    cost(cands * bytesPerFftSample > maxBytes) = Inf;
    [cmin, i] = min(cost);
    if isinf(cmin)
        Nfft = nMin;                      % report below
    else
        Nfft = cands(i);
    end
end
if Nfft * bytesPerFftSample > maxBytes
    error('applyInverseDispersion:memory', ...
        ['The smallest usable FFT (2^%d) needs about %.1f GB, above ' ...
         'MaxMemoryGB = %.1f. Raise the budget, lower fs, or split the band ' ...
         'into channels (smear per channel shrinks with its bandwidth).'], ...
        log2(Nfft), Nfft*bytesPerFftSample/1e9, maxBytes/1e9);
end
end