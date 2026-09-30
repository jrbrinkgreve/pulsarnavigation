function info = applyIQmodulation(inFile, outFile, f_in, fs, fLO, opts)
%APPLYIQMODULATION  Stream a real RF file to complex baseband (I/Q).
%{
Mimics an analog I/Q receiver: mix with a complex local oscillator, low-
pass filter, decimate. Streamed with overlap-add, so memory is bounded
by the FFT size, never by the file size.

  info = applyIQmodulation(inFile, outFile, f_in, fs, fLO)
  info = applyIQmodulation(..., 'FilterOrder', 2048, 'Band', [fLow fHigh])

Inputs:
  inFile   real float32, little-endian RF file (e.g. dispersed signal).
  outFile  complex output: interleaved float32 I,Q pairs (I0,Q0,I1,Q1,...,
           the .cf32 convention), little-endian.
  f_in     [Hz] sample rate of inFile.
  fs       [Hz] requested complex output rate. Decimation factor
           D = round(f_in/fs); the achieved rate f_in/D is reported in
           info.actualFsOut (warning if it differs from fs).
  fLO      [Hz] local oscillator. RF frequency f maps to baseband f - fLO
           (no spectral inversion: above-LO content is at positive
           baseband frequencies).

Name-value options:
  'FilterOrder'  FIR order (default 2048). Rounded up to even, so the
                 filter length is odd and its group delay is an integer
                 number of input samples, which is removed exactly.
  'Cutoff'       [Hz] one-sided low-pass cutoff (default 0.9*actualFsOut/2).
  'StopbandDB'   Kaiser-window stopband attenuation (default 80).
  'Band'         [fLow fHigh] RF band; if given, a warning is issued when
                 it does not fit inside the flat part of the passband.
  'Gain'         'envelope' (default): z = LPF{2*x*exp(-j*2*pi*fLO*t)},
                 the complex envelope, so x(t) = Re{z(t)*exp(j*2*pi*fLO*t)},
                 a tone A*cos(2*pi*f*t+phi) becomes A*exp(j*phi)*
                 exp(j*2*pi*(f-fLO)*t), and mean|z|^2 = 2*mean(x^2) for
                 in-band signals.
                 'unity': no factor 2 (amplitudes halved; the old behaviour).
  'BlockLen'     input samples per block; [] (default) = automatic.
  'MaxMemoryGB'  budget for the block buffers (default 2).
  'SaveInfo'     save info to <outFile>_info.mat (default true).
  'Verbose'      print progress (default true).

Output:
  info  struct: file format, N (complex samples), rates, LO, filter,
        timing and gain conventions, samples read/written.

--- Timing ------------------------------------------------------------------
The low-pass filter is linear-phase with group delay G = (M-1)/2 input
samples. That delay is removed, so output sample k (0-based) corresponds
exactly to input sample k*D, i.e. time k/actualFsOut, and the LO phase is
referenced to input sample 0 (phase 0). Without this, every arrival time
would be late by G/f_in (256 ns for 2049 taps at 4 GHz).
Output length is ceil(N_in / D). The first and last G input samples of
the filter response see zero-padding beyond the file edges.

--- Method --------------------------------------------------------------------
Per block of n input samples:
  1. LO phase in cycles, mod(n0*fLO/f_in, 1) + (0:n-1)*fLO/f_in, computed
     in double from the absolute start index n0 (no accumulated drift).
  2. Mix: z = gain * x .* exp(-j*2*pi*cycles)   (complex single).
  3. Overlap-add FFT convolution with the low-pass kernel; carry M-1
     samples into the next block (blockLen >= M-1 is enforced, so the
     carry never spans more than one block).
  4. Remove the group delay G, keep every D-th sample on the global grid,
     write as interleaved I,Q.
After the last block the carry is flushed the same way.

The low-pass is a Kaiser-windowed sinc (no Signal Processing Toolbox
needed), normalised to unit DC gain. Its transition width is about
(StopbandDB - 8) / (2.285*(M-1)) * f_in/(2*pi); ~9.8 MHz for the defaults
at 4 GHz. The stopband should start below actualFsOut/2 so nothing aliases
into the output band on decimation; a warning is issued if it does not.
%}

arguments
    inFile                 {mustBeTextScalar}
    outFile                {mustBeTextScalar}
    f_in             (1,1) double {mustBePositive, mustBeFinite}
    fs               (1,1) double {mustBePositive, mustBeFinite}
    fLO              (1,1) double {mustBePositive, mustBeFinite}
    opts.FilterOrder (1,1) double {mustBeInteger, mustBePositive} = 2048
    opts.Cutoff            double = []
    opts.StopbandDB  (1,1) double {mustBePositive} = 80
    opts.Band              double = []
    opts.Gain              {mustBeTextScalar} = 'envelope'
    opts.BlockLen          double = []
    opts.MaxMemoryGB (1,1) double {mustBePositive} = 2
    opts.SaveInfo    (1,1) logical = true
    opts.Verbose     (1,1) logical = true
end

tStart  = tic;
inFile  = char(inFile);
outFile = char(outFile);
if strcmp(inFile, outFile)
    error('applyIQmodulation:sameFile', 'inFile and outFile must differ.');
end
if fLO >= f_in/2
    error('applyIQmodulation:lo', 'fLO must be below f_in/2.');
end

switch lower(char(opts.Gain))
    case 'envelope', gainFactor = 2;
    case 'unity',    gainFactor = 1;
    otherwise
        error('applyIQmodulation:gain', 'Gain must be ''envelope'' or ''unity''.');
end

% ---- Decimation ------------------------------------------------------------------
D = max(1, round(f_in / fs));
actualFsOut = f_in / D;
if abs(actualFsOut - fs) > 1e-6*fs
    warning('applyIQmodulation:rate', ...
        ['fs = %.6g Hz is not f_in/integer; using %.6g Hz (D = %d).'], ...
        fs, actualFsOut, D);
end

% ---- Low-pass filter (Kaiser-windowed sinc, odd length) ---------------------------
M = opts.FilterOrder + 1;
if mod(M, 2) == 0
    M = M + 1;
end
G = (M - 1) / 2;                          % group delay, input samples

fc = opts.Cutoff;
if isempty(fc)
    fc = 0.9 * actualFsOut / 2;
end
if fc <= 0 || fc >= f_in/2
    error('applyIQmodulation:cutoff', 'Cutoff must lie in (0, f_in/2).');
end
h = kaiserLowpass(M, fc / f_in, opts.StopbandDB);
transW = (opts.StopbandDB - 8) / (2.285*(M - 1)) * f_in / (2*pi);

if fc + transW/2 > actualFsOut/2
    warning('applyIQmodulation:alias', ...
        ['Stopband starts at %.4g Hz, above the output Nyquist %.4g Hz: ' ...
         'content there aliases on decimation. Lower Cutoff or raise FilterOrder.'], ...
        fc + transW/2, actualFsOut/2);
end
if ~isempty(opts.Band)
    bb = opts.Band - fLO;
    if any(abs(bb) > fc - transW/2)
        warning('applyIQmodulation:band', ...
            ['Band edges map to %.4g..%.4g Hz at baseband, outside the flat ' ...
             'passband (+-%.4g Hz).'], bb(1), bb(end), fc - transW/2);
    end
end

% ---- Input size ------------------------------------------------------------------------
fInfo = dir(inFile);
if isempty(fInfo)
    error('applyIQmodulation:notFound', 'inFile "%s" not found.', inFile);
end
if mod(fInfo.bytes, 4) ~= 0
    warning('applyIQmodulation:partial', ...
        'inFile size is not a multiple of 4 bytes; ignoring trailing bytes.');
end
Nx = floor(fInfo.bytes / 4);
if Nx == 0
    error('applyIQmodulation:empty', 'inFile "%s" is empty.', inFile);
end

% ---- Block size ------------------------------------------------------------------------
% Bytes per FFT sample: LO cycles (double, 8) + LO complex double temp (16)
% + LO single (8) + input (4) + mixed (8) + spectrum (8) + H (8)
% + product temp (8) + ifft output (8) + slices (8) = ~84.
bytesPerFftSample = 84;
maxBytes = opts.MaxMemoryGB * 1e9;
nOne = 2^nextpow2(Nx + M - 1);            % whole file in one block
nMin = 2^nextpow2(max(2*M - 2, 1));       % blockLen >= M-1
if isempty(opts.BlockLen)
    % Filter is short: a few-million-sample FFT is already >99.9% efficient;
    % larger only costs memory.
    Nfft = min(2^nextpow2(max(16*M, 2^22)), nOne);
    while Nfft * bytesPerFftSample > maxBytes && Nfft > nMin
        Nfft = Nfft / 2;
    end
    blockLen = min(Nfft - M + 1, Nx);
    autoBlock = true;
else
    blockLen = min(round(opts.BlockLen), Nx);
    if blockLen < 1
        error('applyIQmodulation:blockLen', 'BlockLen must be positive.');
    end
    if blockLen < Nx && blockLen < M - 1
        error('applyIQmodulation:blockLen', ...
            'BlockLen (%d) must be >= filter length - 1 (%d).', blockLen, M - 1);
    end
    Nfft = 2^nextpow2(blockLen + M - 1);
    autoBlock = false;
end
if Nfft * bytesPerFftSample > maxBytes
    warning('applyIQmodulation:memory', ...
        'Block buffers need about %.2f GB, above MaxMemoryGB = %.2f.', ...
        Nfft*bytesPerFftSample/1e9, opts.MaxMemoryGB);
end
nBlocks = ceil(Nx / blockLen);
nOutExpected = ceil(Nx / D);

if opts.Verbose
    fprintf(['applyIQmodulation: N_in = %d @ %.4g Hz -> N_out = %d @ %.4g Hz ' ...
             '(D = %d), LO %.4g Hz, LPF %d taps, cutoff %.4g Hz, ' ...
             'Nfft = 2^%d, %d block(s)\n'], ...
        Nx, f_in, nOutExpected, actualFsOut, D, fLO, M, fc, log2(Nfft), nBlocks);
end

% ---- Files (explicit little-endian) --------------------------------------------------
[fidIn, msg] = fopen(inFile, 'r', 'ieee-le');
if fidIn == -1
    error('applyIQmodulation:openIn', 'Could not open "%s": %s', inFile, msg);
end
cleanupIn = onCleanup(@() fclose(fidIn)); %#ok<NASGU>

outDir = fileparts(outFile);
if ~isempty(outDir) && ~isfolder(outDir)
    mkdir(outDir);
end
[fidOut, msg] = fopen(outFile, 'w', 'ieee-le');
if fidOut == -1
    error('applyIQmodulation:openOut', 'Could not open "%s": %s', outFile, msg);
end
cleanupOut = onCleanup(@() fclose(fidOut)); %#ok<NASGU>

% ---- Stream ----------------------------------------------------------------------------------
H = fft(single(h), Nfft);
r = fLO / f_in;                            % LO cycles per input sample
carry    = complex(zeros(1, M - 1, 'single'));
fullPos  = 0;                              % input samples processed so far
nWritten = 0;

for b = 1:nBlocks
    x = fread(fidIn, [1 blockLen], 'single=>single');
    n = numel(x);
    if n == 0
        break
    end

    % 1-2. Mix with the LO; phase from the absolute index, in double
    cyc = mod(fullPos * r, 1) + (0:n-1) * r;
    lo  = single(gainFactor * exp(-2i*pi*cyc));
    clear cyc
    z = x .* lo;
    clear x lo

    % 3. Overlap-add low-pass
    Z = fft(z, Nfft);
    clear z
    Z = Z .* H;
    y = ifft(Z);
    clear Z
    y = y(1:n + M - 1);
    y(1:M-1) = y(1:M-1) + carry;
    carry = y(n+1:end);

    % 4. Remove group delay, decimate, write
    nWritten = nWritten + writeDecimated(fidOut, y(1:n), fullPos, G, Nx, D, b);
    fullPos = fullPos + n;

    if opts.Verbose && (b == nBlocks || mod(b, 20) == 0)
        fprintf('  block %d/%d (%.1f%%)\n', b, nBlocks, 100*fullPos/Nx);
    end
end
% Flush: carry holds full-convolution samples Nx+1 .. Nx+M-1
nWritten = nWritten + writeDecimated(fidOut, carry, Nx, G, Nx, D, nBlocks + 1);

if nWritten ~= nOutExpected
    warning('applyIQmodulation:count', ...
        'Wrote %d complex samples, expected %d.', nWritten, nOutExpected);
end

% ---- Info --------------------------------------------------------------------------------------
info = struct();
info.file             = outFile;
info.inFile           = inFile;
info.format           = 'cf32 interleaved (I0,Q0,I1,Q1,...)';
info.precision        = 'single';
info.byteOrder        = 'ieee-le';
info.isComplex        = true;
info.N                = nWritten;           % complex samples
info.fs               = actualFsOut;
info.actualFsOut      = actualFsOut;
info.fsIn             = f_in;
info.samplesRead      = fullPos;
info.samplesWritten   = nWritten;
info.decimationFactor = D;
info.fLO              = fLO;
info.freqMapping      = 'f_baseband = f_RF - fLO (no inversion)';
info.loPhaseRef       = 'LO phase 0 at input sample 0';
info.gainConvention   = lower(char(opts.Gain));
info.gainFactor       = gainFactor;
info.cutoff           = fc;
info.transitionWidth  = transW;
info.stopbandDB       = opts.StopbandDB;
info.filterLen        = M;
info.groupDelayRemoved = G;                 % input samples (G/f_in seconds)
info.t0               = 0;                  % output sample k at t = t0 + k/fs
info.Nfft             = Nfft;
info.blockLen         = blockLen;
info.nBlocks          = nBlocks;
info.autoBlockSize    = autoBlock;
info.elapsed          = toc(tStart);

if opts.SaveInfo
    [d, name] = fileparts(outFile);
    infoFile = fullfile(d, [name '_info.mat']);
    save(infoFile, 'info');
    info.infoFile = infoFile;
end

if opts.Verbose
    fprintf('applyIQmodulation: wrote %d complex samples to %s in %.1f s\n', ...
        nWritten, outFile, info.elapsed);
end
end


% =====================================================================================
function cnt = writeDecimated(fid, y, fullPos, G, Nx, D, blockNum)
%WRITEDECIMATED  Remove group delay G and keep every D-th sample.
% Local sample i (1-based) of y is full-convolution index fullPos + i; after
% removing G it is aligned full-rate index j = fullPos + i - G (1-based).
% Keep 1 <= j <= Nx with mod(j-1, D) == 0.
n  = numel(y);
i1 = max(1, G + 1 - fullPos);
i2 = min(n, Nx + G - fullPos);
cnt = 0;
if i2 < i1
    return
end
j1 = fullPos + i1 - G;
i1 = i1 + mod(-(j1 - 1), D);               % first index on the decimation grid
if i1 > i2
    return
end
v  = y(i1:D:i2);
iq = [real(v); imag(v)];                   % 2 x K -> column-major = I,Q,I,Q,...
c  = fwrite(fid, iq(:), 'single');
if c ~= 2*numel(v)
    error('applyIQmodulation:write', ...
        'Wrote %d of %d values at block %d (disk full?).', c, 2*numel(v), blockNum);
end
cnt = numel(v);
end


% =====================================================================================
function h = kaiserLowpass(M, fcNorm, A)
%KAISERLOWPASS  Odd-length Kaiser-windowed sinc low-pass, unit DC gain.
% fcNorm = cutoff / sample rate. Base MATLAB only (besseli).
if A > 50
    beta = 0.1102 * (A - 8.7);
elseif A >= 21
    beta = 0.5842 * (A - 21)^0.4 + 0.07886 * (A - 21);
else
    beta = 0;
end
G = (M - 1) / 2;
n = (0:M-1) - G;
h = 2 * fcNorm * ones(1, M);
nz = n ~= 0;
h(nz) = sin(2*pi*fcNorm*n(nz)) ./ (pi*n(nz));
w = besseli(0, beta * sqrt(1 - (2*(0:M-1)/(M-1) - 1).^2)) / besseli(0, beta);
h = h .* w;
h = h / sum(h);
end