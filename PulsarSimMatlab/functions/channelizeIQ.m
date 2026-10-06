function info = channelizeIQ(inFile, outBase, fs, fLO, fLow, fHigh, opts)
%CHANNELIZEIQ  Split a complex baseband stream into frequency channels (oversampled PFB).
%{
First stage of the channelized front end: the IQ stream (e.g. 800 MHz wide)
is split into nChan narrow channels that tile [fLow, fHigh] exactly, each a
complex voltage stream at a reduced rate. RFI excision and coherent
dedispersion then work per channel: narrowband RFI stays in one or two
channels, and the dispersion sweep inside one channel is short.

  info = channelizeIQ(inFile, outBase, fs, fLO, fLow, fHigh)
  info = channelizeIQ(..., 'ChanWidth', 3.125e6, 'Oversampling', 4/3)

Inputs:
  inFile   complex IQ, interleaved float32 I,Q (cf32), little-endian
           (output of applyIQmodulation).
  outBase  path without extension, e.g. "data/chan/test_rx_IQ_chan".
           Channel j is written to <outBase>_ch001.dat ... (cf32), the info
           struct to <outBase>_info.mat (so loadInfo(outBase) works).
  fs       [Hz] complex sample rate of inFile.
  fLO      [Hz] LO of inFile (RF = fLO + baseband, no inversion).
  fLow, fHigh  [Hz] RF band to cover; must be a whole number of channels.

Name-value options:
  'ChanWidth'    [Hz] channel spacing = useful width (default 3.125e6).
                 fs/ChanWidth must be an integer K (FFT size).
  'Oversampling' channel rate / ChanWidth (default 4/3). K/Oversampling
                 must be an integer D (decimation).
  'StopbandDB'   prototype filter stopband attenuation (default 80).
  'T0'           [s] time of input sample 0 (default 0).
  'BlockOut'     output samples per channel per block (default 1024).
  'SaveInfo'     save info to <outBase>_info.mat (default true).
  'Verbose'      print a summary (default true).

--- What each channel is ------------------------------------------------------------
Channel j (1-based) has centre fc_j = fLow + (j - 1/2)*ChanWidth (RF). Its
stream is exactly what a separate receiver would give: mix the input down
by fc_j (LO phase 0 at input sample 0), low-pass with the zero-phase
prototype h, keep every D-th sample:

  y_j(m) = sum_n h(m*D - n) * x(n) * exp(-i*2*pi*(fc_j - fLO)*n/fs)

  - sample m (0-based) is at time T0 + m*D/fs; the filter's group delay
    is removed (no TOA offset);
  - RF = fc_j + f_bb (no inversion), so each channel file is a normal IQ
    file for later stages with fLO = fc_j and fs = fs/D;
  - unit DC gain: a tone keeps its complex amplitude, and the passband PSD
    equals the input PSD. A tone A*exp(i*(2*pi*f*n/fs + phi)) gives
    y_j(m) = A*exp(i*phi) * Hc(df) * exp(i*2*pi*df*m*D/fs), df = f - fc_j,
    Hc = frequency response of the prototype (real, ~1 in the passband).

--- Why oversampled ----------------------------------------------------------------------
The channel rate fs/D (4.17 MHz for 3.125 MHz channels at 4/3) is above the
channel width. The prototype is flat to +-ChanWidth/2 and reaches its
stopband at fs/D - ChanWidth/2: content beyond the stopband aliases on
decimation, but only to frequencies outside +-ChanWidth/2. The useful band
of every channel is therefore alias-free; the transition region is
removed later (per-channel dedispersion keeps exactly +-ChanWidth/2), so
neighbouring channels tile the band without gaps or double counting.
(A critically sampled filterbank, D = K, aliases into the channel edges.)

--- Method (polyphase form of the formula above) -------------------------------------------
With K = fs/ChanWidth, channel offsets fc_j - fLO = (k_j + beta)*fs/K
(integer k_j, beta = 1/2 when the band edges sit on the K-grid) and the
prototype h(r), r = 0..L-1, centred at G = (M-1)/2 (zero-padded to
L = P*K taps):

  y_j(m) = exp(-i*2*pi*(k_j+beta)*(m*D+G)/K)
           * sum_p exp(+i*2*pi*k_j*p/K) * u_m(p),
  u_m(p) = sum_q hb(p + q*K) * x(m*D + G - p - q*K),
  hb(r)  = h(r) * exp(+i*2*pi*beta*r/K).

So per output sample: weight a frame of L input samples, fold it into K
bins (sum over q), one K-point IFFT gives all channels, then a phase
correction from the absolute sample index (double precision). This is the
same arithmetic as K separate mix-filter-decimate receivers, at ~1/K of
the cost.

Prototype: Kaiser-windowed sinc (as applyIQmodulation), cutoff fs/(2D),
transition width fs/D - ChanWidth, odd length M from the Kaiser formula.

--- Edges --------------------------------------------------------------------------------------
Output m uses inputs m*D - G .. m*D + G; outside the file they are zero.
info.fullySupported is the 1-based output range that sees only real data.
Output length per channel: ceil(N_in/D).
%}

arguments
    inFile                 {mustBeTextScalar}
    outBase                {mustBeTextScalar}
    fs               (1,1) double {mustBePositive, mustBeFinite}
    fLO              (1,1) double {mustBePositive, mustBeFinite}
    fLow             (1,1) double {mustBePositive, mustBeFinite}
    fHigh            (1,1) double {mustBePositive, mustBeFinite}
    opts.ChanWidth   (1,1) double {mustBePositive, mustBeFinite} = 3.125e6
    opts.Oversampling (1,1) double {mustBeGreaterThanOrEqual(opts.Oversampling, 1)} = 4/3
    opts.StopbandDB  (1,1) double {mustBePositive} = 80
    opts.T0          (1,1) double {mustBeFinite} = 0
    opts.BlockOut    (1,1) double {mustBeInteger, mustBePositive} = 1024
    opts.SaveInfo    (1,1) logical = true
    opts.Verbose     (1,1) logical = true
end

tStart  = tic;
inFile  = char(inFile);
outBase = char(outBase);
dF      = opts.ChanWidth;

% ---- Channel grid -------------------------------------------------------------------------
K = fs / dF;
if abs(K - round(K)) > 1e-9 * K
    error('channelizeIQ:grid', 'fs/ChanWidth = %.6g must be an integer.', K);
end
K = round(K);
D = K / opts.Oversampling;
if abs(D - round(D)) > 1e-9 * D
    error('channelizeIQ:grid', 'K/Oversampling = %d/%.6g must be an integer.', K, opts.Oversampling);
end
D = round(D);
if fLow >= fHigh
    error('channelizeIQ:band', 'fLow must be below fHigh.');
end
nChan = (fHigh - fLow) / dF;
if abs(nChan - round(nChan)) > 1e-6
    error('channelizeIQ:band', ...
        'Band width %.6g Hz is not a whole number of %.6g Hz channels.', fHigh - fLow, dF);
end
nChan = round(nChan);

a     = (fLow - fLO) / dF + 0.5;         % first channel centre, in channel units
kFirst = floor(a + 1e-9);
beta  = a - kFirst;                      % fractional grid offset, in [0, 1)
if abs(beta) < 1e-9, beta = 0; end
kj    = kFirst + (0:nChan-1).';          % signed channel indices (column)
chanFreqs = fLow + ((1:nChan).' - 0.5) * dF;   % RF centres
fsOut = fs / D;

% ---- Prototype filter -----------------------------------------------------------------------
transW = fsOut - dF;                     % flat to dF/2, stopband from fsOut - dF/2
fcut   = fsOut / 2;
M = ceil((opts.StopbandDB - 8) / (2.285 * 2*pi * transW / fs)) + 1;
if mod(M, 2) == 0, M = M + 1; end
G = (M - 1) / 2;                         % group delay, input samples
P = ceil(M / K);                         % taps per polyphase branch
Lf = P * K;                              % frame length
h = zeros(Lf, 1);
h(1:M) = kaiserLowpass(M, fcut / fs, opts.StopbandDB);
r  = (0:Lf-1).';
hb = single(h .* exp(2i*pi*beta*r/K));   % beta modulation, phase in double
clear r

% Coverage: every channel's stopband must stay inside +-fs/2 (no wrap-around)
bbEdge = [kj(1) + beta, kj(end) + beta] * dF + [-1, 1] * (fsOut - dF/2);
if bbEdge(1) < -fs/2 || bbEdge(2) > fs/2
    error('channelizeIQ:band', ...
        'Band %.4g-%.4g Hz plus channel transitions does not fit inside fLO +- fs/2.', ...
        fLow, fHigh);
end

% ---- Input size -----------------------------------------------------------------------------------
fInfo = dir(inFile);
if isempty(fInfo)
    error('channelizeIQ:notFound', 'inFile "%s" not found.', inFile);
end
if mod(fInfo.bytes, 8) ~= 0
    warning('channelizeIQ:partial', 'inFile size is not a multiple of 8 bytes; ignoring trailing bytes.');
end
Nx = floor(fInfo.bytes / 8);
if Nx == 0
    error('channelizeIQ:empty', 'inFile "%s" is empty.', inFile);
end
nOut = ceil(Nx / D);                     % m = 0..nOut-1, m*D <= Nx-1
Mb   = min(opts.BlockOut, nOut);
nBlocks = ceil(nOut / Mb);

if opts.Verbose
    fprintf(['channelizeIQ: N = %d @ %.4g Hz -> %d channels x %d @ %.6g Hz ' ...
             '(K = %d, D = %d, oversampling %.4g), %.4f-%.4f GHz\n' ...
             '  prototype %d taps (%d per branch), flat to %.4g Hz, stop from %.4g Hz, ' ...
             '%d block(s)\n'], ...
        Nx, fs, nChan, nOut, fsOut, K, D, K/D, fLow/1e9, fHigh/1e9, ...
        M, P, dF/2, fsOut - dF/2, nBlocks);
end

% ---- Files -------------------------------------------------------------------------------------------
[fidIn, msg] = fopen(inFile, 'r', 'ieee-le');
if fidIn == -1
    error('channelizeIQ:openIn', 'Could not open "%s": %s', inFile, msg);
end
cleanupIn = onCleanup(@() fclose(fidIn)); %#ok<NASGU>

outDir = fileparts(outBase);
if ~isempty(outDir) && ~isfolder(outDir)
    mkdir(outDir);
end
chanFiles = strings(nChan, 1);
fidOut = zeros(nChan, 1);
for j = 1:nChan
    chanFiles(j) = sprintf('%s_ch%03d.dat', outBase, j);
    [fidOut(j), msg] = fopen(chanFiles(j), 'w', 'ieee-le');
    if fidOut(j) == -1
        arrayfun(@fclose, fidOut(1:j-1));
        error('channelizeIQ:openOut', 'Could not open "%s": %s', chanFiles(j), msg);
    end
end
cleanupOut = onCleanup(@() arrayfun(@fclose, fidOut)); %#ok<NASGU>

% ---- Stream -------------------------------------------------------------------------------------------
fftRows = mod(kj, K) + 1;                % rows of the K-point IFFT that are kept
offs    = (0:Lf-1).';                    % frame tap r
kb      = kj + beta;                     % channel offsets in units of fs/K

for b = 1:nBlocks
    m  = (b-1)*Mb : min(b*Mb, nOut) - 1;            % 0-based output indices
    nb = numel(m);
    c  = m*D + G;                                    % input index at tap r = 0
    i0 = c(1) - (Lf - 1);                            % first input index needed
    i1 = c(end);                                     % last
    buf = readRange(fidIn, i0, i1, Nx);              % column, zeros outside the file

    X = buf((c - i0 + 1) - offs);                    % Lf x nb frames, x(c - r)
    X = X .* hb;
    U = reshape(sum(reshape(X, K, P, nb), 2), K, nb);
    clear X
    Y = ifft(U, [], 1) * K;
    Y = Y(fftRows, :);
    cyc = mod(kb .* c, K) / K;                       % exact in double
    Y = Y .* single(exp(-2i*pi*cyc));

    for j = 1:nChan
        iq = [real(Y(j, :)); imag(Y(j, :))];
        cnt = fwrite(fidOut(j), iq, 'single');
        if cnt ~= 2*nb
            error('channelizeIQ:write', ...
                'Wrote %d of %d values to channel %d at block %d (disk full?).', cnt, 2*nb, j, b);
        end
    end

    if opts.Verbose && (b == nBlocks || mod(b, 100) == 0)
        fprintf('  block %d/%d (%.1f%%)\n', b, nBlocks, 100*(m(end)+1)/nOut);
    end
end

% ---- Info ----------------------------------------------------------------------------------------------
info = struct();
info.file             = string(outBase);
info.chanFiles        = chanFiles;
info.inFile           = inFile;
info.format           = 'cf32 interleaved (I0,Q0,I1,Q1,...), one file per channel';
info.precision        = 'single';
info.byteOrder        = 'ieee-le';
info.isComplex        = true;
info.nChan            = nChan;
info.N                = nOut;                       % complex samples per channel
info.fs               = fsOut;
info.actualFsOut      = fsOut;
info.fsIn             = fs;
info.fLOIn            = fLO;
info.chanFreqs        = chanFreqs;                  % RF centre = fLO of each channel file
info.chanWidth        = dF;
info.fLow             = fLow;
info.fHigh            = fHigh;
info.K                = K;
info.decimation       = D;
info.oversampling     = K / D;
info.freqMapping      = 'channel j: RF = chanFreqs(j) + f_bb (no inversion)';
info.loPhaseRef       = 'mixing phase 0 at input sample 0';
info.gainConvention   = 'unit DC gain (tone amplitude and passband PSD preserved)';
info.filterLen        = M;
info.prototype        = h(1:M).';                   % zero-phase taps, centre at index G+1
info.tapsPerBranch    = P;
info.cutoff           = fcut;
info.passbandEdge     = dF / 2;
info.stopbandEdge     = fsOut - dF/2;
info.stopbandDB       = opts.StopbandDB;
info.groupDelayRemoved = G;                         % input samples (G/fs seconds)
info.t0               = opts.T0;                    % sample m at t0 + m/fs
info.fullySupported   = [ceil(G/D) + 1, floor((Nx - 1 - G)/D) + 1];   % 1-based
info.samplesRead      = Nx;
info.blockOut         = Mb;
info.nBlocks          = nBlocks;
info.elapsed          = toc(tStart);

if opts.SaveInfo
    infoFile = [outBase '_info.mat'];
    save(infoFile, 'info');
    info.infoFile = infoFile;
end
if opts.Verbose
    fprintf('channelizeIQ: wrote %d channels x %d samples (%s_ch###.dat) in %.1f s\n', ...
        nChan, nOut, outBase, info.elapsed);
end
end


% =====================================================================================
function x = readRange(fid, i0, i1, Nx)
%READRANGE  Complex samples i0..i1 (0-based) as a column; zeros outside 0..Nx-1.
n  = i1 - i0 + 1;
x  = complex(zeros(n, 1, 'single'));
a  = max(i0, 0);
b  = min(i1, Nx - 1);
if b >= a
    fseek(fid, a * 8, 'bof');
    raw = fread(fid, [2, b - a + 1], 'single=>single');
    x(a - i0 + 1 : b - i0 + 1) = complex(raw(1, :), raw(2, :)).';
end
end


% =====================================================================================
function h = kaiserLowpass(M, fcNorm, A)
%KAISERLOWPASS  Odd-length Kaiser-windowed sinc low-pass, unit DC gain.
% Same design as the local function in applyIQmodulation (copied, so that
% validated stage stays untouched). fcNorm = cutoff / sample rate.
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
h = (h / sum(h)).';
end
