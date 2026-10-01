function info = generatePulsarSignal(outFile, T, f_in, A, L, dutycycle, opts)
%GENERATEPULSARSIGNAL  Synthetic pulsar signal, streamed to disk.
%{
Gaussian white noise gated by a periodic, unit-peak Gaussian envelope
(one bump per rotation period). Output is always streamed to outFile in
blocks, matching the file-in/file-out style of the rest of the pipeline.

  info = generatePulsarSignal(outFile, T, f_in, A, L)
  info = generatePulsarSignal(outFile, T, f_in, A, L, dutycycle)
  info = generatePulsarSignal(..., 'Seed', 42, 'EnvelopeMode', 'power', ...)

Inputs:
  outFile    path of the binary output file. Written as real-valued
             float32, little-endian ('ieee-le'), no header. The parent
             folder is created if it does not exist.
  T          [s]   pulsar rotation period
  f_in       [Hz]  sampling frequency
  A          [V]   noise std dev at the peak of the envelope
  L          [s]   total signal length
  dutycycle  [%]   pulse FWHM as percent of T (default 5). Which profile
                   the FWHM refers to is set by EnvelopeMode.

Name-value options:
  'Seed'          RNG seed (integer in [0, 2^32-1]). If omitted, a seed
                  is derived from the clock and reported in info.seed, so
                  every run is reproducible after the fact. Uses a private
                  RandStream; the global RNG state is not touched.
  'EnvelopeMode'  'power' (default): the POWER profile (x^2) has
                  FWHM = dutycycle% of T, which is what a power-based
                  envelope detector recovers.
                  'amplitude': the AMPLITUDE envelope has that FWHM
                  (power profile is then narrower by sqrt(2)); this
                  matches the behaviour of the old version.
  'BlockSize'     samples per processed block (default 4e6). Peak RAM is
                  roughly 4 double arrays of this length (~130 MB).
  'SaveInfo'      true (default): also save info to <outFile>_info.mat,
                  so the ground truth travels with the data file.
  'Verbose'       true to print progress (default false).

Output:
  info  struct with file format, sample count, sampling rate, envelope
        parameters, ground-truth pulse centre times, seed, and bytes
        written. info.actualFsOut mirrors the field name used by the
        downstream modules so the info structs can be chained.

Envelope notes:
  Pulse centres are at (k+0.5)*T, so the first pulse peaks at T/2 and
  the signal starts near zero. The envelope is a periodic sum of
  Gaussians, evaluated per sample over all pulses within ~6 sigma, and
  normalised by its exact peak value, so overlapping pulses (large
  dutycycle) keep their true shape instead of being clipped. Time is
  computed in double precision; only the final samples are single.

Reading the file back (small tests only):
  fid = fopen(info.file, 'r', info.byteOrder);
  x = fread(fid, [1 info.N], '*single'); fclose(fid);
  t = (0:info.N-1) / info.fs;
%}

arguments
    outFile                {mustBeTextScalar}
    T                (1,1) double {mustBePositive, mustBeFinite}
    f_in             (1,1) double {mustBePositive, mustBeFinite}
    A                (1,1) double {mustBeNonnegative, mustBeFinite}
    L                (1,1) double {mustBePositive, mustBeFinite}
    dutycycle        (1,1) double {mustBePositive, mustBeFinite} = 5
    opts.Seed              double = []
    opts.EnvelopeMode      {mustBeTextScalar} = 'power'
    opts.BlockSize   (1,1) double {mustBeInteger, mustBePositive} = 4e6
    opts.SaveInfo    (1,1) logical = true
    opts.Verbose     (1,1) logical = false
end

tStart  = tic;
outFile = char(outFile);
mode    = lower(char(opts.EnvelopeMode));
if ~any(strcmp(mode, {'power', 'amplitude'}))
    error('generatePulsarSignal:badMode', ...
        'EnvelopeMode must be ''power'' or ''amplitude'', got ''%s''.', mode);
end
powerMode = strcmp(mode, 'power');

% ---------------------------------------------------------------------
% Envelope geometry
% G(t) is a periodic sum of unit-peak Gaussians with std dev sigma.
%   power mode:     x = A*sqrt(G)*n  ->  power profile ~ G, FWHM = FWHM
%   amplitude mode: x = A*G*n        ->  amplitude envelope FWHM = FWHM
% ---------------------------------------------------------------------
FWHM  = (dutycycle/100) * T;
sigma = FWHM / (2*sqrt(2*log(2)));
if powerMode
    FWHM_power     = FWHM;
    FWHM_amplitude = FWHM * sqrt(2);   % sqrt of a Gaussian is sqrt(2) wider
else
    FWHM_amplitude = FWHM;
    FWHM_power     = FWHM / sqrt(2);
end
if dutycycle > 50
    warning('generatePulsarSignal:overlap', ...
        'dutycycle = %g%%: pulses overlap strongly; envelope never nears zero.', ...
        dutycycle);
end

% Neighbouring pulses to include per sample. Each sample is assigned its
% nearest pulse k0, which is at most T/2 away; pulse k0+j is then at
% least (|j|-0.5)*T away, so |j| <= 6*sigma/T + 0.5 covers everything
% within 6 sigma.
nNeighbors = max(1, ceil(6*sigma/T + 0.5));
kk         = -nNeighbors:nNeighbors;
peakNorm   = sum(exp(-0.5*(kk*T/sigma).^2));   % exact peak of periodic sum

% ---------------------------------------------------------------------
% Sizes and ground truth
% ---------------------------------------------------------------------
N    = round(L*f_in) + 1;             % inclusive of t = L
tEnd = (N-1) / f_in;
if N > flintmax
    error('generatePulsarSignal:tooLong', 'N = %g exceeds exact double range.', N);
end
kMax         = floor(tEnd/T - 0.5);
pulseCenters = ((0:max(kMax, -1)) + 0.5) * T;  % empty if no full centre fits

% ---------------------------------------------------------------------
% Random stream (private, reproducible)
% ---------------------------------------------------------------------
if isempty(opts.Seed)
    seed = mod(floor(posixtime(datetime('now')) * 1e6), 2^32);
else
    seed = opts.Seed;
    if ~isscalar(seed) || seed ~= floor(seed) || seed < 0 || seed >= 2^32
        error('generatePulsarSignal:badSeed', ...
            'Seed must be an integer in [0, 2^32-1].');
    end
end
rngType = 'mt19937ar';
rs = RandStream(rngType, 'Seed', seed);

% ---------------------------------------------------------------------
% Output file (explicit little-endian)
% ---------------------------------------------------------------------
byteOrder = 'ieee-le';
outDir = fileparts(outFile);
if ~isempty(outDir) && ~isfolder(outDir)
    mkdir(outDir);
end
[fid, msg] = fopen(outFile, 'w', byteOrder);
if fid == -1
    error('generatePulsarSignal:openFailed', ...
        'Could not open "%s" for writing: %s', outFile, msg);
end
cleanupObj = onCleanup(@() fclose(fid)); %#ok<NASGU>

% ---------------------------------------------------------------------
% Block-wise generation
% ---------------------------------------------------------------------
blockSize = opts.BlockSize;
nBlocks   = ceil(N / blockSize);
nWritten  = 0;
Asingle   = single(A);

for b = 1:nBlocks
    i0 = (b-1) * blockSize;             % zero-based index of first sample
    n  = min(blockSize, N - i0);
    tt = (i0 : i0+n-1) / f_in;          % double: sub-ns accuracy even at large t

    k0 = round(tt/T - 0.5);             % nearest pulse index, per sample
    G  = zeros(1, n);
    for j = kk
        delta = tt - (k0 + j + 0.5)*T;
        G = G + exp(-0.5*(delta/sigma).^2);
    end
    G = G / peakNorm;                   % peak exactly 1, shape preserved

    if powerMode
        env = single(sqrt(G));
    else
        env = single(G);
    end

    block = Asingle * env .* randn(rs, 1, n, 'single');

    cnt = fwrite(fid, block, 'single');
    if cnt ~= n
        error('generatePulsarSignal:writeFailed', ...
            'Wrote %d of %d samples in block %d (disk full?).', cnt, n, b);
    end
    nWritten = nWritten + cnt;

    if opts.Verbose && (b == nBlocks || mod(b, max(1, round(nBlocks/10))) == 0)
        fprintf('generatePulsarSignal: block %d/%d (%.1f%%)\n', ...
            b, nBlocks, 100*nWritten/N);
    end
end

% ---------------------------------------------------------------------
% Info struct
% ---------------------------------------------------------------------
info = struct();
info.file           = outFile;
info.precision      = 'single';
info.byteOrder      = byteOrder;
info.isComplex      = false;
info.N              = N;
info.fs             = f_in;
info.actualFsOut    = f_in;           % same name as downstream modules
info.duration       = tEnd;
info.T              = T;
info.A              = A;
info.L              = L;
info.dutycycle      = dutycycle;
info.envelopeMode   = mode;
info.sigma          = sigma;          % std dev of G (power profile in 'power' mode)
info.FWHM_power     = FWHM_power;
info.FWHM_amplitude = FWHM_amplitude;
info.pulseCenters   = pulseCenters;   % ground-truth peak times [s]
info.nPulses        = numel(pulseCenters);
info.seed           = seed;
info.rngType        = rngType;
info.blockSize      = blockSize;
info.samplesWritten = nWritten;
info.bytesWritten   = nWritten * 4;
info.elapsed        = toc(tStart);

if opts.SaveInfo
    [d, name] = fileparts(outFile);
    infoFile = fullfile(d, [name '_info.mat']);
    save(infoFile, 'info');
    info.infoFile = infoFile;
end
end