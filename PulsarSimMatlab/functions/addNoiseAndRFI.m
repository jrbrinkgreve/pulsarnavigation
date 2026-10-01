function info = addNoiseAndRFI(inFile, outFile, fs, opts)
%ADDNOISEANDRFI  Add receiver noise and terrestrial interference at RF.
%{
Physical placement: the pulsar signal travels through the interstellar
medium (dispersion), then arrives at the receiver, where thermal noise of
the receiver and man-made interference are ADDED. This stage therefore runs
after applyDispersionStream and before applyIQmodulation: noise and RFI are
not dispersed. (Dedispersion later applies the inverse chirp to them,
which is exactly what happens with real data.)

  info = addNoiseAndRFI(inFile, outFile, fs, 'Band', [fLow fHigh], ...
             'SNRdB', -20, 'SignalInfo', info_gen, 'RFI', rfi, 'Seed', 1)

Inputs:
  inFile   real float32 RF file (the dispersed sky signal), little-endian
  outFile  real float32 RF file: sky signal + receiver noise + RFI
  fs       [Hz] sample rate (f_in)

Name-value options:
  'Band'        [fLow fHigh] analysis band; reference for SNR and INR
  'SNRdB'       receiver-end signal-to-noise ratio (definition below)
  'SignalInfo'  info_gen (gives the pulse amplitude A, and enables the
                SNR / TOA-precision predictions); or give 'SignalAmp' = A
  'NoiseStd'    alternative to SNRdB: receiver-noise std per RF sample
  'RFI'         array of rfiSource(...) structs (default: none)
  'Seed'        base seed for noise and RFI (default 1); independent of the
                generator seed so the pulsar realization stays the same
  'BlockSize'   samples per block (default 4e6)
  'SaveInfo', 'Verbose'

SNR definition (receiver end):
  SNRdB = 10*log10( S_peak / N ), with S_peak the pulsar power spectral
  density at the peak of the pulse and N the receiver-noise power spectral
  density, both in the analysis band. In radio-astronomy terms this is the
  peak flux density over the system-equivalent flux density, S_peak/SEFD.
  It does not depend on bandwidth, time resolution or dispersion (each
  frequency still sees the full pulse height, just at a delayed time).
  With the generator's white signal of std A (at the pulse peak) and white
  receiver noise of std sigma_n, both at f_in:  SNR = A^2 / sigma_n^2.

  What it means in practice is printed and stored in info.prediction:
  the matched-filter SNR of one pulse and of all pulses folded, and the
  best achievable TOA uncertainty per pulse and for the whole file
  (radiometer equation for a Gaussian power profile, including self-noise):
     SNR_pulse^2 = B * int( (rho*p)^2 / (1 + rho*p)^2 ) dt
     1/sigma_TOA^2 = B * int( (rho*p')^2 / (1 + rho*p)^2 ) dt
  with rho = 10^(SNRdB/10), p(t) the unit-peak power profile, B = fHigh - fLow.
  For weak pulses: SNR_pulse ~ rho*sqrt(B*sigma_t*sqrt(pi)), about 383*rho
  for the current 0.5 ms pulses in a 400 MHz band.

RFI power: see rfiSource ('INRdB' relative to the receiver-noise power in
the band). Continuous and chip-level phases are computed from the absolute
sample index, and random parts come from counter-based substreams, so the
result does not depend on the block size.
%}

arguments
    inFile                  {mustBeTextScalar}
    outFile                 {mustBeTextScalar}
    fs                (1,1) double {mustBePositive, mustBeFinite}
    opts.Band         (1,2) double {mustBePositive} = [1.2e9 1.6e9]
    opts.SNRdB              double = []
    opts.SignalAmp          double = []
    opts.SignalInfo         struct = struct([])
    opts.NoiseStd           double = []
    opts.RFI                struct = struct([])
    opts.Seed         (1,1) double {mustBeInteger, mustBeNonnegative} = 1
    opts.BlockSize    (1,1) double {mustBeInteger, mustBePositive} = 4e6
    opts.SaveInfo     (1,1) logical = true
    opts.Verbose      (1,1) logical = true
end

tStart  = tic;
inFile  = char(inFile);
outFile = char(outFile);
if strcmp(inFile, outFile)
    error('addNoiseAndRFI:sameFile', 'inFile and outFile must differ.');
end
fLow = opts.Band(1); fHigh = opts.Band(2);
if fLow >= fHigh || fHigh >= fs/2
    error('addNoiseAndRFI:band', 'Band must satisfy 0 < fLow < fHigh < fs/2.');
end
B = fHigh - fLow;

% ---- Receiver-noise level ----------------------------------------------------------
A = opts.SignalAmp;
if isempty(A) && ~isempty(opts.SignalInfo), A = opts.SignalInfo.A; end
if ~isempty(opts.NoiseStd)
    sigN = opts.NoiseStd;
    snrDB = NaN;
    if ~isempty(A) && sigN > 0, snrDB = 10*log10(A^2 / sigN^2); end
elseif ~isempty(opts.SNRdB)
    if isempty(A)
        error('addNoiseAndRFI:amp', 'SNRdB needs ''SignalInfo'' (info_gen) or ''SignalAmp''.');
    end
    snrDB = opts.SNRdB;
    sigN  = A / sqrt(10^(snrDB/10));               % Inf dB -> 0
else
    sigN = 0; snrDB = Inf;
end
PnBand = sigN^2 * B / (fs/2);                       % receiver-noise power in the band
if sigN > 0 && snrDB < -60
    warning('addNoiseAndRFI:weak', ['SNR %.1f dB: the pulse is below single-precision ' ...
        'resolution of the noise in the RF file; results may be limited by rounding.'], snrDB);
end

% ---- RFI sources ---------------------------------------------------------------------
rfi = opts.RFI(:);
nR = numel(rfi);
if nR > 0 && PnBand == 0
    error('addNoiseAndRFI:ref', 'RFI power is defined relative to receiver noise; add noise.');
end
fInfo = dir(inFile);
if isempty(fInfo)
    error('addNoiseAndRFI:notFound', 'inFile "%s" not found.', inFile);
end
Nx = floor(fInfo.bytes / 4);
tTotal = Nx / fs;

src = struct('type', {}, 'label', {}, 'amp', {}, 'r', {}, 'rs', {}, ...
             'evStart', {}, 'evLen', {}, 'def', {});
for s = 1:nR
    d = rfi(s);
    INR = 10^(d.INRdB/10);
    e = struct('type', d.type, 'label', d.label, 'amp', NaN, 'r', NaN, 'rs', [], ...
               'evStart', [], 'evLen', [], 'def', d);
    if ~strcmp(d.type, 'impulse')
        if d.freq <= 0 || d.freq >= fs/2
            error('addNoiseAndRFI:freq', 'RFI "%s": Freq must lie in (0, fs/2).', d.label);
        end
        e.amp = sqrt(2 * INR * PnBand);             % carrier: power = amp^2/2
        e.r   = d.freq / fs;                        % cycles per sample
    else
        e.amp = sqrt(INR) * sigN;                   % white burst: same in-band ratio
    end
    e.rs = RandStream('mrg32k3a', 'Seed', opts.Seed + 1000*s);
    if strcmp(d.type, 'impulse')                    % Poisson event list, whole file
        e.rs.Substream = 1;
        t = []; tc = 0;
        while true
            tc = tc - log(rand(e.rs)) / d.rate;
            if tc >= tTotal, break; end
            t(end+1) = tc; %#ok<AGROW>
        end
        e.evStart = floor(t * fs);                  % 0-based first sample
        e.evLen   = max(1, round(d.duration * fs)) * ones(size(t));
    end
    src(s) = e;
end

% ---- Predictions (what the SNR means) ------------------------------------------------
pred = struct();
if ~isempty(opts.SignalInfo) && sigN > 0 && isfinite(snrDB)
    pred = predictPerformance(opts.SignalInfo, 10^(snrDB/10), B);
end

if opts.Verbose
    fprintf('addNoiseAndRFI: N = %d @ %.4g Hz, band %.3f-%.3f GHz\n', Nx, fs, fLow/1e9, fHigh/1e9);
    if sigN > 0
        fprintf('  receiver noise: SNR = %.2f dB (peak S/N per Hz), sigma_n = %.4g per sample\n', ...
            snrDB, sigN);
        if isfield(pred, 'snrPulse')
            fprintf(['  predicted: SNR per pulse %.2f, folded over %d pulses %.1f; ' ...
                     'best TOA error per pulse %.3g us, whole file %.3g us\n'], ...
                pred.snrPulse, pred.nPulses, pred.snrFolded, ...
                pred.toaErrPulse*1e6, pred.toaErrFolded*1e6);
        end
    else
        fprintf('  receiver noise: off\n');
    end
    for s = 1:nR
        d = src(s).def;
        fprintf('  RFI %d: %-8s %-24s INR %+.1f dB', s, d.type, d.label, d.INRdB);
        if ~strcmp(d.type, 'impulse'), fprintf(' @ %.3f MHz', d.freq/1e6); end
        if strcmp(d.type, 'impulse'), fprintf(', %d bursts', numel(src(s).evStart)); end
        fprintf('\n');
    end
end

% ---- Files ---------------------------------------------------------------------------
[fidIn, msg] = fopen(inFile, 'r', 'ieee-le');
if fidIn == -1
    error('addNoiseAndRFI:openIn', 'Could not open "%s": %s', inFile, msg);
end
cleanupIn = onCleanup(@() fclose(fidIn)); %#ok<NASGU>
outDir = fileparts(outFile);
if ~isempty(outDir) && ~isfolder(outDir), mkdir(outDir); end
[fidOut, msg] = fopen(outFile, 'w', 'ieee-le');
if fidOut == -1
    error('addNoiseAndRFI:openOut', 'Could not open "%s": %s', outFile, msg);
end
cleanupOut = onCleanup(@() fclose(fidOut)); %#ok<NASGU>

% ---- Stream ----------------------------------------------------------------------------
rsN = RandStream('mt19937ar', 'Seed', opts.Seed);
blk = opts.BlockSize;
nBlocks = ceil(Nx / blk);
CH = 2^16;                                          % chips per random substream
nW = 0;
for b = 1:nBlocks
    n0 = (b - 1) * blk;
    n  = min(blk, Nx - n0);
    y  = fread(fidIn, [1 n], 'single=>single');
    if numel(y) < n
        error('addNoiseAndRFI:read', 'Unexpected end of "%s".', inFile);
    end
    if sigN > 0
        y = y + single(sigN) * randn(rsN, 1, n, 'single');
    end
    if nR > 0
        nAbs = n0 + (0:n-1);                        % absolute sample index (double)
        for s = 1:nR
            y = y + single(rfiBlock(src(s), nAbs, n0, n, fs, CH));
        end
    end
    c = fwrite(fidOut, y, 'single');
    if c ~= n
        error('addNoiseAndRFI:write', 'Wrote %d of %d samples in block %d.', c, n, b);
    end
    nW = nW + n;
    if opts.Verbose && (b == nBlocks || mod(b, max(1, round(nBlocks/10))) == 0)
        fprintf('  %3.0f%%\n', 100*nW/Nx);
    end
end

% ---- Info ------------------------------------------------------------------------------
info = struct();
info.file          = outFile;
info.inFile        = inFile;
info.precision     = 'single';
info.byteOrder     = 'ieee-le';
info.isComplex     = false;
info.N             = nW;
info.fs            = fs;
info.actualFsOut   = fs;
info.band          = [fLow fHigh];
info.noiseAddedAt  = 'RF, after dispersion (receiver input); not dispersed';
info.snrDB         = snrDB;
info.snrDefinition = 'pulse-peak signal PSD / receiver-noise PSD in the band (S_peak/SEFD)';
info.noiseStd      = sigN;
info.noisePSD      = sigN^2 / fs;                   % two-sided, per Hz
info.noisePowerBand = PnBand;
info.rfi           = rmfield(src, {'rs', 'def'});
if nR > 0
    [info.rfi.params] = src.def;
end
info.prediction    = pred;
info.seed          = opts.Seed;
info.elapsed       = toc(tStart);

if opts.SaveInfo
    [d, name] = fileparts(outFile);
    infoFile = fullfile(d, [name '_info.mat']);
    save(infoFile, 'info');
    info.infoFile = infoFile;
end
if opts.Verbose
    fprintf('addNoiseAndRFI: wrote %d samples to %s in %.1f s\n', nW, outFile, info.elapsed);
end
end


% =====================================================================================
function x = rfiBlock(e, nAbs, n0, n, fs, CH)
%RFIBLOCK  One RFI source over samples nAbs (absolute, 0-based). Double precision.
d = e.def;
t = nAbs / fs;
switch d.type
    case 'cw'
        cyc = mod(nAbs * e.r, 1) + 0.5 * d.drift * t.^2 + d.phase/(2*pi);
        x = e.amp * cos(2*pi*cyc);

    case 'bpsk'
        m  = floor(nAbs * (d.chipRate / fs));          % chip index
        ck = floor(m / CH);                             % chip chunk -> substream
        uc = unique(ck);
        vals = zeros(CH, numel(uc));
        for i = 1:numel(uc)
            e.rs.Substream = uc(i) + 1;
            vals(:, i) = 2*randi(e.rs, [0 1], CH, 1) - 1;
        end
        [~, col] = ismember(ck, uc);
        chips = vals(sub2ind([CH numel(uc)], m - ck*CH + 1, col));
        cyc = mod(nAbs * e.r, 1) + d.phase/(2*pi);
        x = e.amp * chips .* cos(2*pi*cyc);

    case 'pulsed'
        pri = 1 / d.prf;
        tin = mod(t - d.startTime, pri);                 % time since pulse start
        on  = find(tin < d.pulseWidth);
        x = zeros(1, n);
        if ~isempty(on)
            k   = d.chirpBW / d.pulseWidth;              % chirp rate [Hz/s]
            cyc = mod(nAbs(on) * e.r, 1) + 0.5 * k * (tin(on) - d.pulseWidth/2).^2;
            x(on) = e.amp * cos(2*pi*cyc);
        end

    case 'impulse'
        x = zeros(1, n);
        hit = find(e.evStart < n0 + n & e.evStart + e.evLen > n0);
        for i = hit(:).'
            e.rs.Substream = i + 1;                      % burst i: its own substream
            w  = randn(e.rs, 1, e.evLen(i));
            a0 = max(n0, e.evStart(i));
            a1 = min(n0 + n, e.evStart(i) + e.evLen(i)) - 1;
            x(a0 - n0 + 1 : a1 - n0 + 1) = x(a0 - n0 + 1 : a1 - n0 + 1) + ...
                e.amp * w(a0 - e.evStart(i) + 1 : a1 - e.evStart(i) + 1);
        end
end
end


% =====================================================================================
function pred = predictPerformance(g, rho, B)
%PREDICTPERFORMANCE  Radiometer-equation predictions for a Gaussian power profile.
if strcmpi(g.envelopeMode, 'power')
    sp = g.sigma;
else
    sp = g.sigma / sqrt(2);
end
t  = linspace(-8*sp, 8*sp, 40001);
p  = exp(-0.5*(t/sp).^2);                            % unit-peak power profile
dp = -t/sp^2 .* p;
den = (1 + rho*p).^2;
snr2  = B * trapz(t, (rho*p).^2 ./ den);
fish  = B * trapz(t, (rho*dp).^2 ./ den);
nP = numel(g.pulseCenters);
pred = struct('rho', rho, 'snrPulse', sqrt(snr2), 'nPulses', nP, ...
    'snrFolded', sqrt(snr2*nP), 'toaErrPulse', 1/sqrt(fish), ...
    'toaErrFolded', 1/sqrt(fish*nP), 'bandwidth', B, ...
    'note', 'optimal estimator, flat band B, Gaussian power profile, incl. self-noise');
end