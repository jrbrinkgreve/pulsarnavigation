function info = blankingWeights(info_chan, info_dc, info_det, mask, outFile, opts)
%BLANKINGWEIGHTS  Exact data weights W, V, X(L) of detected time bins from a blanking mask.
%{
RFI excision sets channel samples x(n) to zero BEFORE per-channel coherent
dedispersion. A dedispersed sample y(n) = sum_l h(l) x(n-l) (h: the
channel's dedispersion filter, including the delay to the reference
frequency) then carries a fraction of valid data, and the detected power
of time bin k (mean of |y|^2 over n samples) has, for Gaussian x (receiver
noise and the noise-like pulsar signal alike):
  mean            W_k * m                     (m = unblanked mean power)
  variance        V_k * sigma^2
  cov with k+L    X_k(L) * sigma^2,  L = 1..Lmax
sigma^2 = m^2 * radiometer, radiometer = 1/(Bnoise*binDt) (detectChannels
info.noise), the normalization foldProfile / estimateTOA / detectPulsar use.
With C(a, a+tau) = sum_l h(l) conj(h(l+tau)) keep(a-l) the covariance of
the blanked dedispersed samples (unit-variance white input):
  w(a)   = C(a, a) / m0,  m0 = sum |h|^2;      W_k = mean of w over bin k
  V_k    = (1/n^2) sum_{a,b in k} |C_ab|^2 / m0^2 / radiometer
  X_k(L) = (1/n^2) sum_{a in k, b in k+L} |C_ab|^2 / m0^2 / radiometer
C(., . + tau) is the mask convolved with q_tau(l) = h(l) conj(h(l+tau)): one
FFT convolution per sample lag tau = 0 .. (Lmax+1)*n - 1 (experiment 3a,
tests/expBlankingVariance.m; generalized here to lags L = 1..Lmax).

The channel samples are not white: the channelizer output is flat over the
useful band but rolls off between its edge and Nyquist (spectrum S_x =
|prototype response|^2 aliased at the channel rate; autocorrelation R_x,
significant up to ~20 lags). H is zero outside the useful band, so without
blanking only the flat part matters, but blanking mixes the roll-off region
into the band. The MEAN is therefore computed with the true R_x:
  E|y(a)|^2 = sum_d R_x(d) sum_l h(l) conj(h(l+d)) keep(a-l) keep(a-l-d)
            = R_x(0) (q_0 * keep)(a) + 2 Re sum_{d>0} R_x(d) (q_d * keep.keep_d)(a),
  w(a) = E|y(a)|^2 / (the same without blanking),
one extra FFT of the mask product per lag d (the q_d are the q_tau above).
White input (R_x = delta) is the simple form above; the true spectrum lowers
the power of nearly empty bins (valid fraction < 1 %) by ~4 % (A3a). V and X
keep the white form (error bars only; ~1-2 % in heavily blanked bins).

  info = blankingWeights(info_chan, info_dc, info_det, mask, outFile)

Inputs:
  info_chan info of channelizeIQ (prototype, decimation, fs, fsIn): the
            channel noise spectrum. Not used with 'InputSpectrum' 'white'.
  info_dc   info of dedisperseChannels (chanFreqs, chanWidth, fs, N, DM,
            refFreq, edgeWidth, nPast, nFuture, Nfft): h of every channel
            is taken from applyInverseDispersion itself (an impulse through
            the same filter, same Nfft).
  info_det  info of detectChannels (N bins, binLen, noise.Lmax,
            noise.radiometer): bin k = channel samples (k-1)*binLen+1 .. k*binLen.
  mask      [M x 3] blanked intervals [channel, firstSample, lastSample]
            (1-based, inclusive, on the channel sample grid of info_dc;
            clipped to the file; overlaps allowed). [] = nothing blanked.
  outFile   data-weight file for foldProfile 'DataWeights': float32
            [nChan x (2+Lmax) x nBins], per bin W of all channels, then V,
            then X(1) .. X(Lmax).

Name-value options:
  'InputSpectrum' 'channelizer' (default): the mean with the true channel
                spectrum from info_chan; 'white': white channel input (the
                simple form; for comparison).
  'SpectrumTol' lags d with |R_x(d)| / R_x(0) below this are dropped
                (default 1e-9; ~20 lags remain).
  'MinWeight'   W below this is set to 0 (with V, X): such a bin holds no
                usable data, and FFT rounding / filter tails would leave tiny
                positive values (default 1e-6).
  'MaxMemoryGB' limit for the weight array held in memory (default 4).
  'SaveInfo'    save info to <outFile>_info.mat (default true).
  'Verbose'     print a summary (default true).

Output info: the fields foldProfile needs (file, nChan, N, Lmax, byteOrder)
  plus binLen, blankedFraction (per channel, samples), validFraction (per
  channel, mean W over info_det.fullySupportedBins), V0 / X0 (per channel:
  the unblanked values from h), nMaskRows, minWeight, inputSpectrum,
  Rx (normalized autocorrelation used for the mean, lags 0..), elapsed.
Unblanked channels get W = 1 and their unblanked V0, X0 in every bin.
Bins inside info_det.fullySupportedBins are exact; outside them blanked
channels include the zero padding at the file ends, unblanked channels
not (the fold uses only fully supported bins).
%}

arguments
    info_chan         struct
    info_dc           struct
    info_det          struct
    mask              double
    outFile                 {mustBeTextScalar}
    opts.InputSpectrum      {mustBeTextScalar} = 'channelizer'
    opts.SpectrumTol  (1,1) double {mustBePositive} = 1e-9
    opts.MinWeight    (1,1) double {mustBeNonnegative} = 1e-6
    opts.MaxMemoryGB  (1,1) double {mustBePositive} = 4
    opts.SaveInfo     (1,1) logical = true
    opts.Verbose      (1,1) logical = true
end

tStart  = tic;
outFile = char(outFile);
nChan = info_dc.nChan;
Nc    = info_dc.N;                                  % samples per channel
n     = info_det.binLen;
nb    = info_det.N;                                 % detected time bins
Lmax  = info_det.noise.Lmax;
rad   = info_det.noise.radiometer;
nQ    = 2 + Lmax;
if info_det.nChan ~= nChan
    error('blankingWeights:chan', 'info_dc and info_det have different channel counts.');
end
if nb * n > Nc
    error('blankingWeights:bins', 'info_det has more bins than the channel files hold.');
end
if isempty(mask), mask = zeros(0, 3); end
if size(mask, 2) ~= 3 || any(mask(:, 1) < 1 | mask(:, 1) > nChan | mask(:, 1) ~= round(mask(:, 1)))
    error('blankingWeights:mask', 'mask must be [channel, first, last] rows with valid channels.');
end
memGB = nChan * nQ * nb * 4 / 1e9;
if memGB > opts.MaxMemoryGB
    error('blankingWeights:memory', ['Weight array needs %.2f GB; use shorter data or raise ' ...
        'MaxMemoryGB.'], memGB);
end

[d0, name] = fileparts(outFile);
if ~isempty(d0) && ~isfolder(d0), mkdir(d0); end
impIn  = fullfile(tempdir, sprintf('%s_imp_in.dat', name));
impOut = fullfile(tempdir, sprintf('%s_imp_out.dat', name));
cleanImp = onCleanup(@() deleteFiles({impIn, impOut}));

Y = zeros(nChan, nQ, nb, 'single');                 % the weight file, in memory
blankedFraction = zeros(nChan, 1);
validFraction   = zeros(nChan, 1);
V0 = zeros(nChan, 1); X0 = zeros(nChan, Lmax);
sb = info_det.fullySupportedBins;
nTau = (Lmax + 1) * n;                              % sample lags 0 .. nTau-1
dF = info_dc.chanWidth;

% ---- Channel noise autocorrelation for the mean ---------------------------------------
spec = lower(char(opts.InputSpectrum));
switch spec
    case 'white'
        Rx = 1;
    case 'channelizer'
        Rx = channelAutocorr(info_chan, info_dc, opts.SpectrumTol);
    otherwise
        error('blankingWeights:spectrum', 'InputSpectrum must be ''channelizer'' or ''white''.');
end
nR = numel(Rx) - 1;                                 % lags d = 0..nR
wR = [Rx(1); 2*Rx(2:end)];                          % d > 0 counted twice (d and -d)

for j = 1:nChan
    % ---- h of this channel: an impulse through the same filter (same Nfft) -----------
    nF = info_dc.nFuture(j); nP = info_dc.nPast(j);
    Nimp = nF + nP + 1;
    imp = zeros(Nimp, 1); imp(nF + 1) = 1;          % y(nF+1+l) = h(l), l = -nF..nP
    writeCF32(impIn, imp);
    fc = info_dc.chanFreqs(j);
    d = applyInverseDispersion(impIn, impOut, info_dc.fs, fc, info_dc.DM, fc - dF/2, fc + dF/2, ...
        'RefFreq', info_dc.refFreq, 'AllowRefOutsideBand', true, ...
        'EdgeFrac', info_dc.edgeWidth / dF, 'Nfft', info_dc.Nfft(j), ...
        'SaveInfo', false, 'Verbose', false);
    if d.nPast ~= nP || d.nFuture ~= nF
        error('blankingWeights:filter', ['Channel %d: filter reach %d/%d differs from ' ...
            'dedisperseChannels (%d/%d); other GuardTime?'], j, d.nPast, d.nFuture, nP, nF);
    end
    h  = double(readCF32(impOut));                  % h(l), l = -nF..nP
    m0 = sum(abs(h).^2);

    % unblanked values: C_tau = sum_l h(l) conj(h(l+tau)), the same in every bin
    cU = zeros(nTau, 1);
    for tau = 0:nTau-1
        cU(tau + 1) = sum(h(1:end-tau) .* conj(h(1+tau:end))) / m0;
    end
    [V0(j), X0(j, :)] = binSums(reshape(repmat(abs(cU.').^2, n, 1), n, 1, nTau), n, Lmax, 1);
    V0(j) = V0(j) / rad; X0(j, :) = X0(j, :) / rad;
    cR = zeros(nR + 1, 1);                          % unblanked mean power with R_x
    for dd = 0:nR
        cR(dd + 1) = sum(h(1:end-dd) .* conj(h(1+dd:end)));
    end
    den = real(sum(wR .* cR));

    % ---- mask of this channel ----------------------------------------------------------
    rows = mask(mask(:, 1) == j, 2:3);
    keep = true(Nc, 1);
    for r = 1:size(rows, 1)
        a = max(1, rows(r, 1)); b = min(Nc, rows(r, 2));
        if a <= b, keep(a:b) = false; end
    end
    blankedFraction(j) = mean(~keep);
    if all(keep)
        Y(j, 1, :) = 1;
        Y(j, 2, :) = V0(j);
        for L = 1:Lmax, Y(j, 2 + L, :) = X0(j, L); end
        validFraction(j) = 1;
        continue
    end

    % ---- exact statistics of every bin: one FFT convolution of the mask per lag ------
    % (q_tau * keep)(a) = sum_l q_tau(l) keep(a-l): keep at index a-1, q at l+nF
    % (0-based) -> result index a-1+nF; zero padding = no data outside the file
    M  = 2^nextpow2(Nc + Nimp);
    Kf = fft(double(keep), M);
    A  = zeros(n, nb, nTau);                        % |C(a, a+tau)|^2 / m0^2, a in bin k
    aIdx = nF + (1:nb*n).';                         % 1-based index of a = 1..nb*n
    Acc = zeros(M, 1);                              % mean: sum_d wR(d) F(keep.keep_d) F(q_d)
    for tau = 0:max(nTau, nR + 1) - 1
        q  = h(1:end-tau) .* conj(h(1+tau:end));    % l = -nF .. nP-tau
        Qf = fft(q, M);
        if tau <= nR
            if tau == 0
                Acc = Acc + wR(1) * (Kf .* Qf);
            else                                    % keep(s) keep(s-d)
                kd = double(keep & [false(tau, 1); keep(1:end-tau)]);
                Acc = Acc + wR(tau + 1) * (fft(kd, M) .* Qf);
            end
        end
        if tau < nTau
            Ct = ifft(Kf .* Qf);
            A(:, :, tau + 1) = reshape(abs(Ct(aIdx)).^2 / m0^2, n, nb);
        end
    end
    num = real(ifft(Acc));
    w = num(aIdx) / den;                            % valid fraction of each sample
    W = mean(reshape(w, n, nb), 1).';
    [V, X] = binSums(A, n, Lmax, nb);
    V = V / rad; X = X / rad;
    for L = 1:Lmax
        X(nb-L+1:nb, L) = 0;                        % no partner bin k+L
    end
    zero = W < opts.MinWeight;                      % no usable data
    W(zero) = 0; V(zero) = 0; X(zero, :) = 0;
    for L = 1:Lmax                                  % and no covariance with an empty bin k+L
        X([zero(L+1:end); false(L, 1)], L) = 0;
    end
    Y(j, 1, :) = W;
    Y(j, 2, :) = V;
    for L = 1:Lmax, Y(j, 2 + L, :) = X(:, L); end
    validFraction(j) = mean(W(sb(1):sb(2)));
    if opts.Verbose && mod(j, 32) == 0
        fprintf('  blankingWeights: %d/%d channels\n', j, nChan);
    end
end

% ---- Write ---------------------------------------------------------------------------------
[fid, msg] = fopen(outFile, 'w', 'ieee-le');
if fid == -1
    error('blankingWeights:open', 'Could not open "%s": %s', outFile, msg);
end
cleanOut = onCleanup(@() fclose(fid));
c = fwrite(fid, Y, 'single');
if c ~= numel(Y)
    error('blankingWeights:write', 'Wrote %d of %d values (disk full?).', c, numel(Y));
end
clear cleanOut

info = struct();
info.file            = outFile;
info.format          = 'float32 [nChan x (2+Lmax) x N]: per bin W (all channels), V, X(1..Lmax)';
info.byteOrder       = 'ieee-le';
info.nChan           = nChan;
info.N               = nb;
info.Lmax            = Lmax;
info.binLen          = n;
info.powerFile       = info_det.file;
info.normalization   = 'V, X relative to m^2 * radiometer (info_det.noise.radiometer)';
info.blankedFraction = blankedFraction;
info.validFraction   = validFraction;
info.V0              = V0;
info.X0              = X0;
info.nMaskRows       = size(mask, 1);
info.minWeight       = opts.MinWeight;
info.inputSpectrum   = spec;
info.Rx              = Rx;                          % lags 0..nR, used for the mean
info.elapsed         = toc(tStart);
if opts.SaveInfo
    infoFile = fullfile(d0, [name '_info.mat']);
    save(infoFile, 'info');
    info.infoFile = infoFile;
end
if opts.Verbose
    nB = nnz(blankedFraction > 0);
    fprintf(['blankingWeights: %d of %d channels blanked (%.3f %% of samples), valid fraction ' ...
             '%.4f (fully supported bins); unblanked V %.4f, X(1) %.4f; %.1f s\n'], nB, nChan, ...
        100*mean(blankedFraction), mean(validFraction), mean(V0), mean(X0(:, 1)), info.elapsed);
end
end


% =========================================================================================
function Rx = channelAutocorr(info_chan, info_dc, tol)
%CHANNELAUTOCORR  Autocorrelation R_x(d), d = 0..nR, of the channel samples for white
% input through the channelizer: S_x = |prototype response|^2 aliased at the channel
% rate (decimation D), normalized to 1 in the useful band; R_x = its inverse FFT,
% cut where |R_x(d)| < tol * R_x(0).
if abs(info_chan.fs - info_dc.fs) > 1e-6 * info_dc.fs || info_chan.nChan ~= info_dc.nChan
    error('blankingWeights:chanInfo', 'info_chan does not belong to info_dc.');
end
p  = double(info_chan.prototype(:));
D  = info_chan.decimation;
Mc = 2^12;                                          % channel-rate frequency grid
Sx = sum(reshape(abs(fft(p, Mc * D)).^2, Mc, D), 2);  % input bin k -> channel bin mod(k, Mc)
f  = (0:Mc-1).' * info_chan.fs / Mc;
f(f >= info_chan.fs/2) = f(f >= info_chan.fs/2) - info_chan.fs;
Sx = Sx / mean(Sx(abs(f) <= 0.45 * info_chan.chanWidth));
R  = real(ifft(Sx));
nR = find(abs(R(1:Mc/2)) > tol * R(1), 1, 'last') - 1;
Rx = R(1:nR + 1);
end

function [V, X] = binSums(A, n, Lmax, nb)
%BINSUMS  Sum |C(a, a+tau)|^2 over the sample pairs of bin k (V) and of bins k, k+L (X).
% A(i+1, k, tau+1) = |C(a, a+tau)|^2 for sample a = (k-1)*n + 1 + i of bin k.
% Pairs inside bin k: b = a + tau with i + tau <= n-1, counted twice for tau > 0.
% Pairs with bin k+L: b in bin k+L, i.e. i + tau in [L*n, L*n + n-1].
V = zeros(nb, 1); X = zeros(nb, Lmax);
for tau = 0:size(A, 3) - 1
    B = A(:, :, tau + 1);
    if tau < n
        V = V + (1 + (tau > 0)) * sum(B(1:n-tau, :), 1).';
    end
    for L = 1:Lmax
        i = max(0, L*n - tau) : min(n-1, L*n + n-1 - tau);
        if ~isempty(i)
            X(:, L) = X(:, L) + sum(B(i + 1, :), 1).';
        end
    end
end
V = V / n^2; X = X / n^2;
end

function writeCF32(file, x)
fid = fopen(file, 'w', 'ieee-le');
fwrite(fid, [real(x(:)).'; imag(x(:)).'], 'single');
fclose(fid);
end

function y = readCF32(file)
fid = fopen(file, 'r', 'ieee-le');
c = onCleanup(@() fclose(fid));
raw = fread(fid, [2 Inf], 'single=>single');
y = complex(raw(1, :), raw(2, :)).';
end

function deleteFiles(files)
for i = 1:numel(files)
    if isfile(files{i}), delete(files{i}); end
end
end
