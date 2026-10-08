function [mask, info] = detectRFI(info_chan, opts)
%DETECTRFI  Find impulsive / pulsed RFI in the channel IQ files: a blanking mask (B2).
%{
RFI excision, detection part. Works on the raw channel files of channelizeIQ,
BEFORE dedispersion, where a radar pulse or a broadband burst is a few
samples long. The mask goes to blankChannels (zero the samples) and to
blankingWeights (exact data weights of the detected bins).

  [mask, info] = detectRFI(info_chan)
  [mask, info] = detectRFI(info_chan, 'Scales', [1 2 4 8 16], 'PFA', 1e-6)

--- Statistic --------------------------------------------------------------------------
Per channel, the power p(m) = |x(m)|^2 of the complex samples is averaged
over windows of n samples (n = each of 'Scales'; half-overlapping, stride
max(1, n/2)). A window is flagged when its mean power exceeds eta_n times
the channel's baseline power m:

  mean_n / m > eta_n,   P(mean_n / m > eta_n | Gaussian noise) = PFA.

Receiver noise (and the noise-like pulsar signal) is complex Gaussian, so a
single sample's power is exponential. The samples of a window are slightly
correlated, with autocorrelation rho(d) for white input through the
channelizer prototype h at decimation D: rho(d) = sum_r h(r) h(r + d*D) /
sum h^2 (the mixing phase drops out; rho(1) = 0.05 for the 3.125 MHz
channels at 4/3 oversampling). Then exactly

  sum of the n powers / m = sum_k lambda_k E_k,

E_k independent unit exponentials, lambda_k the eigenvalues of the n x n
correlation matrix toeplitz(rho(0..n-1)) (sum lambda_k = n; 0.61..1.05 for
n = 16). Its survival function is that of a chain of exponential stages
(phase-type), computed exactly with a matrix exponential:
  P(sum > x) = [1 0 .. 0] * expm(Q x) * ones,  Q = bidiagonal, -1/lambda_k
on the diagonal, +1/lambda_k above it,
and eta_n solves P(sum > n*eta_n) = PFA (13.82, 8.39, 5.39, 3.68, 2.69 for
n = 1, 2, 4, 8, 16 at PFA 1e-6). The Gamma distribution with the matched
variance, n_eff = n / (1 + 2 * sum_{d=1}^{n-1} (1 - d/n) rho(d)^2) (n = 16:
15.8), has the right width but a thinner tail (eigenvalues above 1): it
would give 1.02 / 1.04 / 1.08-1.10 x the intended false flags at PFA 1e-3 /
1e-4 / 1e-6 (n = 4-16). n_eff is kept in info for reference.

Short windows catch strong short RFI (radar ~+38 dB per channel, impulses);
longer windows catch weaker, longer bursts. Windows much longer than ~16
samples would also flag a bright pulsar (at -5 dB the pulse raises the
power by 32 %: a 256-sample window on the pulse peak sits at its
threshold); at the target -25 ... -54 dB the pulsar cannot trigger any
window.

--- Baseline -------------------------------------------------------------------------
m per channel and per block of 'BaselineTime' (default 16 ms; slow gain
changes), from the median: for exponential power median = m * ln 2, so
m = median / ln 2 (robust: RFI in a fraction f of the samples raises it by
~1.4 f). Only fully supported samples are used (info_chan.fullySupported).
Pass 2 ('Passes', default 2) re-estimates m without the samples flagged in
pass 1 (incl. the guard) and flags again. In a channel dominated by a
constant-envelope signal (strong carrier, GNSS) the power is not
exponential and median / ln 2 overestimates the level (up to 1/ln 2): such
channels flag less, which is intended (constant-envelope RFI is not
blanked; the 'optimal' fit gives these channels little weight).

--- Mask ---------------------------------------------------------------------------------
Every flagged window [i, i+n-1] is widened by 'Guard' samples on each side
(default: the prototype half-length in channel samples, (filterLen-1)/2/D,
rounded = 10 for the current channelizer: a short RF pulse reaches that far
into its neighbouring channel samples), then merged into intervals.

Frequency guard (8 Oct 2026, B5b): a strong event also has spectral
sidelobes in the neighbouring channels, at the same time, too weak to be
seen per sample there (a radar at +38 dB leaves ~0.3x the noise per sample
3 channels away). For RFI locked to the pulsar period even such sub-noise
leftovers add up coherently in the fold (B5b: 1 us TOA bias at -24 dB, ~150
us at -54 dB). So every flagged interval whose peak sample power is at
least 'FreqGuardMin' x baseline (default 100, +20 dB) is blanked in the same
time interval in 'FreqGuard' channels on each side. Weak flags (noise false
flags peak at ~14x, a bright pulsar too) do not spread. Default 7 (Jasper,
8 Oct): a locked +38 dB radar with 0.1 us edges then biases TOAs by 1 us
only below -51 dB (2 us at -54 dB; 5 channels: -46 dB / 6.4 us), for ~0.004 %
of the data per guard channel (runLockedRadar.m).

Uses only the channel data and the channelizer settings (no ground truth).

Inputs:
  info_chan  info of channelizeIQ (chanFiles, nChan, N, fs, prototype,
             decimation, filterLen, fullySupported).

Name-value options:
  'Scales'        window lengths n [samples] (default [1 2 4 8 16]).
  'PFA'           false-flag probability per window for Gaussian noise
                  (default 1e-6).
  'Guard'         samples added on each side of a flagged window
                  (default [] = prototype half-length, 10).
  'BaselineTime'  [s] block length of the baseline (default 16e-3).
  'Passes'        1 or more (default 2), see Baseline.
  'FreqGuard'     channels blanked on each side of a strong event, same
                  time interval (default 7; 0 = off).
  'FreqGuardMin'  peak sample power / baseline that makes an event strong
                  (default 100).
  'KeepWindows'   store the start samples of the flagged windows (final
                  pass) per channel and scale in info.flaggedWindows
                  (default false; for tests).
  'MaxMemoryGB'   limit for one channel in memory (default 4).
  'Verbose'       print a summary (default true).

Outputs:
  mask  [M x 3] rows [channel, firstSample, lastSample] (1-based, inclusive,
        channel sample grid; sorted, merged): the format of blankChannels
        and blankingWeights.
  info  scales, stride, nEff (matched-variance Gamma dof, reference), lambda
        (eigenvalues per scale: the exact distribution), threshold (eta per
        scale), pFA, guard,
        passes, baselineTime, blockEdges (0-based edges: block b = samples
        blockEdges(b)+1 .. blockEdges(b+1)), baseline [nChan x nBlocks],
        rho (lags 0..), nWindows (per scale, per channel), nFlaggedWindows
        [nChan x nScales], expectedFalse (per scale, all channels: nChan *
        nWindows * PFA, Gaussian noise), ownFlaggedSamples (per channel,
        own detections), nStrong (strong events per channel), events
        [M x 5]: every own flagged interval (before the frequency guard) as
        [channel, first, last, peakSample, peak power / baseline] (for
        periodicRFI), freqGuard,
        freqGuardMin, flaggedSamples / flaggedFraction (per channel, the
        final mask incl. the frequency guard), nIntervals, flaggedWindows
        (if KeepWindows), elapsed.
%}

arguments
    info_chan          struct
    opts.Scales        (1,:) double {mustBeInteger, mustBePositive} = [1 2 4 8 16]
    opts.PFA           (1,1) double {mustBePositive, mustBeLessThan(opts.PFA, 1)} = 1e-6
    opts.Guard               double = []
    opts.BaselineTime  (1,1) double {mustBePositive} = 16e-3
    opts.Passes        (1,1) double {mustBeInteger, mustBePositive} = 2
    opts.FreqGuard     (1,1) double {mustBeInteger, mustBeNonnegative} = 7
    opts.FreqGuardMin  (1,1) double {mustBePositive} = 100
    opts.KeepWindows   (1,1) logical = false
    opts.MaxMemoryGB   (1,1) double {mustBePositive} = 4
    opts.Verbose       (1,1) logical = true
end

tStart = tic;
nChan  = info_chan.nChan;
Nc     = info_chan.N;                               % complex samples per channel
sup    = info_chan.fullySupported;
memGB  = Nc * 48 / 1e9;                             % samples, power, cumsum, masks
if memGB > opts.MaxMemoryGB
    error('detectRFI:memory', 'One channel needs %.2f GB; use shorter data or raise MaxMemoryGB.', memGB);
end
g = opts.Guard;
if isempty(g), g = round((info_chan.filterLen - 1) / 2 / info_chan.decimation); end
if ~isscalar(g) || g < 0 || g ~= round(g)
    error('detectRFI:guard', 'Guard must be a whole number of samples >= 0.');
end

% ---- Thresholds from the channel noise statistics -----------------------------------------
rho    = channelRho(info_chan);
scales = unique(opts.Scales);
if scales(end) > Nc
    error('detectRFI:scales', 'Scales longer than the file.');
end
nS     = numel(scales);
stride = max(1, floor(scales / 2));
nEff   = zeros(1, nS);
eta    = zeros(1, nS);
nWin   = zeros(1, nS);
lambda = cell(1, nS);
for s = 1:nS
    n = scales(s);
    d = 1:min(n - 1, numel(rho) - 1);
    nEff(s) = n / (1 + 2 * sum((1 - d/n) .* rho(d + 1).'.^2));
    r = zeros(n, 1); r(1:min(n, numel(rho))) = rho(1:min(n, numel(rho)));
    lam = eig(toeplitz(r));
    lambda{s} = lam(lam > 1e-9 * max(lam));         % ~0: no contribution
    % exact threshold on the sum; the Gamma value (slightly low) brackets it
    x0 = gammaincinv(opts.PFA, nEff(s), 'upper') * n / nEff(s);
    f  = @(x) log(sumExpSurvival(lambda{s}, x)) - log(opts.PFA);
    eta(s)  = fzero(f, [0.5 * x0, 2 * x0]) / n;
    nWin(s) = numel(1:stride(s):Nc - n + 1);
end

% ---- Baseline blocks ------------------------------------------------------------------------
nBlk  = max(1, round(Nc / round(opts.BaselineTime * info_chan.fs)));
edges = round(linspace(0, Nc, nBlk + 1));
blkOf = zeros(1, Nc);
for b = 1:nBlk, blkOf(edges(b) + 1 : edges(b + 1)) = b; end
inSup = false(1, Nc); inSup(sup(1):sup(2)) = true;
minUse = 1000;                                      % samples needed for a median

baseline = zeros(nChan, nBlk);
nFlagged = zeros(nChan, nS);
ownFlagged = zeros(nChan, 1);
nStrong  = zeros(nChan, 1);
strong   = cell(nChan, 1);                          % strong intervals [first last] per channel
events   = cell(nChan, 1);                          % [channel, first, last, peakSample, peak]
rows = cell(nChan, 1);
if opts.KeepWindows, flaggedWindows = cell(nChan, nS); end

for j = 1:nChan
    p = readPower(info_chan.chanFiles(j), Nc);
    c = [0, cumsum(p)];
    covered = false(1, Nc);                         % flagged windows + guard
    mb = nan(1, nBlk);
    for pass = 1:opts.Passes
        % baseline per block: median of the (unflagged) supported samples / ln 2
        for b = 1:nBlk
            idx = edges(b) + 1 : edges(b + 1);
            use = inSup(idx) & ~covered(idx);
            if nnz(use) >= minUse
                mb(b) = median(p(idx(use))) / log(2);
            elseif pass == 1                         % short block: all its samples
                mb(b) = median(p(idx)) / log(2);
            end                                      % else: keep the previous pass
        end
        % windows of every scale against eta * baseline of the block they start in
        starts = cell(1, nS);
        for s = 1:nS
            n  = scales(s);
            i0 = 1:stride(s):Nc - n + 1;
            hit = (c(i0 + n) - c(i0)) > eta(s) * n * mb(blkOf(i0));
            starts{s} = i0(hit);
        end
        % coverage: windows widened by the guard, as a difference array
        a = []; e = [];
        for s = 1:nS
            a = [a, starts{s} - g]; %#ok<AGROW>
            e = [e, starts{s} + scales(s) - 1 + g]; %#ok<AGROW>
        end
        a = max(a, 1); e = min(e, Nc);
        dc = accumarray([a, e + 1].', [ones(1, numel(a)), -ones(1, numel(e))].', [Nc + 1, 1]);
        covered = cumsum(dc(1:Nc)).' > 0;
    end
    baseline(j, :) = mb;
    for s = 1:nS
        nFlagged(j, s) = numel(starts{s});
        if opts.KeepWindows, flaggedWindows{j, s} = starts{s}; end
    end
    ownFlagged(j) = nnz(covered);
    dcv  = diff([false, covered, false]);
    first = find(dcv == 1); last = find(dcv == -1) - 1;
    rows{j} = [j * ones(numel(first), 1), first.', last.'];
    % strong events: peak sample power / baseline of their block >= FreqGuardMin
    pk = zeros(numel(first), 1); iPk = pk;
    for i = 1:numel(first)
        [pmax, im] = max(p(first(i):last(i)));
        pk(i) = pmax / mb(blkOf(first(i)));
        iPk(i) = first(i) + im - 1;
    end
    events{j} = [j * ones(numel(first), 1), first.', last.', iPk, pk];
    isS = pk >= opts.FreqGuardMin;
    strong{j} = [first(isS).', last(isS).'];
    nStrong(j) = nnz(isS);
    if opts.Verbose && (j == nChan || mod(j, 32) == 0)
        fprintf('  detectRFI: %d/%d channels\n', j, nChan);
    end
end
% frequency guard: the strong intervals also in FreqGuard channels on each side
gRows = cell(nChan, 1);
for j = 1:nChan
    if isempty(strong{j}) || opts.FreqGuard == 0, continue; end
    nb = setdiff(max(1, j - opts.FreqGuard) : min(nChan, j + opts.FreqGuard), j);
    m  = size(strong{j}, 1);
    gRows{j} = [kron(nb(:), ones(m, 1)), repmat(strong{j}, numel(nb), 1)];
end
mask = normalizeMask([vertcat(rows{:}); vertcat(gRows{:})], nChan, Nc);
flaggedSamples = accumarray(mask(:, 1), mask(:, 3) - mask(:, 2) + 1, [nChan, 1]);

info = struct();
info.scales          = scales;
info.stride          = stride;
info.nEff            = nEff;                        % matched-variance Gamma dof (reference)
info.lambda          = lambda;                      % eigenvalues per scale (exact distribution)
info.threshold       = eta;                         % x baseline, per scale
info.pFA             = opts.PFA;
info.guard           = g;
info.passes          = opts.Passes;
info.baselineTime    = opts.BaselineTime;
info.baselineMethod  = 'median of |x|^2 per block / ln 2 (exponential power), fully supported samples';
info.blockEdges      = edges;
info.baseline        = baseline;
info.rho             = rho;
info.nWindows        = nWin;                        % per channel
info.nFlaggedWindows = nFlagged;
info.expectedFalse   = nChan * nWin * opts.PFA;     % all channels, Gaussian noise
info.ownFlaggedSamples = ownFlagged;
info.nStrong         = nStrong;
info.events          = vertcat(events{:});            % own detections, before the frequency guard
if isempty(info.events), info.events = zeros(0, 5); end
info.freqGuard       = opts.FreqGuard;
info.freqGuardMin    = opts.FreqGuardMin;
info.flaggedSamples  = flaggedSamples;
info.flaggedFraction = flaggedSamples / Nc;
info.nIntervals      = size(mask, 1);
if opts.KeepWindows, info.flaggedWindows = flaggedWindows; end
info.chanFile        = info_chan.file;
info.elapsed         = toc(tStart);

if opts.Verbose
    fprintf(['detectRFI: %d channels, windows %s samples, PFA %.1e (thresholds %s x baseline), ' ...
             'guard %d\n'], nChan, mat2str(scales), opts.PFA, mat2str(round(eta, 3)), g);
    fprintf(['  flagged windows per scale %s (Gaussian noise would give %s); %d strong events ' ...
             '(frequency guard +-%d channels); %d intervals, %.4f %% of samples in %d channels; ' ...
             '%.1f s\n'], mat2str(sum(nFlagged, 1)), mat2str(round(info.expectedFalse, 1)), ...
        sum(nStrong), opts.FreqGuard, info.nIntervals, 100*mean(info.flaggedFraction), ...
        nnz(flaggedSamples), info.elapsed);
end
end


% =====================================================================================
function rho = channelRho(info_chan)
%CHANNELRHO  Autocorrelation of the channel samples, lags d = 0.., for white input
% through the channelizer: rho(d) = sum_r h(r) h(r + d*D) / sum h^2 (real; the
% mixing phase cancels). Lags up to the prototype length.
h  = double(info_chan.prototype(:));
D  = info_chan.decimation;
nL = floor((numel(h) - 1) / D);
rho = zeros(nL + 1, 1);
for d = 0:nL
    rho(d + 1) = sum(h(1:end - d*D) .* h(1 + d*D:end));
end
rho = rho / rho(1);
end

function mask = normalizeMask(rows, nChan, Nc)
%NORMALIZEMASK  Rows [channel, first, last] -> per channel clipped to 1..Nc, sorted,
% merged where they overlap or touch (as blankChannels does), channel by channel.
out = cell(nChan, 1);
for j = unique(rows(:, 1)).'
    r = rows(rows(:, 1) == j, 2:3);
    r = [max(r(:, 1), 1), min(r(:, 2), Nc)];
    r = sortrows(r(r(:, 1) <= r(:, 2), :));
    if isempty(r), continue; end
    endMax = cummax(r(:, 2));
    newRun = [true; r(2:end, 1) > endMax(1:end-1) + 1];
    lastRow = [find(newRun(2:end)); size(r, 1)];
    out{j} = [j * ones(nnz(newRun), 1), r(newRun, 1), endMax(lastRow)];
end
mask = vertcat(out{:});
if isempty(mask), mask = zeros(0, 3); end
end

function P = sumExpSurvival(lam, x)
%SUMEXPSURVIVAL  P(sum_k lam_k E_k > x), E_k independent unit exponentials: the
% survival function of a chain of exponential stages with rates 1/lam_k (the order of
% the stages does not matter), [1 0 .. 0] * expm(Q x) * ones. Exact; no cancellation
% when eigenvalues are close (unlike the partial-fraction formula).
mu = 1 ./ lam(:);
Q  = diag(-mu);
if numel(mu) > 1, Q = Q + diag(mu(1:end-1), 1); end
E  = expm(Q * x);
P  = max(sum(E(1, :)), realmin);                    % far tail: rounding, keep log finite
end

function p = readPower(file, Nc)
%READPOWER  |x|^2 of a cf32 channel file, as a double row.
[fid, msg] = fopen(file, 'r', 'ieee-le');
if fid == -1
    error('detectRFI:open', 'Could not open "%s": %s', file, msg);
end
cleanIn = onCleanup(@() fclose(fid));
x = fread(fid, [2 Nc], 'single=>double');
if size(x, 2) < Nc
    error('detectRFI:read', 'Unexpected end of "%s".', file);
end
p = x(1, :).^2 + x(2, :).^2;
end
