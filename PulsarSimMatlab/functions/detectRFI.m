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
Limit: valid for a pulsar at most ~0.3 x the receiver noise (SNRdB <= -5 dB).
A brighter pulse blanks itself: on the pulse peak 10 % of the 16-sample
windows are flagged at 0 dB, 63 % at +3 dB, almost all from +5 dB, and from
+10 dB single samples pass 100 x the baseline, so the frequency guard
blanks neighbouring channels too. Noise-free data (no receiver noise) is
the extreme case: the baseline is ~0 off the pulse and every on-pulse
window is flagged. The mask then depends on the pulse itself, which
blankingWeights does not model: the TOA stays unbiased by symmetry, but
amplitude, SNR and red. chi^2 are wrong (deep analysis 13 Oct 2026, s. 4).

--- Per-channel thresholds (21 Oct 2026) ----------------------------------------------
The thresholds above assume a channel holds white receiver noise (correlation rho(d) of
the channelizer only). Continuous noise-like RFI narrower than the channel (a 200 kHz
carrier in a 3.125 MHz channel, the edge of an LTE block) is still complex Gaussian - its
power is still exponential per sample - but correlated over ~20 samples, so a window's
mean fluctuates far more and the white thresholds flag its own power peaks (narrow
Inmarsat on the voltage chain: ~1400 false windows per carrier channel per second
against ~2 expected). Those false flags cut the peaks, and the valid samples next to a
cut keep the peak's correlated shoulders: the phase bins with blanking come out 4-9 %
noisier than any model of independent blanking (testRFIChain, 21 Oct). So from pass 2 on
('ChannelThresholds', default true) each channel's own complex autocorrelation
r_c(d) = <x(m) x*(m+d)> / <|x|^2>, d = 0 .. max(Scales)-1, is measured on all supported
samples except the previous pass's strong flagged intervals (peak sample >= 30 x the
baseline: pulsed RFI; Gaussian noise of any spectrum has exponential power and reaches
that with probability e^-30 per sample, so the noise's own peaks stay in and the estimate
is unbiased), and its thresholds come from the eigenvalues of the Hermitian
toeplitz(r_c) exactly as above: the channel's continuous RFI counts as part of its noise,
and only impulsive / pulsed RFI stands out. Such a channel also takes its baseline from
those samples (leaving out its weaker flags, which are its own noise peaks, would bias the
baseline low and double the false flags at PFA 1e-3). A channel keeps the white thresholds unless its own move one by more than
1 % (estimation noise moves them by ~1e-4): clean channels are bit for bit as before.

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
  'ChannelThresholds'  thresholds from each channel's own correlation from pass 2
                  on (default true; false = the white thresholds everywhere, as
                  before 21 Oct 2026), see Per-channel thresholds.
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
        (if KeepWindows), channelThresholds, ownThreshold [nChan x 1] (channel used
        its own thresholds in the final pass), thresholdUsed [nChan x nScales],
        chanRho [nChan x max(Scales)] (r_c, lags 0.., last estimate; NaN if not
        measured), elapsed.
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
    opts.ChannelThresholds (1,1) logical = true
    opts.KeepWindows   (1,1) logical = false
    opts.MaxMemoryGB   (1,1) double {mustBePositive} = 4
    opts.Verbose       (1,1) logical = true
end

tStart = tic;
nChan  = info_chan.nChan;
Nc     = info_chan.N;                               % complex samples per channel
sup    = info_chan.fullySupported;
memGB  = Nc * (48 + 16 * opts.ChannelThresholds) / 1e9;   % samples, power, cumsum, masks
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
[eta, lambda, nEff] = windowThresholds(rho, scales, opts.PFA);
nWin   = zeros(1, nS);
for s = 1:nS
    nWin(s) = numel(1:stride(s):Nc - scales(s) + 1);
end
maxLag = scales(end) - 1;                           % lags of a channel's own correlation
rcPeakMin = 30;                                     % peak / baseline of intervals left out of it

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
ownThr  = false(nChan, 1);
thrUsed = repmat(eta, nChan, 1);
chanRho = nan(nChan, maxLag + 1);

for j = 1:nChan
    [p, x] = readPower(info_chan.chanFiles(j), Nc, opts.ChannelThresholds);
    c = [0, cumsum(p)];
    covered = false(1, Nc);                         % flagged windows + guard
    mb = nan(1, nBlk);
    etaJ = eta;
    for pass = 1:opts.Passes
        % from pass 2: thresholds from the channel's own correlation, measured on the samples
        % outside the previous pass's strong intervals (pulsed RFI); kept white unless one
        % moves by > 1 %. A channel with its own thresholds also takes its baseline from those
        % samples: its weaker flags are its own noise peaks, and leaving them out would bias
        % the baseline low (narrow carrier, PFA 1e-3: ~2x the false flags)
        excl = covered;                             % samples the baseline leaves out
        if opts.ChannelThresholds && pass >= 2
            dcv  = diff([false, covered, false]);    % previous pass's intervals
            f1 = find(dcv == 1); l1 = find(dcv == -1) - 1;
            strongCov = false(1, Nc);
            for i = 1:numel(f1)
                if max(p(f1(i):l1(i))) >= rcPeakMin * mb(blkOf(f1(i))), strongCov(f1(i):l1(i)) = true; end
            end
            use = inSup & ~strongCov;
            if nnz(use) >= 10 * minUse
                rc = channelCorr(x, use, maxLag);
                chanRho(j, :) = rc.';
                etaC = windowThresholds(rc, scales, opts.PFA);
                ownThr(j) = max(abs(etaC ./ eta - 1)) > 0.01;
                if ownThr(j), etaJ = etaC; excl = strongCov; else, etaJ = eta; end
            end
        end
        % baseline per block: median of the (unflagged) supported samples / ln 2
        for b = 1:nBlk
            idx = edges(b) + 1 : edges(b + 1);
            use = inSup(idx) & ~excl(idx);
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
            hit = (c(i0 + n) - c(i0)) > etaJ(s) * n * mb(blkOf(i0));
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
    thrUsed(j, :) = etaJ;
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
info.channelThresholds = opts.ChannelThresholds;
info.ownThreshold    = ownThr;                      % channel used its own thresholds (final pass)
info.thresholdUsed   = thrUsed;                     % x baseline, per channel and scale
info.chanRho         = chanRho;                     % r_c(0..maxLag) per channel (NaN: not measured)
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
    if opts.ChannelThresholds
        fprintf('  own thresholds (correlated continuous RFI) in %d channel(s)%s\n', nnz(ownThr), ...
            sprintf(' %d', find(ownThr)));
    end
end
end


% =====================================================================================
function [eta, lambda, nEff] = windowThresholds(rho, scales, PFA)
%WINDOWTHRESHOLDS  Thresholds eta (x baseline) per window length for complex Gaussian samples
% with correlation rho(d), d = 0.. (real: the channelizer's; complex: a channel's own): the
% window's power sum is sum_k lambda_k E_k, lambda_k the eigenvalues of the Hermitian
% toeplitz(rho(0..n-1)); eta solves P(sum > n*eta) = PFA exactly (sumExpSurvival).
nS = numel(scales);
nEff = zeros(1, nS); eta = zeros(1, nS); lambda = cell(1, nS);
for s = 1:nS
    n = scales(s);
    d = 1:min(n - 1, numel(rho) - 1);
    nEff(s) = n / (1 + 2 * sum((1 - d/n) .* abs(rho(d + 1)).'.^2));
    r = zeros(n, 1); r(1:min(n, numel(rho))) = rho(1:min(n, numel(rho)));
    lam = real(eig(toeplitz(r)));
    lambda{s} = lam(lam > 1e-9 * max(lam));         % ~0: no contribution
    % exact threshold on the sum; the Gamma value (slightly low) brackets it
    x0 = gammaincinv(PFA, nEff(s), 'upper') * n / nEff(s);
    f  = @(x) log(sumExpSurvival(lambda{s}, x)) - log(PFA);
    br = [0.5 * x0, 2 * x0];
    while sign(f(br(1))) == sign(f(br(2))) && br(2) / br(1) < 1e4   % strongly correlated
        br = [br(1) / 2, br(2) * 2];                % channels: widen (white: never needed)
    end
    eta(s) = fzero(f, br) / n;
end
end

function rc = channelCorr(x, use, maxLag)
%CHANNELCORR  Normalized complex autocorrelation r(d) = <x(m) conj(x(m+d))> / <|x|^2> of a
% channel's samples, d = 0..maxLag, over pairs with both samples in use.
rc = zeros(maxLag + 1, 1);
for d = 0:maxLag
    m = use(1:end - d) & use(1 + d:end);
    rc(d + 1) = sum(x([m, false(1, d)]) .* conj(x([false(1, d), m]))) / nnz(m);
end
rc = rc / real(rc(1));
end

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

function [p, z] = readPower(file, Nc, wantIQ)
%READPOWER  |x|^2 of a cf32 channel file, as a double row (and the complex samples if wantIQ).
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
z = [];
if wantIQ, z = complex(x(1, :), x(2, :)); end
end
