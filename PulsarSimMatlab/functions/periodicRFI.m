function [rows, info] = periodicRFI(info_chan, info_rfi, opts)
%PERIODICRFI  Find periodic RFI (radars) in the detections; blank all its predicted pulses (B6).
%{
Why: RFI whose pulses repeat at a fixed period can pile up at fixed pulse
phases in the fold (a radar locked to the pulsar: B5b), and its weak pulses
(antenna sidelobes, B5c: at -35 dB only ~40 % are caught per pulse) leave
far more than the -54 dB pulsar. What makes it dangerous - regular timing -
also makes it predictable: once the period and the time of one pulse are
known, every pulse is known, also those too weak to be seen. RFI with
irregular timing (Poisson impulses, PRF jitter) cannot be predicted, but it
cannot fold coherently either (only extra noise).

  rows = periodicRFI(info_chan, info_rfi)
  mask = [mask; rows]                          % then blankChannels / blankingWeights

Inputs:
  info_chan  info of channelizeIQ (fs, N, nChan, t0).
  info_rfi   info of detectRFI (events, guard, freqGuard): uses only the
             detections, no ground truth.

Method:
  1. Pulses: detected intervals at the same time (centres within
     'GroupTime' samples) in neighbouring channels (gaps of <= 2 channels)
     are one pulse. Its time = the centre of the interval in its strongest
     channel (the time guard is symmetric); its channels = the range.
  2. Emitters: narrowband pulses (spanning < 'MaxSpan' channels) whose
     strongest channels lie within 'ChanGroup' channels of each other belong
     to one candidate emitter; broadband pulses (impulses) form one group of
     their own (their strongest channel is random).
  3. Period: candidates = the gaps between a pulse and its next 3 pulses,
     divided by 1..'MaxMissing' (missed pulses leave gaps of 2P, 3P, ...),
     within 'PeriodRange'. For each candidate the best phase and the m
     pulses (of n) within +-'TimeTol' samples of t0 + k*P (inliers), and
     how likely that is for random pulse times: one pulse sets the phase,
     each other one fits by chance with q = (2 TimeTol + 1)/P, so
     p = n x P(Binomial(n-1, q) >= m-1). The smallest p wins (so P, not P/2:
     same pulses, twice the chance); accepted if p x (number of candidates)
     < 'FalseAlarm' and m >= 'MinPulses'. (Counting inliers alone is not
     enough: for short trial periods random impulses fit 6 pulses somewhere.)
     t0 and P refined by least squares on the inliers (twice); the inliers
     are removed and the group searched again (several emitters in one
     channel range).
  4. Mask: every predicted pulse in the file is blanked unconditionally -
     the mask does not depend on the noise values, so blankingWeights stays
     exact and nothing is truncated - over pulse width + 2 x ('Guard' +
     'TimeTol' + 4 x the prediction uncertainty at that pulse), in the
     emitter's channel range +- 'FreqGuard'.
  Fixed PRF only (version 1): a staggered PRF (a repeating pattern of
  intervals) fits no single period; such a group is reported in
  info.unfitted (its detected pulses are still blanked by detectRFI).

Name-value options:
  'MinPulses'   inliers needed to accept an emitter (default 6).
  'TimeTol'     [samples] fit tolerance of a pulse time (default 4 = 1 us).
  'GroupTime'   [samples] same pulse in several channels (default 8).
  'ChanGroup'   [channels] strongest channels this close = same emitter
                (default 2).
  'MaxSpan'     [channels] pulses spanning this many or more are broadband
                (default 32).
  'PeriodRange' [s] [min max] (default [50e-6 20e-3]: PRF 50 Hz - 20 kHz).
  'MaxMissing'  gaps of up to this many periods tried (default 8).
  'FalseAlarm'  chance of accepting random pulse times as an emitter,
                per group search (default 1e-6).
  'Guard'       [samples] (default info_rfi.guard).
  'FreqGuard'   [channels] (default info_rfi.freqGuard).
  'MaxEmitters' (default 5).
  'Verbose'     print a summary (default true).

Outputs:
  rows  [M x 3] mask rows [channel, first, last] of the predicted pulses
        (per channel sorted and merged; may overlap detectRFI's mask).
  info  emitters (struct array: prf, period, periodErr [s], logChance (log
        of the chance probability incl. the number of candidates), t0 [s] (pulse
        k = 0), nGroup (pulses in the group), nInliers, rmsResid [samples],
        channels [c1 c2] (detected), blankChannels [lo hi], pulseWidth [s],
        nPredicted, maxSigmaPred [samples]), nPulses, unfitted (groups with
        >= MinPulses pulses but no accepted period: [c1 c2 nPulses]),
        blankedSamples (per channel, these rows only), blankedFraction, elapsed.
%}

arguments
    info_chan          struct
    info_rfi           struct
    opts.MinPulses     (1,1) double {mustBeInteger, mustBePositive} = 6
    opts.TimeTol       (1,1) double {mustBePositive} = 4
    opts.GroupTime     (1,1) double {mustBePositive} = 8
    opts.ChanGroup     (1,1) double {mustBeInteger, mustBeNonnegative} = 2
    opts.MaxSpan       (1,1) double {mustBeInteger, mustBePositive} = 32
    opts.PeriodRange   (1,2) double {mustBePositive} = [50e-6 20e-3]
    opts.MaxMissing    (1,1) double {mustBeInteger, mustBePositive} = 8
    opts.FalseAlarm    (1,1) double {mustBePositive} = 1e-6
    opts.Guard               double = []
    opts.FreqGuard           double = []
    opts.MaxEmitters   (1,1) double {mustBeInteger, mustBeNonnegative} = 5
    opts.Verbose       (1,1) logical = true
end

tStart = tic;
fs = info_chan.fs; Nc = info_chan.N; nChan = info_chan.nChan;
g  = opts.Guard;     if isempty(g),  g  = info_rfi.guard;     end
fg = opts.FreqGuard; if isempty(fg), fg = info_rfi.freqGuard; end
tol = opts.TimeTol;
P0 = opts.PeriodRange(1) * fs; P1 = opts.PeriodRange(2) * fs;     % [samples]

% ---- 1. Pulses: events at the same time in neighbouring channels -------------------------
ev = info_rfi.events;                                % [channel, first, last, peakSample, peak]
pul = zeros(0, 6);                                   % [t, strongest ch, c1, c2, width, peak]
if ~isempty(ev)
    tc = (ev(:, 2) + ev(:, 3)) / 2;
    [tc, o] = sort(tc); ev = ev(o, :);
    cid = cumsum([true; diff(tc) > opts.GroupTime]);
    for k = 1:cid(end)
        E = ev(cid == k, :);
        E = sortrows(E, 1);
        sid = cumsum([true; diff(E(:, 1)) > 3]);     % gaps of <= 2 channels: same pulse
        for q = 1:sid(end)
            S = E(sid == q, :);
            [~, m] = max(S(:, 5));
            pul(end+1, :) = [(S(m, 2) + S(m, 3)) / 2, S(m, 1), min(S(:, 1)), max(S(:, 1)), ...
                max(1, S(m, 3) - S(m, 2) + 1 - 2*g), S(m, 5)]; %#ok<AGROW>
        end
    end
end

% ---- 2.-3. Emitters: group by strongest channel, then fit a period ------------------------
em = struct('prf', {}, 'period', {}, 'periodErr', {}, 'logChance', {}, 't0', {}, 'nGroup', {}, 'nInliers', {}, ...
    'rmsResid', {}, 'channels', {}, 'blankChannels', {}, 'pulseWidth', {}, 'nPredicted', {}, ...
    'maxSigmaPred', {});
unfitted = zeros(0, 3);
fits = {};                                           % per emitter: a, b, cov, pulse width, channels
if ~isempty(pul)
    [~, o] = sort(pul(:, 2)); pul = pul(o, :);
    broad = pul(:, 4) - pul(:, 3) + 1 >= opts.MaxSpan;
    gid = zeros(size(pul, 1), 1);
    nb = find(~broad);
    gid(nb) = cumsum([true; diff(pul(nb, 2)) > opts.ChanGroup]);   % narrowband: by channel
    if any(broad), gid(broad) = max([0; gid]) + 1; end              % broadband: one group
    for k = 1:max(gid)
        G = pul(gid == k, :);
        nG = size(G, 1);
        nEm0 = numel(em);
        while size(G, 1) >= opts.MinPulses && numel(em) < opts.MaxEmitters
            [P, inl, lpBest] = bestPeriod(sort(G(:, 1)), P0, P1, tol, opts.MaxMissing);
            if nnz(inl) < opts.MinPulses || lpBest >= log(opts.FalseAlarm), break; end
            [~, o] = sort(G(:, 1)); G = G(o, :);
            [a, b, C, inl, rmsR] = refine(G(:, 1), P, inl, tol);
            if nnz(inl) < opts.MinPulses, break; end
            I = G(inl, :);
            c12 = [min(I(:, 3)), max(I(:, 4))];
            fits{end+1} = struct('a', a, 'b', b, 'C', C, 'w', median(I(:, 5)), 'c12', c12); %#ok<AGROW>
            em(end+1) = struct('prf', fs / b, 'period', b / fs, 'periodErr', sqrt(C(2, 2)) / fs, ...
                'logChance', lpBest, ...
                't0', info_chan.t0 + (a - 1) / fs, 'nGroup', nG, 'nInliers', nnz(inl), ...
                'rmsResid', rmsR, 'channels', c12, ...
                'blankChannels', [max(1, c12(1) - fg), min(nChan, c12(2) + fg)], ...
                'pulseWidth', median(I(:, 5)) / fs, 'nPredicted', 0, 'maxSigmaPred', 0); %#ok<AGROW>
            G = G(~inl, :);
        end
        if numel(em) == nEm0 && nG >= opts.MinPulses       % enough pulses, no period
            unfitted(end+1, :) = [min(pul(gid == k, 3)), max(pul(gid == k, 4)), nG]; %#ok<AGROW>
        end
    end
end

% ---- 4. Mask: every predicted pulse, unconditionally ----------------------------------------
R = cell(numel(em), 1);
for e = 1:numel(em)
    f = fits{e};
    hw0 = f.w / 2 + g + tol;
    kmin = ceil((1 - hw0 - f.a) / f.b) - 1; kmax = floor((Nc + hw0 - f.a) / f.b) + 1;
    kk = (kmin:kmax).';
    sig = sqrt(max(0, f.C(1, 1) + kk.^2 * f.C(2, 2) + 2 * kk * f.C(1, 2)));   % prediction error
    c  = f.a + kk * f.b;
    hw = hw0 + 4 * sig;
    iv = [floor(c - hw), ceil(c + hw)];
    keep = iv(:, 2) >= 1 & iv(:, 1) <= Nc;
    iv = iv(keep, :); sig = sig(keep);
    ch = (em(e).blankChannels(1):em(e).blankChannels(2)).';
    R{e} = [kron(ch, ones(size(iv, 1), 1)), repmat(iv, numel(ch), 1)];
    em(e).nPredicted = size(iv, 1);
    em(e).maxSigmaPred = max([0; sig]);
end
rows = normalizeMask(vertcat(R{:}), nChan, Nc);
blankedSamples = accumarray(rows(:, 1), rows(:, 3) - rows(:, 2) + 1, [nChan, 1]);

info = struct();
info.emitters        = em;
info.nPulses         = size(pul, 1);
info.unfitted        = unfitted;
info.nRows           = size(rows, 1);
info.blankedSamples  = blankedSamples;
info.blankedFraction = blankedSamples / Nc;
info.guard           = g;
info.freqGuard       = fg;
info.timeTol         = tol;
info.minPulses       = opts.MinPulses;
info.elapsed         = toc(tStart);

if opts.Verbose
    fprintf('periodicRFI: %d detected pulses, %d periodic emitter(s)\n', info.nPulses, numel(em));
    for e = 1:numel(em)
        E = em(e);
        fprintf(['  PRF %.4f Hz (period %.4f us +- %.1g ns), %d of %d pulses fit (rms %.1f samples), ' ...
                 'channels %d-%d -> blanked %d-%d, %d predicted pulses (width %.2f us)\n'], E.prf, ...
            E.period * 1e6, E.periodErr * 1e9, E.nInliers, E.nGroup, E.rmsResid, E.channels, ...
            E.blankChannels, E.nPredicted, E.pulseWidth * 1e6);
    end
    for u = 1:size(unfitted, 1)
        fprintf(['  no single period for %d pulses in channels %d-%d (random, e.g. impulses, ' ...
                 'or a staggered PRF)\n'], unfitted(u, 3), unfitted(u, 1:2));
    end
    fprintf('periodicRFI: %.4f %% of samples blanked by prediction; %.2f s\n', ...
        100 * mean(info.blankedFraction), info.elapsed);
end
end


% =====================================================================================
function [P, inl, lpBest] = bestPeriod(t, P0, P1, tol, maxMissing)
%BESTPERIOD  The most significant period for the pulse times t (sorted): the least likely
% to fit this well by chance. lpBest = log(chance probability x number of candidates).
n = numel(t);
cand = [];
for i = 1:n - 1
    for j = i + 1:min(n, i + 3)
        c = (t(j) - t(i)) ./ (1:maxMissing);
        cand = [cand, c(c >= P0 & c <= P1)]; %#ok<AGROW>
    end
end
P = NaN; inl = false(n, 1); lpBest = Inf;
if isempty(cand), return; end
cand = unique(round(cand * 64) / 64);                % 1/64 sample: enough, fewer duplicates
best = Inf;
for k = 1:numel(cand)
    in = inliers(t, cand(k), tol);
    m  = nnz(in);
    if m < 2, continue; end
    q  = min(1, (2*tol + 1) / cand(k));
    lp = log(n) + log(betainc(q, m - 1, n - m + 1));    % n x P(Bin(n-1, q) >= m-1)
    if lp < best || (lp == best && cand(k) > P)
        best = lp; P = cand(k); inl = in;
    end
end
lpBest = best + log(numel(cand));                    % for all candidates tried
end

function in = inliers(t, P, tol)
%INLIERS  Pulses within +-tol of the best phase of period P (each pulse tried as the phase).
r = mod(t, P);
best = 0; in = false(size(t));
for i = 1:numel(t)
    d = abs(mod(r - r(i) + P/2, P) - P/2);
    cur = d <= tol;
    if nnz(cur) > best, best = nnz(cur); in = cur; end
end
end

function [a, b, C, inl, rmsR] = refine(t, P, inl, tol)
%REFINE  Least squares t_i = a + b*k_i on the inliers (k counted from the first inlier,
% then centred), twice; C = covariance of [a b] with the residual scatter (>= 0.5 sample).
a = NaN; b = P; C = zeros(2); rmsR = NaN;
for it = 1:2
    if nnz(inl) < 2, return; end
    ti = t(inl); tRef = ti(1);
    k = round((ti - tRef) / b);
    kc = mean(k);
    A = [ones(size(k)), k - kc];
    x = A \ (ti - tRef);
    b = x(2);
    a = tRef + x(1) - b * kc;                        % time of pulse k = 0 (k from the first inlier)
    res = ti - (a + b * k);
    s2 = max(mean(res.^2), 0.25);
    Cc = s2 * ((A' * A) \ eye(2));                  % covariance of [offset at kc, b]
    % covariance of [a, b] with a = offset - b*kc
    J = [1, -kc; 0, 1];
    C = J * Cc * J';
    rmsR = sqrt(mean(res.^2));
    kall = round((t - a) / b);
    inl = abs(t - (a + b * kall)) <= tol;
end
end

function mask = normalizeMask(rows, nChan, Nc)
%NORMALIZEMASK  Rows [channel, first, last] -> per channel clipped, sorted, merged.
out = cell(nChan, 1);
if isempty(rows), mask = zeros(0, 3); return; end
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
