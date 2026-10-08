function [detection, info] = detectPulsar(fold, info_fold, template, opts)
%DETECTPULSAR  Neyman-Pearson detection of the pulsar in folded profiles.
%{
Is the pulsar there? For every sub-integration (and the total fold) this
tests H0 "baseline + noise only" against H1 "baseline + template-shaped
pulse + noise" at a fixed false-alarm probability P_FA (Neyman-Pearson:
the highest detection probability P_D for that P_FA). Uses only what an
observer knows: the fold, the template, the noise bandwidth and P_FA.

  [detection, info] = detectPulsar(fold, info_fold, template, 'Bnoise', Bn)

Noise normalization (H0 model):
  baseline a = mean of the profile (exact under H0; slightly high when a
  pulse is present). Per time bin var = a^2/(Bnoise*binDt) (radiometer);
  per phase bin j var v_j = that * W2_j/W_j^2, covariances with bins j+d
  from fold.weightX (all lags d), as in estimateTOA. Several channels:
  combined as in estimateTOA (combineChannels: sum of the channel profiles,
  exclusion rule), each channel with its own baseline a_c, and the channel
  variances and covariances add up.
  detection.normProfile = (p - a)./sqrt(v): the profile in units of its noise.

Known phase ('Phase', in turns, default 0 = ephemeris prediction):
  matched filter T0 = c'(p - a) / sqrt(c' Cov c), c = template shifted to
  the phase, minus its mean (the baseline drops out). Under H0 T0 ~ N(0,1),
  so the threshold is etaKnown = Phi^-1(1 - P_FA); under H1 its mean is the
  matched-filter SNR, P_D = 1 - Phi(etaKnown - SNR).

Unknown phase (acquisition): Tmax = max over all bin shifts of T(tau),
  numerator and noise both by FFT correlation. Searching all phases raises
  the false-alarm rate (look-elsewhere effect); Rice's formula for the
  maximum of a stationary Gaussian process over one turn gives
    P_FA(eta) ~ (1 - Phi(eta)) + sqrt(lambda2)/(2*pi) * exp(-eta^2/2),
    lambda2 = sum (2*pi*k)^2 |S_k|^2 / sum |S_k|^2  (template harmonics k >= 1)
  and etaUnknown solves P_FA(eta) = 'PFA'.

T0 and Tmax are calibrated under H0 (for deciding detection). For the
signal-to-noise ratio of a detected pulse use estimateTOA's snr.

Noise check: detection.noiseRatio = measured variance of the off-pulse bins
(template at phaseMax below OnPulseFrac of its peak) / radiometer variance.
About 1 for clean data; > 1 flags RFI or a wrong noise model.

Inputs:
  fold, info_fold  from foldProfile (channels combined as in estimateTOA)
  template         [NBin x 1], phase 0 at bin 1 (as for estimateTOA)

Name-value options:
  'Bnoise'       [Hz] noise-equivalent bandwidth of one channel of the fold
                 (required; scalar or one value per channel), as estimateTOA
  'PFA'          false-alarm probability per profile (default 1e-3)
  'Phase'        [turns] known (predicted) pulse phase (default 0)
  'Weighting'    'equal' (default): as above. 'optimal': the weighted matched
                 filter for several channels (as estimateTOA 'optimal'): every
                 bin of every channel weighted by 1/variance under H0 (channel
                 level a_c from its data-weighted mean), baselines removed per
                 channel, known relative gains s_c ('ChannelGain'):
                   N(tau) = sum_c s_c sum_j w_cj (p_cj - pbar_c) T(phi_j - tau)
                 normalized by its exact H0 standard deviation (all lags of
                 weightX), so T0 and T(tau) are N(0,1) under H0 and the
                 thresholds are unchanged. No channel is left out (empty bins
                 weigh nothing). normProfile = sum_c s_c w_cj (p_cj - pbar_c) /
                 sqrt(sum_c s_c^2 w_cj) (the weighted combination per bin, in
                 noise units); noiseRatio from it. ('baseline' = sum of a_c.)
  'ChannelGain'  s_c for 'optimal' (default 1 for every channel).
  'MinCoverage'  fraction of phase bins with data (default 1, as estimateTOA)
  'OnPulseFrac'  template fraction defining on-pulse for the noise check (0.01)
  'Verbose'      print a summary (default true)

Outputs:
  detection  per sub-integration (column vectors, NaN / false when not tested):
               tested, coverage, baseline, T0, Tmax, phaseMax [turns],
               detectedKnown (T0 > etaKnown), detectedUnknown (Tmax > etaUnknown),
               noiseRatio, nChanUsed; normProfile [NBin x nSub];
             etaKnown, etaUnknown; total (same fields for the whole fold).
             (Not named 'det': that would shadow MATLAB's det().)
  info       PFA, phaseKnown, Bnoise, lambda2, riceFactor (= sqrt(lambda2)/2pi), method.
%}

arguments
    fold      struct
    info_fold struct
    template  double
    opts.Bnoise      double = []
    opts.PFA         (1,1) double {mustBeInRange(opts.PFA, 0, 1, 'exclusive')} = 1e-3
    opts.Phase       (1,1) double = 0
    opts.Weighting   {mustBeTextScalar} = 'equal'
    opts.ChannelGain double = []
    opts.MinCoverage (1,1) double {mustBeInRange(opts.MinCoverage, 0, 1)} = 1
    opts.OnPulseFrac (1,1) double {mustBePositive} = 0.01
    opts.Verbose     (1,1) logical = true
end

N    = info_fold.NBin;
nSub = info_fold.nSub;
if isempty(opts.Bnoise)
    error('detectPulsar:bnoise', 'Bnoise is required (radiometer noise model).');
end
template = template(:);
if numel(template) ~= N
    error('detectPulsar:template', 'Template must have NBin = %d points.', N);
end
template = template / max(template);

% ---- Template quantities ---------------------------------------------------------------
c  = template - mean(template);                       % zero mean: baseline drops out
C  = fft(c);
K  = floor(N/2) - 1;
k  = (1:K).';
lambda2 = sum((2*pi*k).^2 .* abs(C(k + 1)).^2) / sum(abs(C(k + 1)).^2);
riceFactor = sqrt(lambda2) / (2*pi);
ksg = [0:ceil(N/2)-1, -floor(N/2):-1].';              % signed harmonic numbers
c0  = real(ifft(C .* exp(-2i*pi*ksg*opts.Phase)));    % template at the known phase

% ---- Thresholds ---------------------------------------------------------------------------
Q = @(x) 0.5 * erfc(x / sqrt(2));                     % 1 - Phi(x)
etaKnown   = sqrt(2) * erfcinv(2 * opts.PFA);
etaUnknown = fzero(@(e) Q(e) + riceFactor*exp(-e^2/2) - opts.PFA, [0 40]);

% Weights: [NBin x nSub x nW (x nLag)], nW = nChan (per channel) or 1 (shared)
nChan = size(fold.prof, 3);
if isfield(fold, 'weightX'), WXall = fold.weightX; else, WXall = zeros(N, nSub); end
nW   = size(fold.weight, 3);
nLag = size(WXall, 4);
if ~any(numel(opts.Bnoise) == [1 nChan])
    error('detectPulsar:bnoise', 'Bnoise must be a scalar or one value per channel (%d).', nChan);
end
Bc = reshape(opts.Bnoise, 1, []);                     % per channel, or scalar
tmplCC = zeros(N, nLag);                              % template products at lags 1..nLag
for lg = 1:nLag
    tmplCC(:, lg) = fft(c .* circshift(c, -lg));
end

cst = struct('N', N, 'c', c, 'C', C, 'c0', c0, 'tmplC2', fft(c.^2), ...
    'tmplCC', tmplCC, 'shapeOn', template, ...
    'onFrac', opts.OnPulseFrac, 'binDt', info_fold.binDt);
weighting = lower(char(opts.Weighting));
if ~any(strcmp(weighting, {'equal', 'optimal'}))
    error('detectPulsar:weighting', 'Weighting must be ''equal'' or ''optimal''.');
end
if strcmp(weighting, 'optimal')
    sGain = opts.ChannelGain;
    if isempty(sGain), sGain = ones(1, nChan); end
    if numel(sGain) ~= nChan || any(~isfinite(sGain)) || any(sGain < 0) || ~any(sGain > 0)
        error('detectPulsar:gain', 'ChannelGain must have one value >= 0 per channel (%d).', nChan);
    end
    sGain = reshape(sGain, 1, []);
end

% ---- Per sub-integration -------------------------------------------------------------------
nanv = nan(nSub, 1); falv = false(nSub, 1);
detection = struct('tested', falv, 'coverage', nanv, 'baseline', nanv, 'T0', nanv, ...
    'Tmax', nanv, 'phaseMax', nanv, 'detectedKnown', falv, 'detectedUnknown', falv, ...
    'noiseRatio', nanv, 'normProfile', nan(N, nSub), 'nChanUsed', zeros(nSub, 1));
for s = 1:nSub
    Ps  = reshape(fold.prof(:, s, :), N, nChan);
    Ws  = reshape(fold.weight(:, s, :), N, nW);
    W2s = reshape(fold.weight2(:, s, :), N, nW);
    WXs = reshape(WXall(:, s, :, :), N, nW, nLag);
    if strcmp(weighting, 'optimal')
        detection.coverage(s) = mean(any(Ws > 0, 2));
        if detection.coverage(s) < opts.MinCoverage || detection.coverage(s) == 0, continue; end
        r = testOptimal(Ps, Ws, W2s, WXs, Bc, sGain, cst);
        detection.nChanUsed(s) = r.nUsed;
    else
        ch = combineChannels(Ps, Ws, W2s, WXs, Bc);
        detection.coverage(s)  = ch.coverage;
        detection.nChanUsed(s) = ch.nUsed;
        if detection.coverage(s) < opts.MinCoverage || detection.coverage(s) == 0, continue; end
        r = testOne(ch, cst);
    end
    detection.tested(s)          = true;
    detection.baseline(s)        = r.a;
    detection.T0(s)              = r.T0;
    detection.Tmax(s)            = r.Tmax;
    detection.phaseMax(s)        = r.phaseMax;
    detection.detectedKnown(s)   = r.T0 > etaKnown;
    detection.detectedUnknown(s) = r.Tmax > etaUnknown;
    detection.noiseRatio(s)      = r.noiseRatio;
    detection.normProfile(:, s)  = r.z;
end
detection.etaKnown   = etaKnown;
detection.etaUnknown = etaUnknown;

% ---- Total fold -----------------------------------------------------------------------------
Wtot  = reshape(sum(fold.weight, 2), N, nW);
W2tot = reshape(sum(fold.weight2, 2), N, nW);
WXtot = reshape(sum(WXall, 2), N, nW, nLag);
if strcmp(weighting, 'optimal')
    rt = testOptimal(fold.profTotal, Wtot, W2tot, WXtot, Bc, sGain, cst);
    nUsedT = rt.nUsed;
else
    cht = combineChannels(fold.profTotal, Wtot, W2tot, WXtot, Bc);
    rt = testOne(cht, cst);
    nUsedT = cht.nUsed;
end
detection.total = struct('baseline', rt.a, 'T0', rt.T0, 'Tmax', rt.Tmax, ...
    'phaseMax', rt.phaseMax, 'detectedKnown', rt.T0 > etaKnown, ...
    'detectedUnknown', rt.Tmax > etaUnknown, 'noiseRatio', rt.noiseRatio, ...
    'normProfile', rt.z, 'nChanUsed', nUsedT);

info = struct('method', 'Neyman-Pearson matched filter (known phase) and max over phase (Rice threshold)', ...
    'PFA', opts.PFA, 'phaseKnown', opts.Phase, 'Bnoise', opts.Bnoise, ...
    'lambda2', lambda2, 'riceFactor', riceFactor, 'NBin', N, ...
    'noiseModel', 'H0 radiometer: per time bin var = baseline^2/(Bnoise*binDt)', ...
    'weighting', weighting);
if strcmp(weighting, 'optimal')
    info.method = ['weighted matched filter over channels (weights 1/var under H0), exact H0 ' ...
                   'normalization; known phase and max over phase (Rice threshold)'];
    info.channelGain = sGain;
end

if opts.Verbose
    t = detection.tested;
    fprintf(['detectPulsar: P_FA %.3g -> threshold %.2f (known phase), %.2f (unknown phase, ' ...
             'Rice factor %.2f)\n'], opts.PFA, etaKnown, etaUnknown, riceFactor);
    fprintf(['  %d sub-int(s) tested: detected %d (known phase), %d (unknown phase); ' ...
             'median T0 %.1f, median noise ratio %.3f\n'], nnz(t), nnz(detection.detectedKnown), ...
        nnz(detection.detectedUnknown), median(detection.T0(t)), median(detection.noiseRatio(t)));
    fprintf('  total fold: T0 %.1f, Tmax %.1f at phase %+.4f turns, noise ratio %.3f\n', ...
        rt.T0, rt.Tmax, rt.phaseMax, rt.noiseRatio);
end
end


% =====================================================================================
function r = testOne(ch, c)
%TESTONE  Noise normalization and detection statistics of one (channel-combined)
% profile; ch from combineChannels.
N = c.N;
p = ch.p; W = ch.W; W2 = ch.W2; WX = ch.WX; have = ch.have;
nLag = size(WX, 3);
a = mean(p(have));

% H0 noise model: variance per bin and covariances with bins j+d, summed over the
% (independent) channels; per time bin var = level^2/(Bnoise*binDt), level = the
% channel's baseline
if size(ch.Pc, 2) == 1                                % one channel (full band)
    s2 = 1 / (ch.B * c.binDt) * a^2;
else
    s2 = 1 ./ (ch.B * c.binDt) .* mean(ch.Pc(have, :), 1).^2;
end
v  = sum(s2 .* W2 ./ W.^2, 2);
cv = zeros(N, nLag);
for lg = 1:nLag
    cv(:, lg) = sum(s2 .* WX(:, :, lg) ./ (W .* circshift(W, -lg)), 2);
end
v(~have | ~isfinite(v)) = 0;
for lg = 1:nLag
    cl = cv(:, lg);
    cl(~have | ~circshift(have, -lg) | ~isfinite(cl)) = 0;
    cv(:, lg) = cl;
end
d = p - a; d(~have) = 0;

% Known phase: matched filter with the template at the predicted phase
den0 = sum(c.c0.^2 .* v);
for lg = 1:nLag
    den0 = den0 + 2*sum(c.c0 .* circshift(c.c0, -lg) .* cv(:, lg));
end
T0 = (c.c0.' * d) / sqrt(den0);

% Unknown phase: all bin shifts at once; element m+1 = template shifted by m bins
num = real(ifft(fft(d) .* conj(c.C)));
den = real(ifft(fft(v) .* conj(c.tmplC2)));
for lg = 1:nLag
    den = den + 2*real(ifft(fft(cv(:, lg)) .* conj(c.tmplCC(:, lg))));
end
T = num ./ sqrt(max(den, realmin));
[Tmax, im] = max(T);
phaseMax = mod((im - 1)/N + 0.5, 1) - 0.5;

% Noise check on the off-pulse bins (template at phaseMax)
shape = circshift(c.shapeOn, im - 1);
off = have & shape < c.onFrac;
po = p(off);                                          % expected variance with the off-pulse
noiseRatio = mean((po - mean(po)).^2) / (mean(v(off)) * (mean(po)/a)^2);   % baseline

z = nan(N, 1);
z(have) = d(have) ./ sqrt(v(have));
r = struct('a', a, 'T0', T0, 'Tmax', Tmax, 'phaseMax', phaseMax, ...
    'noiseRatio', noiseRatio, 'z', z);
end


% =====================================================================================
function r = testOptimal(Pc, W, W2, WX, B, sGain, c)
%TESTOPTIMAL  Weighted matched filter over channels under H0 (see 'Weighting').
% Pc [N x nChan] channel profiles; W, W2 [N x nW], WX [N x nW x nLag] fold weights
% (nW = 1: shared); B = Bnoise (scalar or per channel); sGain [1 x nChan].
N = c.N; nChan = size(Pc, 2); nLag = size(WX, 3);
r = struct('a', NaN, 'T0', NaN, 'Tmax', NaN, 'phaseMax', NaN, 'noiseRatio', NaN, ...
           'z', nan(N, 1), 'nUsed', 0);
if size(W, 2) == 1
    W = repmat(W, 1, nChan); W2 = repmat(W2, 1, nChan); WX = repmat(WX, 1, nChan, 1);
end
hasC = W > 0;
use  = any(hasC, 1) & sGain > 0;
if ~any(use), return; end
Pc = Pc(:, use); W = W(:, use); W2 = W2(:, use); WX = WX(:, use, :); hasC = hasC(:, use);
s = sGain(use);
if numel(B) > 1, B = B(use); end
Pc(~hasC) = 0;
haveAny = any(hasC, 2);

% H0 noise of every channel: level a_c (data-weighted mean), radiometer, all lags
aC  = sum(W .* Pc, 1) ./ sum(W, 1);
s2  = aC.^2 ./ (B * c.binDt);                       % per time bin variance, per channel
v   = s2 .* W2 ./ W.^2;
v(~hasC) = 0;
cv  = zeros(N, size(Pc, 2), nLag);
for d = 1:nLag
    cl = s2 .* WX(:, :, d) ./ (W .* circshift(W, -d, 1));
    cl(~hasC | ~circshift(hasC, -d, 1) | ~isfinite(cl)) = 0;
    cv(:, :, d) = cl;
end
w = zeros(size(v)); okv = hasC & v > 0; w(okv) = 1 ./ v(okv);
Wc   = sum(w, 1);
pbar = sum(w .* Pc, 1) ./ Wc;
u    = sum(s .* w .* (Pc - pbar), 2);               % weighted, baseline-free combination
Om   = sum(s.^2 .* w, 2);

% all bin shifts m at once (element m+1 = template shifted by m bins)
Num = real(ifft(fft(u) .* conj(c.C)));
WT  = real(ifft(fft(w) .* conj(c.C)));              % sum_j w_cj c(j-m), [N x nU]
Tb  = WT ./ Wc;                                     % weighted template mean per channel
Var = real(ifft(fft(Om) .* conj(c.tmplC2))) - sum(s.^2 .* Wc .* Tb.^2, 2);
for d = 1:nLag
    h  = w .* circshift(w, -d, 1) .* cv(:, :, d);   % weights x lag-d covariance
    Hd = sum(s.^2 .* h, 2);
    E  = real(ifft(fft(h + circshift(h, d, 1)) .* conj(c.C)));
    Ld = real(ifft(fft(Hd) .* conj(c.tmplCC(:, d)))) - sum(s.^2 .* Tb .* E, 2) ...
         + sum(s.^2 .* Tb.^2 .* sum(h, 1), 2);
    Var = Var + 2 * Ld;
end
T = Num ./ sqrt(max(Var, realmin));
[Tmax, im] = max(T);
phaseMax = mod((im - 1)/N + 0.5, 1) - 0.5;

% known phase: the template at the predicted phase, direct sums
g  = s .* w .* (c.c0 - sum(w .* c.c0, 1) ./ Wc);    % N0 = sum_c g_c' p_c
V0 = sum(g.^2 .* v, 'all');
for d = 1:nLag
    V0 = V0 + 2 * sum(g .* circshift(g, -d, 1) .* cv(:, :, d), 'all');
end
T0 = (c.c0.' * u) / sqrt(V0);

% weighted combination per bin in noise units; noise check off-pulse (template at phaseMax)
z = nan(N, 1);
z(haveAny) = u(haveAny) ./ sqrt(Om(haveAny));
shape = circshift(c.shapeOn, im - 1);
off = haveAny & shape < c.onFrac;
aOff = sum(W .* Pc .* off, 1) ./ sum(W .* off, 1);  % off-pulse level per channel
rho2 = (aOff ./ aC).^2;  rho2(~isfinite(rho2)) = 1;
eVar = sum(s.^2 .* w .* rho2, 2) ./ Om;              % expected var of z (level correction)
zo = z(off);
noiseRatio = mean((zo - mean(zo)).^2) / mean(eVar(off));

r.a = sum(aC); r.T0 = T0; r.Tmax = Tmax; r.phaseMax = phaseMax;
r.noiseRatio = noiseRatio; r.z = z; r.nUsed = nnz(use);
end
