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

% ---- Per sub-integration -------------------------------------------------------------------
nanv = nan(nSub, 1); falv = false(nSub, 1);
detection = struct('tested', falv, 'coverage', nanv, 'baseline', nanv, 'T0', nanv, ...
    'Tmax', nanv, 'phaseMax', nanv, 'detectedKnown', falv, 'detectedUnknown', falv, ...
    'noiseRatio', nanv, 'normProfile', nan(N, nSub), 'nChanUsed', zeros(nSub, 1));
for s = 1:nSub
    ch = combineChannels(reshape(fold.prof(:, s, :), N, nChan), ...
        reshape(fold.weight(:, s, :), N, nW), reshape(fold.weight2(:, s, :), N, nW), ...
        reshape(WXall(:, s, :, :), N, nW, nLag), Bc);
    detection.coverage(s)  = ch.coverage;
    detection.nChanUsed(s) = ch.nUsed;
    if detection.coverage(s) < opts.MinCoverage || detection.coverage(s) == 0, continue; end
    r = testOne(ch, cst);
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
cht = combineChannels(fold.profTotal, reshape(sum(fold.weight, 2), N, nW), ...
    reshape(sum(fold.weight2, 2), N, nW), reshape(sum(WXall, 2), N, nW, nLag), Bc);
rt = testOne(cht, cst);
detection.total = struct('baseline', rt.a, 'T0', rt.T0, 'Tmax', rt.Tmax, ...
    'phaseMax', rt.phaseMax, 'detectedKnown', rt.T0 > etaKnown, ...
    'detectedUnknown', rt.Tmax > etaUnknown, 'noiseRatio', rt.noiseRatio, ...
    'normProfile', rt.z, 'nChanUsed', cht.nUsed);

info = struct('method', 'Neyman-Pearson matched filter (known phase) and max over phase (Rice threshold)', ...
    'PFA', opts.PFA, 'phaseKnown', opts.Phase, 'Bnoise', opts.Bnoise, ...
    'lambda2', lambda2, 'riceFactor', riceFactor, 'NBin', N, ...
    'noiseModel', 'H0 radiometer: per time bin var = baseline^2/(Bnoise*binDt)');

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
