function [toa, info] = estimateTOA(fold, info_fold, template, opts)
%ESTIMATETOA  Times of arrival from folded profiles by template matching.
%{
Fourier-domain template matching (Taylor 1992, "FFTFIT") of every
sub-integration profile, and of the total fold, against a known template:

  profile(phi) = a + b * template(phi - tau) + noise

giving the phase offset tau (turns, sub-bin precision), amplitude b and
baseline a. The TOA of sub-integration s is

  TOA_s = tRef_s + tau_s / fRef_s          (fold.subint.tRef, .fRef)

i.e. the arrival time of the template's phase 0 (pulse peak for a
template peaking at bin 1), at the dedispersion reference frequency.

  [toa, info] = estimateTOA(fold, info_fold, template, 'Bnoise', Bn)

Inputs:
  fold, info_fold  from foldProfile. Several channels: see 'Weighting'.
                   'equal' combines them into one profile p = sum of the
                   channel profiles; a channel with no data in some phase
                   bins of a profile, where other channels have data, is
                   left out of that profile (its absence would leave a
                   dip). Channels with partial data (0 < weight) stay in.
  template         [NBin x 1] profile shape, phase 0 at bin 1 (the value at
                   bin j is the template at phase (j-1)/NBin). Normalized
                   here to peak 1, so b is the peak height above baseline.

Name-value options:
  'NoiseModel'  'radiometer' (default): per time bin, var = m^2/(Bnoise*binDt)
                with m = a + b*template the fitted power model, per channel
                (a_c, b_c: the same Fourier projection as a, b, with tau
                fixed); the fold weights carry the rest (time-bin
                correlations, blanking). Exact for Gaussian signal +
                Gaussian receiver noise (self-noise and system noise alike);
                needs 'Bnoise'.
                'offpulse': variance and lag covariances estimated from the
                residuals in off-pulse bins (classic approach; needs
                off-pulse noise, i.e. receiver noise in the data).
  'Bnoise'      [Hz] noise-equivalent bandwidth of ONE channel of the fold,
                (int W^2)^2 / int W^4 for its bandpass W: scalar (same for
                all channels) or one value per channel. Single-channel
                (full-band) fold: the band's noise bandwidth.
  'Weighting'   'equal' (default): FFTFIT of the equal-weight channel sum
                (above). 'optimal': weighted fit of all channels at once,
                  prof_cj = a_c + b * s_c * T(phi_j - tau) + noise,
                with a baseline a_c per channel, one amplitude b times the
                known relative channel gain s_c ('ChannelGain'), one tau;
                every bin of every channel weighted by 1/variance (noise
                model as below, iterated with the fitted model). Empty bins
                weigh nothing, so no channel is left out and there are no
                dips; nearly empty bins count for almost nothing. For fixed
                tau, a_c and b follow by weighted least squares; tau
                maximizes N(tau)/sqrt(Dn(tau)) (the weighted matched filter;
                FFT correlations on the Upsample grid, then a bounded 1-D
                search). Error bars: sandwich Cov = A^-1 B A^-1 of the
                linearized fit (A = J'WJ, B = J'W Sigma W J) with the full
                covariance Sigma of every channel (all lags of weightX).
                'amp' and 'ampErr' are b * sum(s_c) (the summed pulse, as
                for 'equal'); 'baseline' is sum(a_c).
  'ChannelGain' s_c for 'optimal': relative pulse amplitude per channel
                (pulsar spectrum x bandpass, observer knowledge); default
                1 for every channel.
  'MaxIterations' 'optimal': fit -> weights -> fit until tau and b change
                by < 1e-12 (default 10 rounds at most).
  'MinCoverage' fraction of phase bins that must have data (default 1:
                only complete turns). Missing bins are filled by circular
                linear interpolation and get zero weight in the noise.
  'MaxHarmonic' highest harmonic used (default floor(NBin/2) - 1).
  'Upsample'    oversampling of the coarse cross-correlation (default 8).
  'OnPulseFrac' bins with model above this fraction of the peak (above
                baseline) define on-pulse (chi^2 and off-pulse mask;
                default 0.01).
  'MaxRedChi2'  quality flag: toa.flagChi2 = redChi2 > MaxRedChi2
                (default 2). The reduced chi^2 compares the fit residuals
                with the noise model; ~1 when the data are template + noise,
                >> 1 when something else is in the profile (e.g. a radar
                pulse the fit locked onto, which gives a wrong TOA with a
                small error bar).
  'Verbose'     print a summary (default true).

Outputs:
  toa   struct, per sub-integration (column vectors, NaN when invalid):
          valid, coverage, phase [turns], phaseErr, toa [s], toaErr [s],
          amp (b), ampErr, baseline (a), snr (= b/ampErr), redChi2,
          flagChi2 (redChi2 > MaxRedChi2), tRef, fRef, turnRef,
          nChanUsed (channels combined), chanBaseline [nSub x nChan]
          ('optimal': a_c per channel, NaN where unused)
        plus toa.total (the whole fold): phase, phaseErr, timeOffset
        [s] (= phase/f0), timeOffsetErr, amp, ampErr, snr, redChi2, flagChi2,
        nChanUsed.
  info  method, noise model, Bnoise, harmonics, template.

Uncertainties:
  tau solves C'(tau) = 0 for the cross-correlation C(tau). C'(tau) is a
  linear combination sum_j d_j * p_j of the profile bins, so
    var(tau) = d' * Cov(p) * d / C''(tau)^2,
  with Cov(p) from the noise model, INCLUDING the covariance of
  neighbouring phase bins created by linear assignment in the fold
  (fold.weightX). Ignoring it underestimates var by ~1.5x for 'linear'.
  Channels are independent, so Cov(p) is the sum of the channel
  covariances: var(p_j) = sum_c s_cj^2 W2_cj / W_cj^2 and, for every lag
  d = 1..D of weightX, cov(p_j, p_j+d) = sum_c s_cj s_c,j+d WX_cjd /
  (W_cj W_c,j+d), with s_cj = m_c(phi_j) / sqrt(Bnoise_c * binDt).
  Verified by Monte Carlo on the self-noise case (predicted 0.469 us vs
  empirical 0.46-0.47 us).

Note on optimality: FFTFIT weights bins by the template derivative, which
is optimal for white (receiver-dominated) noise. With pure self-noise
(variance proportional to the profile squared, as in the current noise-free
synthetic data) a flatter weighting such as the power centroid is more
precise; once receiver noise dominates (the realistic case) FFTFIT is the
better estimator.
%}

arguments
    fold      struct
    info_fold struct
    template  double
    opts.NoiseModel  {mustBeTextScalar} = 'radiometer'
    opts.Bnoise      double = []
    opts.Weighting   {mustBeTextScalar} = 'equal'
    opts.ChannelGain double = []
    opts.MaxIterations (1,1) double {mustBeInteger, mustBePositive} = 10
    opts.MinCoverage (1,1) double {mustBeInRange(opts.MinCoverage, 0, 1)} = 1
    opts.MaxHarmonic double = []
    opts.Upsample    (1,1) double {mustBeInteger, mustBePositive} = 8
    opts.OnPulseFrac (1,1) double {mustBePositive} = 0.01
    opts.MaxRedChi2  (1,1) double {mustBePositive} = 2
    opts.Verbose     (1,1) logical = true
end

N     = info_fold.NBin;
nSub  = info_fold.nSub;
binDt = info_fold.binDt;
model = lower(char(opts.NoiseModel));
if ~any(strcmp(model, {'radiometer', 'offpulse'}))
    error('estimateTOA:noise', 'NoiseModel must be ''radiometer'' or ''offpulse''.');
end
if strcmp(model, 'radiometer') && isempty(opts.Bnoise)
    error('estimateTOA:bnoise', 'Bnoise is required for the radiometer noise model.');
end
weighting = lower(char(opts.Weighting));
if ~any(strcmp(weighting, {'equal', 'optimal'}))
    error('estimateTOA:weighting', 'Weighting must be ''equal'' or ''optimal''.');
end
template = template(:);
if numel(template) ~= N
    error('estimateTOA:template', 'Template must have NBin = %d points.', N);
end
template = template / max(template);

K = floor(N/2) - 1;
if ~isempty(opts.MaxHarmonic), K = min(K, opts.MaxHarmonic); end
k  = (1:K).';
S  = fft(template);
Sk = S(k + 1);
sumS2 = sum(abs(Sk).^2);
ksg = [0:ceil(N/2)-1, -floor(N/2):-1].';               % signed harmonic numbers

% Weights: [NBin x nSub x nW (x nLag)], nW = nChan (per channel) or 1 (shared)
nChan = size(fold.prof, 3);
if isfield(fold, 'weightX'), WXall = fold.weightX; else, WXall = zeros(N, nSub); end
nW   = size(fold.weight, 3);
nLag = size(WXall, 4);
if ~isempty(opts.Bnoise) && ~any(numel(opts.Bnoise) == [1 nChan])
    error('estimateTOA:bnoise', 'Bnoise must be a scalar or one value per channel (%d).', nChan);
end
Bc = reshape(opts.Bnoise, 1, []);                       % per channel, or scalar
cst = struct('N', N, 'k', k, 'S', S, 'Sk', Sk, 'sumS2', sumS2, 'ksg', ksg, ...
    'model', model, 'binDt', binDt, 'M', opts.Upsample * N, 'onFrac', opts.OnPulseFrac);
if strcmp(weighting, 'optimal')
    sGain = opts.ChannelGain;
    if isempty(sGain), sGain = ones(1, nChan); end
    if numel(sGain) ~= nChan || any(~isfinite(sGain)) || any(sGain < 0) || ~any(sGain > 0)
        error('estimateTOA:gain', 'ChannelGain must have one value >= 0 per channel (%d).', nChan);
    end
    sGain = reshape(sGain, 1, []);
    % band-limited template (harmonics |k| <= K and DC) on the bin grid and on the fine grid
    Sbl = zeros(N, 1); Sbl([1; k + 1; N - k + 1]) = S([1; k + 1; N - k + 1]);
    M = cst.M;
    Sf = zeros(M, 1); Sf([1; k + 1; M - k + 1]) = Sbl([1; k + 1; N - k + 1]) * (M / N);
    Tf = real(ifft(Sf));
    cst.Sbl = Sbl; cst.TfF = fft(Tf); cst.T2fF = fft(Tf.^2);
    cst.maxIt = opts.MaxIterations;
end

% ---- Per sub-integration -------------------------------------------------------------
nanv = nan(nSub, 1);
toa = struct('valid', false(nSub, 1), 'coverage', nanv, 'phase', nanv, ...
    'phaseErr', nanv, 'toa', nanv, 'toaErr', nanv, 'amp', nanv, 'ampErr', nanv, ...
    'baseline', nanv, 'snr', nanv, 'redChi2', nanv, 'flagChi2', false(nSub, 1), ...
    'tRef', fold.subint.tRef, 'fRef', fold.subint.fRef, 'turnRef', fold.subint.turnRef, ...
    'nChanUsed', zeros(nSub, 1), 'chanBaseline', nan(nSub, nChan));

for s = 1:nSub
    Ps  = reshape(fold.prof(:, s, :), N, nChan);
    Ws  = reshape(fold.weight(:, s, :), N, nW);
    W2s = reshape(fold.weight2(:, s, :), N, nW);
    WXs = reshape(WXall(:, s, :, :), N, nW, nLag);
    if strcmp(weighting, 'optimal')
        toa.coverage(s) = mean(any(Ws > 0, 2));
        if toa.coverage(s) < opts.MinCoverage || toa.coverage(s) == 0
            continue
        end
        r = fitOptimal(Ps, Ws, W2s, WXs, Bc, sGain, cst);
        toa.nChanUsed(s) = r.nUsed;
        toa.chanBaseline(s, :) = r.aC;
    else
        ch = combineChannels(Ps, Ws, W2s, WXs, Bc);
        toa.coverage(s)  = ch.coverage;
        toa.nChanUsed(s) = ch.nUsed;
        if toa.coverage(s) < opts.MinCoverage || toa.coverage(s) == 0
            continue
        end
        r = fitOne(ch, cst);
    end
    if ~r.ok, continue; end
    toa.valid(s)    = true;
    toa.phase(s)    = r.tau;
    toa.phaseErr(s) = r.tauErr;
    toa.toa(s)      = toa.tRef(s) + r.tau / toa.fRef(s);
    toa.toaErr(s)   = r.tauErr / toa.fRef(s);
    toa.amp(s)      = r.b;
    toa.ampErr(s)   = r.bErr;
    toa.baseline(s) = r.a;
    toa.snr(s)      = r.b / r.bErr;
    toa.redChi2(s)  = r.redChi2;
    toa.flagChi2(s) = r.redChi2 > opts.MaxRedChi2;
end

% ---- Total fold -------------------------------------------------------------------------
Wtot  = reshape(sum(fold.weight, 2), N, nW);
W2tot = reshape(sum(fold.weight2, 2), N, nW);
WXtot = reshape(sum(WXall, 2), N, nW, nLag);
if strcmp(weighting, 'optimal')
    rt = fitOptimal(fold.profTotal, Wtot, W2tot, WXtot, Bc, sGain, cst);
    nUsedT = rt.nUsed;
else
    cht = combineChannels(fold.profTotal, Wtot, W2tot, WXtot, Bc);
    rt = fitOne(cht, cst);
    nUsedT = cht.nUsed;
end
toa.total = struct('phase', rt.tau, 'phaseErr', rt.tauErr, ...
    'timeOffset', rt.tau / info_fold.f0, 'timeOffsetErr', rt.tauErr / info_fold.f0, ...
    'amp', rt.b, 'ampErr', rt.bErr, 'snr', rt.b / rt.bErr, 'redChi2', rt.redChi2, ...
    'flagChi2', rt.redChi2 > opts.MaxRedChi2, 'nChanUsed', nUsedT);

info = struct('method', 'FFTFIT (Fourier-domain template matching), Newton-refined', ...
    'noiseModel', model, 'Bnoise', opts.Bnoise, 'binDt', binDt, 'harmonics', K, ...
    'upsample', opts.Upsample, 'template', template, 'NBin', N, 'maxRedChi2', opts.MaxRedChi2, ...
    'toaConvention', 'TOA = tRef + phase/fRef; time of template phase 0 at the dedispersion reference frequency', ...
    'weighting', weighting);
if strcmp(weighting, 'optimal')
    info.method = ['weighted multi-channel fit (a_c + b*s_c*T, weights 1/var, iterated), ' ...
                   'FFT correlation + bounded search; sandwich error bars'];
    info.channelGain = sGain;
end

if opts.Verbose
    v = toa.valid;
    fprintf(['estimateTOA: %d/%d sub-int(s) fitted (%s noise, %s weighting), median SNR %.1f, ' ...
             'median TOA error %.3f us, median red. chi2 %.3f\n'], nnz(v), nSub, model, weighting, ...
        median(toa.snr(v)), median(toa.toaErr(v))*1e6, median(toa.redChi2(v)));
    if any(toa.flagChi2)
        fprintf('  %d sub-int(s) flagged: red. chi2 > %g (profile is not template + noise)\n', ...
            nnz(toa.flagChi2), opts.MaxRedChi2);
    end
    fprintf('  total fold: phase offset %+.3e turns = %+.4f +- %.4f us, SNR %.1f\n', ...
        rt.tau, toa.total.timeOffset*1e6, toa.total.timeOffsetErr*1e6, toa.total.snr);
end

end


% =====================================================================================
function r = fitOne(ch, c)
%FITONE  FFTFIT of one (channel-combined) profile; see the header of estimateTOA.
% ch from combineChannels: p (sum), Pc (channel profiles), W, W2, WX, B, have.
N = c.N; k = c.k; S = c.S; Sk = c.Sk; sumS2 = c.sumS2; ksg = c.ksg;
p = ch.p; Pc = ch.Pc; W = ch.W; W2 = ch.W2; WX = ch.WX; have = ch.have;
nLag = size(WX, 3);
        r = struct('ok', false, 'tau', NaN, 'tauErr', NaN, 'b', NaN, 'bErr', NaN, ...
                   'a', NaN, 'redChi2', NaN);
        if ~any(have), return; end
        if ~all(have)                                       % fill gaps circularly
            idx = find(have);
            xi  = [idx - N; idx; idx + N];
            yi  = repmat(p(idx), 3, 1);
            p(~have) = interp1(xi, yi, find(~have), 'linear');
            if size(Pc, 2) > 1                              % the channels likewise
                Pc(~have, :) = interp1(xi, repmat(Pc(idx, :), 3, 1), find(~have), 'linear');
            end
        end
        P = fft(p);
        X = P(k + 1) .* conj(Sk);

        % Coarse: band-limited upsampled cross-correlation
        M = c.M;
        spec = zeros(M, 1);
        spec(k + 1) = X;
        cc = real(ifft(spec));
        [~, im] = max(cc);
        tau = (im - 1) / M;
        if tau >= 0.5, tau = tau - 1; end

        % Refine: Newton on C'(tau) = 0
        w1 = 2*pi*k;
        for it = 1:50
            e  = exp(1i * w1 * tau);
            c1 = real(sum(1i * w1 .* X .* e));
            c2 = real(sum(-(w1.^2) .* X .* e));
            if c2 >= 0, return; end                         % not a maximum
            d = -c1 / c2;
            d = max(min(d, 1/M), -1/M);                     % stay near the coarse peak
            tau = tau + d;
            if abs(d) < 1e-14, break; end
        end
        e  = exp(1i * w1 * tau);
        c2 = real(sum(-(w1.^2) .* X .* e));
        b  = real(sum(X .* e)) / sumS2;
        a  = (real(P(1)) - b * real(S(1))) / N;
        tau = mod(tau + 0.5, 1) - 0.5;

        % Fitted model on the phase grid (sum and, for the noise, per channel:
        % a_c, b_c by the same projection, tau fixed; they add up to a, b)
        shape = real(ifft(S .* exp(-2i*pi*ksg*tau)));
        m = a + b * shape;
        onP = (m - a) > c.onFrac * max(m - a);
        if size(Pc, 2) == 1
            mC = m;
        else
            PC = fft(Pc);
            bC = real(sum(PC(k + 1, :) .* conj(Sk) .* e, 1)) / sumS2;
            aC = (real(PC(1, :)) - bC * real(S(1))) / N;
            mC = aC + bC .* shape;                          % [N x nChan used]
        end

        % Noise model: variance of each profile bin and covariance with bin j+d,
        % summed over the (independent) channels, lags d = 1..nLag
        if strcmp(c.model, 'radiometer')
            sTB = abs(mC) ./ sqrt(ch.B * c.binDt);          % per time bin, per channel
            v   = sum(sTB.^2 .* W2 ./ W.^2, 2);
            cv  = zeros(N, nLag);
            for lg = 1:nLag
                cv(:, lg) = sum(sTB .* circshift(sTB, -lg) .* WX(:, :, lg) ./ ...
                    (W .* circshift(W, -lg)), 2);
            end
        else
            res = p - m;
            off = ~onP & have;
            if nnz(off) < 10, return; end
            v0 = mean(res(off).^2);
            v  = v0 * ones(N, 1);
            cv = zeros(N, nLag);
            for lg = 1:nLag
                offN = off & circshift(off, -lg);
                resN = circshift(res, -lg);
                cv(:, lg) = mean(res(offN) .* resN(offN)) * ones(N, 1);
            end
        end
        v(~have | ~isfinite(v)) = 0;
        for lg = 1:nLag
            cl = cv(:, lg);
            cl(~have | ~circshift(have, -lg) | ~isfinite(cl)) = 0;
            cv(:, lg) = cl;
        end

        % tau: C'(tau) = sum_j d_j p_j
        D = zeros(N, 1);
        D(k + 1) = 1i * w1 .* conj(Sk) .* exp(1i * w1 * tau);
        dvec = real(fft(D));
        varC1 = sum(dvec.^2 .* v);
        % b: C(tau) = sum_j c_j p_j
        Cv = zeros(N, 1);
        Cv(k + 1) = conj(Sk) .* exp(1i * w1 * tau);
        cvec = real(fft(Cv));
        varC0 = sum(cvec.^2 .* v);
        for lg = 1:nLag
            varC1 = varC1 + 2*sum(dvec .* circshift(dvec, -lg) .* cv(:, lg));
            varC0 = varC0 + 2*sum(cvec .* circshift(cvec, -lg) .* cv(:, lg));
        end

        use = onP & have & v > 0;
        r.ok      = true;
        r.tau     = tau;
        r.tauErr  = sqrt(max(varC1, 0)) / abs(c2);
        r.b       = b;
        r.bErr    = sqrt(max(varC0, 0)) / sumS2;
        r.a       = a;
        r.redChi2 = sum((p(use) - m(use)).^2 ./ v(use)) / max(nnz(use) - 3, 1);
    end


% =====================================================================================
function r = fitOptimal(Pc, W, W2, WX, B, sGain, c)
%FITOPTIMAL  Weighted multi-channel fit prof_cj = a_c + b*s_c*T(phi_j - tau) + noise.
% Pc [N x nChan] channel profiles; W, W2 [N x nW], WX [N x nW x nLag] fold weights
% (nW = 1: shared by all channels); B = Bnoise (scalar or per channel). Weights 1/var
% from the noise model with the current fit, iterated to convergence; see the header.
N = c.N; nChan = size(Pc, 2);
r = struct('ok', false, 'tau', NaN, 'tauErr', NaN, 'b', NaN, 'bErr', NaN, 'a', NaN, ...
           'redChi2', NaN, 'nUsed', 0, 'aC', nan(1, nChan));
if size(W, 2) == 1                                  % shared weights
    W = repmat(W, 1, nChan); W2 = repmat(W2, 1, nChan); WX = repmat(WX, 1, nChan, 1);
end
hasC = W > 0;
use  = any(hasC, 1) & sGain > 0;
if ~any(use), return; end
Pc = Pc(:, use); W = W(:, use); W2 = W2(:, use); WX = WX(:, use, :); hasC = hasC(:, use);
s = sGain(use);
if numel(B) > 1, B = B(use); end
Pc(~hasC) = 0;
shapeAt = @(t) real(ifft(c.Sbl .* exp(-2i*pi*c.ksg*t)));

% iterate: noise model of the current fit -> weights 1/var -> fit
aC = sum(W .* Pc, 1) ./ sum(W, 1);                  % start: channel level, no pulse
b = 0; tau = 0; conv = false;
for it = 1:c.maxIt
    mC = repmat(aC, N, 1);
    if b ~= 0, mC = mC + b * s .* shapeAt(tau); end
    [v, ~] = chanNoise(Pc, mC, W, W2, WX, B, hasC, c, it == 1);
    w = zeros(size(v)); okv = hasC & v > 0; w(okv) = 1 ./ v(okv);
    [tauN, bN, aCN, ok] = wlsFit(Pc, w, s, c, shapeAt);
    if ~ok, return; end
    dTau = abs(mod(tauN - tau + 0.5, 1) - 0.5);
    conv = it > 1 && dTau < 1e-12 && abs(bN - b) <= 1e-12 * abs(bN);
    tau = tauN; b = bN; aC = aCN;
    if conv, break; end
end

% final model, its noise (all lags) and the weights it implies
Tt = shapeAt(tau);
dT = real(ifft(c.Sbl .* (-2i*pi*c.ksg) .* exp(-2i*pi*c.ksg*tau)));   % dT/dtau
mC = aC + b * s .* Tt;
[v, cv] = chanNoise(Pc, mC, W, W2, WX, B, hasC, c, false);
w = zeros(size(v)); okv = hasC & v > 0; w(okv) = 1 ./ v(okv);

% sandwich error bars: parameters (tau, b), baselines a_c profiled out per channel
Wc  = sum(w, 1);
G1  = b * s .* (dT - sum(w .* dT, 1) ./ Wc);        % d model / d tau, centred per channel
G2  = s .* (Tt - sum(w .* Tt, 1) ./ Wc);            % d model / d b
A   = [sum(w .* G1.^2, 'all'), sum(w .* G1 .* G2, 'all'); 0, sum(w .* G2.^2, 'all')];
A(2, 1) = A(1, 2);
x1 = w .* G1; x2 = w .* G2;
Bm = [quadForm(x1, x1, v, cv), quadForm(x1, x2, v, cv); 0, quadForm(x2, x2, v, cv)];
Bm(2, 1) = Bm(1, 2);
Cov = (A \ Bm) / A;                                % A^-1 B A^-1 (A symmetric)

% reduced chi^2 over the on-pulse bins (diagonal variances, as for 'equal')
onP = Tt - min(Tt) > c.onFrac * (max(Tt) - min(Tt));
sel = repmat(onP, 1, size(Pc, 2)) & okv;
nT  = nnz(sel);
res = Pc - mC;

sumS = sum(s);
r.ok      = true;
r.tau     = mod(tau + 0.5, 1) - 0.5;
r.tauErr  = sqrt(max(Cov(1, 1), 0));
r.b       = b * sumS;
r.bErr    = sqrt(max(Cov(2, 2), 0)) * sumS;
r.a       = sum(aC);
r.redChi2 = sum(res(sel).^2 .* w(sel)) / max(nT - 3, 1);
r.nUsed   = nnz(use);
r.aC(use) = aC;
r.converged = conv;
end


% =====================================================================================
function [tau, b, aC, ok] = wlsFit(Pc, w, s, c, shapeAt)
%WLSFIT  For weights w [N x nU]: tau maximizing N(tau)/sqrt(Dn(tau)), then b and a_c.
% With channel-weighted means pbar_c, Tbar_c(tau): b = N/Dn,
%   N(tau)  = sum_j T(phi_j - tau) u_j,   u_j = sum_c s_c w_cj (p_cj - pbar_c)
%   Dn(tau) = sum_j T^2 Om_j - sum_c s_c^2 (sum_j w_cj T)^2 / W_c,  Om_j = sum_c s_c^2 w_cj
N = c.N; M = c.M; nU = size(Pc, 2);
tau = NaN; b = NaN; aC = nan(1, nU);
Wc = sum(w, 1);
ok = all(Wc > 0);
if ~ok, return; end
pbar = sum(w .* Pc, 1) ./ Wc;
u  = sum(s .* w .* (Pc - pbar), 2);
Om = sum(s.^2 .* w, 2);
% coarse: all shifts tau = m/M at once (circular correlations on the fine grid)
idx = (0:N-1).' * (M / N) + 1;
uf = zeros(M, 1); uf(idx) = u;
Of = zeros(M, 1); Of(idx) = Om;
wf = zeros(M, nU); wf(idx, :) = w;
Nn = real(ifft(fft(uf) .* conj(c.TfF)));
D1 = real(ifft(fft(Of) .* conj(c.T2fF)));
WT = real(ifft(fft(wf) .* conj(c.TfF)));
Dn = D1 - sum((s.^2 ./ Wc) .* WT.^2, 2);
F  = Nn ./ sqrt(max(Dn, realmin));
[~, im] = max(F);
t0 = (im - 1) / M;
% refine: bounded search of the exact objective around the best grid point
negF = @(t) -objective(t, u, Om, w, s, Wc, shapeAt);
tau = fminbnd(negF, t0 - 1/M, t0 + 1/M, optimset('TolX', 1e-13));
% polish: Newton on dF/dtau = 0 (analytic slope, numerical curvature); a search on F
% values alone stops at ~1e-10 turns
dS = @(t) real(ifft(c.Sbl .* (-2i*pi*c.ksg) .* exp(-2i*pi*c.ksg*t)));
gF = @(t) slope(t, u, Om, w, s, Wc, shapeAt, dS);
hs = 1e-6;
for it = 1:20
    g  = gF(tau);
    cF = (gF(tau + hs) - gF(tau - hs)) / (2*hs);
    if ~(cF < 0), break; end                        % not at a maximum: keep the search value
    step = -g / cF;
    if abs(step) > 1/M, break; end
    tau = tau + step;
    if abs(step) < 1e-15, break; end
end
Tt  = shapeAt(tau);
WTt = sum(w .* Tt, 1);
Dt  = Om.' * Tt.^2 - sum(s.^2 ./ Wc .* WTt.^2);
b   = (u.' * Tt) / Dt;
aC  = pbar - b * s .* WTt ./ Wc;
end

function F = objective(t, u, Om, w, s, Wc, shapeAt)
Tt = shapeAt(t);
WTt = sum(w .* Tt, 1);
Dn = Om.' * Tt.^2 - sum(s.^2 ./ Wc .* WTt.^2);
F = (u.' * Tt) / sqrt(max(Dn, realmin));
end

function g = slope(t, u, Om, w, s, Wc, shapeAt, dS)
% dF/dtau of F = N / sqrt(Dn) (see wlsFit)
Tt = shapeAt(t); dT = dS(t);
WTt = sum(w .* Tt, 1); WdT = sum(w .* dT, 1);
Nn  = u.' * Tt;            dN  = u.' * dT;
Dn  = Om.' * Tt.^2 - sum(s.^2 ./ Wc .* WTt.^2);
dDn = 2 * Om.' * (Tt .* dT) - 2 * sum(s.^2 ./ Wc .* WTt .* WdT);
g = dN / sqrt(Dn) - 0.5 * Nn * dDn / Dn^1.5;
end


% =====================================================================================
function [v, cv] = chanNoise(Pc, mC, W, W2, WX, B, hasC, c, first)
%CHANNOISE  Variance v [N x nU] and lag covariances cv [N x nU x nLag] of every
% channel's profile bins: radiometer model with the model power mC, or (offpulse)
% the residual statistics of each channel (uniform before the first fit).
[N, nU] = size(Pc); nLag = size(WX, 3);
cv = zeros(N, nU, nLag);
if strcmp(c.model, 'radiometer')
    sTB = abs(mC) ./ sqrt(B * c.binDt);
    v = sTB.^2 .* W2 ./ W.^2;
    for d = 1:nLag
        cv(:, :, d) = sTB .* circshift(sTB, -d, 1) .* WX(:, :, d) ./ (W .* circshift(W, -d, 1));
    end
elseif first
    v = ones(N, nU);
else
    res = Pc - mC;
    sh  = mC - min(mC, [], 1);
    off = hasC & sh <= c.onFrac * max(sh, [], 1);
    v = zeros(N, nU);
    for ch = 1:nU
        o = off(:, ch);
        if nnz(o) < 10, continue; end
        v(:, ch) = mean(res(o, ch).^2);
        for d = 1:nLag
            oN = o & circshift(o, -d);
            rN = circshift(res(:, ch), -d);
            cv(:, ch, d) = mean(res(oN, ch) .* rN(oN));
        end
    end
end
v(~hasC | ~isfinite(v)) = 0;
for d = 1:nLag
    cl = cv(:, :, d);
    cl(~hasC | ~circshift(hasC, -d, 1) | ~isfinite(cl)) = 0;
    cv(:, :, d) = cl;
end
end

function q = quadForm(x, y, v, cv)
%QUADFORM  sum_c x_c' Sigma_c y_c for per-channel covariances given by the diagonal v
% and the lag-d covariances cv(:, :, d) = cov(bin j, bin j+d) (circular).
q = sum(x .* y .* v, 'all');
for d = 1:size(cv, 3)
    q = q + sum((x .* circshift(y, -d, 1) + circshift(x, -d, 1) .* y) .* cv(:, :, d), 'all');
end
end
