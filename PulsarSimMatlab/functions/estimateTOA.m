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
  fold, info_fold  from foldProfile. Several channels are combined into one
                   profile p = sum of the channel profiles (equal weights);
                   a channel with no data in some phase bins of a profile,
                   where other channels have data, is left out of that
                   profile (its absence would leave a dip). Channels with
                   partial data (0 < weight) stay in.
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
          nChanUsed (channels combined)
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

% ---- Per sub-integration -------------------------------------------------------------
nanv = nan(nSub, 1);
toa = struct('valid', false(nSub, 1), 'coverage', nanv, 'phase', nanv, ...
    'phaseErr', nanv, 'toa', nanv, 'toaErr', nanv, 'amp', nanv, 'ampErr', nanv, ...
    'baseline', nanv, 'snr', nanv, 'redChi2', nanv, 'flagChi2', false(nSub, 1), ...
    'tRef', fold.subint.tRef, 'fRef', fold.subint.fRef, 'turnRef', fold.subint.turnRef, ...
    'nChanUsed', zeros(nSub, 1));

for s = 1:nSub
    ch = combineChannels(reshape(fold.prof(:, s, :), N, nChan), ...
        reshape(fold.weight(:, s, :), N, nW), reshape(fold.weight2(:, s, :), N, nW), ...
        reshape(WXall(:, s, :, :), N, nW, nLag), Bc);
    toa.coverage(s)  = ch.coverage;
    toa.nChanUsed(s) = ch.nUsed;
    if toa.coverage(s) < opts.MinCoverage || toa.coverage(s) == 0
        continue
    end
    r = fitOne(ch, cst);
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
cht = combineChannels(fold.profTotal, reshape(sum(fold.weight, 2), N, nW), ...
    reshape(sum(fold.weight2, 2), N, nW), reshape(sum(WXall, 2), N, nW, nLag), Bc);
rt = fitOne(cht, cst);
toa.total = struct('phase', rt.tau, 'phaseErr', rt.tauErr, ...
    'timeOffset', rt.tau / info_fold.f0, 'timeOffsetErr', rt.tauErr / info_fold.f0, ...
    'amp', rt.b, 'ampErr', rt.bErr, 'snr', rt.b / rt.bErr, 'redChi2', rt.redChi2, ...
    'flagChi2', rt.redChi2 > opts.MaxRedChi2, 'nChanUsed', cht.nUsed);

info = struct('method', 'FFTFIT (Fourier-domain template matching), Newton-refined', ...
    'noiseModel', model, 'Bnoise', opts.Bnoise, 'binDt', binDt, 'harmonics', K, ...
    'upsample', opts.Upsample, 'template', template, 'NBin', N, 'maxRedChi2', opts.MaxRedChi2, ...
    'toaConvention', 'TOA = tRef + phase/fRef; time of template phase 0 at the dedispersion reference frequency');

if opts.Verbose
    v = toa.valid;
    fprintf(['estimateTOA: %d/%d sub-int(s) fitted (%s noise), median SNR %.1f, ' ...
             'median TOA error %.3f us, median red. chi2 %.3f\n'], nnz(v), nSub, model, ...
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
function ch = combineChannels(prof, W, W2, WX, B)
%COMBINECHANNELS  Channels used for one profile, their sum and their weights.
% prof [N x nChan] channel profiles; W, W2 [N x nW], WX [N x nW x nLag] their
% fold weights (nW = nChan, or 1 when all channels share them); B = Bnoise,
% scalar or [1 x nChan]. Exclusion rule: a channel is used only if it has data
% (W > 0) in every phase bin where any channel has data, so every bin of the
% sum holds the same channels (a missing channel would leave a dip).
nChan = size(prof, 2);
have  = any(W > 0, 2);
if size(W, 2) == 1
    use = true(1, nChan);
else
    use = all(W(have, :) > 0, 1);
    W = W(:, use); W2 = W2(:, use); WX = WX(:, use, :);
end
if numel(B) > 1, B = B(use); end
Pc = prof(:, use);
ok = any(use) && any(have);
ch = struct('p', sum(Pc, 2), 'Pc', Pc, 'W', W, 'W2', W2, 'WX', WX, 'B', B, ...
    'have', have & ok, 'nUsed', nnz(use) * ok, 'coverage', mean(have) * ok);
end
