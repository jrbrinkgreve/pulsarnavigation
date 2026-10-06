function expBlankingVariance()
%EXPBLANKINGVARIANCE  Experiment 3a: detected power after blanking + per-channel dedispersion.
%{
Run from the PulsarSimMatlab folder: run('tests/expBlankingVariance.m'). ~15 s.

Question (unit 3 design): RFI excision sets input samples of a channel to 0
before coherent dedispersion. Each dedispersed output sample y(n) =
sum_l h(l) x(n-l) then has a valid fraction
    w(n) = sum_l |h(l)|^2 keep(n-l) / sum_l |h(l)|^2,
and E|y|^2 = w * m exactly. What is the VARIANCE of a detected time bin
(mean of |y|^2 over 4 samples) as a function of its weight W (mean of w)?
  - "w^2 model": var = var0 * W^2  (the gap only scales the amplitude)
  - "w model":   var = var0 * W    (whole samples lost)
  - exact: for Gaussian x the bin samples have covariance C = H*K*H'
    (K = diag(keep)), so var = (1/n^2) * sum_ab |C_ab|^2, and the
    covariance with the next bin is (1/n^2) * sum_{a in k, b in k+1} |C_ab|^2.
var0 = variance of an unblanked bin (includes the narrow-channel factor nu0).

Setup: noise only (unit variance complex white), bottom channel (centre
1.2015625 GHz, the longest intra-channel sweep), channel rate 4.1667 MHz,
applyInverseDispersion with the channel's band and reference = channel top,
DM 5 (sweep 75 us) and DM 100 (1.5 ms). Masks:
  radar   4 us (17 samples) blanked at PRF 373 Hz
  random  50 us (208 samples) blanks at random (Poisson) starts, ~10 % blanked
  gaps    1 ms (4167 samples) blanked every 5 ms (20 %)
The exact covariance is cheap: C_{a,a+tau} = sum_l h(l) conj(h(l+tau)) keep(a-l)
is a convolution of the mask with q_tau(l) = h(l)*conj(h(l+tau)), one FFT per
lag tau = 0..2n-1, so it is computed for EVERY bin (cross-checked against the
direct matrix on 300 bins).
Per mask: time bins grouped by W; measured mean P/W (expect 1) and variance
of P - W*m against the exact prediction and the two simple models, all on
the same bins. Checks: mean = W, exact = measured (all bins and per group,
4 sigma), conv = matrix method. The w^2 / w comparison is the result of the
experiment. Saves data/mc/expBlankingVariance.mat.
%}

root = fullfile(fileparts(mfilename('fullpath')), '..');
addpath(fullfile(root, 'functions'));
tmp = fullfile(tempdir, 'expBlankingVariance');
if ~isfolder(tmp), mkdir(tmp); end

fsC = 800e6/192; dF = 3.125e6; fc = 1.2e9 + dF/2;
band = fc + [-1 1]*dF/2;
nPerBin = 4;
N = 2^22;                                       % ~1 s per case
DMs = [5 100];
rs = RandStream('mt19937ar', 'Seed', 2026);
x = (randn(rs, N, 1) + 1i*randn(rs, N, 1)) / sqrt(2);

% Masks (keep = true for valid samples), 0-based sample n
prf = round(fsC / 373);
keepRadar = true(N, 1);
for s = 1000 : prf : N - 17, keepRadar(s+1 : s+17) = false; end
keepRandom = true(N, 1);
s = 0;
while true
    s = s + round(-log(rand(rs)) * 208 / 0.1 * 0.9);    % mean gap so that ~10 % is blanked
    if s + 208 > N, break; end
    keepRandom(s+1 : s+208) = false;
    s = s + 208;
end
keepGaps = true(N, 1);
for s = 2000 : 20833 : N - 4167, keepGaps(s+1 : s+4167) = false; end
masks = struct('name', {'radar 4 us', 'random 50 us', 'gaps 1 ms'}, ...
               'keep', {keepRadar, keepRandom, keepGaps});

fails = {};
results = struct([]);
for DM = DMs
    % ---- Impulse response of this channel's dedispersion filter ----------------------
    Nimp = 2^16; n0 = Nimp/2;
    imp = zeros(Nimp, 1); imp(n0 + 1) = 1;
    writeCF32(fullfile(tmp, 'imp.dat'), imp);
    d = applyInverseDispersion(fullfile(tmp, 'imp.dat'), fullfile(tmp, 'imp_out.dat'), fsC, ...
        fc, DM, band(1), band(2), 'RefFreq', band(2), 'SaveInfo', false, 'Verbose', false);
    yi = readCF32(fullfile(tmp, 'imp_out.dat'));
    nF = d.nFuture; nP = d.nPast;
    hv = yi(n0 + 1 + (-nF:nP));                    % h(l), l = -nF..nP (y(n) = sum h(l) x(n-l))
    m0 = sum(abs(hv).^2);                           % unblanked mean power (unit-variance input)
    g  = abs(hv).^2 / m0;
    gC = zeros(N, 1); gC(mod((-nF:nP).', N) + 1) = g; % circular kernel at lag l
    Gf = fft(gC);
    edge = nP + nF + 8;                             % stay clear of the file edges

    for im = 1:numel(masks)
        keep = masks(im).keep;
        writeCF32(fullfile(tmp, 'in.dat'), x .* keep);
        applyInverseDispersion(fullfile(tmp, 'in.dat'), fullfile(tmp, 'out.dat'), fsC, fc, DM, ...
            band(1), band(2), 'RefFreq', band(2), 'SaveInfo', false, 'Verbose', false);
        y = readCF32(fullfile(tmp, 'out.dat'));
        w = real(ifft(fft(double(keep)) .* Gf));      % w(n) = sum_l g(l) keep(n-l)

        % Detected time bins (normalized to the unblanked mean power)
        nb = floor(N / nPerBin);
        P  = mean(reshape(abs(double(y(1:nb*nPerBin))).^2, nPerBin, nb), 1).' / m0;
        W  = mean(reshape(w(1:nb*nPerBin), nPerBin, nb), 1).';
        kOK = (ceil(edge/nPerBin) + 1 : nb - ceil(edge/nPerBin) - 2).';
        r  = P - W;                                     % E[P] = W exactly

        % Exact second moments for ALL bins: C_{a,a+tau} = (q_tau * keep)(a),
        % q_tau(l) = h(l)*conj(h(l+tau)): one FFT convolution of the mask per lag
        [vEx, cEx, var0Ex, cov0Ex] = exactBinStatsConv(hv, keep, nPerBin, nF, nP, m0, nb);
        % Cross-check against the direct 2n x 2n covariance matrix on 300 bins
        chk = kOK(1) + 1000 + (0:299).';
        [vG, cG] = exactBinStats(hv, keep, chk, nPerBin, nF, nP, m0);
        dChk = max(abs([vG - vEx(chk); cG - cEx(chk)])) / var0Ex;

        unb  = kOK(W(kOK) > 0.99999);
        var0 = mean(r(unb).^2);                         % measured unblanked bin variance
        allR = mean(r(kOK).^2) / mean(vEx(kOK));        % all bins: measured / exact
        sAll = sqrt(2 / numel(kOK)) * 1.4;
        fprintf(['\nDM %g, mask "%s": %.1f %% blanked; filter span %.0f us; unblanked bin var ' ...
                 '%.4f (exact %.4f); conv vs matrix method max diff %.1e of var0\n'], DM, ...
            masks(im).name, 100*mean(~keep), (nF + nP)/fsC*1e6, var0, var0Ex, dChk);
        fprintf('  all %d bins: measured / exact variance %.4f +- %.4f\n', numel(kOK), allR, sAll);
        fprintf('  %-13s %7s %9s %9s %9s %9s %9s %8s\n', 'W range', 'bins', 'mean P/W', ...
            'var meas', 'exact', 'w^2 mod', 'w model', '+-1sig');
        if dChk > 1e-6, fails{end+1} = sprintf('conv method DM%g %s', DM, masks(im).name); end %#ok<SAGROW>
        if abs(allR - 1) > 4*sAll, fails{end+1} = sprintf('exact all DM%g %s', DM, masks(im).name); end %#ok<SAGROW>

        % Groups by W: same bins for measurement, exact and the two simple models
        edgesW = [0 0.02 0.25 0.5 0.75 0.9 0.97 0.995 0.99999 1.0001];
        for gI = 1:numel(edgesW) - 1
            sel = kOK(W(kOK) >= edgesW(gI) & W(kOK) < edgesW(gI+1));
            if numel(sel) < 200, continue; end
            vm  = mean(r(sel).^2);
            ex  = mean(vEx(sel));
            mPW = mean(P(sel)) / mean(W(sel));
            % scatter of a variance estimate: sqrt(2/n) for Gaussian r; P is a mean of 4
            % |y|^2 (excess kurtosis ~1.5 -> x1.3) and neighbouring bins correlate -> x1.4
            sv  = ex * sqrt(2 / numel(sel)) * 1.4;
            m2  = var0Ex * mean(W(sel).^2);
            m1  = var0Ex * mean(W(sel));
            fprintf('  %5.3f-%5.3f %7d %9.4f %9.5f %9.5f %9.5f %9.5f %8.5f\n', edgesW(gI), ...
                min(edgesW(gI+1), 1), numel(sel), mPW, vm, ex, m2, m1, sv);
            sm = sqrt(ex / numel(sel)) * 1.2 / max(mean(W(sel)), eps);
            if abs(mPW - 1) > 4*sm + 1e-6, fails{end+1} = sprintf('mean DM%g %s', DM, masks(im).name); end %#ok<SAGROW>
            if abs(vm - ex) > 4*sv, fails{end+1} = sprintf('exact DM%g %s', DM, masks(im).name); end %#ok<SAGROW>
            results(end+1).DM = DM; %#ok<SAGROW>
            results(end).mask = masks(im).name; results(end).Wrange = edgesW(gI:gI+1);
            results(end).nBins = numel(sel); results(end).meanPW = mPW;
            results(end).varMeas = vm; results(end).varExact = ex;
            results(end).varW2 = m2; results(end).varW1 = m1;
        end
        % Next-bin covariance: measured vs exact (all bins, and bins at blanking edges)
        kk = kOK(1:end-1); bl = kk(W(kk) < 0.97 & W(kk) > 0.02);
        cm  = mean(r(kk) .* r(kk+1));  cmB = mean(r(bl) .* r(bl+1));
        fprintf(['  next-bin covariance / var0: all bins %.4f (exact %.4f, unblanked %.4f); ' ...
                 'blanked-edge bins %.4f (exact %.4f; w-scaled unblanked %.4f)\n'], ...
            cm/var0Ex, mean(cEx(kk))/var0Ex, cov0Ex/var0Ex, cmB/var0Ex, mean(cEx(bl))/var0Ex, ...
            cov0Ex/var0Ex * mean(W(bl).*W(bl+1)));
    end
end

save(fullfile(root, 'data', 'mc', 'expBlankingVariance.mat'), 'results');
if isempty(fails)
    fprintf('\nexpBlankingVariance: mean = W and exact variance = measured in all groups: PASS\n');
else
    error('expBlankingVariance:failed', 'FAILED: %s', strjoin(unique(fails), ', '));
end
end


% =================================================================================
function [vEx, cEx, var0, cov0] = exactBinStatsConv(hv, keep, n, nF, nP, m0, nb)
%EXACTBINSTATSCONV  Exact variance of every time bin and covariance with the next bin
% (normalized by m0^2), for unit-variance complex Gaussian input blanked by keep:
%   C_{a,a+tau} = sum_l h(l) conj(h(l+tau)) keep(a-l) = (q_tau * keep)(a),
%   var_k  = (1/n^2) sum_{a,b in k} |C_ab|^2,
%   cov_k  = (1/n^2) sum_{a in k, b in k+1} |C_ab|^2,
% one FFT convolution per lag tau = 0..2n-1. var0/cov0: the same without blanking.
N  = numel(keep);
Kf = fft(double(keep));
lags = (-nF:nP).';
A  = zeros(nb*n, 2*n);                 % |C_tau(a)|^2 / m0^2, a = 1..nb*n (1-based)
cU = zeros(1, 2*n);                    % unblanked C_tau / m0
for tau = 0:2*n-1
    q = zeros(N, 1);
    q(mod(lags(1:end-tau), N) + 1) = hv(1:end-tau) .* conj(hv(1+tau:end));
    Ct = ifft(Kf .* fft(q));
    A(:, tau+1) = abs(Ct(1:nb*n)).^2 / m0^2;
    cU(tau+1) = sum(q) / m0;
end
vEx = zeros(nb, 1); cEx = zeros(nb, 1);
var0 = 0; cov0 = 0;
for tau = 0:n-1                        % pairs inside a bin: i = 0..n-1-tau, b = a + tau
    B = reshape(A(:, tau+1), n, nb);
    f = 1 + (tau > 0);                 % (a,b) and (b,a)
    vEx = vEx + f * sum(B(1:n-tau, :), 1).';
    var0 = var0 + f * (n - tau) * abs(cU(tau+1))^2;
end
for tau = 1:2*n-1                      % a in bin k (offset i), b in bin k+1: tau = n + j - i
    B = reshape(A(:, tau+1), n, nb);
    i = max(0, n - tau) : min(n-1, 2*n-1-tau);
    cEx = cEx + sum(B(i + 1, :), 1).';
    cov0 = cov0 + numel(i) * abs(cU(tau+1))^2;
end
vEx = vEx / n^2; cEx = cEx / n^2; var0 = var0 / n^2; cov0 = cov0 / n^2;
end


% =================================================================================
function [vEx, cEx] = exactBinStats(hv, keep, bins, n, nF, nP, m0)
%EXACTBINSTATS  Same quantities from the direct 2n x 2n sample covariance matrix
% C = H*K*H' (slow; cross-check of exactBinStatsConv on a few bins).
vEx = zeros(numel(bins), 1); cEx = zeros(numel(bins), 1);
for i = 1:numel(bins)
    a  = (bins(i) - 1)*n + (0:2*n-1).';             % 0-based samples of bins k and k+1
    j  = (a(1) - nP : a(end) + nF);                  % inputs that reach them
    lag = a - j;                                     % 2n x numel(j)
    in  = lag >= -nF & lag <= nP;
    Hm  = zeros(size(lag));
    Hm(in) = hv(lag(in) + nF + 1);
    Hm  = Hm .* sqrt(double(keep(j + 1))).';
    G   = Hm * Hm';                                  % 2n x 2n sample covariance
    A2  = abs(G).^2;
    vEx(i) = sum(sum(A2(1:n, 1:n))) / n^2 / m0^2;
    cEx(i) = sum(sum(A2(1:n, n+1:2*n))) / n^2 / m0^2;
end
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
