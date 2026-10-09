%DEMOBINNINGMATCHEDFILTER  Demo: matched filter at full rate vs after binning vs after decimation.
%{
Run from the PulsarSimMatlab folder: run('tests/demoBinningMatchedFilter.m'), or
section by section (Ctrl+Enter). ~15 s, ~1 GB of memory. A demo for
explaining, not a unit test (no PASS/FAIL); it leaves its variables in the
workspace.

The signal: like the pulsar, noise that gets stronger during the pulse.
Complex baseband (I/Q), fs = 100 MHz, independent samples (white over the
band), receiver noise power 1, pulsar power a*g(t) on top: a slow Gaussian
envelope g, period 10 ms, FWHM 0.5 ms, peak a (a few % of the noise).
Detected power p = |x|^2. Matched filter = correlation of the detected
power with the zero-mean pulse shape, at every lag (periodic: 10 whole turns).

Three ways:
  1. FULL RATE: matched filter on p at 100 MHz (1e7 samples).
  2. BINNED:    average p over 10 kHz bins (10,000 samples per bin; what
                detectPower / detectChannels do), matched filter at 10 kHz.
  3. DECIMATED: keep one p sample per 10 kHz bin (no averaging), matched
                filter at 10 kHz.

Expected (radiometer / detection theory, weak signal):
  1 and 2 the same SNR: way 2 IS way 1 with a template that is constant
  within each bin (sum over bins of mean(p)*T = (1/L) * sum over samples of
  p*T). The only loss is the pulse smearing over a 100 us bin (bin = FWHM/5
  here: ~0.5 %). Way 2 is ~100x faster (binning included; printed with
  operation counts and memory).
  3 loses a factor sqrt(10,000) = 100 in SNR: it uses 1 sample in 10,000.
  Binning = boxcar anti-alias filter + decimation; way 3 is decimation
  WITHOUT the filter, so all the noise above 5 kHz folds into the band.
  SNR = a * sqrt(sum of independent samples used * <(g - <g>)^2>).
%}

%% Parameters
fs   = 100e6;          % [Hz] complex sample rate (= bandwidth: independent samples)
P    = 10e-3;          % [s] pulse period
fwhm = 0.5e-3;         % [s] FWHM of the pulse (power envelope)
Tobs = 0.1;            % [s] observation: 10 turns
fBin = 10e3;           % [Hz] bin rate for ways 2 and 3
dTarget = 20;          % SNR of way 1 we aim for (sets the pulse amplitude a)
nTr  = 100;            % Monte Carlo trials for the measured SNR
rng(1);

N  = round(fs*Tobs);   % 1e7 samples
L  = round(fs/fBin);   % 10,000 samples per bin
nb = N/L;              % 1000 bins
t  = ((0:N-1).' + 0.5) / fs;
sg = fwhm / (2*sqrt(2*log(2)));
g  = exp(-0.5*((mod(t, P) - P/2)/sg).^2);   % pulse shape, peak 1, centred in each turn

% Templates (zero mean: the constant noise level drops out of the correlation)
t1 = g - mean(g);                           % 1. full rate
G2 = mean(reshape(g, L, nb), 1).';          % 2. pulse as seen by the bin average
t2 = G2 - mean(G2);
g3 = g(L/2 : L : end);                      % 3. pulse at the kept samples
t3 = g3 - mean(g3);

% Pulse amplitude for SNR dTarget at full rate, and the predicted SNR of each way.
% Per sample var(|x|^2) = (mean power)^2 = 1 (complex Gaussian), so the
% correlation at the true lag has signal a*sum(g.*t) = a*sum(t.^2) and noise
% std sqrt(sum(t.^2) * var per value): var 1 (ways 1, 3) or 1/L (way 2).
a  = dTarget / sqrt(sum(t1.^2));
dPred = [a*sqrt(sum(t1.^2)), a*sqrt(L*sum(t2.^2)), a*sqrt(sum(t3.^2))];
fprintf('pulse peak a = %.4f of the noise power; predicted SNR: full %.2f, binned %.2f, decimated %.3f\n', ...
    a, dPred);

%% One realization: signal, the three detected-power series
pw = 1 + a*g;                                       % power: noise 1 + pulsar
x  = sqrt(pw/2) .* (randn(N, 1) + 1i*randn(N, 1)); % complex Gaussian with that power
p  = abs(x).^2;                                     % detected power, 100 MHz
p2 = mean(reshape(p, L, nb), 1).';                  % 2. binned, 10 kHz
p3 = p(L/2 : L : end);                              % 3. decimated, 10 kHz

%% Matched filter (correlation over all lags, via FFT; periodic: whole turns)
corrF = @(y, tm) real(ifft(fft(y) .* conj(fft(tm))));
c1 = corrF(p, t1);
c2 = corrF(p2, t2);
c3 = corrF(p3, t3);
% in units of the noise std of each output (so the peak height = SNR)
c1 = c1 / sqrt(sum(t1.^2));
c2 = c2 / sqrt(sum(t2.^2) / L);
c3 = c3 / sqrt(sum(t3.^2));
fprintf('this realization, output at lag 0: full %.2f, binned %.2f, decimated %.2f (noise std = 1)\n', ...
    c1(1), c2(1), c3(1));

%% Computational cost (tic/toc of one run each)
tPre = zeros(1, 3); tMF = zeros(1, 3);           % binning / decimation, matched filter
tic; corrF(p, t1);                           tMF(1)  = toc;
tic; pb = mean(reshape(p, L, nb), 1).';      tPre(2) = toc;
tic; corrF(pb, t2);                          tMF(2)  = toc;
tic; pd = p(L/2 : L : end);                  tPre(3) = toc;
tic; corrF(pd, t3);                          tMF(3)  = toc;
nPts = [N, nb, nb];                               % points in the matched filter
fftOps    = 3 * 5 * nPts .* log2(nPts);           % ~3 FFTs of n points, ~5 n log2(n) flops each
preOps    = [0, N, 0];                            % binning: one addition per sample
directOps = nPts.^2;                              % direct (time-domain) correlation, all lags
arrayMB   = nPts * 16 / 1e6;                      % one complex double array of n points
tTot = tPre + tMF;
fprintf('Computational cost of the matched filter over all lags:\n');
fprintf('  %-28s %9s %10s %10s %11s %12s %10s\n', 'way', 'points', 'prep [s]', 'MF [s]', ...
    'ops (FFT)', 'ops (direct)', 'array [MB]');
names = {'1. full rate (100 MHz)', '2. binned (10 kHz bins)', '3. decimated (1 in 10,000)'};
for k = 1:3
    fprintf('  %-28s %9d %10.2g %10.2g %11.2g %12.2g %10.3g\n', names{k}, nPts(k), tPre(k), ...
        tMF(k), fftOps(k) + preOps(k), directOps(k), arrayMB(k));
end
fprintf('  binned is %.0fx faster than full rate (binning included)\n', tTot(1)/tTot(2));

%% Monte Carlo: measured SNR = mean / std of the output at the true lag
s = zeros(nTr, 3);
fprintf('Monte Carlo, %d trials ', nTr);
for r = 1:nTr
    xr = sqrt(pw/2) .* (randn(N, 1) + 1i*randn(N, 1));
    pr = abs(xr).^2;
    s(r, :) = [pr.'*t1, mean(reshape(pr, L, nb), 1)*t2, pr(L/2 : L : end).'*t3];
    if mod(r, 10) == 0, fprintf('.'); end
end
dMeas = mean(s) ./ std(s);
dErr  = sqrt((1 + dMeas.^2/2) / nTr);              % 1-sigma of an SNR estimate
fprintf('\n');
nUsed = [N, N, nb];
for k = 1:3
    fprintf('  %-28s samples used %9d: SNR measured %6.2f +- %.2f, predicted %6.2f\n', ...
        names{k}, nUsed(k), dMeas(k), dErr(k), dPred(k));
end
fprintf('  binned / full: %.3f (predicted %.3f: smearing over a %.0f us bin, FWHM %.0f us)\n', ...
    dMeas(2)/dMeas(1), dPred(2)/dPred(1), 1e6/fBin, fwhm*1e6);
fprintf('  decimated / full: %.4f (predicted %.4f ~ 1/sqrt(L) = %.4f)\n', ...
    dMeas(3)/dMeas(1), dPred(3)/dPred(1), 1/sqrt(L));

%% Figure
lagOf = @(n, rate) mod((0:n-1).' / rate + P/2, P) - P/2;   % lag folded into one turn
fig = figure('Name', 'Matched filter: full rate vs binned vs decimated', 'Color', 'w', ...
    'Position', [80 80 1250 780]);
tl = tiledlayout(fig, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
show = 30e-3;                                       % first 30 ms of the power series
k1 = 1:500:round(show*fs);                          % display only: every 500th sample
ax = nexttile(tl); plot(ax, t(k1)*1e3, p(k1), '.', 'MarkerSize', 2); hold(ax, 'on');
plot(ax, t(k1)*1e3, pw(k1), 'r-', 'LineWidth', 1.2);
title(ax, '1. detected power, 100 MHz (every 500th sample shown)'); xlabel(ax, 't [ms]'); ylabel(ax, '|x|^2');
legend(ax, 'samples', 'true power 1 + a g', 'Location', 'northeast');
tb = ((0:nb-1).' + 0.5) / fBin; kb = tb < show;
ax = nexttile(tl); plot(ax, tb(kb)*1e3, p2(kb), 'o-', 'MarkerSize', 3); hold(ax, 'on');
plot(ax, t(k1)*1e3, pw(k1), 'r-', 'LineWidth', 1.2);
title(ax, '2. binned: mean of 10,000 samples per 100 \mus'); xlabel(ax, 't [ms]'); ylabel(ax, 'mean |x|^2');
t3s = t(L/2 : L : end); k3 = t3s < show;
ax = nexttile(tl); plot(ax, t3s(k3)*1e3, p3(k3), 'o-', 'MarkerSize', 3); hold(ax, 'on');
plot(ax, t(k1)*1e3, pw(k1), 'r-', 'LineWidth', 1.2);
title(ax, '3. decimated: 1 sample per 100 \mus'); xlabel(ax, 't [ms]'); ylabel(ax, '|x|^2');
lag1 = lagOf(N, fs); [lag1s, o1] = sort(lag1); o1 = o1(1:100:end); lag1s = lag1s(1:100:end);
lag2 = lagOf(nb, fBin); [lag2s, o2] = sort(lag2);
yl = [-4, dTarget + 6];
ax = nexttile(tl); plot(ax, lag1s*1e3, c1(o1), '-');
title(ax, sprintf('1. matched filter, full rate: peak %.1f', c1(1))); ylim(ax, yl);
xlabel(ax, 'lag [ms]'); ylabel(ax, 'MF output normalised by noise level'); grid(ax, 'on');
ax = nexttile(tl); plot(ax, lag2s*1e3, c2(o2), 'o-', 'MarkerSize', 3);
title(ax, sprintf('2. matched filter, binned: peak %.1f', c2(1))); ylim(ax, yl);
xlabel(ax, 'lag [ms]'); ylabel(ax, 'MF output normalised by noise level'); grid(ax, 'on');
ax = nexttile(tl); plot(ax, lag2s*1e3, c3(o2), 'o-', 'MarkerSize', 3);
title(ax, sprintf('3. matched filter, decimated: peak %.1f', c3(1))); ylim(ax, yl);
xlabel(ax, 'lag [ms]'); ylabel(ax, 'MF output normalised by noise level'); grid(ax, 'on');
title(tl, sprintf(['Same data, three ways. SNR full / binned / decimated: %.1f / %.1f / %.2f (%d trials); ' ...
    'matched-filter time %.3g / %.2g / %.2g s'], dMeas, nTr, tTot));
