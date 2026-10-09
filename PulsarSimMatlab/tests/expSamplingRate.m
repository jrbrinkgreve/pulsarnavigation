function expSamplingRate()
%EXPSAMPLINGRATE  Does a higher sample rate help a matched filter? (9 Oct 2026)
%{
Run from the PulsarSimMatlab folder: run('tests/expSamplingRate.m'). ~20 s.

Question (meeting 9 Oct): a matched filter on a noisy signal sampled at
100 kHz looks less noisy than at 1 kHz, so should a higher f_out improve
detection? And (professor) the criterion is signal energy over noise PSD,
E/N0, and the square law might give a gain beyond a linear system.
Answer in short: the matched-filter SNR is sqrt(2E/N0) with the N0 the
filter sees; with the physical N0 fixed (T_sys) and all data used, it does
not depend on the rate. Details: logs/2026-10-09_fout-study.md.

Experiment 1 - known pulse in additive white noise (linear).
  "Analog" rate fa = 1 MHz, white noise with a FIXED one-sided PSD
  N0 = 2*sigA^2/fa; pulse train P = 100 ms, FWHM 5 ms, 1 s. Signal at rate
  fs, matched filter, three set-ups:
    A  the same sigma per sample at every fs (noise generated at fs; the
       usual quick test): the PSD the filter sees, 2*sigma^2/fs, DROPS as fs
       rises -> SNR ~ sqrt(fs)
    B  fixed PSD, integrate-and-dump to fs (= detectPower / detectChannels)
       -> SNR constant
    C  fixed PSD, one analog sample kept per output sample (no anti-alias
       filter): all noise aliases into fs/2 -> SNR ~ sqrt(fs) up to B's
  Plus one realization of the matched-filter output versus lag at 1 kHz and
  100 kHz (case B: the same analog noise at both rates -> the same curve).
Experiment 2 - the pulsar case: noise-like signal, square law.
  Complex noise at B = 1 MHz (independent samples), power 1 + a*g(t). At
  bin rate f_out:
    i   square every sample, then average over the bin (= the pipeline)
    ii  square every sample, keep one per bin (detector + ADC, no averaging)
    iii average the voltages over the bin, then square
  Envelope-domain E/N0 (E_env ~ (pulse power)^2, N_env = 2*var/rate of the
  squared samples that enter): i = radiometer bound a*sqrt(B*int(g-<g>)^2 dt)
  at every f_out; ii, iii ~ sqrt(f_out), reaching it only at f_out = B.
Checks: every measured SNR within 4 sigma of its prediction (Monte Carlo
sigma of mean/std with nTr trials: sqrt((1 + d^2/2)/nTr)). Errors at the end
if any check fails. Saves data/mc/expSamplingRate.mat; shows one figure.
%}

root = fullfile(fileparts(mfilename('fullpath')), '..');
rngOld = rng; rng(1); restoreRng = onCleanup(@() rng(rngOld));
nTr = 300;                                    % Monte Carlo trials per point
mcSig = @(d) sqrt((1 + d.^2/2) / nTr);        % 1-sigma of an SNR estimate mean/std

% ---------------- common pulse ------------------------------------------------------
fa = 1e6; Tobs = 1; P = 0.1; fwhm = 5e-3;
Na = round(fa*Tobs);
t  = ((0:Na-1).' + 0.5) / fa;
sg = fwhm / (2*sqrt(2*log(2)));
g  = exp(-0.5*((mod(t, P) - P/2)/sg).^2);
rates = [1e3 2e3 5e3 1e4 5e4 1e5 5e5 1e6];   % divisors of fa
nR = numel(rates);

% ================= Experiment 1 =====================================================
sigA = 1;  N0 = 2*sigA^2/fa;                  % physical noise PSD (one-sided)
A    = 5 / sqrt(2*(sum(g.^2)/fa)/N0);         % sqrt(2E/N0) = 5 for the analog pulse
sigFixed = sigA*sqrt(rates(1)/fa);            % case A: equals case B at 1 kHz
snr1 = zeros(nR, 3); pred1 = zeros(nR, 3); n0rel = zeros(nR, 3);
for i = 1:nR
    fs = rates(i); L = round(fa/fs); nb = Na/L;
    gB = mean(reshape(g, L, nb), 1).';        % pulse seen by a bin average
    gC = g(ceil(L/2):L:end);                  % pulse at the kept samples
    st = zeros(nTr, 3);
    for r = 1:nTr
        n  = sigA*randn(Na, 1);
        xA = A*gB + sigFixed*randn(nb, 1);
        xB = A*gB + mean(reshape(n, L, nb), 1).';
        xC = A*gC + n(ceil(L/2):L:end);
        st(r, :) = [xA.'*gB, xB.'*gB, xC.'*gC];
    end
    snr1(i, :) = mean(st) ./ std(st);
    sig2 = [sigFixed^2, sigA^2/L, sigA^2];    % per-sample noise variance per case
    N0c  = 2*sig2/fs;                         % PSD the filter sees
    Ec   = A^2*[sum(gB.^2), sum(gB.^2), sum(gC.^2)]/fs;   % pulse energy
    pred1(i, :) = sqrt(2*Ec./N0c);
    n0rel(i, :) = N0c/N0;
end
z1 = (snr1 - pred1) ./ mcSig(pred1);
fprintf('Experiment 1: known pulse in white noise, matched-filter SNR (measured / sqrt(2E/N0))\n');
fprintf('   fs [Hz] |  A same sigma/sample    |  B fixed PSD, integrate  |  C fixed PSD, sample only\n');
fprintf('           |  SNR   pred  N0/N0phys  |  SNR   pred  N0/N0phys  |  SNR   pred  N0/N0phys\n');
for i = 1:nR
    fprintf('  %8.0f | %6.2f %6.2f %8.3g  | %5.2f %5.2f %8.3g   | %5.2f %5.2f %8.3g\n', rates(i), ...
        snr1(i,1), pred1(i,1), n0rel(i,1), snr1(i,2), pred1(i,2), n0rel(i,2), ...
        snr1(i,3), pred1(i,3), n0rel(i,3));
end
pass1 = all(abs(z1(:)) < 4);
fprintf('   measured vs sqrt(2E/N0): max |z| = %.1f sigma: %s\n', max(abs(z1(:))), passStr(pass1));

% one realization: matched-filter output versus lag at 1 kHz and 100 kHz
mf = struct();
n  = sigA*randn(Na, 1);                       % case B: the SAME analog noise at both rates
for fs = [1e3 1e5]
    L = round(fa/fs); nb = Na/L;
    gB = mean(reshape(g, L, nb), 1).';
    xB = A*gB + mean(reshape(n, L, nb), 1).';
    xA = A*gB + sigFixed*randn(nb, 1);
    cc = @(x, s) real(ifft(fft(x) .* conj(fft(gB)))) / (s*norm(gB));   % in output noise std
    lag = (0:nb-1).' / fs; lag(lag >= Tobs/2) = lag(lag >= Tobs/2) - Tobs;
    [lag, o] = sort(lag);
    cA = cc(xA, sigFixed); cB = cc(xB, sigA/sqrt(L));
    mf.(sprintf('f%d', fs)) = struct('lag', lag, 'A', cA(o), 'B', cB(o));
end

% ================= Experiment 2 =====================================================
B  = fa;                                      % complex samples at the bandwidth
gz = g - mean(g);
a  = 5 / sqrt(sum(gz.^2));                    % radiometer d = a*sqrt(B*int gz^2 dt) = 5
pw = 1 + a*g;
st2 = zeros(nTr, nR, 3);
for r = 1:nTr
    x = sqrt(pw/2) .* (randn(Na, 1) + 1i*randn(Na, 1));
    p = abs(x).^2;
    for i = 1:nR
        L = round(B/rates(i)); nb = Na/L;
        gB = mean(reshape(g, L, nb), 1).'; tB = gB - mean(gB);
        gC = g(ceil(L/2):L:end);           tC = gC - mean(gC);
        pI   = mean(reshape(p, L, nb), 1).';
        pII  = p(ceil(L/2):L:end);
        pIII = abs(mean(reshape(x, L, nb), 1).').^2;
        st2(r, i, :) = [pI.'*tB, pII.'*tC, pIII.'*tB];
    end
end
snr2  = squeeze(mean(st2, 1) ./ std(st2, 0, 1));
pred2 = zeros(nR, 3);
for i = 1:nR
    L = round(B/rates(i)); nb = Na/L;
    gB = mean(reshape(g, L, nb), 1).'; tB = gB - mean(gB);
    gC = g(ceil(L/2):L:end);           tC = gC - mean(gC);
    % i: bin variance 1/L; ii: one sample, variance 1; iii: |mean of L|^2, mean and std 1/L
    pred2(i, :) = [a*sqrt(L*sum(tB.^2)), a*sqrt(sum(tC.^2)), a*sqrt(sum(tB.^2))];
end
z2 = (snr2 - pred2) ./ mcSig(pred2);
fprintf('Experiment 2: noise-like pulse (power 1 + a g), square law, matched filter on the binned power\n');
fprintf('   f_out [Hz] | i square->average (pipeline) | ii square->keep one | iii average->square\n');
fprintf('              |   SNR   pred                 |   SNR   pred        |   SNR   pred\n');
for i = 1:nR
    fprintf('   %9.0f  |  %5.2f  %5.2f                 |  %5.2f  %5.2f       |  %5.2f  %5.2f\n', ...
        rates(i), snr2(i,1), pred2(i,1), snr2(i,2), pred2(i,2), snr2(i,3), pred2(i,3));
end
fprintf('   radiometer bound a*sqrt(B*int(g - <g>)^2 dt) = %.2f, a = %.4f (pulse peak / noise power)\n', ...
    a*sqrt(sum(gz.^2)), a);
pass2 = all(abs(z2(:)) < 4);
fprintf('   measured vs envelope-domain prediction: max |z| = %.1f sigma: %s\n', max(abs(z2(:))), passStr(pass2));

% ================= Figure ===========================================================
fig = figure('Name', 'expSamplingRate', 'Color', 'w', 'Position', [100 100 1150 800]);
tl  = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
ax = nexttile(tl);
loglog(ax, rates, snr1(:,1), 'o-', 'LineWidth', 1.4, 'DisplayName', 'A: same \sigma per sample at every rate'); hold(ax, 'on');
loglog(ax, rates, snr1(:,2), 's-', 'LineWidth', 1.4, 'DisplayName', 'B: fixed PSD, integrate (= pipeline)');
loglog(ax, rates, snr1(:,3), 'd-', 'LineWidth', 1.4, 'DisplayName', 'C: fixed PSD, keep 1 sample (aliasing)');
loglog(ax, rates, pred1, 'k:', 'HandleVisibility', 'off');
yline(ax, 5, 'k--', '\surd(2E/N_0), physical N_0', 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off');
grid(ax, 'on'); xlabel(ax, 'sample rate [Hz]'); ylabel(ax, 'matched-filter SNR');
title(ax, 'Exp. 1: known pulse in white noise (dotted: \surd(2E/N_0) with the N_0 the filter sees)');
ylim(ax, [0.1 300]); legend(ax, 'Location', 'northwest');
ax = nexttile(tl);
loglog(ax, rates, snr2(:,1), 's-', 'LineWidth', 1.4, 'DisplayName', 'i: square \rightarrow average (pipeline)'); hold(ax, 'on');
loglog(ax, rates, snr2(:,2), 'd-', 'LineWidth', 1.4, 'DisplayName', 'ii: square \rightarrow keep 1 per bin');
loglog(ax, rates, snr2(:,3), '^-', 'LineWidth', 1.4, 'DisplayName', 'iii: average voltages \rightarrow square');
loglog(ax, rates, pred2, 'k:', 'HandleVisibility', 'off');
yline(ax, 5, 'k--', 'radiometer: a\surd(B \int(g-\langle g \rangle)^2dt)', 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off');
grid(ax, 'on'); xlabel(ax, 'f_{out} [Hz]'); ylabel(ax, 'matched-filter SNR');
title(ax, 'Exp. 2: noise-like pulse, square law, binned power (B = 1 MHz)');
ylim(ax, [0.1 300]); legend(ax, 'Location', 'northwest');
for c = 'AB'
    ax = nexttile(tl);
    plot(ax, mf.f100000.lag*1e3, mf.f100000.(c), '-', 'Color', [0.85 0.33 0.1], 'DisplayName', '100 kHz'); hold(ax, 'on');
    plot(ax, mf.f1000.lag*1e3, mf.f1000.(c), 'o-', 'Color', [0 0.45 0.74], 'MarkerSize', 3, 'DisplayName', '1 kHz');
    xlim(ax, [-60 60]); grid(ax, 'on'); xlabel(ax, 'lag [ms]'); ylabel(ax, 'output / output noise std');
    if c == 'A'
        title(ax, 'MF output, case A (same \sigma per sample): 100 kHz looks far cleaner');
    else
        title(ax, 'MF output, case B (fixed PSD, integrate; same noise): same curve');
    end
    legend(ax, 'Location', 'northeast');
end

save(fullfile(root, 'data', 'mc', 'expSamplingRate.mat'), 'rates', 'snr1', 'pred1', 'n0rel', ...
    'snr2', 'pred2', 'mf', 'nTr', 'A', 'a', 'sigFixed');
if pass1 && pass2
    fprintf('expSamplingRate: ALL PASSED\n');
else
    error('expSamplingRate:failed', 'FAILED: measured SNR differs from the prediction (> 4 sigma).');
end
end


% =================================================================================
function s = passStr(p)
if p, s = 'PASS'; else, s = 'FAIL'; end
end
