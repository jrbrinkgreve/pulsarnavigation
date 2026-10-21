function clk = makeClock(tSpan, opts)
%MAKECLOCK  A receiver clock's error over time: offset, rate, ageing and power-law noise (block D3b; truth side).
%{
Every time stamp of the receiver comes from its oscillator, so the clock's error goes straight into
every TOA. This draws one realization of that error, eps(t) = clock reading - true time, for the
simulator's truth (trueResiduals 'Clock'). The observer never sees it; the navigation filter
estimates the clock's offset and rate as part of its state.

  clk = makeClock([t1 t2])                                   % a perfect clock
  clk = makeClock([t1 t2], 'Preset', "crystal", 'Rate', 1.2e-6, 'Seed', 3)
  e = clk.offset(t); y = clk.rate(t);                        % s and s/s at t (s since J2000, TT)

The model: eps(t) = Offset + Rate (t - t1) + Drift (t - t1)^2 / 2 + x(t), with x(t) the integral of
the oscillator's random fractional-frequency noise y(t). y is described, as usual for oscillators,
by the power-law coefficients h = [h0 h_1 h_2] of its one-sided spectrum
  S_y(f) = h0 + h_1 / f + h_2 / f^2,
white, flicker and random-walk frequency noise. Their Allan variances (Van Dierendonck, McGraw &
Brown 1984, "Relationship between Allan variances and Kalman filter parameters", eqs. 6-8):
  sigma_y^2(tau) = h0 / (2 tau) + 2 ln2 h_1 + (2 pi)^2 / 6 h_2 tau   (clk.sigmaY(tau)).
Drawing on a grid of 'Step' seconds over the span:
  white      x gets independent steps of variance h0 Step / 2 (exact).
  random walk the frequency and its integral as a two-state process with the exact covariance
             q2 [Step^3/3 Step^2/2; Step^2/2 Step], q2 = 2 pi^2 h_2 (exact).
  flicker    no finite-state model is exact (the paper's remark), so the step-averaged frequency
             is synthesized by FFT with the two-sided spectrum h_1 / (2 |f|), over twice the span
             (the FFT's wrap-around then never joins the two ends), and integrated. Its Allan
             variance is right within 3 % from 5 steps up to a twentieth of the span; at 2 steps
             it is 9 % high (the step averages get the full 1/f spectrum up to the grid's Nyquist
             frequency).
Between the grid nodes x is linear (rate piecewise constant). The noise this leaves out is a
white-FM Brownian bridge of sqrt(h0 Step / 8): 0.12 ns for the crystal at 60 s.

Presets ('H' overrides):
  "perfect"  h = 0 (default).
  "crystal"  the paper's example, a crystal oscillator stable to 3 parts in 1e11:
             h = [1.8e-21 6.492e-22 1.37e-24]. Its Figure 3 shows ~11 ns of time error 200 s after
             a start with known time and frequency; that number depends on how flicker is started,
             so here only the Allan deviation is comparable (white and random-walk noise alone
             reproduce the paper's growth exactly: 8.5 ns, testMakeClock).
A commercial SDR's reference (e.g. a TCXO specified to +-2 ppm) enters through 'Rate'; its noise
coefficients are not published and are best measured on the device.

Inputs:
  tSpan   [t1 t2] s since J2000.0 (TT): the span the clock must cover.

Name-value options:
  'Preset'  "perfect" | "crystal" (default "perfect").
  'H'       [h0 h_1 h_2] (1/Hz, 1, Hz) instead of the preset's.
  'Offset'  [s] clock error at t1 (default 0).
  'Rate'    [s/s] fractional frequency offset (default 0).
  'Drift'   [1/s] ageing: the rate's change per second (default 0).
  'Step'    [s] grid step (default 60).
  'Seed'    random seed (default 1; own RandStream, the global one is untouched).

Output clk: offset (@(t) eps [s]), rate (@(t) d eps/dt: Rate + Drift (t - t1) + the step-averaged
noise frequency), sigmaY (@(tau) the model's Allan deviation), tSpan, step, h, preset, deterministic
[Offset Rate Drift], seed, grid (struct t, x: the noise x at the nodes, x(1) = 0), made. Times
outside the grid: error 'makeClock:range'.
%}

arguments
    tSpan (1,2) double
    opts.Preset (1,1) string {mustBeMember(opts.Preset, ["perfect" "crystal"])} = "perfect"
    opts.H (1,3) double = [NaN NaN NaN]
    opts.Offset (1,1) double = 0
    opts.Rate (1,1) double = 0
    opts.Drift (1,1) double = 0
    opts.Step (1,1) double {mustBePositive} = 60
    opts.Seed (1,1) double {mustBeInteger, mustBeNonnegative} = 1
end

if tSpan(2) <= tSpan(1)
    error('makeClock:span', 'tSpan must be [start end] with end > start.');
end
switch opts.Preset
    case "perfect", h = [0 0 0];
    case "crystal", h = [1.8e-21 6.492e-22 1.37e-24];
end
if all(isfinite(opts.H))
    if any(opts.H < 0), error('makeClock:h', 'H = [h0 h_1 h_2] must be nonnegative.'); end
    h = opts.H;
end
t1 = tSpan(1); D = opts.Step;
K = ceil((tSpan(2) - t1) / D);                         % steps; nodes 0..K
tg = t1 + (0:K) * D;
rs = RandStream('mt19937ar', 'Seed', opts.Seed);
x = zeros(1, K + 1);

% white frequency noise: independent steps of x
if h(1) > 0
    x = x + [0, cumsum(sqrt(h(1) * D / 2) * randn(rs, 1, K))];
end

% flicker frequency noise: the step-averaged frequency by FFT synthesis over twice the span
if h(2) > 0
    N = 2^nextpow2(2 * K);
    f = (1:N/2) / (N * D);                             % positive frequencies [Hz]
    S2 = h(2) ./ (2 * f);                              % two-sided spectrum of y
    A = sqrt(N * S2 / D);                              % E|Y_m|^2 = N S2(f_m) / D
    Y = zeros(1, N);
    Y(2:N/2) = A(1:end-1) .* (randn(rs, 1, N/2 - 1) + 1i * randn(rs, 1, N/2 - 1)) / sqrt(2);
    Y(N/2 + 1) = A(end) * randn(rs);
    Y(N/2 + 2:N) = conj(Y(N/2:-1:2));
    yF = real(ifft(Y));
    x = x + [0, cumsum(D * yF(1:K))];
end

% random-walk frequency noise: frequency and its integral, exact two-state steps
if h(3) > 0
    q2 = 2 * pi^2 * h(3);
    L = chol(q2 * [D^3/3, D^2/2; D^2/2, D], 'lower');
    w = L * randn(rs, 2, K);
    xr = zeros(1, K + 1); yr = 0;
    for k = 1:K
        xr(k + 1) = xr(k) + D * yr + w(1, k);
        yr = yr + w(2, k);
    end
    x = x + xr;
end

off = opts.Offset; rate = opts.Rate; drift = opts.Drift;
yStep = diff(x) / D;                                   % noise frequency, averaged per step
clk = struct();
clk.offset = @(t) off + rate * (t - t1) + 0.5 * drift * (t - t1).^2 + gridInterp(tg, x, t);
clk.rate   = @(t) rate + drift * (t - t1) + gridStep(tg, yStep, t);
clk.sigmaY = @(tau) sqrt(h(1) ./ (2 * tau) + 2 * log(2) * h(2) + (2*pi)^2 / 6 * h(3) * tau);
clk.tSpan  = tSpan;
clk.step   = D;
clk.h      = h;
clk.preset = char(opts.Preset);
clk.deterministic = [off rate drift];
clk.seed   = opts.Seed;
clk.grid   = struct('t', tg, 'x', x);
clk.made   = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm'));
end


% =====================================================================================
function v = gridInterp(tg, x, t)
%GRIDINTERP  The noise x at times t, linear between the grid nodes.
checkRange(tg, t);
v = reshape(interp1(tg, x, t(:), 'linear'), size(t));
end

function v = gridStep(tg, yStep, t)
%GRIDSTEP  The step-averaged noise frequency at times t (the last node belongs to the last step).
checkRange(tg, t);
k = min(discretize(t, tg), numel(yStep));
v = reshape(yStep(k), size(t));
end

function checkRange(tg, t)
if any(t(:) < tg(1) | t(:) > tg(end))
    error('makeClock:range', 'Time outside the clock''s grid [%.10g, %.10g] s.', tg(1), tg(end));
end
end
