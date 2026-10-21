function pred = makePredictor(psr, traj, tSpan, opts)
%MAKEPREDICTOR  Predicted pulse phase at the receiver: Chebyshev segments from the timing model and
%delayModel (block D2a).
%{
The pulsar's rotation is a clock at the solar-system barycentre (SSB):
  phi(t_SSB) = Phi0 + f0 (t_SSB - TRef) + F1 (t_SSB - TRef)^2 / 2,
and a pulse received at time t on the trajectory arrives at the SSB at t_SSB = t + delta(t)
(delayModel: Roemer, Einstein, Shapiro). So the phase the receiver should see is
phi(t + delta(t)): Earth's orbit makes its rate change by up to 1e-4 (Doppler), the Earth's
rotation by 1e-6, a ship by 3e-8. Evaluating delayModel for every time bin of a fold would put all
that physics into the stream; instead, as at every radio observatory (TEMPO's polycos), the phase
is precomputed once per observation as short polynomials, and the fold only evaluates them. The
trajectory (a fixed site, a ship track, a probe, or a Kalman filter's state: makeTrajectory) is
the only input about the receiver: the entry point for a motion predictor (decision 2 of block
D, 21 Oct 2026).

  pred = makePredictor(psr, traj, [t1 t2])
  [phi, f] = deal(pred.phase(tau), pred.freq(tau))   % tau: local time, s after pred.tOrigin

Inputs:
  psr    timing model: f0 [Hz], F1 [Hz/s] (default 0), TRef [s since J2000, at the SSB], Phi0
         [turns] (default 0), raDeg, decDeg, distance [m] (optional; delayModel).
  traj   makeTrajectory struct (the receiver).
  tSpan  [t1 t2] s since J2000 (TT): the span to predict; t1 is the local time origin. With
         'Clock' these are clock readings (the data's own time stamps).

Name-value options:
  'SegLength'  [s] segment length before refinement (default 3600).
  'Degree'     Chebyshev degree per segment (default 12).
  'Tol'        [s] largest timing error allowed (default 1e-9): each segment is checked against
               the direct computation halfway between its nodes and split in two until the
               phase error is below Tol x f0 (a ship's turn on a linear track is a kink in the
               delay; the split isolates it). Shortest segment 'MinSeg' (default 1 s).
  'DelayArgs'  name-value pairs for delayModel (default {}), e.g. {'Terms', "roemer"}.
  'Clock'      [] (default: the time stamps are true time) or the observer's ESTIMATE of its clock
               error, as a function of the clock reading c: @(c) epsHat [s], or a struct with such
               a handle in the field 'offset'. The receiver's samples are timed by its own
               oscillator, so the data's times are clock readings; the predictor turns each into
               its best estimate of true time, t = c - epsHat(c), before the delays. An SDR's
               +-2 ppm rate would otherwise smear the fold like an uncorrected Doppler shift
               (240 us in a 120 s sub-int). Local time tau and freq are then in clock seconds.
               Block D3c, option (a): the observer folds with its own clock estimate.

Precision: absolute times (~8e8 s) and absolute phases (up to ~1e11 turns) never enter the
polynomials. The phase is kept relative to N0 = the integer turn at t1:
  phi_rel(tau) = frac0 + f0 (tau + delta(t) - delta0) + F1 (u - u0)(u + u0) / 2,
u = t + delta(t) - TRef; only frac0 = phi(t1) - N0 carries the double rounding of the big
product (~1e-6 turns), and it is one constant for the whole span (no steps between segments).
With 'Clock', tau is replaced by tau - (epsHat(t1 + tau) - epsHat(t1)), also a small difference.

Output pred: tOrigin (= t1), tSpanLocal [0, t2 - t1], N0, frac0, seg [nSeg x 2] (local start and
end), coef [nSeg x (Degree+1)] (Chebyshev coefficients of phi_rel), phase (@(tau) phi_rel, turns
since N0), freq (@(tau) d phi / d tau [Hz]: the observed spin frequency), fRange [min max],
maxErr [s] (largest timing error found by the checks), psr, trajType, delayArgs, clock (the
estimate's handle, or []), tOriginTrue (the true time estimated at t1), made.
Evaluation outside tSpanLocal: error 'makePredictor:range'.
%}

arguments
    psr struct
    traj struct
    tSpan (1,2) double
    opts.SegLength (1,1) double {mustBePositive} = 3600
    opts.Degree    (1,1) double {mustBeInteger, mustBePositive} = 12
    opts.Tol       (1,1) double {mustBePositive} = 1e-9
    opts.MinSeg    (1,1) double {mustBePositive} = 1
    opts.DelayArgs cell = {}
    opts.Clock = []
end

if tSpan(2) <= tSpan(1)
    error('makePredictor:span', 'tSpan must be [start end] with end > start.');
end
F1 = 0; if isfield(psr, 'F1') && ~isempty(psr.F1), F1 = psr.F1; end
Phi0 = 0; if isfield(psr, 'Phi0') && ~isempty(psr.Phi0), Phi0 = psr.Phi0; end
f0 = psr.f0; t0 = tSpan(1); L = tSpan(2) - t0; K = opts.Degree;
epsHat = clockHandle(opts.Clock);

% reference at t0: the only place where absolute time and phase meet
e0 = epsHat(t0);
tt0 = t0 - e0;                                         % the true time estimated at the first reading
d0 = delayModel(tt0, traj, psr, opts.DelayArgs{:});
u0 = (tt0 - psr.TRef) + d0.total;
phAbs0 = Phi0 + f0 * u0 + 0.5 * F1 * u0^2;
N0 = round(phAbs0); frac0 = phAbs0 - N0;
phiRel = @(tau) phiRelative(tau, t0, tt0, e0, epsHat, traj, psr, opts.DelayArgs, d0.total, u0, f0, F1, frac0);

% segments: Chebyshev fits, split until the checks pass
edges = unique([0:opts.SegLength:L, L]);
todo = [edges(1:end-1).', edges(2:end).'];
seg = zeros(0, 2); coef = zeros(0, K + 1); maxErr = 0;
xn = cos(pi * ((0:K) + 0.5) / (K + 1));               % Chebyshev nodes on [-1, 1]
xc = cos(pi * (0:K) / K);                             % checks: the extrema, between the nodes
Tn = cos((0:K).' * acos(xn));                         % T_k(x_j), [K+1 x K+1]
while ~isempty(todo)
    a = todo(1, 1); b = todo(1, 2); todo(1, :) = [];
    m = (a + b) / 2; h = (b - a) / 2;
    c = (2 / (K + 1)) * (Tn * phiRel(m + h * xn).').';
    c(1) = c(1) / 2;
    err = max(abs(chebEval(c, xc) - phiRel(m + h * xc))) / f0;
    if err > opts.Tol && (b - a) / 2 >= opts.MinSeg
        todo = [todo; a, m; m, b]; %#ok<AGROW>
        continue
    end
    seg = [seg; a, b]; coef = [coef; c]; maxErr = max(maxErr, err); %#ok<AGROW>
end
[seg, ord] = sortrows(seg); coef = coef(ord, :);
dcoef = zeros(size(coef));
for i = 1:size(coef, 1)
    dcoef(i, :) = chebDeriv(coef(i, :)) / ((seg(i, 2) - seg(i, 1)) / 2);
end

pred = struct();
pred.tOrigin    = t0;
pred.tSpanLocal = [0 L];
pred.N0         = N0;
pred.frac0      = frac0;
pred.seg        = seg;
pred.coef       = coef;
pred.phase      = @(tau) evalSeg(seg, coef, tau);
pred.freq       = @(tau) evalSeg(seg, dcoef, tau);
fs = pred.freq(linspace(0, L, max(3, ceil(L / 600) + 1)));
pred.fRange     = [min(fs), max(fs)];
pred.maxErr     = maxErr;
pred.psr        = psr;
pred.trajType   = traj.type;
pred.delayArgs  = opts.DelayArgs;
pred.clock      = [];
if ~isempty(opts.Clock), pred.clock = epsHat; end
pred.tOriginTrue = tt0;
pred.made       = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm'));
end


% =====================================================================================
function ph = phiRelative(tau, t0, tt0, e0, epsHat, traj, psr, dArgs, delta0, u0, f0, F1, frac0)
%PHIRELATIVE  Phase since N0 at local (clock) times tau, from small differences only.
dt = tau - (epsHat(t0 + tau) - e0);                    % true time since tt0 (= tau for a perfect clock)
d = delayModel(tt0 + dt, traj, psr, dArgs{:});
du = dt + (d.total - delta0);                          % u - u0
ph = frac0 + f0 * du + 0.5 * F1 * du .* (du + 2 * u0);
end

function h = clockHandle(clk)
%CLOCKHANDLE  The observer's clock estimate as a handle @(c) (zero for a perfect clock).
if isempty(clk)
    h = @(c) zeros(size(c));
elseif isa(clk, 'function_handle')
    h = clk;
elseif isstruct(clk) && isfield(clk, 'offset') && isa(clk.offset, 'function_handle')
    h = clk.offset;
else
    error('makePredictor:clock', 'Clock must be [], a function handle @(c) or a struct with a field offset = @(c).');
end
end

function y = chebEval(c, x)
%CHEBEVAL  sum_k c(k+1) T_k(x) by Clenshaw's recurrence (c(1) already halved).
b1 = zeros(size(x)); b2 = b1;
for k = numel(c):-1:2
    b0 = 2 * x .* b1 - b2 + c(k);
    b2 = b1; b1 = b0;
end
y = x .* b1 - b2 + c(1);
end

function d = chebDeriv(c)
%CHEBDERIV  Coefficients (same convention) of the x-derivative of a Chebyshev series.
n = numel(c) - 1;
cc = c; cc(1) = 2 * c(1);                              % full convention for the recurrence
dd = zeros(1, n + 2);
for k = n:-1:1
    dd(k) = dd(k + 2) + 2 * k * cc(k + 1);
end
d = dd(1:n + 1); d(1) = d(1) / 2;
end

function y = evalSeg(seg, coef, tau)
%EVALSEG  Evaluate the segment series at local times tau (the segment that contains each).
if any(tau < seg(1, 1) | tau > seg(end, 2))
    error('makePredictor:range', 'Time outside the predictor [%g, %g] s (local).', seg(1, 1), seg(end, 2));
end
k = discretize(tau, [seg(:, 1); seg(end, 2)]);
y = zeros(size(tau));
for i = unique(k(:)).'
    sel = k == i;
    m = (seg(i, 1) + seg(i, 2)) / 2; h = (seg(i, 2) - seg(i, 1)) / 2;
    y(sel) = chebEval(coef(i, :), (tau(sel) - m) / h);
end
end
