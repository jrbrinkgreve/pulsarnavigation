function res = trueResiduals(t, psr, trajTrue, trajAssumed, opts)
%TRUERESIDUALS  The pulse's true arrival against the observer's prediction (block D3a; truth side).
%{
The observer folds with a predictor made from the trajectory it ASSUMES (makePredictor), and it
time-stamps the data with its own clock. The pulses arrive where the TRUE trajectory puts them,
at true time. This gives, for pulses arriving at true times t, how late each one is against the
prediction - the timing residual that navigation measures - and where it therefore sits in the
observer's fold. Truth side, like the generator: only the simulator calls it (the ground-truth
rule); the observer sees only the shifted pulse.

  res = trueResiduals(t, psr, trajTrue, trajAssumed)
  res = trueResiduals(t, psr, makeTrajectory('fixed', 'Lat', 52, 'Lon', 5, 'OffsetENU', [100e3 0 0]), ...
                      makeTrajectory('fixed', 'Lat', 52, 'Lon', 5), 'Clock', @(t) 1e-6 + 1e-9 * (t - t(1)))

The physics. A pulse that passes the solar-system barycentre (SSB) at time S arrives at the
receiver at the true time t, with t + delta_T(t) = S (delayModel on the true trajectory). The
observer expects that pulse at t_A, with t_A + delta_A(t_A) = S (the assumed trajectory). Its
clock reads c = t + eps(t) when the pulse arrives (eps: clock reading - true time), and the
observer corrects that reading with its own estimate epsHat(c) (makePredictor 'Clock'). The
residual is
  r = (c - epsHat(c)) - t_A = x + eps(t) - epsHat(c),
  x = t - t_A = (delta_A - delta_T) / (1 + d delta_A/dt)
to first order (the factor is 1 + Doppler, |.| < 1e-4 on the Earth). x is solved exactly here,
by Newton steps on delta_A(t - x) - delta_T(t) - x = 0 in small numbers only (absolute times of
~8e8 s carry 1.2e-7 s of rounding; through the Doppler rate that is < 1e-10 s in x).

In the fold the pulse then sits at phase offset f r turns, f the observed spin frequency of the
assumed prediction (f = f_spin(S) (1 + d delta_A/dt), as makePredictor's freq): positive = late.
A true receiver closer to the pulsar than assumed gives r < 0 (the pulse arrives early).
Example: a site 100 km east of the assumed one, B0329+54 at 52 N: r is a daily sine of
193.33 us (2.7e-4 turns). A TOA is unambiguous only while |f r| < 0.5 turns: +-107,000 km for
B0329+54, +-863 km for J0437 (warning 'trueResiduals:wrap').

Inputs:
  t            [1 x N] true arrival times, s since J2000.0 (TT), as delayModel. For a sub-int,
               the observer's reference time is close enough: for a 100 km offset r changes by
               1.4e-8 s per second, for a clock by its rate error.
               Each t - x must lie inside both trajectories' time ranges.
  psr          timing model: f0 [Hz], F1 [Hz/s] (default 0), TRef [s since J2000, SSB], raDeg,
               decDeg, distance [m] (optional), as makePredictor.
  trajTrue     makeTrajectory struct: where the receiver really is.
  trajAssumed  makeTrajectory struct: where the observer thinks it is (its predictor's input).

Name-value options:
  'Clock'      [] (default: a perfect clock), a function handle @(t) -> clock error eps [s]
               (clock reading - true time, at true times t), or a struct with such a handle in
               the field 'offset' (makeClock, block D3b).
  'ClockEstimate'  [] (default: none, epsHat = 0) or the observer's estimate of its clock error
               as a function of the clock reading, @(c) epsHat [s] or a struct with field
               offset: the same estimate its predictor uses (makePredictor 'Clock'). Only the
               part of the clock error it has not estimated shifts the pulse.
  'DelayArgs'  name-value pairs for delayModel, used for both trajectories (default {}).

Output res: t; residual [1 x N] (s: observed - predicted arrival, the clock included);
geometric (x), clock (eps) and clockEstimate (epsHat(c)) parts [1 x N] s; fObs [1 x N] Hz; phaseOffset [1 x N] turns (where
the pulse sits in the observer's fold); deltaTrue (delta_T(t)) and deltaAssumed (delta_A(t - x))
[1 x N] s; solveErr (largest |delta_A(t - x) - delta_T(t) - x| left, s); trajTypes; made.
%}

arguments
    t (1,:) double
    psr struct
    trajTrue struct
    trajAssumed struct
    opts.Clock = []
    opts.ClockEstimate = []
    opts.DelayArgs cell = {}
end

c = 299792458;
F1 = 0; if isfield(psr, 'F1') && ~isempty(psr.F1), F1 = psr.F1; end
dArgs = opts.DelayArgs;

% the clock: its true error at t, and the observer's estimate at the reading c = t + eps(t)
epsC = clockValues(opts.Clock, t, 'Clock');
epsHat = clockValues(opts.ClockEstimate, t + epsC, 'ClockEstimate');

% the true arrival: delta_T at t
dT = delayModel(t, trajTrue, psr, dArgs{:});

% the predicted arrival t_A = t - x: Newton on g(x) = delta_A(t - x) - delta_T(t) - x
dA = delayModel(t, trajAssumed, psr, dArgs{:});
x = dA.total - dT.total;                              % first guess (exact up to the Doppler factor)
for it = 1:3
    dA = delayModel(t - x, trajAssumed, psr, dArgs{:});
    rate = (dA.nHat.' * dA.vObs) / c;                  % d delta_A / dt (Roemer rate; the rest < 5e-10)
    g = dA.total - dT.total - x;
    x = x + g ./ (1 + rate);
end
dA = delayModel(t - x, trajAssumed, psr, dArgs{:});
solveErr = max(abs(dA.total - dT.total - x));
rate = (dA.nHat.' * dA.vObs) / c;

% where the pulse sits in the observer's fold
fSpin = psr.f0 + F1 * ((t - psr.TRef) + dT.total);   % spin frequency at the SSB time S
fObs = fSpin .* (1 + rate);
r = x + epsC - epsHat;
ph = fObs .* r;
if any(abs(ph) >= 0.5)
    warning('trueResiduals:wrap', ['The pulse is %.3g turns from its predicted phase: beyond half a turn ' ...
        'its TOA is ambiguous by whole turns.'], max(abs(ph)));
end

res = struct();
res.t            = t;
res.residual     = r;
res.geometric    = x;
res.clock        = epsC;
res.clockEstimate = epsHat;
res.fObs         = fObs;
res.phaseOffset  = ph;
res.deltaTrue    = dT.total;
res.deltaAssumed = dA.total;
res.solveErr     = solveErr;
res.trajTypes    = {trajTrue.type, trajAssumed.type};
res.made         = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm'));
end


% =====================================================================================
function v = clockValues(clk, t, name)
%CLOCKVALUES  A clock (or clock estimate) evaluated at t: [] -> 0, a handle, or a struct's offset.
if isempty(clk)
    v = zeros(size(t));
elseif isa(clk, 'function_handle')
    v = clk(t);
elseif isstruct(clk) && isfield(clk, 'offset') && isa(clk.offset, 'function_handle')
    v = clk.offset(t);
else
    error('trueResiduals:clock', '%s must be [], a function handle or a struct with a field offset = @(t).', name);
end
v = reshape(v, size(t));
end
