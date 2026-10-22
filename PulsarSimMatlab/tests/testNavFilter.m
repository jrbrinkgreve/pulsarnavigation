function testNavFilter()
%TESTNAVFILTER  Unit tests for navFilter (block D4d: the Kalman filter and smoother for a moving receiver).
%{
Run from the PulsarSimMatlab folder: run('tests/testNavFilter.m'), or runAllTests("testNavFilter").
~1 min; no files. The truth comes from trueResiduals (delayModel on the true and the assumed track)
and makeClock, independent of the filter's own model.

  1. A fixed site is the batch fit: with no velocity error (prior and noise 0) and no height
     constraint, the smoother's position equals navFit's (B0329+54, 300 TOAs over 30 days, 100 km
     east): east and radial within 1e-3 of their errors, the errors within 1e-3 relative -
     "self" without a clock, "self" with the crystal clock (also the clock estimate), and
     "calibrated" with the clock (also the polar component).
  2. A ship: a dead-reckoning track along 50 N from 5 W at 2.1 m/s for 14 days; the true track
     drifts from it by an unknown current (0.3 m/s per axis, drawn) plus a random walk (as the
     filter's model), added in the local tangent plane of the DR track (the filter's height model;
     on the real surface the drop d^2/2R is left to the closed loop, which keeps d small);
     B0329+54 (10 TOAs a day) and B0950+08 (0.75 a day), 734 us each; the crystal clock with an
     unknown rate. 40 realizations, the smoothed east and north errors (estimate - truth) / stated
     error at four fixed times each (days 2, 5.5, 9 and 14: a smoothed track's errors are strongly
     correlated in time, so the pulls of all TOAs are far from independent): std 1 (+- 0.2), mean 0
     (+- 4 / sqrt(40)); shown: the dead-reckoning error after 14 days against the filter's.
  3. The same with 2 % outliers at 5-10 sigma (20 realizations, pulls at the same four times): the
     robust filter keeps pulls of std 1 (+- 0.3) and flags at least 70 % of all outliers (an
     outlier early on, against a still vague prediction, cannot be recognised); without 'Robust'
     the pulls are wider.
  4. 'Predict': the covariances equal those of the run with data ('Robust' false), 1e-9 relative.
  5. A trajectory in the 'SSB' frame, missing residuals and a bad clock are refused
     ('navFilter:traj', 'navFilter:toa', 'navFilter:clock').
  6. 'ClockPrior': the ships of check 2 with a calibrated start whose port clock is off by a drawn
     offset (1 ms) and rate (1e-9) - the observer's residuals against its port clock estimate. With
     'ClockPrior' [1e-3 1e-9] the pulls of east, north and the clock offset at the four times have
     std 1 (+- 0.2) and mean 0 (+- 4 / sqrt(40)); with the prior taken as exact ([0 0]) the rms of
     the north or clock pulls exceeds 2 (the option is needed). 'ClockPrior' without 'Clock' is
     refused ('navFilter:clock'). Its own random stream: checks 2 and 3 keep their numbers.
Errors at the end if any check fails. (22 Oct 2026, block D4d; check 6 for D4c-1.)
%}

root = fullfile(fileparts(mfilename('fullpath')), '..');
addpath(fullfile(root, 'functions'));
fails = {};
rs = RandStream('mt19937ar', 'Seed', 52);
c = 299792458; t1 = 9555.5 * 86400; sig = 734e-6;
b0329 = struct('raDeg', 53.2474, 'decDeg', 54.5787);
b0950 = struct('raDeg', 15 * (9 + 53/60 + 9.3/3600), 'decDeg', 7 + 55/60 + 36/3600);
psrT = @(p) struct('f0', 1, 'F1', 0, 'TRef', t1, 'raDeg', p.raDeg, 'decDeg', p.decDeg);
hCr = makeClock([0 1], 'Preset', "crystal").h;

% 1. A fixed site = the batch fit
N = 300; t = t1 + sort(30 * 86400 * rand(rs, 1, N));
siteA = makeTrajectory('fixed', 'Lat', 52, 'Lon', 5);
siteE = makeTrajectory('fixed', 'Lat', 52, 'Lon', 5, 'OffsetENU', [100e3 0 0]);
rE = trueResiduals(t, psrT(b0329), siteE, siteA).residual;
clk = makeClock([t(1), t(end)], 'Preset', "crystal", 'Rate', 1e-9, 'Seed', 7);
rN = rE + sig * randn(rs, 1, N);
fixedArgs = {'VelocityPrior', 0, 'VelocityNoise', 0, 'Height', false, 'Robust', false};
eq = [-sind(5) cosd(5) 0; cosd(5) sind(5) 0];         % east, radial in terms of X, Y, Z
dev = zeros(1, 3); ok1 = true;
cases = {"self", [], rN; "self", hCr, rN + clk.offset(t); "calibrated", hCr, rN + (clk.offset(t) - clk.offset(t(1)) - 1e-9 * (t - t(1)))};
for ic = 1:3
    toa = struct('t', t, 'r', cases{ic, 3}, 'sigma', sig * ones(1, N));
    fB = navFit(toa, b0329, siteA, 'Start', cases{ic, 1}, 'Clock', cases{ic, 2}, 'ClockStart', t(1), 'Robust', false);
    fK = navFilter(toa, b0329, siteA, 'Start', cases{ic, 1}, 'Clock', cases{ic, 2}, fixedArgs{:});
    Ms = siteA.enu;                                     % local east, north, up -> Earth-fixed
    pK = Ms * fK.x(1:3, 1); PK = Ms * fK.P(1:3, 1:3, 1) * Ms';
    estK = eq * pK; errK = sqrt(diag(eq * PK * eq'));
    estB = [fB.east; fB.radial]; errB = [fB.eastErr; fB.radialErr];
    d = max([abs(estK - estB) ./ errB; abs(errK ./ errB - 1)]);
    if ~isempty(cases{ic, 2})                           % the clock estimate too
        d = max(d, max(abs(fK.clock - fB.clock) ./ fB.clockErr));
    end
    if cases{ic, 1} == "calibrated"                     % and the polar component
        d = max([d, abs(pK(3) - fB.polar) / fB.polarErr, abs(sqrt(PK(3, 3)) / fB.polarErr - 1)]);
    end
    dev(ic) = d; ok1 = ok1 && d < 1e-3;
end
fprintf(['1. fixed site = navFit: largest difference %.1e ("self"), %.1e ("self" + clock), %.1e ("calibrated" + clock) ' ...
         '(of the errors): %s\n'], dev, passStr(ok1));
if ~ok1, fails{end+1} = 'batch equivalence'; end

% 2. A ship with dead-reckoning drift, two pulsars, the crystal clock
L = 14 * 86400; tg = t1 + (0:3600:L); vShip = 2.1;
lonDR = -5 - rad2deg(vShip * (tg - t1) / (6371e3 * cosd(50)));
trajDR = makeTrajectory('track', 'Time', tg, 'Lat', 50 * ones(size(tg)), 'Lon', lonDR, 'Interp', 'pchip');
nR = 40;
[pE, pN, dr, kfErr, pER, pNR, hitR, pEW, pC, pX] = deal([]);
rs6 = RandStream('mt19937ar', 'Seed', 53);             % check 6's draws, apart from checks 2 and 3
ws = warning('off', 'trueResiduals:wrap');             % seconds of SDR clock drift: no fold here, only residuals
for m = 1:nR
    [trajT, offT] = shipTruth(trajDR, tg, rs);
    [toa, ~] = twoPulsars(t1, L, rs);
    clk = makeClock([t1, t1 + L + 86400], 'Preset', "crystal", 'Rate', 2e-6 * (2 * rand(rs) - 1), 'Seed', 100 + m);
    r = zeros(1, numel(toa.t));
    for k = 1:2
        s = toa.k == k; p = {b0329, b0950};
        r(s) = trueResiduals(toa.t(s), psrT(p{k}), trajT, trajDR, 'Clock', clk).residual;
    end
    toa.r = r + sig * randn(rs, 1, numel(r));
    f = navFilter(toa, [b0329, b0950], trajDR, 'Clock', hCr, 'PositionPrior', 1e3);
    truthENU = offT(f.t);                               % true - DR in the local frame, [3 x N]
    at = arrayfun(@(dd) find(f.t >= t1 + dd * 86400, 1), [2 5.5 9]); at(end+1) = numel(f.t);   % four times
    pE = [pE, (f.enu(1, at) - truthENU(1, at)) ./ f.enuErr(1, at)]; %#ok<AGROW>
    pN = [pN, (f.enu(2, at) - truthENU(2, at)) ./ f.enuErr(2, at)]; %#ok<AGROW>
    dr = [dr, norm(truthENU(1:2, end))]; kfErr = [kfErr, norm(f.enu(1:2, end) - truthENU(1:2, end))]; %#ok<AGROW>
    if m <= 20                                          % 3. with outliers
        bad = rand(rs, 1, numel(r)) < 0.02;
        toaO = toa; toaO.r(bad) = toaO.r(bad) + sig * (5 + 5 * rand(rs, 1, nnz(bad))) .* sign(randn(rs, 1, nnz(bad)));
        fR = navFilter(toaO, [b0329, b0950], trajDR, 'Clock', hCr, 'PositionPrior', 1e3);
        fW = navFilter(toaO, [b0329, b0950], trajDR, 'Clock', hCr, 'PositionPrior', 1e3, 'Robust', false);
        pER = [pER, (fR.enu(1, at) - truthENU(1, at)) ./ fR.enuErr(1, at)]; %#ok<AGROW>
        pNR = [pNR, (fR.enu(2, at) - truthENU(2, at)) ./ fR.enuErr(2, at)]; %#ok<AGROW>
        pEW = [pEW, (fW.enu(1, at) - truthENU(1, at)) ./ fW.enuErr(1, at)]; %#ok<AGROW>
        badSorted = bad(f.order);
        hitR = [hitR, [nnz(fR.outlier(badSorted)); nnz(badSorted)]]; %#ok<AGROW>
    end
    % 6. a calibrated start, the port clock off by a drawn offset and rate
    t0 = toa.t(1); dB = 1e-3 * randn(rs6); dD = 1e-9 * randn(rs6);
    cRead = toa.t + clk.offset(toa.t);                  % the clock's readings
    eHat = clk.offset(t0) + dB + (clk.rate(t0) + dD) * (cRead - t0);   % the port estimate there
    toaC = toa; toaC.r = toa.r - eHat;
    bTrue = clk.offset(toa.t) - eHat;                   % the clock error left in the residuals
    calArgs = {[b0329, b0950], trajDR, 'Clock', hCr, 'PositionPrior', 1e3, 'Start', "calibrated"};
    fC = navFilter(toaC, calArgs{:}, 'ClockPrior', [1e-3 1e-9]);
    fX = navFilter(toaC, calArgs{:});                   % the prior taken as exact
    bT = bTrue(fC.order);
    pC = [pC, [(fC.enu(1:2, at) - truthENU(1:2, at)) ./ fC.enuErr(1:2, at); (fC.clock(at) - bT(at)) ./ fC.clockErr(at)]]; %#ok<AGROW>
    pX = [pX, [(fX.enu(1:2, at) - truthENU(1:2, at)) ./ fX.enuErr(1:2, at); (fX.clock(at) - bT(at)) ./ fX.clockErr(at)]]; %#ok<AGROW>
end
warning(ws);
ok2 = abs(std(pE) - 1) < 0.2 && abs(std(pN) - 1) < 0.2 && abs(mean(pE)) < 4 / sqrt(nR) && abs(mean(pN)) < 4 / sqrt(nR);
fprintf(['2. ship, 14 days, %d realizations: pulls east mean %+.2f std %.3f, north mean %+.2f std %.3f; after 14 days ' ...
         'the dead-reckoning error %.0f km (median), the filter''s %.0f km: %s\n'], nR, mean(pE), std(pE), mean(pN), std(pN), ...
    median(dr) / 1e3, median(kfErr) / 1e3, passStr(ok2));
if ~ok2, fails{end+1} = 'ship'; end
fHit = sum(hitR(1, :)) / sum(hitR(2, :));
ok3 = abs(std(pER) - 1) < 0.3 && abs(std(pNR) - 1) < 0.3 && fHit >= 0.7 && std(pEW) > std(pER);
fprintf('3. 2 %% outliers: robust pulls east %.3f, north %.3f, outliers flagged %.0f %% (%d of %d); without ''Robust'' east %.3f: %s\n', ...
    std(pER), std(pNR), 100 * fHit, sum(hitR(1, :)), sum(hitR(2, :)), std(pEW), passStr(ok3));
if ~ok3, fails{end+1} = 'outliers'; end

% 4. Predict
ws = warning('off', 'trueResiduals:wrap');
fD = navFilter(toa, [b0329, b0950], trajDR, 'Clock', hCr, 'PositionPrior', 1e3, 'Robust', false);
fP = navFilter(rmfield(toa, 'r'), [b0329, b0950], trajDR, 'Clock', hCr, 'PositionPrior', 1e3, 'Predict', true);
e4 = max(abs(fP.P(:) - fD.P(:))) / max(abs(fD.P(:)));
warning(ws);
ok4 = e4 < 1e-9;
fprintf('4. Predict = the covariances of the run with data: %.1e relative: %s\n', e4, passStr(ok4));
if ~ok4, fails{end+1} = 'Predict'; end

% 5. Input checks
probe = makeTrajectory('custom', 'Frame', 'SSB', 'PosFcn', @(tt) repmat([1e11; 0; 0], 1, numel(tt)));
ok5 = strcmp(errId(@() navFilter(toa, [b0329, b0950], probe)), 'navFilter:traj') && ...
      strcmp(errId(@() navFilter(rmfield(toa, 'r'), [b0329, b0950], trajDR)), 'navFilter:toa') && ...
      strcmp(errId(@() navFilter(toa, [b0329, b0950], trajDR, 'Clock', [1 -1 0])), 'navFilter:clock');
fprintf('5. SSB trajectory, missing residuals and a bad clock refused %d: %s\n', ok5, passStr(ok5));
if ~ok5, fails{end+1} = 'input checks'; end

% 6. ClockPrior (the realizations ran with check 2)
mC = mean(pC, 2); sC = std(pC, 0, 2); rX = rms(pX, 2);
ok6 = all(abs(sC - 1) < 0.2) && all(abs(mC) < 4 / sqrt(nR)) && max(rX(2:3)) > 2 && ...
      strcmp(errId(@() navFilter(toa, [b0329, b0950], trajDR, 'Start', "calibrated", 'ClockPrior', [1e-3 1e-9])), 'navFilter:clock');
fprintf(['6. calibrated start, the port clock off by 1 ms and 1e-9: pulls east %+.2f / %.3f, north %+.2f / %.3f, clock %+.2f / %.3f ' ...
         '(mean / std); the prior taken as exact: rms east %.1f, north %.1f, clock %.1f: %s\n'], [mC, sC].', rX, passStr(ok6));
if ~ok6, fails{end+1} = 'ClockPrior'; end

if isempty(fails)
    fprintf('testNavFilter: ALL PASSED\n');
else
    error('testNavFilter:failed', 'FAILED: %s', strjoin(fails, ', '));
end
end


% =================================================================================
function [trajT, offT] = shipTruth(trajDR, tg, rs)
% The true track: the dead-reckoning track plus a drift in the local horizontal frame - an unknown current
% (0.3 m/s per axis) plus a velocity random walk (1.2e-7 m^2/s^3, as the filter's default), integrated
% on the hourly grid. offT(t) gives true - DR in the local east-north-up frame of the DR track.
dt = diff(tg(1:2)); K = numel(tg);
v = 0.3 * randn(rs, 2, 1) + [zeros(2, 1), cumsum(sqrt(1.2e-7 * dt) * randn(rs, 2, K - 1), 2)];
d = [zeros(2, 1), cumsum(dt * (v(:, 1:end-1) + v(:, 2:end)) / 2, 2)];   % east, north [m]
rDR = trajDR.pos(tg);
[E, Nn] = frameAt(rDR);
rT = rDR + E .* d(1, :) + Nn .* d(2, :);
trajT = makeTrajectory('custom', 'Frame', 'ECEF', 'PosFcn', @(t) interp1(tg, rT.', t, 'pchip').');
offT = @(t) [interp1(tg, d(1, :), t, 'pchip'); interp1(tg, d(2, :), t, 'pchip'); zeros(size(t))];
end

function [E, N] = frameAt(r)
% Local east and north unit vectors at Earth-fixed positions r [3 x M], WGS84 normal (Bowring's latitude):
% the filter's frame, so that the drift is exactly horizontal for its height model.
a = 6378137; f = 1 / 298.257223563; e2 = f * (2 - f); b = a * (1 - f); ep2 = (a^2 - b^2) / b^2;
p = hypot(r(1, :), r(2, :)); lon = atan2(r(2, :), r(1, :));
th = atan2(r(3, :) * a, p * b);
lat = atan2(r(3, :) + ep2 * b * sin(th).^3, p - e2 * a * cos(th).^3);
E = [-sin(lon); cos(lon); zeros(size(lon))];
N = [-sin(lat) .* cos(lon); -sin(lat) .* sin(lon); cos(lat)];
end

function [toa, n] = twoPulsars(t1, L, rs)
% TOA times: B0329+54 10 a day, B0950+08 0.75 a day, at random times; sigma 734 us.
n = [round(10 * L / 86400), round(0.75 * L / 86400)];
t = [t1 + L * rand(rs, 1, n(1)), t1 + L * rand(rs, 1, n(2))];
k = [ones(1, n(1)), 2 * ones(1, n(2))];
[t, o] = sort(t);
toa = struct('t', t, 'sigma', 734e-6 * ones(size(t)), 'k', k(o));
end

function id = errId(f)
id = '';
try
    f();
catch err
    id = err.identifier;
end
end

function s = passStr(p)
if p, s = 'PASS'; else, s = 'FAIL'; end
end
