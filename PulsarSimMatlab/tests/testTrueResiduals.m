function testTrueResiduals()
%TESTTRUERESIDUALS  Unit tests for trueResiduals (block D3a: the true arrival against the prediction).
%{
Run from the PulsarSimMatlab folder: run('tests/testTrueResiduals.m'), or
runAllTests("testTrueResiduals"). ~2 s; small files in its own temporary folder. The truth for
each check is computed independently of the function's Newton solution: closed forms, or
delayModel on the true trajectory directly.

  1. The same trajectory as truth and assumption, a perfect clock: residual exactly 0.
  2. A site 100 km east of the assumed one (52 N 5 E, B0329+54, one day, Roemer only): the
     residual times (1 + d delta_A/dt) equals -(R dr . n)/c, the closed form with the Earth
     rotation angle, a daily sine of 193.33 us (1e-10 s: the analytic model's own rounding - the
     Earth's mean longitude, ~9500 deg in 2026, carries 1.8e-12 deg in double, 5 mm of orbit,
     1.6e-11 s, at t and t - x independently; and t - x itself is rounded to 6e-8 s, times the
     Doppler 1e-4); negative where the true site is closer to the pulsar. The residual alone
     is that sine divided by 1 + Doppler (-8.3e-5 on 1 March). With the full model the residual
     differs by the Einstein site term only (< 100 ns).
  3. A probe at 300 km/s along the line of sight (Doppler 1e-3, 'SSB' frame), the true one 1000 km
     closer to the pulsar (Roemer only): the residual is exactly -(1000 km / c) / (1 + v/c)
     (1e-10 s); dropping the Doppler factor would be off by 3.3 us.
  4. Clock: a handle and a struct with field offset give the same result; the clock part equals
     the handle's values, the geometric part is unchanged, residual = geometric + clock and
     phaseOffset = fObs x residual.
  5. End to end with the fold. Synthetic noise-free data (50 Hz, 20 s, 20 us bins) whose pulses
     arrive where the TRUE probe of check 3 puts them (delayModel on the true trajectory, all
     terms), folded with the predictor of the ASSUMED probe (foldProfile 'Predictor'):
     a. at 200 times the predicted phase minus the data's phase equals phaseOffset (1e-7 turns,
        the predictor's 1 ns tolerance at 50 Hz);
     b. the folded pulse peak (centroid) sits at phaseOffset, -0.167 turns: early (1e-4 turns).
  6. A bad clock is refused ('trueResiduals:clock'); a pulse more than half a turn off (4000 km
     at 50 Hz) gives the warning 'trueResiduals:wrap'.
  7. End to end with an SDR-like clock (D3c): the samples are taken at clock readings of a clock
     with eps(t) = 3 ms + 2e-6 (t - t1); the pulses come from the true probe of check 3. The
     observer folds with the assumed probe's predictor and its clock ESTIMATE (offset 10 us and
     rate 2e-7 off, makePredictor 'Clock'):
     a. at 200 times the predicted phase minus the data's phase equals phaseOffset with
        'ClockEstimate' (1e-7 turns);
     b. the folded pulse sits at phaseOffset (1e-4 turns);
     c. with a perfect estimate the clock drops out: residual = geometric (1e-12 s: the clock
        reading t + eps is an absolute time, rounded to 6e-8 s, times the rate 2e-6).
Errors at the end if any check fails. (21 Oct 2026, block D3a; check 7 block D3c.)
%}

root = fullfile(fileparts(mfilename('fullpath')), '..');
addpath(fullfile(root, 'functions'));
fails = {};
tmp = tempname(tempdir); mkdir(tmp);
cleanTmp = onCleanup(@() rmdir(tmp, 's'));
c = 299792458; AU = 149597870700;
t1 = 9555.5 * 86400;                                  % 2026-03-01 00:00 TT [s since J2000]
rs = RandStream('mt19937ar', 'Seed', 5);
b0329 = struct('f0', 1.399541538720, 'F1', 0, 'TRef', t1 - 1e7, 'Phi0', 0, 'raDeg', 53.24738, 'decDeg', 54.57878);
siteA = makeTrajectory('fixed', 'Lat', 52, 'Lon', 5);
siteT = makeTrajectory('fixed', 'Lat', 52, 'Lon', 5, 'OffsetENU', [100e3 0 0]);
t = t1 + sort(86400 * rand(rs, 1, 2000));

% 1. Same trajectory
r1 = trueResiduals(t, b0329, siteA, siteA);
ok1 = all(r1.residual == 0) && all(r1.phaseOffset == 0);
fprintf('1. same trajectory, perfect clock: largest |residual| %.1e s: %s\n', max(abs(r1.residual)), passStr(ok1));
if ~ok1, fails{end+1} = 'same trajectory'; end

% 2. The site 100 km east
r2 = trueResiduals(t, b0329, siteT, siteA, 'DelayArgs', {'Terms', "roemer"});
dr = siteT.ecef - siteA.ecef; th = era(t);
ra = deg2rad(b0329.raDeg); de = deg2rad(b0329.decDeg);
drn = ((cos(th) * dr(1) - sin(th) * dr(2)) * cos(de) * cos(ra) + (sin(th) * dr(1) + cos(th) * dr(2)) ...
       * cos(de) * sin(ra) + dr(3) * sin(de)) / c;    % (delta_T - delta_A) at t: the true site's R dr . n / c
rate = (delayModel(t + 1, siteA, b0329, 'Terms', "roemer").total - ...
        delayModel(t - 1, siteA, b0329, 'Terms', "roemer").total) / 2;
x2 = r2.geometric .* (1 + rate);                       % the residual without the Doppler factor
e2 = max(abs(x2 + drn));
amp = (max(x2) - min(x2)) / 2;
closer = drn > 0.1 * max(drn);                        % the true site clearly closer to the pulsar
r2f = trueResiduals(t, b0329, siteT, siteA);
eF = max(abs(r2f.residual - r2.residual));
ok2 = e2 < 1e-10 && abs(amp / (100e3 * cos(de) / c) - 1) < 1e-5 && all(r2.residual(closer) < 0) && ...
      eF < 1e-7 && r2.solveErr < 1e-12;
fprintf(['2. site 100 km east: residual (1 + rate) vs -(R dr . n)/c %.1e s, amplitude %.3f us (100 km cos(dec)/c ' ...
         '%.3f us); Doppler factor %.1e; early where closer %d; full model - Roemer %.1f ns; solve %.0e s: %s\n'], ...
    e2, amp * 1e6, 100e3 * cos(de) / c * 1e6, mean(rate), all(r2.residual(closer) < 0), eF * 1e9, ...
    r2.solveErr, passStr(ok2));
if ~ok2, fails{end+1} = 'site offset'; end

% 3. A probe along the line of sight, 1000 km closer
psr5 = struct('f0', 50, 'F1', 0, 'TRef', t1, 'Phi0', 0, 'raDeg', 30, 'decDeg', 10);
n = [cosd(10) * cosd(30); cosd(10) * sind(30); sind(10)];
v = 3e5 * n; r0 = 1.3 * AU * [0.6; 0.8; 0];
probeA = makeTrajectory('custom', 'Frame', 'SSB', 'PosFcn', @(tt) r0 + v * (tt - t1), ...
    'VelFcn', @(tt) repmat(v, 1, numel(tt)));
probeT = makeTrajectory('custom', 'Frame', 'SSB', 'PosFcn', @(tt) r0 + 1e6 * n + v * (tt - t1), ...
    'VelFcn', @(tt) repmat(v, 1, numel(tt)));
r3 = trueResiduals(t, psr5, probeT, probeA, 'DelayArgs', {'Terms', "roemer"});
x3 = -(1e6 / c) / (1 + 3e5 / c);
e3 = max(abs(r3.geometric - x3));
ok3 = e3 < 1e-10;
fprintf(['3. probe 1000 km closer, Doppler 1e-3: residual vs -(1000 km/c)/(1 + v/c) %.1e s (%.4f ms); ' ...
         'without the factor off by %.2f us: %s\n'], e3, x3 * 1e3, (1e6 / c + x3) * 1e6, passStr(ok3));
if ~ok3, fails{end+1} = 'probe'; end

% 4. Clock
clkF = @(tt) 1e-6 + 1e-9 * (tt - t1) + 2e-9 * sin(2 * pi * (tt - t1) / 3600);
r4a = trueResiduals(t, b0329, siteT, siteA, 'Clock', clkF);
r4b = trueResiduals(t, b0329, siteT, siteA, 'Clock', struct('offset', clkF));
ok4 = isequal(r4a.residual, r4b.residual) && isequal(r4a.clock, clkF(t)) && ...
      isequal(r4a.geometric, r2f.geometric) && isequal(r4a.residual, r4a.geometric + r4a.clock) && ...
      isequal(r4a.phaseOffset, r4a.fObs .* r4a.residual);
fprintf('4. clock: handle = struct, clock part = the handle, geometry unchanged, sums exact %d: %s\n', ok4, passStr(ok4));
if ~ok4, fails{end+1} = 'clock'; end

% 5. End to end: the data's pulses from the true probe, folded with the assumed probe's predictor
f0 = psr5.f0; dt = 20e-6; Ld = 20; NB = 1024; nPer = 100;
pred = makePredictor(psr5, probeA, [t1, t1 + Ld]);
tf = (0:round(Ld / dt) - 1) * dt;                     % file time = predictor local time
tg = 0:0.01:Ld;
dTg = delayModel(t1 + tg, probeT, psr5).total;        % the true delays (all terms), 10 ms grid
phiD = f0 * (tf + spline(tg, dTg, tf)) - pred.N0;      % the data's phase (Phi0 = 0, TRef = t1)
tau = sort(0.5 + (Ld - 1) * rand(rs, 1, 200));
r5 = trueResiduals(t1 + tau, psr5, probeT, probeA);
e5a = max(abs(pred.phase(tau) - (f0 * (tau + spline(tg, dTg, tau)) - pred.N0) - r5.phaseOffset));
sg = 0.02 / (2 * sqrt(2 * log(2)));
pulse = @(ph) 1 + 5 * exp(-0.5 * ((ph - round(ph)) / sg).^2);
det = writePower(fullfile(tmp, 'true.dat'), pulse(phiD), dt);
[~, fP] = foldProfile(det, 'x', f0, 'Predictor', pred, 'NBin', NB, 'SubintPeriods', nPer, ...
    'SaveFile', false, 'Verbose', false);
p = fP.profTotal(:) - 1; ph = (0:NB-1).' / NB;
[~, jm] = max(p);
d = mod(ph - ph(jm) + 0.5, 1) - 0.5; w = abs(d) < 0.1;
cen = ph(jm) + sum(d(w) .* p(w)) / sum(p(w));
off = mean(r5.phaseOffset);
e5b = abs(mod(cen - off + 0.5, 1) - 0.5);
ok5 = e5a < 1e-7 && e5b < 1e-4 && off < 0;
fprintf(['5. end to end: predicted - data phase vs phaseOffset %.1e turns; folded peak %.5f turns vs ' ...
         'phaseOffset %.5f (diff %.1e): %s\n'], e5a, mod(cen + 0.5, 1) - 0.5, off, e5b, passStr(ok5));
if ~ok5, fails{end+1} = 'end to end'; end

% 6. Input checks
probeW = makeTrajectory('custom', 'Frame', 'SSB', 'PosFcn', @(tt) r0 + 4e6 * n + v * (tt - t1), ...
    'VelFcn', @(tt) repmat(v, 1, numel(tt)));
okC = strcmp(errId(@() trueResiduals(t(1:5), b0329, siteT, siteA, 'Clock', 5)), 'trueResiduals:clock');
fprintf('6. (the next warning is expected)\n');
lastwarn('');
ws = warning('off', 'backtrace');
trueResiduals(t1 + tau(1:5), psr5, probeW, probeA);
warning(ws);
[~, wid] = lastwarn;
ok6 = okC && strcmp(wid, 'trueResiduals:wrap');
fprintf('6. bad clock refused %d; 4000 km at 50 Hz warns %d: %s\n', okC, strcmp(wid, 'trueResiduals:wrap'), passStr(ok6));
if ~ok6, fails{end+1} = 'input checks'; end

% 7. End to end with an SDR-like clock and the observer's estimate of it
b = 3e-3; d = 2e-6;
epsT = @(tt) b + d * (tt - t1);                       % the true clock error at true time
ehat = @(cc) (b - 1e-5) + (d - 2e-7) * (cc - t1);     % the observer's estimate at clock readings
pred7 = makePredictor(psr5, probeA, [t1, t1 + Ld], 'Clock', ehat);
trueOf = @(tc) (tc - b) / (1 + d);                    % true time since t1 of clock reading t1 + tc
tgT = -0.01:0.01:Ld;                                  % true-time grid of the true delays
dTg7 = delayModel(t1 + tgT, probeT, psr5).total;
phiOf = @(tc) f0 * (trueOf(tc) + spline(tgT, dTg7, trueOf(tc))) - pred7.N0;   % the data's phase
r7 = trueResiduals(t1 + trueOf(tau), psr5, probeT, probeA, 'Clock', epsT, 'ClockEstimate', ehat);
e7a = max(abs(pred7.phase(tau) - phiOf(tau) - r7.phaseOffset));
det7 = writePower(fullfile(tmp, 'clock.dat'), pulse(phiOf(tf)), dt);
[~, f7] = foldProfile(det7, 'x', f0, 'Predictor', pred7, 'NBin', NB, 'SubintPeriods', nPer, ...
    'SaveFile', false, 'Verbose', false);
p7 = f7.profTotal(:) - 1;
[~, jm] = max(p7);
d7 = mod(ph - ph(jm) + 0.5, 1) - 0.5; w7 = abs(d7) < 0.1;
cen7 = ph(jm) + sum(d7(w7) .* p7(w7)) / sum(p7(w7));
off7 = mean(r7.phaseOffset);
e7b = abs(mod(cen7 - off7 + 0.5, 1) - 0.5);
ehatP = @(cc) b + d * trueOf(cc - t1);                % a perfect estimate: eps at the reading's true time
r7p = trueResiduals(t1 + trueOf(tau), psr5, probeT, probeA, 'Clock', epsT, 'ClockEstimate', ehatP);
e7c = max(abs(r7p.residual - r7p.geometric));
ok7 = e7a < 1e-7 && e7b < 1e-4 && e7c < 1e-12;
fprintf(['7. SDR clock (3 ms, 2e-6) with the estimate 10 us / 2e-7 off: predicted - data phase vs phaseOffset %.1e turns; ' ...
         'folded peak %.5f vs %.5f (diff %.1e); perfect estimate: residual - geometric %.1e s: %s\n'], ...
    e7a, mod(cen7 + 0.5, 1) - 0.5, off7, e7b, e7c, passStr(ok7));
if ~ok7, fails{end+1} = 'clock end to end'; end

if isempty(fails)
    fprintf('testTrueResiduals: ALL PASSED\n');
else
    error('testTrueResiduals:failed', 'FAILED: %s', strjoin(fails, ', '));
end
end


% =================================================================================
function th = era(t)
% Earth rotation angle (IERS 2003), UT1 = TT - 69.184 s (as testDelayModel, independent of delayModel).
th = 2*pi * mod(0.7790572732640 + 1.00273781191135448 * (t - 69.184) / 86400, 1);
end

function det = writePower(file, p, dt)
% One-channel float32 power file and the detectChannels-like info foldProfile needs.
fid = fopen(file, 'w', 'ieee-le'); fwrite(fid, single(p), 'float32'); fclose(fid);
det = struct('file', file, 'nChan', 1, 'N', numel(p), 'binTime0', 0, 'binDt', dt, ...
    'fullySupportedBins', [1 numel(p)], 'byteOrder', 'ieee-le', 'chanFreqs', 0);
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
