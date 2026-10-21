function testMakeClock()
%TESTMAKECLOCK  Unit tests for makeClock (block D3b: the receiver clock's error).
%{
Run from the PulsarSimMatlab folder: run('tests/testMakeClock.m'), or runAllTests("testMakeClock").
~2 s; no files. The truths are the power-law formulas of Van Dierendonck, McGraw & Brown (1984):
the Allan variance sigma_y^2(tau) = h0 / (2 tau) + 2 ln2 h_1 + (2 pi)^2 / 6 h_2 tau (eqs. 6-8), and
the time error's growth from a known start (eq. 60, without flicker):
Q11 = h0 / 2 Dt + 2 pi^2 / 3 h_2 Dt^3.

  1. Perfect clock: offset(t) = Offset + Rate (t - t1) + Drift (t - t1)^2 / 2 and rate(t) = Rate +
     Drift (t - t1) (1e-15 relative); as trueResiduals' 'Clock' with one trajectory for truth and
     assumption, the residual is exactly the clock error.
  2. Seeds and evaluation: the same seed gives the same clock, another seed another; at the grid
     nodes offset equals the drawn x exactly; between them it is linear, and rate is its slope
     (1e-9 relative).
  3. Allan variance: 40 clocks of 20,000 steps of 60 s per noise type (white, flicker, random walk,
     and the paper's crystal with all three), overlapping Allan variance at 5 ... 1000 steps
     (5 min ... 17 h) against sigma_y^2: within 3 % + 3 x the scatter of the mean. Shown for
     information: flicker at 2 steps, ~9 % high (help text).
  4. The h-to-Kalman relation without flicker: the crystal's h0 and h_2, 2000 clocks on a 1 s
     grid, the variance of x 200 s after the start against Q11 (8.5 ns; 3 x the scatter of a
     variance from 2000 samples, ~10 %).
  5. Times outside the grid, a reversed span and a negative h are refused with their ids.
Errors at the end if any check fails. (21 Oct 2026, block D3b.)
%}

root = fullfile(fileparts(mfilename('fullpath')), '..');
addpath(fullfile(root, 'functions'));
fails = {};
t1 = 9555.5 * 86400;                                  % 2026-03-01 00:00 TT [s since J2000]
hCr = [1.8e-21 6.492e-22 1.37e-24];                   % the paper's crystal (its Figure 3)

% 1. Perfect clock
off = 2.5e-3; rate = 1.7e-6; drift = 3e-13;           % an SDR-like rate (+-2 ppm class), ageing
c1 = makeClock([t1, t1 + 86400], 'Offset', off, 'Rate', rate, 'Drift', drift);
tq = t1 + sort(86400 * rand(RandStream('mt19937ar', 'Seed', 2), 1, 500));
e1 = max(abs(c1.offset(tq) - (off + rate * (tq - t1) + 0.5 * drift * (tq - t1).^2))) / max(abs(c1.offset(tq)));
e1r = max(abs(c1.rate(tq) - (rate + drift * (tq - t1)))) / rate;
site = makeTrajectory('fixed', 'Lat', 52, 'Lon', 5);
psr = struct('f0', 1.399541538720, 'F1', 0, 'TRef', t1, 'raDeg', 53.24738, 'decDeg', 54.57878);
r1 = trueResiduals(tq, psr, site, site, 'Clock', c1);
ok1 = e1 < 1e-15 && e1r < 1e-15 && isequal(r1.residual, c1.offset(tq)) && all(c1.grid.x == 0);
fprintf(['1. perfect clock: offset vs formula %.1e, rate %.1e (relative); trueResiduals residual = clock ' ...
         'error %d: %s\n'], e1, e1r, isequal(r1.residual, c1.offset(tq)), passStr(ok1));
if ~ok1, fails{end+1} = 'perfect clock'; end

% 2. Seeds and evaluation
cA = makeClock([t1, t1 + 86400], 'Preset', "crystal", 'Seed', 7);
cB = makeClock([t1, t1 + 86400], 'Preset', "crystal", 'Seed', 7);
cC = makeClock([t1, t1 + 86400], 'Preset', "crystal", 'Seed', 8);
tg = cA.grid.t; xg = cA.grid.x;
mid = (tg(1:end-1) + tg(2:end)) / 2;
eLin = max(abs(cA.offset(mid) - (xg(1:end-1) + xg(2:end)) / 2)) / max(abs(xg));
eRate = max(abs(cA.rate(mid) - diff(xg) / cA.step)) / max(abs(diff(xg) / cA.step));
ok2 = isequal(cA.grid.x, cB.grid.x) && ~isequal(cA.grid.x, cC.grid.x) && isequal(cA.offset(tg), xg) && ...
      eLin < 1e-9 && eRate < 1e-9;
fprintf(['2. same seed same clock %d, other seed differs %d; nodes exact %d; linear between %.1e, rate = ' ...
         'slope %.1e: %s\n'], isequal(cA.grid.x, cB.grid.x), ~isequal(cA.grid.x, cC.grid.x), ...
    isequal(cA.offset(tg), xg), eLin, eRate, passStr(ok2));
if ~ok2, fails{end+1} = 'seeds and evaluation'; end

% 3. Allan variance
D = 60; K = 20000; ms = [2 5 10 30 100 300 1000]; nR = 40;
types = {[1e-20 0 0], [0 1e-21 0], [0 0 1e-24], hCr};
names = {'white', 'flicker', 'random walk', 'crystal'};
ok3 = true;
for i = 1:numel(types)
    s2 = zeros(nR, numel(ms));
    for r = 1:nR
        c = makeClock([t1, t1 + K * D], 'H', types{i}, 'Step', D, 'Seed', r);
        s2(r, :) = allanVar(c.grid.x, ms, D);
    end
    ratio = mean(s2) ./ c.sigmaY(ms * D).^2;
    sc = std(s2) ./ sqrt(nR) ./ mean(s2);
    use = ms >= 5;
    okI = all(abs(ratio(use) - 1) < 0.03 + 3 * sc(use));
    ok3 = ok3 && okI;
    fprintf('3. %-11s sim / sigma_y^2 at %s steps: %s (scatter %s): %s\n', names{i}, mat2str(ms), ...
        sprintf('%.3f ', ratio), sprintf('%.3f ', sc), passStr(okI));
end
if ~ok3, fails{end+1} = 'Allan variance'; end

% 4. The h-to-Kalman relation (eq. 60 without flicker) from a known start
nR4 = 2000; Dt = 200; x200 = zeros(nR4, 1);
for r = 1:nR4
    c = makeClock([t1, t1 + Dt], 'H', [hCr(1) 0 hCr(3)], 'Step', 1, 'Seed', 100 + r);
    x200(r) = c.offset(t1 + Dt);
end
q11 = hCr(1) / 2 * Dt + 2 * pi^2 / 3 * hCr(3) * Dt^3;
rat4 = mean(x200.^2) / q11;
ok4 = abs(rat4 - 1) < 3 * sqrt(2 / nR4);
fprintf(['4. crystal without flicker, 200 s after the start: rms %.2f ns vs sqrt(Q11) %.2f ns (variance ratio %.3f ' ...
         '+- %.3f): %s\n'], rms(x200) * 1e9, sqrt(q11) * 1e9, rat4, sqrt(2 / nR4), passStr(ok4));
if ~ok4, fails{end+1} = 'Kalman relation'; end

% 5. Input checks
ok5 = strcmp(errId(@() cA.offset(t1 - 1)), 'makeClock:range') && ...
      strcmp(errId(@() cA.rate(t1 + 2 * 86400)), 'makeClock:range') && ...
      strcmp(errId(@() makeClock([t1, t1 - 1])), 'makeClock:span') && ...
      strcmp(errId(@() makeClock([t1, t1 + 10], 'H', [1 -1 0])), 'makeClock:h');
fprintf('5. outside the grid, reversed span, negative h refused %d: %s\n', ok5, passStr(ok5));
if ~ok5, fails{end+1} = 'input checks'; end

if isempty(fails)
    fprintf('testMakeClock: ALL PASSED\n');
else
    error('testMakeClock:failed', 'FAILED: %s', strjoin(fails, ', '));
end
end


% =================================================================================
function s2 = allanVar(x, ms, D)
% Overlapping Allan variance from the time error x on a grid of step D, at tau = m D.
s2 = zeros(1, numel(ms));
for j = 1:numel(ms)
    m = ms(j);
    d = x(1 + 2*m:end) - 2 * x(1 + m:end - m) + x(1:end - 2*m);
    s2(j) = mean(d.^2) / (2 * (m * D)^2);
end
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
