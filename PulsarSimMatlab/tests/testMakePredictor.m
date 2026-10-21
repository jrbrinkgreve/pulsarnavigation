function testMakePredictor()
%TESTMAKEPREDICTOR  Unit tests for makePredictor (block D2a: the predicted pulse phase).
%{
Run from the PulsarSimMatlab folder: run('tests/testMakePredictor.m'), or
runAllTests("testMakePredictor"). ~5 s; no files. The truth is the timing model evaluated
directly with delayModel at the test times (phase differences phi(tau) - phi(0), so the
comparison does not depend on the predictor's reference N0 / frac0). A 100 Hz pulsar in B0329+54's
direction (fast spin stresses the precision), 2026.

  1. Fixed site, 52 N 5 E, two days, F1 = 0: phi(tau) - phi(0) = f0 (tau + delta(tau) - delta(0))
     at 2000 random times within Tol x f0 (1e-9 s); N0 + frac0 = the absolute phase at t1 within
     1e-5 turns (doubles of ~1e10 turns).
  2. Frequency: pred.freq = f0 (1 + d delta / dt) (central difference of delayModel, 1 s; 1e-9
     relative); the Doppler factor f / f0 - 1 within +-1.1e-4 (Earth's orbit + rotation).
  3. Segments: no step at any boundary beyond 2 Tol f0 (each side may be off by Tol f0), every
     segment's check within Tol.
  4. Spin-down: F1 = -1e-10 Hz/s (exaggerated): the F1 (u^2 - u0^2) / 2 = F1 du (du + 2 u0) / 2 term
     reproduced (Tol); the truth from differences (du = tau + delta - delta0), since absolute times
     (~8e8 s) carry 1e-7 s of rounding in double.
  5. A ship's 90 deg turn on a linear track, off-centre in a segment (3.3 h; a kink in the delay:
     the turn at a midpoint would be resolved by one split): the predictor splits segments around
     it repeatedly (more than the unrefined count) and stays within Tol next to the turn.
  6. A probe ('SSB' frame) moving at 30 km/s: phase within Tol.
  7. Evaluation outside the span and a reversed span are refused with their ids.
  8. The observer's clock estimate ('Clock', block D3c): epsHat(c) = 2.5 ms + 1.7e-6 (c - t1) +
     1.5e-13 (c - t1)^2 (an SDR-class rate and ageing). The phase at clock readings equals the
     direct computation at the estimated true times t = c - epsHat(c) (Tol); N0 + frac0 the
     absolute phase there (1e-5 turns); freq the phase's slope per clock second (central
     difference, 1e-9 relative); a zero estimate gives bit-identical coefficients; a bad clock is
     refused ('makePredictor:clock').
Errors at the end if any check fails. (21 Oct 2026, block D2a; test 8 block D3c.)
%}

root = fullfile(fileparts(mfilename('fullpath')), '..');
addpath(fullfile(root, 'functions'));
fails = {};
t1 = 9555.5 * 86400;                                  % 2026-03-01 00:00 TT [s since J2000]
psr = struct('f0', 100, 'F1', 0, 'TRef', t1 - 1e7, 'Phi0', 0.3, 'raDeg', 53.24738, 'decDeg', 54.57878);
site = makeTrajectory('fixed', 'Lat', 52, 'Lon', 5);
L = 2 * 86400; tol = 1e-9;
rs = RandStream('mt19937ar', 'Seed', 3);

% ---------------------------------------------------------------------------------
% 1. Fixed site, two days
% ---------------------------------------------------------------------------------
pr = makePredictor(psr, site, [t1, t1 + L]);
tau = sort(L * rand(rs, 1, 2000));
dl = delayModel(t1 + [0, tau], site, psr).total;
truth = psr.f0 * (tau + dl(2:end) - dl(1));
e1 = max(abs((pr.phase(tau) - pr.phase(0)) - truth)) / psr.f0;
abs0 = psr.Phi0 + psr.f0 * ((t1 - psr.TRef) + dl(1));
eAbs = abs((pr.N0 + pr.frac0) - abs0);
ok1 = e1 < tol && eAbs < 1e-5;
fprintf(['1. fixed site, 2 days, %d segment(s): phase vs direct %.2e s (Tol %.0e); N0 + frac0 vs absolute ' ...
         '%.1e turns: %s\n'], size(pr.seg, 1), e1, tol, eAbs, passStr(ok1));
if ~ok1, fails{end+1} = 'fixed site'; end

% ---------------------------------------------------------------------------------
% 2. Frequency = f0 (1 + d delta/dt)
% ---------------------------------------------------------------------------------
tq = tau(1:20:end);
dd = (delayModel(t1 + tq + 1, site, psr).total - delayModel(t1 + tq - 1, site, psr).total) / 2;
e2 = max(abs(pr.freq(tq) ./ (psr.f0 * (1 + dd)) - 1));
dop = pr.fRange / psr.f0 - 1;
ok2 = e2 < 1e-9 && all(abs(dop) < 1.1e-4);
fprintf('2. frequency vs f0 (1 + d delta/dt) %.1e relative; Doppler factor %.2e ... %.2e: %s\n', e2, dop, passStr(ok2));
if ~ok2, fails{end+1} = 'frequency'; end

% ---------------------------------------------------------------------------------
% 3. Segments
% ---------------------------------------------------------------------------------
b = pr.seg(1:end-1, 2).';
jump = max(abs(pr.phase(b - 1e-6) - pr.phase(b + 1e-6) + 2e-6 * pr.freq(b)));
ok3 = jump < 2 * tol * psr.f0 && pr.maxErr <= tol;
fprintf('3. segment boundaries: largest step %.1e turns; largest check error %.1e s: %s\n', jump, pr.maxErr, ...
    passStr(ok3));
if ~ok3, fails{end+1} = 'segments'; end

% ---------------------------------------------------------------------------------
% 4. Spin-down
% ---------------------------------------------------------------------------------
psF = psr; psF.F1 = -1e-10;
pF = makePredictor(psF, site, [t1, t1 + L]);
du = tau + (dl(2:end) - dl(1));                      % (t_SSB - TRef) - u0, from differences
u0 = (t1 - psF.TRef) + dl(1);
truthF = psF.f0 * du + 0.5 * psF.F1 * du .* (du + 2 * u0);
e4 = max(abs((pF.phase(tau) - pF.phase(0)) - truthF)) / psF.f0;
ok4 = e4 < tol;
fprintf('4. F1 = %.0e Hz/s: phase vs direct %.2e s (F1 term up to %.1f turns): %s\n', psF.F1, e4, ...
    max(abs(0.5 * psF.F1 * du .* (du + 2 * u0))), passStr(ok4));
if ~ok4, fails{end+1} = 'spin-down'; end

% ---------------------------------------------------------------------------------
% 5. A ship's 90 deg turn on a linear track
% ---------------------------------------------------------------------------------
dlon = 10 * 3 * 3600 / (111.32e3 * cosd(52)); dlat = 10 * 3 * 3600 / 111.2e3;
ship = makeTrajectory('track', 'Time', t1 + [0 3.3 6] * 3600, 'Lat', [52 52 52 + dlat], 'Lon', [5 5 + dlon 5 + dlon]);
ps5 = makePredictor(psr, ship, [t1, t1 + 6 * 3600]);
tt = 3.3 * 3600 + [-200:10:-10, -1, 1, 10:10:200];    % next to the turn
d5 = delayModel(t1 + [0, tt], ship, psr).total;
e5 = max(abs((ps5.phase(tt) - ps5.phase(0)) - psr.f0 * (tt + d5(2:end) - d5(1)))) / psr.f0;
ok5 = size(ps5.seg, 1) > 6 && e5 < tol;
fprintf('5. ship turning 90 deg: %d segments (6 unrefined, shortest %.0f s); phase next to the turn vs direct %.2e s: %s\n', ...
    size(ps5.seg, 1), min(diff(ps5.seg, 1, 2)), e5, passStr(ok5));
if ~ok5, fails{end+1} = 'ship turn'; end

% ---------------------------------------------------------------------------------
% 6. A probe in the barycentric frame
% ---------------------------------------------------------------------------------
r0 = [1.2; 0.3; -0.1] * 149597870700; v0 = [-5e3; 3e4; 1e3];
probe = makeTrajectory('custom', 'Frame', 'SSB', 'PosFcn', @(t) r0 + v0 * (t - t1), 'VelFcn', @(t) repmat(v0, 1, numel(t)));
pp = makePredictor(psr, probe, [t1, t1 + L]);
d6 = delayModel(t1 + [0, tau], probe, psr).total;
e6 = max(abs((pp.phase(tau) - pp.phase(0)) - psr.f0 * (tau + d6(2:end) - d6(1)))) / psr.f0;
ok6 = e6 < tol;
fprintf('6. probe (SSB frame, 30 km/s): phase vs direct %.2e s, Doppler %.2e ... %.2e: %s\n', e6, ...
    pp.fRange / psr.f0 - 1, passStr(ok6));
if ~ok6, fails{end+1} = 'probe'; end

% ---------------------------------------------------------------------------------
% 7. Input checks
% ---------------------------------------------------------------------------------
ok7 = strcmp(errId(@() pr.phase(-1)), 'makePredictor:range') && ...
      strcmp(errId(@() makePredictor(psr, site, [t1, t1 - 1])), 'makePredictor:span');
fprintf('7. outside the span and a reversed span refused %d: %s\n', ok7, passStr(ok7));
if ~ok7, fails{end+1} = 'input checks'; end

% ---------------------------------------------------------------------------------
% 8. The observer's clock estimate
% ---------------------------------------------------------------------------------
eh = @(c) 2.5e-3 + 1.7e-6 * (c - t1) + 1.5e-13 * (c - t1).^2;
pc = makePredictor(psr, site, [t1, t1 + L], 'Clock', eh);
phiC = @(tv) psr.f0 * ((tv - (eh(t1 + tv) - eh(t1))) + delayModel(t1 + tv - eh(t1 + tv), site, psr).total);
tv = [0, tau];
ph8 = phiC(tv);
e8 = max(abs((pc.phase(tau) - pc.phase(0)) - (ph8(2:end) - ph8(1)))) / psr.f0;
tt0 = t1 - eh(t1);
abs8 = psr.Phi0 + psr.f0 * ((tt0 - psr.TRef) + delayModel(tt0, site, psr).total);
eAbs8 = abs((pc.N0 + pc.frac0) - abs8);
e8f = max(abs(pc.freq(tq) ./ ((phiC(tq + 1) - phiC(tq - 1)) / 2) - 1));
pz = makePredictor(psr, site, [t1, t1 + L], 'Clock', @(c) zeros(size(c)));
ok8 = e8 < tol && eAbs8 < 1e-5 && e8f < 1e-9 && isequal(pz.coef, pr.coef) && isequal(pz.seg, pr.seg) && ...
      strcmp(errId(@() makePredictor(psr, site, [t1, t1 + 10], 'Clock', 5)), 'makePredictor:clock');
fprintf(['8. clock estimate (2.5 ms, 1.7e-6, ageing): phase vs direct %.2e s; N0 + frac0 %.1e turns; freq vs ' ...
         'slope %.1e; zero estimate bit-identical %d; Doppler + rate %.3e ... %.3e: %s\n'], e8, eAbs8, e8f, ...
    isequal(pz.coef, pr.coef), pc.fRange / psr.f0 - 1, passStr(ok8));
if ~ok8, fails{end+1} = 'clock estimate'; end

if isempty(fails)
    fprintf('testMakePredictor: ALL PASSED\n');
else
    error('testMakePredictor:failed', 'FAILED: %s', strjoin(fails, ', '));
end
end


% =================================================================================
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
