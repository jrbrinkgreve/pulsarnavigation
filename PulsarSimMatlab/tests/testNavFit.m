function testNavFit()
%TESTNAVFIT  Unit tests for navFit (block D4a: the position of a fixed site from TOA residuals).
%{
Run from the PulsarSimMatlab folder: run('tests/testNavFit.m'), or runAllTests("testNavFit").
~20 s; no files. The truth comes from trueResiduals (delayModel on the true and the assumed site),
independent of navFit's own design matrix. Schedule: B0329+54 from 52 N 5 E, 300 TOAs at random
times over 30 days (2026-03-01 on), sigma 734 us each - the scenario's nominal case.

  1. Noise-free residuals of a true site [30 -20 10] km ENU from the assumed one:
     "calibrated" recovers all three ENU components, "self" the east and the radial (equatorial)
     components and gives no polar one; both within 50 m (the residual's 1/(1 + Doppler) factor
     and the Einstein site term, ~2e-4 of the offset, are not in navFit's model).
  2. Monte Carlo, 400 realizations of Gaussian noise on the residuals of a site 100 km east,
     "self": the mean estimate within 4 x its scatter / sqrt(400) of the truth; the pulls
     (estimate - truth) / stated error have std 1 (+- 4 / sqrt(800)) for east and radial; the
     stated error matches the expected +-31 km.
  3. The same with 2 % outliers at 5-10 sigma (random sign): the robust fit stays unbiased with
     pulls of std 1 (+- 0.15), flags at least 90 % of the outliers and at most 1 % of the good
     TOAs; plain weighted least squares ('Robust' false) has pulls of std > 1.3.
  4. "calibrated", 400 realizations: the polar component is recovered (pulls std 1 +- 0.15) with
     an error of about +-16 km, better than the equatorial ones (the constant uses every TOA).
  5. 'Predict': the covariance equals that of the plain weighted fit of the same schedule
     (1e-12 relative).
  6. Two pulsars (B0329+54 and one at RA 120, Dec +20), noise-free: both starts recover the truth
     within 50 m, with one constant per pulsar in "self".
  7. A site that is not a fixed trajectory, missing residuals and a bad clock are refused
     ('navFit:site', 'navFit:toa', 'navFit:clock').
  The clock (block D4b):
  8. 'Clock' [0 0 0] (no wander), "calibrated" (no rate column), gives the fit without a clock (1e-9
     relative).
  9. A crystal clock (makeClock "crystal", plus a 1 ms offset and a 1e-9 rate the fit does not
     know) on the residuals of the site 100 km east, 200 realizations, "self": with 'Clock' the
     east estimate is unbiased with pulls of std 1 (+- 0.2) and an error close to +-31 km; without
     it the clock's milliseconds of wander make the robust fit flag a large part of the TOAs
     (shown: the fraction and the core chi2).
 10. The clock estimate, "calibrated", a clock without flicker (h_1 = 0: the fit's covariance model is
     then exact), 200 realizations: the fitted wander at the TOAs and the prediction one day after
     the last TOA match the true clock within their stated errors (pulls of std 1 +- 0.2). Shown
     for information: the same with the full crystal (flicker left out of the fit's model).
Errors at the end if any check fails. (22 Oct 2026, block D4a; tests 8-10 block D4b.)
%}

root = fullfile(fileparts(mfilename('fullpath')), '..');
addpath(fullfile(root, 'functions'));
fails = {};
rs = RandStream('mt19937ar', 'Seed', 41);
t1 = 9555.5 * 86400;                                  % 2026-03-01 00:00 TT
N = 300; sig = 734e-6;
t = t1 + sort(30 * 86400 * rand(rs, 1, N));
b0329 = struct('raDeg', 53.2474, 'decDeg', 54.5787);
psrT = struct('f0', 1.399541538720, 'F1', 0, 'TRef', t1, 'raDeg', b0329.raDeg, 'decDeg', b0329.decDeg);
siteA = makeTrajectory('fixed', 'Lat', 52, 'Lon', 5);
lon = 5;
eRad = [cosd(lon); sind(lon); 0];                      % radial: equatorial plane, away from the axis

% 1. Noise-free, a general offset
off = [30e3 -20e3 10e3];
siteT = makeTrajectory('fixed', 'Lat', 52, 'Lon', 5, 'OffsetENU', off);
r1 = trueResiduals(t, psrT, siteT, siteA).residual;
toa = struct('t', t, 'r', r1, 'sigma', sig * ones(1, N));
fc = navFit(toa, b0329, siteA, 'Start', "calibrated");
fs = navFit(toa, b0329, siteA);
dx = siteT.ecef - siteA.ecef;
e1c = max(abs(fc.enu - off(:)));
e1s = max(abs([fs.east - off(1), fs.radial - eRad' * dx]));
ok1 = e1c < 50 && e1s < 50 && isnan(fs.polar) && abs(fc.polar - dx(3)) < 50;
fprintf(['1. noise-free [30 -20 10] km ENU: "calibrated" ENU %s km (max error %.1f m), polar %.3f km; "self" east %.3f, ' ...
         'radial %.3f km (truth %.3f; max error %.1f m), polar NaN %d: %s\n'], mat2str(round(fc.enu.' / 1e3, 3)), e1c, ...
    fc.polar / 1e3, fs.east / 1e3, fs.radial / 1e3, eRad' * dx / 1e3, e1s, isnan(fs.polar), passStr(ok1));
if ~ok1, fails{end+1} = 'noise-free'; end

% 2. Monte Carlo, Gaussian noise, 100 km east
siteE = makeTrajectory('fixed', 'Lat', 52, 'Lon', 5, 'OffsetENU', [100e3 0 0]);
rE = trueResiduals(t, psrT, siteE, siteA).residual;
truthE = 100e3; truthR = eRad' * (siteE.ecef - siteA.ecef);
nMC = 400;
[est, err] = deal(zeros(nMC, 2));
for m = 1:nMC
    f = navFit(struct('t', t, 'r', rE + sig * randn(rs, 1, N), 'sigma', sig * ones(1, N)), b0329, siteA);
    est(m, :) = [f.east, f.radial]; err(m, :) = [f.eastErr, f.radialErr];
end
pull = (est - [truthE, truthR]) ./ err;
bias = mean(est) - [truthE, truthR]; sc = std(est);
ok2 = all(abs(bias) < 4 * sc / sqrt(nMC)) && all(abs(std(pull) - 1) < 4 / sqrt(2 * nMC)) && ...
      abs(median(err(:, 1)) / 31e3 - 1) < 0.1;
fprintf(['2. Gaussian noise, %d realizations, 100 km east: east %.1f +- %.1f km (stated %.1f), radial %+.1f +- %.1f km; ' ...
         'bias %+.2f / %+.2f km; pull std %.3f / %.3f: %s\n'], nMC, mean(est(:, 1)) / 1e3, sc(1) / 1e3, ...
    median(err(:, 1)) / 1e3, mean(est(:, 2)) / 1e3, sc(2) / 1e3, bias / 1e3, std(pull), passStr(ok2));
if ~ok2, fails{end+1} = 'Monte Carlo'; end

% 3. 2 % outliers at 5-10 sigma
[estR, errR, estW, errW] = deal(zeros(nMC, 1)); [hit, falseFlag, nBad] = deal(zeros(nMC, 1));
for m = 1:nMC
    rr = rE + sig * randn(rs, 1, N);
    bad = rand(rs, 1, N) < 0.02;
    rr(bad) = rr(bad) + sig * (5 + 5 * rand(rs, 1, nnz(bad))) .* sign(randn(rs, 1, nnz(bad)));
    tq = struct('t', t, 'r', rr, 'sigma', sig * ones(1, N));
    fR = navFit(tq, b0329, siteA);
    fW = navFit(tq, b0329, siteA, 'Robust', false);
    estR(m) = fR.east; errR(m) = fR.eastErr; estW(m) = fW.east; errW(m) = fW.eastErr;
    hit(m) = sum(fR.outlier & bad); falseFlag(m) = sum(fR.outlier & ~bad);
    nBad(m) = nnz(bad);
end
pR = (estR - truthE) ./ errR; pW = (estW - truthE) ./ errW;
fHit = sum(hit) / sum(nBad); fFalse = sum(falseFlag) / (nMC * N - sum(nBad));
ok3 = abs(mean(estR) - truthE) < 4 * std(estR) / sqrt(nMC) && abs(std(pR) - 1) < 0.15 && std(pW) > 1.3 && ...
      fHit >= 0.9 && fFalse <= 0.01;
fprintf(['3. 2 %% outliers at 5-10 sigma: robust east %.1f km, pull std %.3f; least squares pull std %.3f; outliers flagged ' ...
         '%.1f %%, good TOAs flagged %.2f %%: %s\n'], mean(estR) / 1e3, std(pR), std(pW), 100 * fHit, 100 * fFalse, passStr(ok3));
if ~ok3, fails{end+1} = 'outliers'; end

% 4. "calibrated": the polar component
truthP = siteE.ecef(3) - siteA.ecef(3);
[estP, errP] = deal(zeros(nMC, 1));
for m = 1:nMC
    f = navFit(struct('t', t, 'r', rE + sig * randn(rs, 1, N), 'sigma', sig * ones(1, N)), b0329, siteA, 'Start', "calibrated");
    estP(m) = f.polar; errP(m) = f.polarErr;
end
pP = (estP - truthP) ./ errP;
ok4 = abs(std(pP) - 1) < 0.15 && abs(mean(pP)) < 4 / sqrt(nMC) * std(pP) && median(errP) < median(err(:, 1));
fprintf(['4. "calibrated": polar %.1f +- %.1f km (stated %.1f; east and radial %.1f), pulls mean %+.3f std %.3f: %s\n'], ...
    mean(estP) / 1e3, std(estP) / 1e3, median(errP) / 1e3, median(err(:, 1)) / 1e3, mean(pP), std(pP), passStr(ok4));
if ~ok4, fails{end+1} = 'calibrated'; end

% 5. Predict
fPr = navFit(struct('t', t, 'sigma', sig * ones(1, N)), b0329, siteA, 'Predict', true);
fWl = navFit(toa, b0329, siteA, 'Robust', false);
e5 = max(abs(fPr.cov(:) - fWl.cov(:))) / max(abs(fWl.cov(:)));
ok5 = e5 < 1e-12 && all(isnan(fPr.p));
fprintf('5. Predict = covariance of the weighted fit: %.1e relative; east +-%.1f km: %s\n', e5, fPr.eastErr / 1e3, passStr(ok5));
if ~ok5, fails{end+1} = 'Predict'; end

% 6. Two pulsars
p2 = struct('raDeg', 120, 'decDeg', 20);
psr2 = psrT; psr2.raDeg = p2.raDeg; psr2.decDeg = p2.decDeg;
k = 1 + (rand(rs, 1, N) < 0.5);
r6 = r1;
r6(k == 2) = trueResiduals(t(k == 2), psr2, siteT, siteA).residual;
toa6 = struct('t', t, 'r', r6, 'sigma', sig * ones(1, N), 'k', k);
f6c = navFit(toa6, [b0329, p2], siteA, 'Start', "calibrated");
f6s = navFit(toa6, [b0329, p2], siteA);
e6 = max([abs(f6c.enu - off(:)); abs(f6s.east - off(1)); abs(f6s.radial - eRad' * dx)]);
ok6 = e6 < 50 && numel(f6s.p) == 4 && isnan(f6s.polar);
fprintf('6. two pulsars, noise-free: max error %.1f m; "self" has %d parameters (dX dY c1 c2): %s\n', e6, numel(f6s.p), passStr(ok6));
if ~ok6, fails{end+1} = 'two pulsars'; end

% 7. Input checks
ok7 = strcmp(errId(@() navFit(toa, b0329, struct('type', 'track'))), 'navFit:site') && ...
      strcmp(errId(@() navFit(struct('t', t, 'sigma', sig * ones(1, N)), b0329, siteA)), 'navFit:toa') && ...
      strcmp(errId(@() navFit(toa, b0329, siteA, 'Clock', [1 -1 0])), 'navFit:clock');
fprintf('7. bad site, missing residuals and a bad clock refused %d: %s\n', ok7, passStr(ok7));
if ~ok7, fails{end+1} = 'input checks'; end

% 8. A clock without wander = no clock
f8a = navFit(toa, b0329, siteA, 'Start', "calibrated"); f8b = navFit(toa, b0329, siteA, 'Start', "calibrated", 'Clock', [0 0 0]);
e8 = max(abs(f8b.p - f8a.p) ./ sqrt(diag(f8a.cov)));
ok8 = e8 < 1e-9 && max(abs(f8b.clock)) == 0;
fprintf('8. ''Clock'' [0 0 0] = no clock: %.1e of the errors, wander 0 %d: %s\n', e8, max(abs(f8b.clock)) == 0, passStr(ok8));
if ~ok8, fails{end+1} = 'clock off'; end

% 9. A crystal clock with an unknown offset and rate, "self"
hCr = [1.8e-21 6.492e-22 1.37e-24];
nC = 200; tEnd = t(end) + 2 * 86400;
[e9, s9, flagged, chi2No, noConv] = deal(zeros(nC, 1));
ws = warning('off', 'navFit:converge');                % the fit without the clock model may not settle
for m = 1:nC
    clk = makeClock([t(1), tEnd], 'Preset', "crystal", 'Offset', 1e-3, 'Rate', 1e-9, 'Seed', 300 + m);
    rr = rE + clk.offset(t) + sig * randn(rs, 1, N);
    tq = struct('t', t, 'r', rr, 'sigma', sig * ones(1, N));
    f = navFit(tq, b0329, siteA, 'Clock', hCr);
    e9(m) = f.east; s9(m) = f.eastErr;
    fN = navFit(tq, b0329, siteA);
    flagged(m) = mean(fN.outlier); chi2No(m) = fN.chi2Core; noConv(m) = fN.nIter == 50;
end
warning(ws);
p9 = (e9 - truthE) ./ s9;
ok9 = abs(mean(p9)) < 4 / sqrt(nC) * std(p9) && abs(std(p9) - 1) < 0.2 && median(s9) < 1.2 * median(err(:, 1));
fprintf(['9. crystal clock + 1 ms + 1e-9, "self", %d realizations: east %.1f km, pulls mean %+.2f std %.3f, error %.1f km ' ...
         '(no clock %.1f); without ''Clock'': %.0f %% of the TOAs flagged, core chi2 %.1f, not settled in %.0f %%: %s\n'], nC, mean(e9) / 1e3, ...
    mean(p9), std(p9), median(s9) / 1e3, median(err(:, 1)) / 1e3, 100 * median(flagged), median(chi2No), 100 * mean(noConv), passStr(ok9));
if ~ok9, fails{end+1} = 'clock in the fit'; end

% 10. The clock estimate and its prediction, "calibrated"
[zc, zp, zcF, zpF] = deal([]);
for m = 1:nC
    for full = [false true]
        hT = hCr; if ~full, hT(2) = 0; end
        clk = makeClock([t(1), tEnd], 'H', hT, 'Seed', 600 + m);
        rr = rE + clk.offset(t) + sig * randn(rs, 1, N);
        f = navFit(struct('t', t, 'r', rr, 'sigma', sig * ones(1, N)), b0329, siteA, 'Start', "calibrated", 'Clock', hCr);
        tp = t(end) + 86400;
        known = f.clockErr == 0;                        % the start: the clock is known exactly there
        z1 = (clk.offset(t(~known)) - f.clock(~known)) ./ f.clockErr(~known);
        z2 = (clk.offset(tp) - f.clockAt(tp)) / f.clockAtErr(tp);
        if full, zcF = [zcF, z1]; zpF = [zpF, z2]; else, zc = [zc, z1]; zp = [zp, z2]; end %#ok<AGROW>
    end
end
ok10 = abs(std(zc) - 1) < 0.2 && abs(std(zp) - 1) < 0.2;
fprintf(['10. clock estimate, "calibrated", %d realizations: at the TOAs pulls std %.3f (error median %.0f us), one day ' ...
         'ahead std %.3f (error %.0f us); with flicker in the true clock: %.3f, %.3f: %s\n'], nC, std(zc), ...
    1e6 * median(f.clockErr), std(zp), 1e6 * f.clockAtErr(t(end) + 86400), std(zcF), std(zpF), passStr(ok10));
if ~ok10, fails{end+1} = 'clock estimate'; end

if isempty(fails)
    fprintf('testNavFit: ALL PASSED\n');
else
    error('testNavFit:failed', 'FAILED: %s', strjoin(fails, ', '));
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
