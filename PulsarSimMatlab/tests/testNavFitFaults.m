function testNavFitFaults()
%TESTNAVFITFAULTS  navFit under the faults of the risk register: each case puts one fault into the truth.
%{
Run from the PulsarSimMatlab folder: run('tests/testNavFitFaults.m'), or runAllTests("testNavFitFaults").
~10 s; no files. The risk register (logs/2026-11-02_risk-register.md) names faults that the fit
would read as position, without a warning. Each case puts one of them into the truth (trueResiduals,
truth side) and checks the fit's response against an analytic prediction, so that the register's
numbers rest on a test. Schedule as testNavFit: B0329+54 from 52 N 5 E, sigma 734 us per TOA, from
2026-03-01 on; here the TOAs are evenly spaced, one every 2.4 h (10 per day), so that the
predictions hold to a fraction of a percent.

  1. FM-11, the oscillator's daily temperature term. A crystal follows its temperature, which
     cycles with the solar day. A daily rate term of amplitude a gives the clock error
       eps(t) = A sin(2 pi (t - t1) / 86400 s + phi),   A = a 86400 s / (2 pi)
     (137.5 us at a = 1e-8). The position signal is a sine at the sidereal day, 4 min shorter.
     A least-squares sidereal sine over a span T keeps sinc(dw T / 2) of the clock term's
     amplitude, dw = 2 pi (1/86164.0905 - 1/86400) rad/s (0.989 over 30 days, 0.881 over 100),
     and an equatorial offset d gives the amplitude d cos(dec) / c (193.33 us per 100 km). So
     the fit moves the site in the equatorial plane by
       d = c A sinc(dw T / 2) / cos(dec):   70.3 km at 30 days, 62.7 km at 100, for a = 1e-8.
     Checked: noise-free, "self", no clock model, the true site = the assumed one:
     sqrt(east^2 + radial^2) within 2 % of d, for three phases phi (0, 120, 240 deg), at 30 and
     at 100 days.
     Shown, not checked: the bias in units of the fit's own error, sqrt(b' covEq^-1 b), at
     a = 1e-9, 1e-8 and 1e-7; one noisy run at 1e-7 per span, with its core chi2 and the
     fraction of TOAs flagged, against the same noise without the clock term (does the fit see
     the bias?); the bias when the fit has the crystal clock model ('Clock'), which expects a
     random walk, not a daily line.
Errors at the end if any check fails. (2 Nov 2026, risk register FM-11.)
%}

root = fullfile(fileparts(mfilename('fullpath')), '..');
addpath(fullfile(root, 'functions'));
fails = {};
rs = RandStream('mt19937ar', 'Seed', 47);
t1 = 9555.5 * 86400;                                  % 2026-03-01 00:00 TT
sig = 734e-6; c = 299792458;
b0329 = struct('raDeg', 53.2474, 'decDeg', 54.5787);
psrT = struct('f0', 1.399541538720, 'F1', 0, 'TRef', t1, 'raDeg', b0329.raDeg, 'decDeg', b0329.decDeg);
siteA = makeTrajectory('fixed', 'Lat', 52, 'Lon', 5);
hCr = [1.8e-21 6.492e-22 1.37e-24];                   % makeClock's "crystal"
tSid = 86164.0905; tSol = 86400;
dw = 2*pi * (1/tSid - 1/tSol);                        % rad/s: the sidereal phase gains on the solar one
sinc1 = @(x) sin(x) ./ x;
dt = 8640;                                            % one TOA every 2.4 h
toaAt = @(days) t1 + dt/2 + (0:days * tSol / dt - 1) * dt;
daily = @(a, phi) @(tt) a * tSol / (2*pi) * sin(2*pi * (tt - t1) / tSol + phi);   % eps(t) for rate amplitude a

% 1. FM-11: a daily (solar-day) clock term
a = 1e-8; A = a * tSol / (2*pi);
spans = [30 100]; phis = [0 2*pi/3 4*pi/3];
[ratio, dPred] = deal(zeros(numel(spans), numel(phis)), zeros(1, numel(spans)));
for i = 1:numel(spans)
    t = toaAt(spans(i)); N = numel(t);
    dPred(i) = c * A * sinc1(dw * spans(i) * tSol / 2) / cosd(b0329.decDeg);
    for j = 1:numel(phis)
        r = trueResiduals(t, psrT, siteA, siteA, 'Clock', daily(a, phis(j))).residual;
        f = navFit(struct('t', t, 'r', r, 'sigma', sig * ones(1, N)), b0329, siteA);
        ratio(i, j) = hypot(f.east, f.radial) / dPred(i);
    end
end
ok1 = all(abs(ratio(:) - 1) < 0.02);
fprintf(['1. FM-11 daily clock term, a = 1e-8 (A = %.1f us), noise-free: predicted %.1f km (30 d), %.1f km (100 d); ' ...
         'fit / prediction at phases 0, 120, 240 deg %s (30 d), %s (100 d): %s\n'], 1e6 * A, dPred / 1e3, ...
    mat2str(ratio(1, :), 4), mat2str(ratio(2, :), 4), passStr(ok1));
if ~ok1, fails{end+1} = 'FM-11 daily clock term'; end

% Shown: the bias in units of the fit's error, the fit's own diagnostics, the clock model
amps = [1e-9 1e-8 1e-7];
ws = warning('off', 'navFit:converge');
for i = 1:numel(spans)
    t = toaAt(spans(i)); N = numel(t); s = sig * ones(1, N);
    nSig = zeros(1, numel(amps));
    for k = 1:numel(amps)
        r = trueResiduals(t, psrT, siteA, siteA, 'Clock', daily(amps(k), 0)).residual;
        f = navFit(struct('t', t, 'r', r, 'sigma', s), b0329, siteA);
        b = [f.east; f.radial];
        nSig(k) = sqrt(b' / f.covEq * b);
    end
    r7 = trueResiduals(t, psrT, siteA, siteA, 'Clock', daily(1e-7, 0)).residual;
    n = sig * randn(rs, 1, N);
    fN = navFit(struct('t', t, 'r', r7 + n, 'sigma', s), b0329, siteA);        % noisy, with the clock term
    fR = navFit(struct('t', t, 'r', n, 'sigma', s), b0329, siteA);             % the same noise, without it
    f0 = navFit(struct('t', t, 'r', r7, 'sigma', s), b0329, siteA);            % noise-free, no clock model
    fC0 = navFit(struct('t', t, 'r', r7, 'sigma', s), b0329, siteA, 'Clock', hCr);   % noise-free, crystal model
    fC = navFit(struct('t', t, 'r', r7 + n, 'sigma', s), b0329, siteA, 'Clock', hCr);
    fprintf(['   shown, %d d: bias in units of the fit''s error at a = 1e-9, 1e-8, 1e-7: %s (east error +-%.1f km). ' ...
             'a = 1e-7, one noisy run: bias %.0f km, core chi2 %.3f, flagged %.1f %% (same noise without the term: %.3f, %.1f %%). ' ...
             'With the crystal clock model: bias %.0f km of %.0f km (%.0f %%), core chi2 %.3f, flagged %.1f %%\n'], ...
        spans(i), mat2str(nSig, 3), f0.eastErr / 1e3, hypot(fN.east, fN.radial) / 1e3, fN.chi2Core, 100 * mean(fN.outlier), ...
        fR.chi2Core, 100 * mean(fR.outlier), hypot(fC0.east, fC0.radial) / 1e3, hypot(f0.east, f0.radial) / 1e3, ...
        100 * hypot(fC0.east, fC0.radial) / hypot(f0.east, f0.radial), fC.chi2Core, 100 * mean(fC.outlier));
end
warning(ws);

if isempty(fails)
    fprintf('testNavFitFaults: ALL PASSED\n');
else
    error('testNavFitFaults:failed', 'FAILED: %s', strjoin(fails, ', '));
end
end


% =================================================================================
function s = passStr(p)
if p, s = 'PASS'; else, s = 'FAIL'; end
end
