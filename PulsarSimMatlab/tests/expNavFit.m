%{
EXPNAVFIT - Block D's first deliverable: B0329+54's site offset from simulated TOAs (block D4a)

Experiment, not a unit test. Loads the 30-day motion runs of runScenario (motion.on: the true site
100 km east of the assumed one at 52 N 5 E, a perfect clock; nominal and worst-case flux) and fits
the good TOAs with navFit:
  "self"        east and radial (equatorial) components plus a constant - what the observer gets
                without a calibrated start;
  "calibrated"  the full 3-D position (the constant known: here 0, as the simulator has no clock
                or pulsar offset).
Prints the estimates with their errors against the truth, the outliers found, and two checks of
the input: the fit to the TRUE residuals (no noise) gives the offset itself, and the fit's
'Predict' error equals the one of the fit.

Needs data/scenario/scenario_B0329_<flux>_clean_30d_motion_perfect.mat (runScenario with motion.on
= true, pulsars "B0329"). Run from the PulsarSimMatlab folder: run('tests/expNavFit.m'). ~1 s.
(22 Oct 2026, block D4a.)
%}

scriptDir = fileparts(mfilename('fullpath'));
addpath(fullfile(scriptDir, '..', 'functions'));
dataDir = fullfile(scriptDir, '..', 'data', 'scenario');
for fl = ["nominal", "worst"]
    fIn = fullfile(dataDir, sprintf('scenario_B0329_%s_clean_30d_motion_perfect.mat', fl));
    R = load(fIn).R;
    G = R.groups([R.groups.good]);
    M = R.motionTruth;
    toa = struct('t', [G.tAbs], 'r', [G.resid], 'sigma', [G.sigma]);
    psr = struct('raDeg', M.psr.raDeg, 'decDeg', M.psr.decDeg);
    siteA = M.trajAssumed; siteT = M.trajTrue;
    dx = siteT.ecef - siteA.ecef;
    lon = siteA.params.Lon;
    truth = [[-sind(lon) cosd(lon) 0] * dx, [cosd(lon) sind(lon) 0] * dx, dx(3)];   % east, radial, polar
    fS = navFit(toa, psr, siteA);
    fC = navFit(toa, psr, siteA, 'Start', "calibrated");
    fT = navFit(setfield(toa, 'r', [G.trueResid]), psr, siteA, 'Robust', false); %#ok<SFLD>
    fP = navFit(rmfield(toa, 'r'), psr, siteA, 'Predict', true);
    fprintf('\nexpNavFit: B0329+54, %s flux, 30 days, %d good TOAs (sigma median %.0f us); truth 100 km east\n', ...
        fl, numel(G), 1e6 * median(toa.sigma));
    fprintf('  "self"       east %6.1f +- %4.1f km (%+.1f sigma), radial %+6.1f +- %4.1f km (%+.1f sigma); outliers %d; chi2 core %.3f\n', ...
        fS.east / 1e3, fS.eastErr / 1e3, (fS.east - truth(1)) / fS.eastErr, fS.radial / 1e3, fS.radialErr / 1e3, ...
        (fS.radial - truth(2)) / fS.radialErr, nnz(fS.outlier), fS.chi2Core);
    fprintf('  "calibrated" ENU %s km +- %s km; polar %+.1f +- %.1f km (%+.1f sigma)\n', ...
        mat2str(round(fC.enu.' / 1e3, 1)), mat2str(round(sqrt(diag(fC.covENU)).' / 1e3, 1)), fC.polar / 1e3, ...
        fC.polarErr / 1e3, (fC.polar - truth(3)) / fC.polarErr);
    fprintf('  checks: the true residuals give east %.3f km, radial %+.3f km; Predict error %.1f km (fit %.1f km)\n', ...
        fT.east / 1e3, fT.radial / 1e3, fP.eastErr / 1e3, fS.eastErr / 1e3);
    fprintf('  detection of the 100 km offset: %.1f sigma\n', fS.east / fS.eastErr);
end
