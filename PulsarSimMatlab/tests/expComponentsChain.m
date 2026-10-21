%{
EXPCOMPONENTSCHAIN - multi-component profiles: the voltage chain against the fast simulator

Experiment (21 Oct 2026, no PASS/FAIL): does simulateFold('Components') give the TOA statistics
the voltage chain gives for a multi-component pulse? So far components were checked against
foldProfile on a synthetic power series (testSimulateFold 14) and in the generator alone
(testGeneratePulsarSignal 7); here the whole chain.

Pulse: B0329+54's five Gaussian components at 1410 MHz (arXiv:1703.03582 Table 2, in turns, as
runScenario) on the pipeline's 10 ms pulsar - a shape test; the narrowest component (0.0055 turns
= 55 us) is 11 phase bins at NBin 2048. Chain (as main, channel path with excision, no RFI): the
generator with 'Components' -> applyDispersionStream at DM 2.6 (J0437's; the forward filter needs
~1.5 GB instead of 12 GB at DM 25, so the experiment runs next to other jobs; the profile comparison
does not depend on DM) -> addNoiseAndRFI (pipelineParams' snrDB, -5 dB) -> applyIQmodulation ->
processIQ -> foldToTOA with the component template -> errors against the generator's reference
times (validateTOA). nSeed runs of L = 0.2 s with their own pulsar and noise seeds, in their own
folder data/compchain (large files deleted on the way).
Simulator: simulateFold with the same components, the chain's channel noise coefficients, Bnoise
and time bins (detectChannels' info of the first run), the true pulse fraction per channel g_c
(forward band-edge taper; as expSimVsChain), NBin and 1-turn sub-ints, observer 'optimal' with
ChannelGain 1 (as main).
Printed side by side: TOAs, median SNR per turn, median sigma, rms error / rms predicted (with its
+-1/sqrt(2n)), mean red. chi^2, fraction of pChi2 < 0.05. Expected: agreement within the chain's
statistics (~90 TOAs: rms/pred +-7 %), as for the single Gaussian (C2 / C4a).
Run from the PulsarSimMatlab folder: run('tests/expComponentsChain.m'). ~45-60 min.
%}

scriptDir = fileparts(mfilename('fullpath'));
root = fileparts(scriptDir);
addpath(fullfile(root, 'functions'));
oldDir = cd(root); restoreDir = onCleanup(@() cd(oldDir));
pipelineParams;
nSeed = 6; Lc = 0.2; DMc = 2.6; nSimRun = 6; nSimSub = 300;   % simulator: 6 x 300 sub-ints (memory)
comps = [0.159 -11.427 1.479; 0.439 -1.123 3.281; 0.585 0 1.055; 0.152 7.719 2.751; 0.170 9.482 0.847];
comps = [comps(:, 1), comps(:, 2) / 360, comps(:, 3) * 2*sqrt(2*log(2)) / 360];
ephC = ephem; ephC.DM = DMc; ephC.profileFWHM = comps;
procC = proc; procC.frontEnd = 'channels';
dDir = fullfile(dataDir, "compchain");
if ~isfolder(dDir), mkdir(dDir); end
pF = struct('chanIQ', fullfile(dDir, "c_chan"), 'chanBlank', fullfile(dDir, "c_chan_blanked"), ...
    'chanDedisp', fullfile(dDir, "c_dedisp_chan"), 'chanPower', fullfile(dDir, "c_dedisp_chan_power.dat"), ...
    'chanDedispB', fullfile(dDir, "c_dedispB_chan"), 'chanPowerB', fullfile(dDir, "c_dedispB_chan_power.dat"), ...
    'chanWeight', fullfile(dDir, "c_dedispB_chan_weights.dat"), 'foldChan', fullfile(dDir, "c_fold_chan.mat"));

% ---- the voltage chain, nSeed runs --------------------------------------------------------
C = struct('err', [], 'sig', [], 'snr', [], 'chi2', [], 'pchi', []);
for s = 1:nSeed
    t0s = tic;
    info_gen = generatePulsarSignal(fullfile(dDir, "c_raw.dat"), T, f_in, A, Lc, dutycycle, ...
        'Seed', 100 + s, 'EnvelopeMode', 'power', 'Components', comps);
    info_disp = applyDispersionStream(info_gen.file, fullfile(dDir, "c_disp.dat"), DMc, ...
        info_gen.actualFsOut, fLow, fHigh, 'RefFreq', refFreq);
    delete(info_gen.file);
    info_rx = addNoiseAndRFI(info_disp.file, fullfile(dDir, "c_rx.dat"), info_disp.actualFsOut, ...
        'Band', [fLow fHigh], 'SNRdB', snrDB, 'SignalInfo', info_gen, 'Seed', 200 + s, 'Verbose', false);
    delete(info_disp.file);
    info_IQ = applyIQmodulation(info_rx.file, fullfile(dDir, "c_iq.dat"), info_rx.actualFsOut, fs, fLO, ...
        'FilterOrder', filterOrder, 'Band', [fLow fHigh], 'Verbose', false);
    delete(info_rx.file);
    out = processIQ(info_IQ, rx, ephC, procC, pF, 'Verbose', false);
    delete(info_IQ.file);
    if s == 1, det1 = out.info_det; end
    toa = foldToTOA(out.fold, out.info_fold, comps, out.Bnoise, 'Weighting', out.weighting, 'Verbose', false);
    val = validateTOA(toa, info_gen, 'Plot', false, 'Label', sprintf('chain seed %d', s));
    v = toa.valid;
    C.err = [C.err; val.err(:)]; C.sig = [C.sig; toa.toaErr(v)];
    C.snr = [C.snr; toa.snr(v)]; C.chi2 = [C.chi2; toa.redChi2(v)]; C.pchi = [C.pchi; toa.pChi2(v)];
    fprintf('expComponentsChain: chain run %d of %d, %d TOAs, %.0f s\n', s, nSeed, nnz(v), toc(t0s));
end

% ---- the fast simulator at the chain's channel settings -----------------------------------
ncC = [det1.noise.V, det1.noise.X]; BnC = det1.noise.Bnoise; dtC = det1.binDt;
nCh = round((fHigh - fLow) / chanWidth);
gc  = chanGain(fLow, fHigh, chanWidth, nCh);
S = struct('err', [], 'sig', [], 'snr', [], 'chi2', [], 'pchi', []);
for k = 1:nSimRun
    [iS, fS, tS] = simulateFold(ephem.f0, snrDB, nSimSub, 'NChan', nCh, 'ChannelGain', gc, 'BinDt', dtC, ...
        'Bnoise', BnC, 'NoiseCoeffs', ncC, 'Components', comps, 'NBin', procC.nBin, 'Seed', 30 + k, ...
        'Verbose', false);
    tt = foldToTOA(fS, iS, comps, BnC, 'Weighting', 'optimal', 'Verbose', false);
    v  = tt.valid;
    S.err = [S.err; tt.toa(v) - tS.pulseCenters(v)]; S.sig = [S.sig; tt.toaErr(v)];
    S.snr = [S.snr; tt.snr(v)]; S.chi2 = [S.chi2; tt.redChi2(v)]; S.pchi = [S.pchi; tt.pChi2(v)];
end

% ---- side by side -----------------------------------------------------------------------------
fprintf(['\nexpComponentsChain: B0329+54 five components on the %g ms pulsar, %g dB, DM %g, NBin %d, ' ...
         '1 turn per sub-int\n'], T*1e3, snrDB, DMc, procC.nBin);
rows = {'voltage chain', C; 'simulateFold', S};
for r = 1:2
    X = rows{r, 2}; n = numel(X.err);
    rp = sqrt(mean(X.sig.^2));
    fprintf(['  %-14s %5d TOAs: median SNR %.1f, median sigma %.4g us, rms/pred %.3f (+- %.3f), ' ...
             'mean red. chi2 %.4f, pChi2 < 0.05 %.3f\n'], rows{r, 1}, n, median(X.snr), ...
        median(X.sig)*1e6, sqrt(mean(X.err.^2)) / rp, 1/sqrt(2*n), mean(X.chi2), mean(X.pchi < 0.05));
end


function g = chanGain(fLow, fHigh, chanWidth, nCh)
%CHANGAIN  Fraction of the pulsar's power per channel relative to the channel's noise (as
% expSimVsChain): the forward 8 MHz band-edge taper Wf weighted by the channel's own
% dedispersion taper Wi_c: g_c = int Wf^2 Wi_c^2 df / int Wi_c^2 df.
eF = 0.02 * (fHigh - fLow); eC = 0.02 * chanWidth;
g = zeros(1, nCh);
for c = 1:nCh
    lo = fLow + (c - 1) * chanWidth; hi = lo + chanWidth;
    f  = linspace(lo, hi, 20001);
    Wi = ones(size(f));
    a = f < lo + eC; Wi(a) = sin(pi/2 * (f(a) - lo) / eC).^2;
    b = f > hi - eC; Wi(b) = sin(pi/2 * (hi - f(b)) / eC).^2;
    Wf = ones(size(f));
    a = f < fLow + eF; Wf(a) = sin(pi/2 * (f(a) - fLow) / eF).^2;
    b = f > fHigh - eF; Wf(b) = sin(pi/2 * (fHigh - f(b)) / eF).^2;
    g(c) = trapz(f, Wf.^2 .* Wi.^2) / trapz(f, Wi.^2);
end
end
