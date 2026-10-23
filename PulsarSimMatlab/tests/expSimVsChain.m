%{
EXPSIMVSCHAIN - block C2: the fast simulator (simulateFold) against the voltage chain

Experiment (no PASS/FAIL): simulated folds at the settings of the saved full-band
voltage-chain sweeps, processed by foldToTOA exactly as runSNRSweep does, with the
statistics side by side.

References (voltage chain, full band, L = 0.1 s, 1 turn per sub-int):
  data/mc/snrSweep_phaseDE.mat        -25 ... +20 dB, 10 realizations (90 TOAs per point)
  data/mc/snrSweep_threshold_6oct.mat -20 ... -15 dB in 1 dB steps, 20 realizations (180)
  each with its H0 pass (zero sky + noise); runH0 full band (6 Oct, notes §7): T0 std
  1.008 +- 0.024, normProfile lag-1 correlation 0.2484 (model 0.2501).
The simulator gets the chain's signal gain g (expectedPowerModel on main's full-band files:
peak power above baseline = g*rho, the forward taper costs a little) and the observer's
Bnoise (noiseBandwidth of the dedispersion taper), and the sweeps' NBin (read from their
normProfile). Per SNR point nSim one-turn sub-ints (true offset 0, as the chain).

Expected: agreement in the weak regime (rho << 1), where the observer's radiometer model is
the whole story. At high rho the chain's noise has the exact taper coefficients (css, csn,
cnn of expectedPowerModel) and the pulsar's own noise, while the simulator draws exactly
the observer's model (m^2/(Bnoise*binDt)): there its rms/pred is 1 by construction, the
chain's 1.04-1.06.

Channel path (C4a, 20 Oct): 128 simulated channels with the chain's noise coefficients
[V X(1..7)], BnoiseChan and 0.96 us bins (detectChannels' info.noise of main's channel run)
and the true pulse fraction per channel g_c = int Wf^2 Wi_c^2 / int Wi_c^2 (Wf the 8 MHz
forward band-edge taper of the simulated pulsar, Wi_c channel c's dedispersion taper;
item 3), observer fit 'optimal' with ChannelGain 1 (as main), against the chain's
channel-path Monte Carlo (runMonteCarlo, 19 Oct): -5 dB 50 seeds with excision (rms/pred
0.956 +- 0.034, rms pred 2.950 us, mean red. chi2 1.0063, pChi2 mean 0.292, fraction < 0.05
0.205), -15 dB 10 seeds (rms/pred 0.971, mean red. chi2 0.9987, pChi2 mean 0.537); also
with g_c = 1 (no forward taper) for contrast. The simulator has no excision (0.0124 % of
the samples blanked at -5 dB, mostly pulse peaks: seed 43 2.943 -> 2.956 us).

Needs main's full-band data files (data/test_IQ_dedispersed_info.mat etc.: main.m with
frontEnd = 'fullband' once), main's channel power file info (data/chan/
test_dedispB_chan_power_info.mat) and the two sweep files. ~2-3 min.
%}

root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root, 'functions'));
clear restoreDir                                    % another script's (shared workspace): it would cd back after cd(root)
oldDir = cd(root); restoreDir = onCleanup(@() cd(oldDir));
pipelineParams;
nSim = 2000;                                        % simulated sub-ints per SNR point

% ---- Chain's signal gain and the observer's Bnoise ------------------------------------------
info_gen    = loadInfo(fileRaw);
info_disp   = loadInfo(fileDispersed);
info_rx     = loadInfo(fileRx);
info_IQ     = loadInfo(fileIQ);
info_dedisp = loadInfo(fileDedisp);
M   = expectedPowerModel(info_gen, info_disp, info_IQ, info_dedisp, info_rx);
rho = 10^(info_rx.snrDB / 10);
g   = M.sigScale / (M.Pn * rho);                     % peak power above baseline / rho
Bn  = noiseBandwidth(fLow, fHigh, info_dedisp.edgeWidth);
dt  = 1 / f_out;
f0  = ephem.f0;
fprintf(['expSimVsChain: chain at %g dB: signal gain g = %.4f, Bnoise %.2f MHz (simulator ' ...
         'variance 1/(Bn*dt); chain cnn %.4g, css %.4g per MHz); optimum SNR per turn %.2f\n'], ...
    info_rx.snrDB, g, Bn/1e6, M.cnn*1e6, M.css*1e6, M.snrPulse);

% ---- Sweep points: simulated vs chain --------------------------------------------------------
files = ["data/mc/snrSweep_phaseDE.mat", "data/mc/snrSweep_threshold_6oct.mat"];
Q = @(x) 0.5 * erfc(x / sqrt(2));
for fi = 1:numel(files)
    S  = load(files(fi));
    r  = S.res;
    Nb = size(S.detC{1, 1}.normProfile, 1);         % NBin of that sweep
    fprintf('\n%s (%d realizations, NBin %d): chain | simulated (%d sub-ints per point)\n', ...
        files(fi), S.nReal, Nb, nSim);
    fprintf(['  snrDB | median sigma [us]   | rms/pred        | |z|>5           | P_D known       ' ...
             '| P_D unknown     | good fraction\n']);
    for s = 1:numel(r.snrDB)
        [i1, f1, t1] = simulateFold(f0, r.snrDB(s), nSim, 'NBin', Nb, 'SignalGain', g, ...
            'Bnoise', Bn, 'BinDt', dt, 'FWHM', ephem.profileFWHM, 'Seed', 1000*fi + s, ...
            'Verbose', false);
        [toa, det] = foldToTOA(f1, i1, ephem.profileFWHM, Bn, 'Verbose', false);
        v   = toa.valid;
        err = toa.toa(v) - t1.pulseCenters(v);
        err = mod(err + 0.5/f0, 1/f0) - 0.5/f0;      % nearest true pulse
        z   = err ./ toa.toaErr(v);
        sim = [median(toa.toaErr(v))*1e6, sqrt(mean(err.^2)) / sqrt(mean(toa.toaErr(v).^2)), ...
               mean(abs(z) > 5), mean(det.detectedKnown(det.tested)), ...
               mean(det.detectedUnknown(det.tested)), mean(toa.good)];
        % chain: median sigma from its TOAs
        sc  = [S.toaC{:, s}];
        sc  = vertcat(sc.toaErr); sc = median(sc(isfinite(sc)))*1e6;
        ch  = [sc, r.ratio(s), r.fracOut(s), r.pdKnown(s), r.pdUnknown(s), r.nGood(s)/r.n(s)];
        eB  = @(p, n) sqrt(max(p.*(1-p), 1/n)/n);   % binomial 1-sigma (chain's n)
        fprintf(['  %+5g | %8.4g %8.4g | %6.3f %6.3f   | %5.3f %5.3f     | %5.3f %5.3f +-%4.3f ' ...
                 '| %5.3f %5.3f +-%4.3f | %5.3f %5.3f\n'], r.snrDB(s), ch(1), sim(1), ch(2), sim(2), ...
            ch(3), sim(3), ch(4), sim(4), eB(ch(4), r.nTested(s)), ch(5), sim(5), ...
            eB(ch(5), r.nTested(s)), ch(6), sim(6));
    end
    fprintf('  H0 (chain, %d sub-ints): T0 mean %+.3f, std %.3f\n', r.h0.n, r.h0.T0mean, r.h0.T0std);
end

% ---- H0: T0 and the normalized-profile correlation -----------------------------------------
[i0, f0f] = simulateFold(f0, -Inf, 5000, 'NBin', nBin, 'Bnoise', Bn, 'BinDt', dt, 'Seed', 7, ...
    'Verbose', false);
d0 = detectPulsar(f0f, i0, gaussianTemplate(nBin, ephem.profileFWHM), 'Bnoise', Bn, 'Verbose', false);
T0 = d0.T0(d0.tested); zp = d0.normProfile(:, d0.tested);
r1 = mean(zp .* circshift(zp, -1), 'all') / mean(zp.^2, 'all');
fprintf(['\nH0, simulated (%d sub-ints, NBin %d): T0 mean %+.4f, std %.4f (+- %.4f); normProfile ' ...
         'lag-1 %.4f; false alarms %d known / %d unknown (P_FA 1e-3: %.1f)\n'], numel(T0), nBin, ...
    mean(T0), std(T0), 1/sqrt(2*numel(T0)), r1, nnz(d0.detectedKnown), nnz(d0.detectedUnknown), ...
    1e-3*numel(T0));
fprintf('   chain (runH0 full band, 6 Oct): T0 std 1.008 +- 0.024, lag-1 0.2484 (model 0.2501)\n');

% ---- -5 dB with main's settings: asymptotic and copied weights -----------------------------
SF = load(fileFold, 'fold', 'info');                % main's full-band fold (NBin of main)
cS = find(all(SF.fold.weight > 0, 1));              % complete sub-ints (first/last are partial)
wC = struct('weight', SF.fold.weight(:, cS), 'weight2', SF.fold.weight2(:, cS), ...
    'weightX', SF.fold.weightX(:, cS));
nC = numel(cS);
sg = zeros(2, 1); snr = zeros(2, 1);
for w = 1:2
    if w == 1, wArg = {}; else, wArg = {'Weights', wC}; end
    sgW = []; snW = [];
    for rep = 1:200                                 % copied weights: nC sub-ints per draw
        [ia, fa] = simulateFold(f0, info_rx.snrDB, nC, 'NBin', SF.info.NBin, 'SignalGain', g, ...
            'Bnoise', Bn, 'BinDt', dt, 'FWHM', ephem.profileFWHM, 'Seed', 5000 + rep, ...
            wArg{:}, 'Verbose', false);
        ta = foldToTOA(fa, ia, ephem.profileFWHM, Bn, 'Verbose', false);
        sgW = [sgW; ta.toaErr(ta.valid)]; snW = [snW; ta.snr(ta.valid)]; %#ok<AGROW>
    end
    sg(w) = median(sgW) * 1e6; snr(w) = median(snW);
end
fprintf(['\n%g dB, NBin %d, 1 turn per sub-int: median sigma / SNR: simulated asymptotic %.4f us / ' ...
         '%.2f, copied weights (main''s fold) %.4f us / %.2f; chain (main, seed 43) 2.880 us / 93.4; ' ...
         'optimum (expectedPowerModel) %.4f us / %.2f\n'], info_rx.snrDB, SF.info.NBin, sg(1), snr(1), ...
    sg(2), snr(2), M.toaErrPulse*1e6, M.snrPulse);

% ---- Channel path (C4a): 128 channels, forward taper, observer ChannelGain 1 ---------------
detC = loadInfo(fileChanPowerB);                    % main's channel run (noise of the time bins)
ncC  = [detC.noise.V, detC.noise.X];
BnC  = detC.noise.Bnoise;
dtC  = detC.binDt;
nCh  = round((fHigh - fLow) / chanWidth);
gc   = chanGain(fLow, fHigh, chanWidth, nCh);       % true pulse fraction per channel
fprintf(['\nChannel path: %d channels, V %.4f, X(1) %.4f (Lmax %d), BnoiseChan %.4f MHz, binDt %.3f us; ' ...
         'g_c outer three channels %s / %s, mean %.4f\n'], nCh, ncC(1), ncC(2), numel(ncC) - 1, ...
    BnC/1e6, dtC*1e6, mat2str(round(gc(1:3), 3)), mat2str(round(gc(end-2:end), 3)), mean(gc));
cases = {-5, gc, 50, 'g_c (forward taper)'; -5, ones(1, nCh), 50, 'g_c = 1'; ...
         -15, gc, 50, 'g_c (forward taper)'};
for ci = 1:size(cases, 1)
    [snrC, gC, nSeed, lab] = cases{ci, :};
    err = []; sgv = []; snv = []; rc = []; pc = [];
    for sd = 1:nSeed
        [ic, fc, tc] = simulateFold(f0, snrC, 8, 'NChan', nCh, 'ChannelGain', gC, 'BinDt', dtC, ...
            'Bnoise', BnC, 'NoiseCoeffs', ncC, 'FWHM', ephem.profileFWHM, 'Seed', 7000 + 100*ci + sd, ...
            'Verbose', false);
        tt = foldToTOA(fc, ic, ephem.profileFWHM, BnC, 'Weighting', 'optimal', 'Verbose', false);
        v  = tt.valid;
        err = [err; tt.toa(v) - tc.pulseCenters(v)];                    %#ok<AGROW>
        sgv = [sgv; tt.toaErr(v)]; snv = [snv; tt.snr(v)];             %#ok<AGROW>
        rc  = [rc; tt.redChi2(v)]; pc = [pc; tt.pChi2(v)];             %#ok<AGROW>
    end
    nE = numel(err);
    fprintf(['  %+g dB, %s, %d TOAs: rms pred %.4g us, median sigma %.4g us, median SNR %.1f, rms/pred ' ...
             '%.3f (+- %.3f), mean red. chi2 %.4f, pChi2 mean %.3f, fraction < 0.05 %.3f\n'], snrC, lab, ...
        nE, sqrt(mean(sgv.^2))*1e6, median(sgv)*1e6, median(snv), sqrt(mean(err.^2))/sqrt(mean(sgv.^2)), ...
        1/sqrt(2*nE), mean(rc), mean(pc), mean(pc < 0.05));
end
fprintf(['  chain (runMonteCarlo, channels, excision): -5 dB 400 TOAs rms pred 2.950 us, rms/pred 0.956, ' ...
         'mean red. chi2 1.0063, pChi2 mean 0.292, < 0.05 0.205; -15 dB 80 TOAs rms pred 26.36 us, ' ...
         'rms/pred 0.971, mean red. chi2 0.9987, pChi2 mean 0.537\n']);


function g = chanGain(fLow, fHigh, chanWidth, nCh)
%CHANGAIN  Fraction of the pulsar's power per channel relative to the channel's noise: the
% forward 8 MHz band-edge taper Wf (applyDispersionStream, EdgeFrac 0.02, amplitude sin^2)
% weighted by the channel's own dedispersion taper Wi_c (dedisperseChannels, EdgeFrac 0.02 of
% the channel, centre +- chanWidth/2): g_c = int Wf^2 Wi_c^2 df / int Wi_c^2 df.
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
