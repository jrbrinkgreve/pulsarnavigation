%{
EXPSIMRFI - block C4b: the fast simulator with RFI profiles against the voltage chain

Experiment (no PASS/FAIL). simulateFold with 'RFI' (profiles from runRFIProfile) at -5 dB,
128 channels, the true pulse fraction per channel (forward band-edge taper, item 3), the
observer fitting with 'optimal' weighting and ChannelGain 1 (as main):
  1. 0.1 s profiles (runRFITest's exact settings, beam passage inside the run) against
     runRFITest (12 Oct, data/rfi/runRFITest_results.mat; one noise seed, 8 TOAs per case,
     with excision): median SNR, median TOA error, median red. chi^2.
  2. 1 s profiles (steady part, beam passages outside the run): the production numbers
     over nSeed x 8 TOAs - SNR and error relative to clean, error-bar honesty (rms/pred),
     red. chi^2 and pChi2 - for noise only, the realistic scenario and the narrow-Inmarsat
     worst case. (The pulsar is resampled every run, so rms/pred is ~1 for honest error
     bars; runMonteCarlo's 0.956 kept the pulsar fixed.)
Needs data/rfi/profile/rfiProfile_*.mat (runRFIProfile) and runRFITest's results. ~2 min.
%}

root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root, 'functions'));
oldDir = cd(root); restoreDir = onCleanup(@() cd(oldDir));
pipelineParams;
nSeed = 50;
nCh = round((fHigh - fLow) / chanWidth);
gc  = chanGain(fLow, fHigh, chanWidth, nCh);
names  = ["noiseOnly", "realistic", "inmarsatNarrow"];
labels = ["noise only", "realistic scenario (B5e)", "realistic, Inmarsat narrow"];

% ---- 1. 0.1 s profiles vs runRFITest -------------------------------------------------------
R = load(fullfile(dataDir, 'rfi', 'runRFITest_results.mat'));
fprintf('expSimRFI 1: 0.1 s profiles vs runRFITest (-5 dB, excision; chain: one seed, 8 TOAs)\n');
fprintf('  %-28s | SNR chain / sim | median err [us] chain / sim | red. chi2 chain / sim\n', 'case');
for i = 1:numel(names)
    f = fullfile(dataDir, 'rfi', 'profile', sprintf('rfiProfile_%s_0.1s.mat', names(i)));
    if ~isfile(f), fprintf('  %s: no profile\n', names(i)); continue; end
    P = load(f);
    st = simStats(P.prof, gc, ephem, 20, 100*i);
    k = find(strcmp({R.res.label}, labels(i)) & [R.res.excision] == 1, 1);
    fprintf('  %-28s | %5.1f / %5.1f   | %6.3f / %6.3f              | %6.3f / %6.3f\n', labels(i), ...
        R.res(k).snr, st.snr, R.res(k).err, st.err*1e6, R.res(k).redChi2, st.chi2);
end

% ---- 2. 1 s profiles: production numbers ------------------------------------------------------
fprintf('\nexpSimRFI 2: 1 s profiles (steady part), -5 dB, %d x 8 TOAs each\n', nSeed);
fprintf(['  %-28s | SNR (/clean) | median err [us] (/clean) | rms/pred (+-%.3f) | median red. chi2 | ' ...
         'pChi2 mean, < 0.05 | channels: level > 1.05, noise ratio > 1.02, min valid\n'], 'case', ...
    1/sqrt(2*8*nSeed));
ref = [];
for i = 1:numel(names)
    f = fullfile(dataDir, 'rfi', 'profile', sprintf('rfiProfile_%s_1s.mat', names(i)));
    if ~isfile(f), fprintf('  %s: no profile\n', names(i)); continue; end
    P = load(f);
    st = simStats(P.prof, gc, ephem, nSeed, 1000*i);
    if isempty(ref), ref = st; end
    nr = (P.prof.noiseCoeffs(:, 1) + 2*sum(P.prof.noiseCoeffs(:, 2:end), 2)) / ...
         (P.prof.modelCoeffs(1) + 2*sum(P.prof.modelCoeffs(2:end)));
    fprintf('  %-28s | %5.1f (%.3f) | %6.3f (%.3f)           | %6.3f            | %7.4f         | %.3f, %.3f        | %d, %d, %.4f\n', ...
        labels(i), st.snr, st.snr/ref.snr, st.err*1e6, st.err/ref.err, st.ratio, st.chi2, st.pMean, ...
        st.pLow, nnz(P.prof.level > 1.05), nnz(nr > 1.02), min(P.prof.validFraction));
end


function st = simStats(prof, gc, ephem, nSeed, seed0)
%SIMSTATS  nSeed x 8 one-turn sub-ints at -5 dB with the profile; pooled statistics.
sn = []; er = []; rc = []; pc = []; zz = [];
for sd = 1:nSeed
    [ic, fc, tc] = simulateFold(ephem.f0, -5, 8, 'NChan', prof.nChan, 'ChannelGain', gc, ...
        'RFI', prof, 'FWHM', ephem.profileFWHM, 'Seed', seed0 + sd, 'Verbose', false);
    tt = foldToTOA(fc, ic, ephem.profileFWHM, prof.BnoiseChan, 'Weighting', 'optimal', 'Verbose', false);
    v  = tt.valid;
    sn = [sn; tt.snr(v)]; er = [er; tt.toaErr(v)]; rc = [rc; tt.redChi2(v)]; pc = [pc; tt.pChi2(v)]; %#ok<AGROW>
    zz = [zz; (tt.toa(v) - tc.pulseCenters(v)) ./ tt.toaErr(v)];                                    %#ok<AGROW>
end
st = struct('snr', median(sn), 'err', median(er), 'chi2', median(rc), 'ratio', sqrt(mean(zz.^2)), ...
    'pMean', mean(pc), 'pLow', mean(pc < 0.05), 'n', numel(zz));
end

function g = chanGain(fLow, fHigh, chanWidth, nCh)
%CHANGAIN  Pulse fraction per channel from the forward band-edge taper (as expSimVsChain).
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
