%{
RUNSCENARIO - TOA budget of the reference scenario with the fast simulator (block C)

How many hours of data does one TOA take, how many good TOAs come per day, how precise are they
and are their error bars honest - for the 3x3 array (dual polarization, 1.2-1.6 GHz, 128
channels) observing B0329+54 from the Netherlands (52 N) and J0437-4715 from 30 S (a ship that
sailed south), each at the nominal and the worst-case flux. Folds are drawn by simulateFold;
the observer forms TOAs with adaptiveSubints (merge short sub-ints until SNR ~7) and foldToTOA.

Assumptions (notes §11, simulator.html):
  - passes: elevation >= 20 deg, one per sidereal day (B0329+54 from 52 N: 19.6 h, never below
    16.6 deg; J0437 from 30 S: 12.3 h); constant gain over a pass (steered array); site fixed;
  - without motion (motion.on false, the default): TOA budgets, not navigation; the true pulse
    sits 0.1 FWHM from the predicted phase (a small ephemeris/position error) - adaptiveSubints
    closes groups on T0 at the predicted phase, which assumes the prediction is well inside the
    pulse (an error of a sizeable fraction of the FWHM would make the groups longer than needed,
    not the TOAs wrong);
  - with motion (motion.on true, block D3c): the truth of a receiver that is not where the
    observer thinks. The observer folds with the assumed site (52 N 5 E for B0329+54) and its
    clock estimate; the true site is offset by motion.offsetENU (default 100 km east), the true
    clock is perfect or an SDR's ('sdr': the paper's crystal noise, makeClock, and a rate drawn
    within +-2 ppm; the observer knows the rate to clockRateErr = 1e-9 and the offset to
    clockOffsetErr = 1 ms, a calibration in port). Each sub-int's pulse sits at trueResiduals'
    phase offset at its reference time. The scenario's t = 0 is set so that its first transit
    (t = tSid/2) is a real transit at the site after motion.epoch (2026-03-01 00:00 TT). Each TOA
    records its absolute time, the measured residual (TOA - tRef) and the true residual with its
    geometric and clock parts: the input of the navigation fit (D4). The measured residual is in
    the simulator's frame (f0, not the Doppler-shifted frequency): 1e-4 relative, 19 ns on the
    193 us of a 100 km offset. A warning if the true offset drifts by more than 5 % of the FWHM
    within a sub-int (one offset per sub-int would then be too coarse);
  - with an SDR clock the observer tracks it (motion.track "daily", block D4b-2): after each day
    it fits all its TOAs so far with navFit and the clock model (the residuals with its own clock
    estimate of the time added back), predicts its clock for the next day (constant + clockAt) and
    folds that day with the prediction. Until 20 TOAs (~2 days) it keeps the calibration in port
    (a check with 10: day 2 extrapolated a rate fitted from one day and was 1.3 ms off). "none": the
    calibration only (the clock then drifts away by milliseconds within weeks). Groups that span a
    day boundary mix two estimates (their difference is the one-day prediction error, ~0.5 ms);
  - rho in the simulator's convention (peak of the Gaussian pulse / receiver noise per
    polarization, rho = S_peak/SEFD), 'NPol' 2; the observer passes Bnoise x 2;
  - pulse and template: B0329+54's five Gaussian components at 1410 MHz (since 20 Oct; a single
    Gaussian with W50 overstated its peak by 3 dB), J0437 eight components fitted to its measured
    profile (since 21 Oct; the single Gaussian with W50 left out its wings: W_eq 2.4x, ~5x the data
    per TOA);
    rho from S1400, the SEFD and the pulse area of the profile;
  - 128 channels with the forward band-edge taper on the truth side (item 3); channel noise
    settings (Bnoise, binDt, clean noise coefficients) of the 1 s noise-only RFI profile;
  - variants: 'clean' (receiver noise only); 'rfi' (realistic RFI profile of the 1 s chain run,
    stationary; the observer's flat-channel noise model, as the pipeline by default); 'rfiCal'
    (the same, the observer calibrates each channel's noise: 'ObserverNoise' 'measured', the
    pipeline's proc.noiseCal); 'rfiCalScint' (+ scintillation, adaptive groups); 'rfiCalScintFixed'
    (fixed groups of the median length of rfiCalScint, for comparison: what the adaptive length buys).
Per day, the pass is cut into chunks (memory) on sub-int boundaries, so the sub-int grid is
the one a single call would give. The sub-ints of the group still open at the end of a chunk
are carried into the next chunk (also over the gap to the next day: 'AcrossGaps'), which is
exact because fold sums add. Each closed group gives one TOA. What is still open after the
last day is counted as leftover.

Outputs per case (printed table; data/scenario/scenario_<pulsar>_<flux>_<variant>_<n>d.mat, with
motion ..._<n>d_motion_<clock>.mat, and ..._sdr_notrack.mat without tracking):
good TOAs per day, hours of data per TOA (median, 10-90 %), calendar span per TOA, median
sigma (us and c*sigma in km), rms of (TOA - truth)/sigma and the fraction beyond 5, mean
reduced chi^2 and pChi2 < 0.05, leftover fraction, figure of merit sigma*sqrt(hours) [us
sqrt(h)]. (20 Oct 2026, block C; design agreed with Jasper the same day. Motion: 21 Oct 2026,
block D3c, choices agreed with Jasper: perfect clock first, rate known to 1e-9, 1 March.)
%}

pulsars  = ["B0329", "J0437"];
fluxes   = ["nominal", "worst"];
variants = "clean";                  % "clean", "rfi", "rfiCal", "rfiCalScint", "rfiCalScintFixed"
                                     % (the fixed variant needs the rfiCalScint result first)
nDays    = 30;                       % simulated days
overwrite = false;                   % false: skip cases whose result file exists
% the truth of a moving receiver (block D3c; on = false: the TOA budgets, unchanged)
motion = struct('on', false, ...
    'epoch', 9555.5 * 86400, ...     % 2026-03-01 00:00 TT [s since J2000]; the first transit after it
    'offsetENU', [100e3 0 0], ...    % true site minus assumed site, east-north-up [m]
    'clock', "perfect", ...          % "perfect" | "sdr" (crystal noise + a rate within +-2 ppm)
    'clockOffsetErr', 1e-3, ...      % "sdr": the observer's estimate is off by this offset [s]
    'clockRateErr', 1e-9, ...        %   ... and this rate [s/s]
    'track', "daily");               % "sdr": "daily" (refit and predict the clock each day) | "none"

scriptDir = fileparts(mfilename('fullpath'));
addpath(fullfile(scriptDir, 'functions'));
pipelineParams;
outDir = fullfile(dataDir, "scenario");
if ~isfolder(outDir), mkdir(outDir); end
P0 = load(fullfile(dataDir, "rfi", "profile", "rfiProfile_noiseOnly_1s.mat")).prof;   % channel settings
PR = load(fullfile(dataDir, "rfi", "profile", "rfiProfile_realistic_1s.mat")).prof;   % realistic RFI
nCh = P0.nChan; BnC = P0.BnoiseChan; dtC = P0.binDt; ncC = P0.modelCoeffs;
gc  = chanGain(fLow, fHigh, chanWidth, nCh);
nPol = 2; elMin = 20; tSid = 86164.0905;

% pulsar settings (ATNF; scintillation: notes §11). Profile: B0329+54 five Gaussian components at
% 1410 MHz (arXiv:1703.03582 Table 2: centre and sigma in degrees, amplitude) -> [amplitude, centre,
% FWHM] in turns; J0437 eight Gaussian components fitted (21 Oct 2026) to the Parkes 1369 MHz profile
% of Dai et al. 2015 (MNRAS 449, 3223; EPN database dhm+15, CC BY 4.0; phase 0 = the peak): W50, W10,
% W_eq, spectral width and timing merit within 1 % of the measured profile, residual 1.6x its noise
% (its core is sharper than a Gaussian: NBin 1024 for estimateTOA's resolution rule). rho (simulator
% convention: pulse peak / receiver noise per polarization) = S1400 / (SEFD * W_eq), W_eq the pulse
% area of the unit-peak profile in turns; SEFD 3.2e5 Jy nominal, 7.4e5 Jy worst case (notes §11).
cB = [0.159 -11.427 1.479; 0.439 -1.123 3.281; 0.585 0 1.055; 0.152 7.719 2.751; 0.170 9.482 0.847];
cB = [cB(:, 1), cB(:, 2) / 360, cB(:, 3) * 2*sqrt(2*log(2)) / 360];
cJ = [0.045599 -0.205164 0.029901; 0.043891 -0.113820 0.028824; 0.103938 -0.038174 0.261016;
      0.215798 -0.000947 0.053271; 0.462737 -0.000766 0.010044; 0.262274  0.007027 0.027132;
      0.310520  0.009091 0.008531; 0.028718  0.147582 0.073914];   % [amplitude, centre, FWHM] (turns)
sefd = struct('nominal', 3.2e5, 'worst', 7.4e5);
psr.B0329 = struct('name', 'B0329+54', 'f0', 1/0.714519699726, 'ra', 53.2474, 'dec', 54.5787, 'lat', 52.0, 'lon', 5.0, ...
    'fwhm', 0.0092, 'comps', cB, 'S1400', 0.203, 'nBin', 2048, ...   % narrowest component 11 bins
    'scint', struct('bandwidth', 7e6, 'timescale', 780, 'band', [fLow fHigh]), ...   % (1024: 5.6 bins,
    'subintTime', 120, 'chunkHours', 2);                       % error bars 18 % small at SNR 7)
psr.J0437 = struct('name', 'J0437-4715', 'f0', 173.6879458, 'ra', 69.3162, 'dec', -47.2525, 'lat', -30.0, 'lon', 5.0, ...
    'fwhm', 0.0245, 'comps', cJ, 'S1400', 0.150, 'nBin', 1024, ...   % fwhm: W50 (truth offset)
    'scint', struct('bandwidth', 740e6, 'timescale', 2290, 'band', [fLow fHigh]), ...
    'subintTime', 300, 'chunkHours', 24);
for pn = string(fieldnames(psr)).'
    ps = psr.(pn);
    prof = gaussianTemplate(2^16, profileOf(ps));
    wEq = sum(prof) / 2^16;                         % pulse area [turns] of the unit-peak profile
    for fl = ["nominal", "worst"]
        psr.(pn).rho.(fl) = 10*log10(ps.S1400 / (sefd.(fl) * wEq));
    end
    fprintf('runScenario: %s, W_eq %.2f ms, rho %.2f / %.2f dB (nominal / worst)\n', ps.name, ...
        wEq / ps.f0 * 1e3, psr.(pn).rho.nominal, psr.(pn).rho.worst);
end

summary = struct([]);
for ip = 1:numel(pulsars)
    ps = psr.(pulsars(ip));
    passes = passTimes(ps.lat, ps.dec, elMin, nDays, tSid);
    nPer = max(1, round(ps.subintTime * ps.f0)); tSub = nPer / ps.f0;
    TRef = 0.5 / ps.f0;
    tmpl = gaussianTemplate(ps.nBin, profileOf(ps));
    dPhi = 0.1 * ps.fwhm;              % true pulse phase against the ephemeris [turns]
    tag = ""; M = [];
    if motion.on                       % the truth of a moving receiver (the offset is set per day below)
        M = motionTruth(motion, ps, tSid, nDays, ip);
        tag = "_motion_" + motion.clock;
        if motion.clock == "sdr" && motion.track == "none", tag = tag + "_notrack"; end
    end
    for fl = fluxes
        for va = variants
            fOut = fullfile(outDir, sprintf('scenario_%s_%s_%s_%dd%s.mat', pulsars(ip), fl, va, nDays, tag));
            if ~overwrite && isfile(fOut)
                fprintf('runScenario: %s exists, skipped\n', fOut);
                R = load(fOut).R; summary = [summary, R.summary]; %#ok<AGROW>
                continue
            end
            tCase = tic;
            simArgs = {'NChan', nCh, 'ChannelGain', gc, 'NPol', nPol, 'NBin', ps.nBin, ...
                'SubintPeriods', nPer, 'Verbose', false};
            if isempty(ps.comps), simArgs = [simArgs, {'FWHM', ps.fwhm}]; %#ok<AGROW>
            else, simArgs = [simArgs, {'Components', ps.comps}]; end %#ok<AGROW>
            adArgs = {'TargetSNR', 7, 'Weighting', 'optimal', 'AcrossGaps', true, 'Verbose', false};
            maxS = Inf;
            switch va
                case "clean"                        % receiver noise only
                    simArgs = [simArgs, {'Bnoise', BnC, 'BinDt', dtC, 'NoiseCoeffs', ncC}]; %#ok<AGROW>
                case {"rfi", "rfiCal", "rfiCalScint", "rfiCalScintFixed"}
                    simArgs = [simArgs, {'RFI', PR}]; %#ok<AGROW>     % realistic RFI after excision
                    if contains(va, "Cal")          % the observer calibrates each channel's noise
                        simArgs = [simArgs, {'ObserverNoise', 'measured'}]; %#ok<AGROW>
                    end
                    if contains(va, "Scint")
                        simArgs = [simArgs, {'Scintillation', ps.scint}]; %#ok<AGROW>
                    end
                    if endsWith(va, "Fixed")        % fixed groups: the adaptive case's median length
                        fA = fullfile(outDir, sprintf('scenario_%s_%s_rfiCalScint_%dd.mat', pulsars(ip), fl, nDays));
                        if ~isfile(fA)
                            error('runScenario:fixed', 'Run the rfiCalScint case first (%s).', fA);
                        end
                        RA = load(fA).R;
                        maxS = round(median([RA.groups.nMerged]));
                        adArgs = {'TargetSNR', Inf, 'Weighting', 'optimal', 'AcrossGaps', true, 'Verbose', false};
                    end
                otherwise
                    error('runScenario:variant', 'Unknown variant %s.', va);
            end
            G = struct('day', {}, 'nMerged', {}, 'hours', {}, 'spanHours', {}, 'T0', {}, ...
                'valid', {}, 'good', {}, 'err', {}, 'sigma', {}, 'snr', {}, 'chi2', {}, 'pChi2', {}, ...
                'tAbs', {}, 'resid', {}, 'trueResid', {}, 'trueGeo', {}, 'trueClock', {}, 'clockUsed', {});
            carry = []; nSubTotal = 0; maxDrift = 0;
            ehat = [];                                          % the observer's clock estimate
            if motion.on, ehat = M.ehat0; end                   % (day 1: the calibration in port)
            for d = 1:nDays
                if motion.on                                     % the truth with today's clock estimate
                    ca = clockArgsOf(M, ehat);
                    truthOf = @(s) trueResiduals(M.t0 + s, M.psr, M.trajTrue, M.trajAssumed, ca{:});
                    dPhi = @(s) offsetAt(truthOf, s);
                end
                % sub-int grid of the day's pass: complete turns inside it, nPer per sub-int
                k1 = ceil(ps.f0 * (passes(d, 1) - TRef) + 0.5);
                k2 = floor(ps.f0 * (passes(d, 2) - TRef) - 0.5);
                nSubDay = floor((k2 - k1 + 1) / nPer);
                mChunk = max(1, floor(ps.chunkHours * 3600 / tSub));    % sub-ints per chunk
                for c0 = 0:mChunk:nSubDay - 1
                    kA = k1 + c0 * nPer;                                 % first turn of the chunk
                    kB = k1 + min(c0 + mChunk, nSubDay) * nPer - 1;     % last turn
                    win = TRef + ([kA - 0.5 - 1e-6, kB + 0.5 + 1e-6]) / ps.f0;
                    [iC, fC, tC] = simulateFold(ps.f0, ps.rho.(fl), [], 'Passes', win, simArgs{:}, ...
                        'PhaseOffset', dPhi, 'Seed', 1e5 * ip + 1000 * d + c0);
                    nSubTotal = nSubTotal + iC.nSub;
                    if iC.nSub > 1                                       % drift of the true offset per sub-int
                        maxDrift = max(maxDrift, max(abs(diff(tC.phaseOffset))));
                    end
                    if ~isempty(carry)
                        [fC, iC] = catFolds(carry.fold, carry.info, fC, iC);
                    end
                    [fM, iM, sel] = adaptiveSubints(fC, iC, tmpl, BnC * nPol, adArgs{:}, 'MaxSubints', maxS);
                    nG = numel(sel.nMerged);
                    closed = sel.complete | sel.nMerged >= maxS;
                    nClosed = nG - ~closed(end);                         % the last group may be open
                    if ~closed(end)
                        carry = struct();
                        [carry.fold, carry.info] = subFold(fC, iC, find(sel.group == nG));
                    else
                        carry = [];
                    end
                    if nClosed == 0, continue, end
                    [fT, iT] = subFold(fM, iM, 1:nClosed);
                    toa = foldToTOA(fT, iT, profileOf(ps), BnC * nPol, 'Weighting', 'optimal', 'Verbose', false);
                    if motion.on                                         % the truth at each TOA's reference time
                        rT = truthOf(fT.subint.tRef.');
                        truth = fT.subint.tRef + rT.phaseOffset.' / ps.f0;
                    else
                        rT = struct('residual', NaN(1, nClosed), 'geometric', NaN(1, nClosed), ...
                            'clock', NaN(1, nClosed), 'clockEstimate', NaN(1, nClosed));
                        truth = fT.subint.tRef + dPhi / ps.f0;
                    end
                    for g = 1:nClosed
                        G(end+1) = struct('day', d, 'nMerged', sel.nMerged(g), ...
                            'hours', sel.nMerged(g) * tSub / 3600, ...
                            'spanHours', (fT.subint.turnLast(g) - fT.subint.turnFirst(g) + 1) / ps.f0 / 3600, ...
                            'T0', sel.T0(g), 'valid', toa.valid(g), 'good', toa.good(g), ...
                            'err', toa.toa(g) - truth(g), 'sigma', toa.toaErr(g), 'snr', toa.snr(g), ...
                            'chi2', toa.redChi2(g), 'pChi2', toa.pChi2(g), ...
                            'tAbs', tAbsOf(M, fT.subint.tRef(g)), 'resid', toa.toa(g) - fT.subint.tRef(g), ...
                            'trueResid', rT.residual(g), 'trueGeo', rT.geometric(g), ...
                            'trueClock', rT.clock(g) - rT.clockEstimate(g), ...
                            'clockUsed', clockUsedAt(M, ehat, fT.subint.tRef(g))); %#ok<SAGROW>
                    end
                end
                if motion.on && motion.clock == "sdr" && motion.track == "daily"
                    ehat = trackClock(G, M, ehat);               % refit, predict tomorrow's clock
                end
            end
            nLeft = 0;
            if ~isempty(carry), nLeft = carry.info.nSub; end
            S = scenarioSummary(G, nDays, nLeft, nSubTotal, ps, fl, va);
            summary = [summary, S]; %#ok<AGROW>
            R = struct('groups', G, 'summary', S, 'pulsar', ps, 'flux', fl, 'variant', va, ...
                'nDays', nDays, 'passes', passes, 'nPer', nPer, 'nPol', nPol, 'dPhi', dPhi, ...
                'motion', motion, 'motionTruth', motionRecord(M), ...
                'maxDriftPerSubint', maxDrift, ...
                'made', char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm')));
            if isa(dPhi, 'function_handle'), R.dPhi = func2str(dPhi); end
            save(fOut, 'R');
            fprintf('runScenario: %s %s %s, %d days: %d good TOAs, %.1f min\n', ps.name, fl, va, ...
                nDays, S.nGood, toc(tCase) / 60);
            if motion.on
                gd = [G.good];
                fprintf(['runScenario:   motion (%s clock): true residual %.1f ... %.1f us (geometric %.1f ... %.1f, ' ...
                         'clock %.1f ... %.1f); largest drift per sub-int %.2e turns (%.1f %% of the FWHM)\n'], ...
                    motion.clock, 1e6 * min([G(gd).trueResid]), 1e6 * max([G(gd).trueResid]), ...
                    1e6 * min([G(gd).trueGeo]), 1e6 * max([G(gd).trueGeo]), 1e6 * min([G(gd).trueClock]), ...
                    1e6 * max([G(gd).trueClock]), maxDrift, 100 * maxDrift / ps.fwhm);
                if maxDrift > 0.05 * ps.fwhm
                    warning('runScenario:drift', ['The true offset drifts by %.1f %% of the FWHM within a sub-int: ' ...
                        'one offset per sub-int is too coarse.'], 100 * maxDrift / ps.fwhm);
                end
            end
        end
    end
end

fprintf(['\nrunScenario summary (%d simulated days; hours = data per TOA, span = calendar time; ' ...
         'sigma median; z = (TOA - truth)/sigma of good TOAs)\n'], nDays);
fprintf(['  %-10s %-7s %-13s | TOAs/day | hours med (10-90 %%) | span med | sigma [us] | c*sigma [km] | ' ...
         'std z  | |z|>5 | chi2   | pChi2<.05 | leftover | sigma*sqrt(h)\n'], 'pulsar', 'flux', 'variant');
for S = summary
    fprintf(['  %-10s %-7s %-13s | %8.2f | %5.2f (%5.2f-%5.2f)  | %8.2f | %10.3f | %12.2f | %.3f  | ' ...
             '%5.3f | %.4f | %9.3f | %8.3f | %.2f\n'], S.pulsar, S.flux, S.variant, S.perDay, ...
        S.hoursMed, S.hours10, S.hours90, S.spanMed, S.sigmaMed, S.cSigmaKm, S.zStd, S.zOut, ...
        S.chi2Mean, S.pLow, S.leftover, S.merit);
end


% =================================================================================
function P = passTimes(lat, dec, elMin, nDays, tSid)
%PASSTIMES  [nDays x 2] start and end [s] of the daily passes above elMin [deg]: centred on
% the transits at (d - 1/2)*tSid, d = 1..nDays (the first pass starts after t = 0).
c = (sind(elMin) - sind(lat) * sind(dec)) / (cosd(lat) * cosd(dec));
if c >= 1
    error('runScenario:visible', 'The pulsar never rises above %g deg at latitude %g.', elMin, lat);
end
H = acos(max(c, -1));                               % hour angle at the mask [rad]
half = H / (2*pi) * tSid;
if c <= -1, half = 0.5 * tSid - 60; end             % always up: one long pass per day
tr = ((1:nDays).' - 0.5) * tSid;
P = [tr - half, tr + half];
end

function [f, info] = catFolds(fA, iA, fB, iB)
%CATFOLDS  Sub-ints of fold A followed by those of fold B (same NBin, channels, weights).
f = fB;
for fn = ["sum", "prof", "weight", "weight2", "weightX"]
    f.(fn) = cat(2, fA.(fn), fB.(fn));
end
for fn = string(fieldnames(fB.subint)).'
    f.subint.(fn) = [fA.subint.(fn); fB.subint.(fn)];
end
nChan = size(f.sum, 3);
f.profTotal = reshape(sum(f.sum, 2), [], nChan) ./ reshape(sum(f.weight, 2), size(f.weight, 1), []);
info = iB;
info.nSub = iA.nSub + iB.nSub;
end

function [f, info] = subFold(f0, i0, idx)
%SUBFOLD  The sub-ints idx of a fold.
f = f0;
for fn = ["sum", "prof", "weight", "weight2", "weightX"]
    f.(fn) = f0.(fn)(:, idx, :, :);
end
for fn = string(fieldnames(f0.subint)).'
    f.subint.(fn) = f0.subint.(fn)(idx);
end
nChan = size(f.sum, 3);
f.profTotal = reshape(sum(f.sum, 2), [], nChan) ./ reshape(sum(f.weight, 2), size(f.weight, 1), []);
info = i0;
info.nSub = numel(idx);
end

function S = scenarioSummary(G, nDays, nLeft, nSubTotal, ps, fl, va)
%SCENARIOSUMMARY  The numbers of one case.
good = [G.good].';
h = [G.hours].'; hg = h(good);
z = [G.err].' ./ [G.sigma].'; zg = z(good);
sg = [G.sigma].'; sg = sg(good);
nTOA = numel(G);
S = struct('pulsar', ps.name, 'flux', fl, 'variant', va, 'nTOA', nTOA, 'nGood', nnz(good), ...
    'perDay', nnz(good) / nDays, 'hoursMed', median(hg), 'hours10', pct(hg, 0.1), ...
    'hours90', pct(hg, 0.9), 'spanMed', median([G(good).spanHours]), ...
    'sigmaMed', median(sg) * 1e6, 'cSigmaKm', median(sg) * 299792.458, ...
    'zMean', mean(zg), 'zStd', std(zg), 'zOut', mean(abs(zg) > 5), ...
    'chi2Mean', mean([G(good).chi2]), 'pLow', mean([G(good).pChi2] < 0.05), ...
    'leftover', nLeft / max(nSubTotal, 1), 'merit', median(sg) * 1e6 * sqrt(median(hg)));
end

function M = motionTruth(motion, ps, tSid, nDays, ip)
%MOTIONTRUTH  The truth of a moving receiver for one pulsar: the scenario's absolute t = 0, the true
% and assumed sites, and the true clock with the observer's estimate (trueResiduals arguments).
% t0: the first transit at the site after motion.epoch + tSid/2, minus tSid/2, so that passTimes'
% transits ((d - 1/2) tSid) are real ones; hour angle = Earth rotation angle (as delayModel) + lon - RA.
ts = motion.epoch + tSid / 2;
H = 2*pi * (0.7790572732640 + 1.00273781191135448 * (ts - 69.184) / 86400) + deg2rad(ps.lon - ps.ra);
M.t0 = ts + mod(-H, 2*pi) / (2*pi) * tSid - tSid / 2;
M.psr = struct('f0', ps.f0, 'F1', 0, 'TRef', M.t0, 'raDeg', ps.ra, 'decDeg', ps.dec);
M.trajAssumed = makeTrajectory('fixed', 'Lat', ps.lat, 'Lon', ps.lon);
M.trajTrue = makeTrajectory('fixed', 'Lat', ps.lat, 'Lon', ps.lon, 'OffsetENU', motion.offsetENU);
switch motion.clock
    case "perfect"
        M.clk = []; M.ehat0 = []; M.clockH = [];
        M.clockRate = 0;
    case "sdr"                                      % crystal noise and a rate within +-2 ppm
        rs = RandStream('mt19937ar', 'Seed', 500 + ip);
        M.clockRate = 2e-6 * (2 * rand(rs) - 1);
        M.clk = makeClock([M.t0, M.t0 + (nDays + 1) * tSid], 'Preset', "crystal", 'Rate', M.clockRate, 'Seed', 600 + ip);
        r0 = M.clockRate + motion.clockRateErr; b0 = motion.clockOffsetErr; t0 = M.t0;
        M.ehat0 = @(c) b0 + r0 * (c - t0);          % the calibration in port
        M.clockH = M.clk.h;                         % the observer knows its oscillator's noise class
    otherwise
        error('runScenario:clock', 'motion.clock must be "perfect" or "sdr".');
end
fprintf(['runScenario: motion for %s: t = 0 at %.1f s after J2000 (first transit %.1f h after the epoch), ' ...
         'true site %s m ENU from the assumed one, %s clock (rate %.2e)\n'], ps.name, M.t0, ...
    (M.t0 + tSid / 2 - motion.epoch) / 3600, mat2str(motion.offsetENU), motion.clock, M.clockRate);
end

function t = tAbsOf(M, s)
%TABSOF  Absolute time (s since J2000, TT) of scenario time s; NaN without motion.
if isempty(M), t = NaN; else, t = M.t0 + s; end
end

function a = clockArgsOf(M, ehat)
%CLOCKARGSOF  trueResiduals' clock arguments: the true clock and the observer's estimate ({} if perfect).
if isempty(M.clk), a = {}; else, a = {'Clock', M.clk, 'ClockEstimate', ehat}; end
end

function v = clockUsedAt(M, ehat, s)
%CLOCKUSEDAT  The observer's clock estimate at scenario time s (0 without a clock).
if isempty(M) || isempty(M.clk), v = 0; else, v = ehat(M.t0 + s); end
end

function ehat = trackClock(G, M, ehat)
%TRACKCLOCK  The observer's daily clock update (block D4b-2): fit all good TOAs so far with navFit and the
% clock model, its own estimate added back (the residuals then hold the uncorrected clock), and predict
% the clock: the fitted constant (offset at t0; it also holds the polar position term) + clockAt.
g = [G.good];
if nnz(g) < 20, return, end                         % ~2 days: a rate fitted from one day extrapolates poorly
toa = struct('t', [G(g).tAbs], 'r', [G(g).resid] + [G(g).clockUsed], 'sigma', [G(g).sigma]);
ws = warning('off', 'navFit:converge');
f = navFit(toa, struct('raDeg', M.psr.raDeg, 'decDeg', M.psr.decDeg), M.trajAssumed, 'Clock', M.clockH, ...
    'ClockStart', M.t0);
warning(ws);
c1 = f.p(strcmp(f.names, 'c1'));
ehat = @(c) c1 + f.clockAt(c);
end

function v = offsetAt(truthOf, s)
%OFFSETAT  The true phase offset [turns] at scenario times s (for simulateFold's 'PhaseOffset').
r = truthOf(s);
v = r.phaseOffset;
end

function S = motionRecord(M)
%MOTIONRECORD  What R keeps of the motion truth: t0, the pulsar, the sites, the clock rate.
S = struct();
if isempty(M), return, end
S = struct('t0', M.t0, 'psr', M.psr, 'trajTrue', M.trajTrue, 'trajAssumed', M.trajAssumed, ...
    'clockRate', M.clockRate);
end

function prof = profileOf(ps)
%PROFILEOF  The pulse profile for gaussianTemplate / foldToTOA: components, or the single FWHM.
prof = ps.comps;
if isempty(prof), prof = ps.fwhm; end
end

function v = pct(x, p)
%PCT  Percentile p (0..1) of x by linear interpolation of the sorted values (NaN if empty).
x = sort(x(:));
if isempty(x), v = NaN; return, end
if isscalar(x), v = x; return, end
v = interp1(linspace(0, 1, numel(x)), x, p);
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
