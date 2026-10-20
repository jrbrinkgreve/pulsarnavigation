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
  - no motion (block D): TOA budgets, not navigation yet; the true pulse sits 0.1 FWHM from the
    predicted phase (a small ephemeris/position error) - adaptiveSubints closes groups on T0 at
    the predicted phase, which assumes the prediction is well inside the pulse (an error of a
    sizeable fraction of the FWHM would make the groups longer than needed, not the TOAs wrong);
  - rho in the simulator's convention (peak of the Gaussian pulse / receiver noise per
    polarization, rho = S_peak/SEFD), 'NPol' 2; the observer passes Bnoise x 2;
  - pulse and template: B0329+54's five Gaussian components at 1410 MHz (since 20 Oct; a single
    Gaussian with W50 overstated its peak by 3 dB), J0437 a single Gaussian (wings not modelled);
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

Outputs per case (printed table; data/scenario/scenario_<pulsar>_<flux>_<variant>_<n>d.mat):
good TOAs per day, hours of data per TOA (median, 10-90 %), calendar span per TOA, median
sigma (us and c*sigma in km), rms of (TOA - truth)/sigma and the fraction beyond 5, mean
reduced chi^2 and pChi2 < 0.05, leftover fraction, figure of merit sigma*sqrt(hours) [us
sqrt(h)]. (20 Oct 2026, block C; design agreed with Jasper the same day.)
%}

pulsars  = ["B0329", "J0437"];
fluxes   = ["nominal", "worst"];
variants = "clean";                  % "clean", "rfi", "rfiCal", "rfiCalScint", "rfiCalScintFixed"
                                     % (the fixed variant needs the rfiCalScint result first)
nDays    = 30;                       % simulated days
overwrite = false;                   % false: skip cases whose result file exists

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
% FWHM] in turns; J0437 a single Gaussian with W50 (its wings are not modelled yet). rho (simulator
% convention: pulse peak / receiver noise per polarization) = S1400 / (SEFD * W_eq), W_eq the pulse
% area of the unit-peak profile in turns; SEFD 3.2e5 Jy nominal, 7.4e5 Jy worst case (notes §11).
cB = [0.159 -11.427 1.479; 0.439 -1.123 3.281; 0.585 0 1.055; 0.152 7.719 2.751; 0.170 9.482 0.847];
cB = [cB(:, 1), cB(:, 2) / 360, cB(:, 3) * 2*sqrt(2*log(2)) / 360];
sefd = struct('nominal', 3.2e5, 'worst', 7.4e5);
psr.B0329 = struct('name', 'B0329+54', 'f0', 1/0.714519699726, 'dec', 54.5787, 'lat', 52.0, ...
    'fwhm', 0.0092, 'comps', cB, 'S1400', 0.203, 'nBin', 1024, ...
    'scint', struct('bandwidth', 7e6, 'timescale', 780, 'band', [fLow fHigh]), ...
    'subintTime', 120, 'chunkHours', 4);
psr.J0437 = struct('name', 'J0437-4715', 'f0', 173.6879458, 'dec', -47.2525, 'lat', -30.0, ...
    'fwhm', 0.0245, 'comps', [], 'S1400', 0.150, 'nBin', 512, ...
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
    for fl = fluxes
        for va = variants
            fOut = fullfile(outDir, sprintf('scenario_%s_%s_%s_%dd.mat', pulsars(ip), fl, va, nDays));
            if ~overwrite && isfile(fOut)
                fprintf('runScenario: %s exists, skipped\n', fOut);
                R = load(fOut).R; summary = [summary, R.summary]; %#ok<AGROW>
                continue
            end
            tCase = tic;
            simArgs = {'NChan', nCh, 'ChannelGain', gc, 'NPol', nPol, 'NBin', ps.nBin, ...
                'SubintPeriods', nPer, 'PhaseOffset', dPhi, 'Verbose', false};
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
                'valid', {}, 'good', {}, 'err', {}, 'sigma', {}, 'snr', {}, 'chi2', {}, 'pChi2', {});
            carry = []; nSubTotal = 0;
            for d = 1:nDays
                % sub-int grid of the day's pass: complete turns inside it, nPer per sub-int
                k1 = ceil(ps.f0 * (passes(d, 1) - TRef) + 0.5);
                k2 = floor(ps.f0 * (passes(d, 2) - TRef) - 0.5);
                nSubDay = floor((k2 - k1 + 1) / nPer);
                mChunk = max(1, floor(ps.chunkHours * 3600 / tSub));    % sub-ints per chunk
                for c0 = 0:mChunk:nSubDay - 1
                    kA = k1 + c0 * nPer;                                 % first turn of the chunk
                    kB = k1 + min(c0 + mChunk, nSubDay) * nPer - 1;     % last turn
                    win = TRef + ([kA - 0.5 - 1e-6, kB + 0.5 + 1e-6]) / ps.f0;
                    [iC, fC] = simulateFold(ps.f0, ps.rho.(fl), [], 'Passes', win, simArgs{:}, ...
                        'Seed', 1e5 * ip + 1000 * d + c0);
                    nSubTotal = nSubTotal + iC.nSub;
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
                    truth = fT.subint.tRef + dPhi / ps.f0;
                    for g = 1:nClosed
                        G(end+1) = struct('day', d, 'nMerged', sel.nMerged(g), ...
                            'hours', sel.nMerged(g) * tSub / 3600, ...
                            'spanHours', (fT.subint.turnLast(g) - fT.subint.turnFirst(g) + 1) / ps.f0 / 3600, ...
                            'T0', sel.T0(g), 'valid', toa.valid(g), 'good', toa.good(g), ...
                            'err', toa.toa(g) - truth(g), 'sigma', toa.toaErr(g), 'snr', toa.snr(g), ...
                            'chi2', toa.redChi2(g), 'pChi2', toa.pChi2(g)); %#ok<SAGROW>
                    end
                end
            end
            nLeft = 0;
            if ~isempty(carry), nLeft = carry.info.nSub; end
            S = scenarioSummary(G, nDays, nLeft, nSubTotal, ps, fl, va);
            summary = [summary, S]; %#ok<AGROW>
            R = struct('groups', G, 'summary', S, 'pulsar', ps, 'flux', fl, 'variant', va, ...
                'nDays', nDays, 'passes', passes, 'nPer', nPer, 'nPol', nPol, 'dPhi', dPhi, ...
                'made', char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm')));
            save(fOut, 'R');
            fprintf('runScenario: %s %s %s, %d days: %d good TOAs, %.1f min\n', ps.name, fl, va, ...
                nDays, S.nGood, toc(tCase) / 60);
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
