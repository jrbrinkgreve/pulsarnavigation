function testRFIChain(cases, opts)
%TESTRFICHAIN  Slow test: the fast simulator's RFI model against the voltage chain (block C4b).
%{
Run from the PulsarSimMatlab folder: runAllTests("testRFIChain") (not in the default run: slow),
or testRFIChain(["realistic" "inmarsatNarrow" "noiseOnly"]) with tests/ on the path;
testRFIChain("noiseOnly", 'Length', 0.3) is a few-minute smoke run of the code (too short for
the 10-turn statistics).
~30 min per case (the 1 s voltage-chain run dominates), ~5-8 GB memory, ~40 GB of temporary
files in a folder of its own in tempdir (removed at the end, also on error; until 20 Oct one
fixed folder, so two runs at once destroyed each other's files). Uses pipelineParams (the RFI scenario lists,
receiver and processing settings, ephem; DM 25 since 20 Oct) - testParams has no RFI lists.

The RFI model (simulator.html, C4b): rfiProfile measures per channel what RFI leaves in the
detected power after excision (level, variance V and lag covariances X(1..Lmax) of the time
bins, valid fraction) and simulateFold('RFI', prof) draws folds with those statistics. That
assumes (a) stationary RFI, (b) time bins correlated only over the measured lags, so second
moments fix the folded noise, (c) Gaussian folded profiles, (d) blanking not locked to the pulse
phase. testRFIModel checks the machinery on synthetic data with closed-form truth; this test
checks the assumptions on the chain's real RFI and excision, with the same comparison:

Per case, one pulsar-free voltage-chain run of 1 s as runRFIProfile makes it (zero sky +
receiver noise + the case's RFI with beam passages outside the run, I/Q, processIQ on the
channel path with excision; noise seed noiseSeed + 1000, independent of the production
profiles), its detected power split in two halves (info_det.fullySupportedBins):
  half A -> rfiProfile (the calibration; its error from the scatter of 8 blocks);
  half B -> foldProfile with the excision weights at the ephemeris period, NBin 1024 (~10 time
            bins per phase bin per turn), 1 and 10 turns per sub-int (complete sub-ints only);
  60 simulated folds of the same layout from simulateFold with the half-A profile (no pulsar).
Per channel and sub-int length, real (half B) against simulated:
  1. level       mean of the fold / the profile's reference level - stationarity of the level;
  2. noise/model mean u^2, u = (prof - level)/s_j with s_j from the fold's own weight2 (the
                 observer's clean model) - THE second-moment check: correlations beyond Lmax or
                 slow variation make the real value exceed the simulated one, more at 10 turns;
  3. nb. corr    correlation of u in neighbouring phase bins;
  4. kurtosis    (information) excess kurtosis of u against the median over the channels of the
                 same run, sigma = 1.4826 x the median absolute deviation over the channels;
  5. tails       (information) count of |u| > 4 sigma against the mean count of the lowest 90 %
                 of channels (Poisson tail per channel). A phase bin of a short sub-int averages
                 few complex samples, so the detected power is gamma-like (excess kurtosis
                 ~6/(Bn x phase-bin time): ~0.18 at 1 turn, NBin 1024), which simulateFold's
                 Gaussian draw does not reproduce; RFI adds more in some channels. What matters
                 for the simulator (sub-ints of minutes) is whether it averages out: check 8;
  6. flatness    variance of the sub-int-averaged profile over phase / the scatter between
                 sub-ints, against the same statistic with every sub-int rotated by a random
                 number of phase bins (100 rotations of the real data: the null keeps all its
                 correlations and non-Gaussianity) - catches RFI that builds up in the fold;
  7. blanking    the same flatness statistic for the blanked fraction per phase bin (1 - the
                 fold weight / the weight of an unblanked fold of the same data) in channels
                 with blanking - catches blanking locked to the period;
  8. kurt. decay (1 turn) the 1-turn fold merged in groups of 4 turns (mergeSubints): with
                 independent turns the excess kurtosis falls to k1/4; z = (k4 - k1/4) / jackknife
                 error (over the sub-ints, valid for non-Gaussian data) - a slower decay is
                 non-Gaussianity that does not average out (bursts longer than a turn), which
                 the Gaussian draw would miss also for long sub-ints.
z = (real - simulated mean) / sqrt(simulated variance + calibration variance) for 1-3, against
the clean-channel reference and spread for 4, from the rotation null for 6-7, jackknife for 8.
A channel is flagged at |z| > 5 (6-8 one-sided: only excess structure or slow decay counts;
less structure than random is harmless; tails: Poisson tail p < 2.9e-7, the same level); with
~1800 tests per case the Gaussian false-flag rate is 0.1 %, ~1 % with the 60-realization
estimate of sigma (t with 59 dof). Pooled per check: sum(z)/sqrt(n) over the channels for 1-3,
|z| > 5 fails (a bias common to all channels; channels sharing RFI make this conservative);
for 6-7 the statistic of the channel average with shared rotations (broadband impulses blank all
channels at once: their per-channel z values are not independent, 20 Oct); none for 4, 5, 8.
A check FAILs if any channel is flagged or the pooled z is out; 4 and 5 never fail (information).
In the realistic case a FAIL can be a finding about where the second-moment model breaks
(radar sidelobe pulses, bursty emitters) rather than a bug: the table of flagged channels
(number, frequency, check, real, simulated or expected, z) makes it traceable. noiseOnly should
pass everything. Checks 4-8 revised 20 Oct after the first full run (Jasper's OK): realistic
case, 1 turn: kurtosis/tails flagged the partly covered channels and the carrier channel, and
the blanking pooled z was -6.5 from correlated channels.
Since 21 Oct 2026 the half-A profile has 30 lags (as runRFIProfile: narrow carriers correlate
beyond 7) and the chain's detectRFI takes per-channel thresholds (no false flags on continuous
narrowband RFI) - the two causes of the narrow-Inmarsat FAILs of 20 Oct.
Errors at the end if any check fails. (20 Oct 2026, block C4b.)
%}

arguments
    cases (1,:) string {mustBeMember(cases, ["realistic", "inmarsatNarrow", "noiseOnly"])} = ...
        ["realistic", "inmarsatNarrow"]
    opts.Length (1,1) double {mustBePositive} = 1   % [s] chain run per case
end

root = fullfile(fileparts(mfilename('fullpath')), '..');
addpath(fullfile(root, 'functions'));
oldDir = cd(root); restoreDir = onCleanup(@() cd(oldDir)); %#ok<NASGU>
P = loadParams(fullfile(root, 'pipelineParams.m'));   % RFI lists, receiver/processing settings,
% ephem - as a struct: in a function, names like dutycycle would resolve to toolbox functions
ephem = P.ephem; rx = P.rx; fLow = P.fLow; fHigh = P.fHigh; chanWidth = P.chanWidth;
rfiRealistic = P.rfiRealistic;
tmp = tempname(tempdir);                         % own folder per run: two runs at once (e.g.
mkdir(tmp);                                      %   Jasper's and Claude's) cannot collide
cleanTmp = onCleanup(@() rmdir(tmp, 's')); %#ok<NASGU>

Lc = opts.Length; NB = 1024; pers = [1 10]; nSim = 60; nBlk = 8; nRot = 100; zThr = 5;
nLag = 30;                                       % profile lags, as runRFIProfile (21 Oct 2026)
pTail = 0.5 * erfc(zThr / sqrt(2));              % one-sided level of zThr
pF = struct('chanIQ', fullfile(tmp, "p_chan"), 'chanBlank', fullfile(tmp, "p_chan_blanked"), ...
    'chanDedisp', fullfile(tmp, "p_dedisp_chan"), 'chanPower', fullfile(tmp, "p_dedisp_chan_power.dat"), ...
    'chanDedispB', fullfile(tmp, "p_dedispB_chan"), 'chanPowerB', fullfile(tmp, "p_dedispB_chan_power.dat"), ...
    'chanWeight', fullfile(tmp, "p_dedispB_chan_weights.dat"), 'foldChan', fullfile(tmp, "p_fold_chan.mat"));
procP = P.proc; procP.frontEnd = 'channels'; procP.excision = true;
procP.subintPeriods = 10;                        % processIQ's own fold is not used
procP.weightsArgs = {'MaxMemoryGB', 8};          % blankingWeights holds ~4.8 GB at 1 s
% Inmarsat as ten narrow carriers (as runRFIProfile / runRFITest case 15)
inmF = fLow + ((105:114) - 0.5) * chanWidth + [-1.1 0.4 -0.6 1.2 -0.2 0.9 -1.3 0.1 0.7 -0.9] * 1e6;
inmNarrow = arrayfun(@(f) rfiSource('noise', 'Freq', f, 'Bandwidth', 200e3, 'INRdB', -21, ...
    'Label', sprintf('Inmarsat carrier %.2f MHz', f / 1e6)), inmF, 'UniformOutput', false);
inmNarrow = [inmNarrow{:}];
isInm = startsWith({rfiRealistic.label}, 'Inmarsat');
checks = ["level", "noise/model", "nb. corr", "kurtosis", "tails", "flatness", "blanking", "kurt. decay"];
infoOnly = [4 5]; oneSided = [6 7 8];         % kurtosis, tails: information; structure: excess only
summary = strings(0, 1); fails = {};

for ci = 1:numel(cases)
    cs = cases(ci);
    tCase = tic;
    switch cs
        case "realistic",      src = rfiRealistic;
        case "inmarsatNarrow", src = [rfiRealistic(~isInm), inmNarrow];
        case "noiseOnly",      src = struct([]);
    end
    for k = 1:numel(src)                         % steady part: beam passages outside the run
        if isfinite(src(k).scanPeriod) && src(k).scanPeriod > 0 && ...
                abs(mod(src(k).beamTime - Lc/2 + src(k).scanPeriod/2, src(k).scanPeriod) ...
                    - src(k).scanPeriod/2) < Lc
            src(k).beamTime = Lc + 5;
        end
    end

    % ---- the voltage chain without pulsar (as runRFIProfile) --------------------------------
    info_zero = generatePulsarSignal(fullfile(tmp, "p_zero.dat"), P.T, P.f_in, 0, Lc, P.dutycycle, ...
        'Seed', P.seed, 'EnvelopeMode', P.genEnvMode, 'Verbose', false);
    info_rx = addNoiseAndRFI(info_zero.file, fullfile(tmp, "p_rx.dat"), info_zero.actualFsOut, ...
        'Band', [fLow fHigh], 'NoiseStd', 1, 'RFI', src, 'Seed', P.noiseSeed + 1000, 'Verbose', false);
    delete(info_zero.file);
    info_IQ = applyIQmodulation(info_rx.file, fullfile(tmp, "p_rx_IQ.dat"), info_rx.actualFsOut, ...
        P.fs, P.fLO, 'FilterOrder', P.filterOrder, 'Band', [fLow fHigh], 'Verbose', false);
    delete(info_rx.file);
    out = processIQ(info_IQ, rx, ephem, procP, pF, 'Verbose', false);
    delete(info_IQ.file);
    det = out.info_det; iw = out.info_w;
    Bn = det.noise.Bnoise; dt = det.binDt; nChan = det.nChan;
    fprintf('testRFIChain, case %s: %g s chain run, DM %g, %d channels, %.0f s\n', cs, Lc, ephem.DM, ...
        nChan, toc(tCase));

    % ---- half A: the profile and its calibration error ---------------------------------------
    k = det.fullySupportedBins; mid = floor((k(1) + k(2)) / 2);
    detA = det; detA.fullySupportedBins = [k(1), mid];
    detB = det; detB.fullySupportedBins = [mid + 1, k(2)];
    fprintf('  half A: ');
    prof = rfiProfile(detA, iw, 'Label', sprintf('%s, half A (testRFIChain)', cs), 'Lmax', nLag);
    e = round(linspace(k(1), mid + 1, nBlk + 1));
    bLev = zeros(nChan, nBlk); bSum = zeros(nChan, nBlk);
    for b = 1:nBlk
        d = detA; d.fullySupportedBins = [e(b), e(b+1) - 1];
        pb = rfiProfile(d, iw, 'RefLevel', prof.refLevel, 'Lmax', nLag, 'Verbose', false);
        bLev(:, b) = pb.level.';
        bSum(:, b) = pb.noiseCoeffs(:, 1) + 2*sum(pb.noiseCoeffs(:, 2:end), 2);
    end
    sumP   = prof.noiseCoeffs(:, 1) + 2*sum(prof.noiseCoeffs(:, 2:end), 2);
    calLev = std(bLev, 0, 2) / sqrt(nBlk);
    calSum = std(bSum, 0, 2) / sqrt(nBlk) ./ sumP;   % relative error of the sum rule
    blanked = prof.validFraction < 1 - 1e-6;

    % ---- half B: real folds vs simulated folds, 1 and 10 turns per sub-int -------------------
    rs = RandStream('mt19937ar', 'Seed', 7);
    for per = pers
        fa = {'TRef', ephem.TRef, 'F1', ephem.F1, 'NBin', NB, 'SubintPeriods', per, ...
              'SaveFile', false, 'Verbose', false};
        [iR, fR] = foldProfile(detB, 'x', ephem.f0, fa{:}, 'DataWeights', iw);
        [~, fU] = foldProfile(detB, 'x', ephem.f0, fa{:});     % unblanked: the time grid's weights
        use = fR.subint.nTimeBins(:) >= 0.999 * per / (ephem.f0 * dt);
        nS = nnz(use);
        if nS < 2
            fprintf('  half B, %2d turn(s) per sub-int: %d complete sub-int(s) - skipped (run too short)\n', ...
                per, nS);
            continue
        end
        sR = foldStats(fR, use, Bn, dt, prof.refLevel, fU.weight, nRot, rs);
        sim = cell(nSim, 1);
        for r = 1:nSim
            [~, fS] = simulateFold(ephem.f0, -Inf, nS * per, 'NChan', nChan, 'RFI', prof, ...
                'NBin', NB, 'SubintPeriods', per, 'Seed', r, 'Verbose', false);
            sim{r} = foldStats(fS, true(nS, 1), Bn, dt, 1, [], 0, rs);
        end
        S = [sim{:}];
        mS = @(f) mean([S.(f)], 2); vS = @(f) var([S.(f)], 0, 2);
        z = NaN(nChan, numel(checks));
        z(:, 1) = (sR.lev - mS('lev')) ./ sqrt(vS('lev') + calLev.^2);
        z(:, 2) = (sR.ratio - mS('ratio')) ./ sqrt(vS('ratio') + (mS('ratio') .* calSum).^2);
        z(:, 3) = (sR.rho1 - mS('rho1')) ./ sqrt(vS('rho1'));
        % kurtosis and tails against the clean-channel reference of the same run: a phase bin of a
        % short sub-int averages few complex samples, so the detected power is gamma-like (excess
        % kurtosis ~6/(Bn*phase-bin time)), which the Gaussian simulation does not reproduce;
        % the median over channels (mostly clean) holds that, so only RFI's extra shows
        kRef = median(sR.kurt);
        kSig = 1.4826 * median(abs(sR.kurt - kRef));   % robust spread over the channels (the
        z(:, 4) = (sR.kurt - kRef) / kSig;              %   estimator scatters more than for a Gaussian)
        tc = sort(sR.tail);
        mT = max(mean(tc(1:ceil(0.9*nChan))), 1e-3);  % reference count: the lowest 90 % of channels
        pT = ones(nChan, 1);
        pT(sR.tail > 0) = gammainc(mT, sR.tail(sR.tail > 0));    % P(N >= count), Poisson
        z(:, 5) = (sR.tail - mT) ./ sqrt(mT);      % shown only
        z(:, 6) = sR.flatZ;
        z(blanked, 7) = sR.blankZ(blanked);
        % kurtosis decay (1 turn only): merged in groups of 4 turns, independent turns divide the
        % excess kurtosis by 4; slower decay = non-Gaussianity that does not average out (bursts)
        k4 = NaN(nChan, 1); kExp = sR.kurt / 4;
        if per == 1 && nS >= 8
            iu = find(use); nG4 = floor(nS / 4);
            grp = zeros(numel(use), 1); grp(iu(1:4*nG4)) = repelem(1:nG4, 4);
            [fR4, iR4] = mergeSubints(fR, iR, grp);
            s4 = foldStats(fR4, true(nG4, 1), Bn, dt, prof.refLevel, [], 0, rs);
            k4 = s4.kurt;
            z(:, 8) = (k4 - kExp) ./ sqrt(s4.kurtSig.^2 + (sR.kurtSig / 4).^2);
        end
        realV = [sR.lev, sR.ratio, sR.rho1, sR.kurt, sR.tail, sR.flat, sR.blank, k4];
        simV  = [mS('lev'), mS('ratio'), mS('rho1'), repmat(kRef, nChan, 1), repmat(mT, nChan, 1), ...
                 sR.flatNull, sR.blankNull, kExp]; % kurtosis, tails: clean-channel reference; decay: k1/4
        zPoolOne = [NaN(1, 5), sR.poolFlatZ, sR.poolBlankZ, NaN];   % channel averages, shared rotations

        fprintf(['  half B, %2d turn(s) per sub-int: %d complete sub-ints; real vs %d simulated folds ' ...
                 '(mean over channels); radiometer non-Gaussianity: median excess kurtosis %.3f ' ...
                 '(Gaussian simulation %.3f)\n'], per, nS, nSim, kRef, mean(mS('kurt')));
        for q = 1:numel(checks)
            ok = ~isnan(z(:, q));
            if ~any(ok)
                if checks(q) == "kurt. decay", why = 'checked at 1 turn only';
                else, why = 'no channel with blanking'; end
                fprintf('    %-11s (%s)\n', checks(q), why);
                continue
            end
            if q == 5
                flag  = pT < pTail;
                zPool = 0;                         % relative to the channels' own reference
            elseif q == 4
                flag  = ok & abs(z(:, q)) > zThr;
                zPool = 0;
            elseif ismember(q, oneSided)
                flag  = ok & z(:, q) > zThr;
                zPool = zPoolOne(q);               % NaN for the decay: per channel only
                if isnan(zPool), zPool = 0; end
            else
                flag  = ok & abs(z(:, q)) > zThr;
                zPool = sum(z(ok, q)) / sqrt(nnz(ok));
            end
            if ismember(q, oneSided)
                pass = ~any(flag) && zPool <= zThr;
            else
                pass = ~any(flag) && abs(zPool) <= zThr;
            end
            res = passStr(pass);
            if ismember(q, infoOnly), res = 'info'; pass = true; end
            fprintf('    %-11s real %9.4f, sim %9.4f, pooled z %+6.2f, flagged %3d of %3d: %s\n', ...
                checks(q), mean(realV(ok, q)), mean(simV(ok, q)), zPool, nnz(flag), nnz(ok), res);
            for c = find(flag).'
                fprintf('      ch %3d  %8.3f MHz  real %9.4f  sim %9.4f  z %+7.2f\n', c, ...
                    det.chanFreqs(c) / 1e6, realV(c, q), simV(c, q), z(c, q));
            end
            summary(end+1, 1) = sprintf('%-15s %2d turn(s)  %-11s %s', cs, per, checks(q), res); %#ok<AGROW>
            if ~pass, fails{end+1} = sprintf('%s %d turn(s) %s', cs, per, checks(q)); end %#ok<AGROW>
        end
    end
    delete(fullfile(tmp, '*'));                  % the case's chain files (~40 GB)
    fprintf('  case %s done in %.1f min\n', cs, toc(tCase) / 60);
end

fprintf('\ntestRFIChain summary:\n');
fprintf('  %s\n', summary);
if isempty(fails)
    fprintf('testRFIChain: all passed\n');
else
    error('testRFIChain:fail', 'testRFIChain: failed: %s', strjoin(fails, ', '));
end
end


% =================================================================================
function P = loadParams(file)
%LOADPARAMS  The variables of a parameter script as a struct.
run(file);
names = setdiff(who, 'file');
for i = 1:numel(names)
    P.(names{i}) = eval(names{i});
end
end

function st = foldStats(fold, use, Bn, dt, ref, W0, nRot, rs)
%FOLDSTATS  Per channel over the sub-ints in use: level (prof/ref), noise relative to the
% observer's model (mean u^2, u = (prof - level)/s_j, s_j = level*sqrt(weight2/(Bn*dt))/weight),
% neighbour-bin correlation and excess kurtosis of u (with its jackknife error over the
% sub-ints: right also for non-Gaussian data), count of |u| > 4 sigma; with nRot > 0 also the
% flatness z of the mean profile and (W0: the unblanked fold's weights) of the blanked
% fraction, against random rotations of the sub-ints - per channel, and for the averages over
% the channels with one shared rotation per sub-int (channels share RFI, e.g. broadband
% impulses blank all of them at once, so their z values are not independent).
P  = fold.prof(:, use, :) / ref;
W  = fold.weight(:, use, :); W2 = fold.weight2(:, use, :);
[N, nS, nChan] = size(P); nW = size(W, 3);
z0 = zeros(nChan, 1); n0 = NaN(nChan, 1);
st = struct('lev', z0, 'ratio', z0, 'rho1', z0, 'kurt', z0, 'kurtSig', z0, 'tail', z0, ...
    'flat', n0, 'flatNull', n0, 'flatZ', n0, 'blank', n0, 'blankNull', n0, 'blankZ', n0, ...
    'poolFlatZ', NaN, 'poolBlankZ', NaN);
if nRot > 0
    Dall = NaN(N, nS, nChan); Ball = NaN(N, nS, nChan);
end
for c = 1:nChan
    Wc = W(:, :, min(c, nW)); W2c = W2(:, :, min(c, nW)); Pc = P(:, :, c);
    ok = Wc > 0 & isfinite(Pc);
    L  = mean(Pc(ok));
    u  = (Pc - L) ./ (L * sqrt(W2c / (Bn*dt)) ./ Wc);
    u(~ok) = NaN;
    r  = mean(u(ok).^2);
    nb = u .* circshift(u, -1, 1);
    v  = u(ok) / sqrt(r);
    S2 = sum(u.^2, 1, 'omitnan'); S4 = sum(u.^4, 1, 'omitnan'); Nn = sum(ok, 1);
    kj = ((sum(S4) - S4) ./ (sum(Nn) - Nn)) ./ ((sum(S2) - S2) ./ (sum(Nn) - Nn)).^2 - 3;
    st.lev(c)     = L;
    st.ratio(c)   = r;
    st.rho1(c)    = mean(nb(~isnan(nb))) / r;
    st.kurt(c)    = mean(v.^4) - 3;
    st.kurtSig(c) = sqrt((nS - 1) / nS * sum((kj - mean(kj)).^2));   % jackknife over sub-ints
    st.tail(c)    = nnz(abs(v) > 4);
    if nRot > 0
        D = Pc - L; D(~ok) = NaN;
        [st.flat(c), st.flatNull(c), st.flatZ(c)] = flatness(D, nRot, rs);
        Dall(:, :, c) = D / L;
        if ~isempty(W0)
            B = 1 - Wc ./ W0(:, use);
            B(~(W0(:, use) > 0)) = NaN;
            [st.blank(c), st.blankNull(c), st.blankZ(c)] = flatness(B, nRot, rs);
            Ball(:, :, c) = B;
        end
    end
end
if nRot > 0
    [~, ~, st.poolFlatZ] = flatness(mean(Dall, 3, 'omitnan'), nRot, rs);
    if ~isempty(W0)
        [~, ~, st.poolBlankZ] = flatness(mean(Ball, 3, 'omitnan'), nRot, rs);
    end
end
end

function [F, Fnull, z] = flatness(D, nRot, rs)
%FLATNESS  Variance over phase of the sub-int-averaged D [NBin x nSub] divided by the expected
% part from the scatter between sub-ints (~1 without phase structure), and its z against nRot
% random rotations of the sub-ints (each column shifted by its own random number of bins).
[N, nS] = size(D);
F  = flatStat(D);
Fn = zeros(nRot, 1);
for m = 1:nRot
    sh  = randi(rs, N, 1, nS) - 1;
    idx = mod((0:N-1).' + sh, N) + 1 + N * (0:nS-1);
    Fn(m) = flatStat(D(idx));
end
Fnull = mean(Fn);
z = (F - Fnull) / std(Fn);
end

function F = flatStat(D)
m = mean(D, 2, 'omitnan');
n = sum(~isnan(D), 2);
v = var(D, 0, 2, 'omitnan') ./ n;
F = sum((m - mean(m, 'omitnan')).^2, 'omitnan') / sum(v, 'omitnan');
end

function s = passStr(p)
if p, s = 'PASS'; else, s = 'FAIL'; end
end
