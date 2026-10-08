%{
RUNLOCKEDRADAR - a radar locked to the pulsar period: what excision leaves in the fold (B5b)

A radar whose PRF is a multiple of the pulsar frequency (here 400 Hz = 4 f0)
puts its pulses at the same pulse phases every turn, so whatever excision
leaves of it folds coherently into the profile. One of its four pulses per
turn is placed on the flank of the pulse (+1 sigma), the worst place for a
TOA. (In reality the Earth's motion makes the pulsar phase drift against any
ground-based RFI over long folds; that belongs to block D.)

Because the pulsar is ~80,000x weaker at -54 dB than at -5 dB, a TOA test at
-5 dB alone says little. So the leftover radar is measured in units that
scale: its folded power in units of the noise baseline (independent of the
pulsar), compared with the pulsar profile height rho at any SNR.

Method (main's noise seed and -5 dB pulsar in every run):
  1. Receiver runs: reference (no radar), and with the locked radar (realistic
     0.1 us edges, and instant edges). Channel by channel, radar run minus
     reference = the radar alone (the chain is linear).
  2. Per case (edges x excision off/on): mask from the full data (detectRFI);
     the reference, the radar-only signal and the full data all go through
     blankChannels -> dedisperseChannels -> detectChannels -> foldProfile
     with the same mask and data weights (blankingWeights).
  3. Leftover profile L(phi) = folded radar-only power / channel baseline,
     averaged over channels (in baseline units the channels have equal noise).
  4. TOA bias at -5 dB: TOAs of the full data minus those of the reference
     (same noise), and as a cross-check the reference + leftover.
  5. Bias at weak SNR: noise-free folds rho*template + leftover (same weights)
     through estimateTOA ('optimal'), rho from -5 to -60 dB; the rho where the
     bias reaches 1 us.
  Cases: edges x excision off/on, and with excision the frequency guard of
  detectRFI (strong events blanked in +-FreqGuard channels) swept; and (B6-3)
  the locked radar on a rotating antenna seen only through -35 dB sidelobes,
  with excision without / with the periodic mask (periodicRFI). runCases
  picks a subset (only the radars it needs are simulated).

Needs: the sky files (main.m with runStage.sky once). Files in data/locked
(overwritten). ~9 min for all cases.
%}

% Paths and parameters
scriptDir = fileparts(mfilename('fullpath'));
addpath(fullfile(scriptDir, 'functions'));
pipelineParams;
doPlot = true;

lkDir = fullfile(dataDir, "locked");
info_gen  = loadInfo(fileRaw);
info_disp = loadInfo(fileDispersed);
template  = gaussianTemplate(nBin, ephem.profileFWHM);
rxArgs = {'Band', [fLow fHigh], 'SNRdB', snrDB, 'SignalInfo', info_gen, 'Seed', noiseSeed, 'Verbose', false};
iqArgs = {'FilterOrder', filterOrder, 'Band', [fLow fHigh], 'Verbose', false};
chArgs = {'ChanWidth', chanWidth, 'Verbose', false};

% The locked radar: PRF 4 f0; one pulse per turn on the flank after dedispersion
% (its channels are moved earlier by the dispersion delay between 1.3 GHz and refFreq)
fR    = 1300e6; pw = 2e-6;
prf   = 4 * ephem.f0;
Kdm   = 4.148808e3 * 1e12 * ephem.DM;                % [s Hz^2]
delta = Kdm * (1/fR^2 - 1/refFreq^2);                % 4.17 ms
sigT  = ephem.profileFWHM / ephem.f0 / (2*sqrt(2*log(2)));
tStart = mod(ephem.TRef + sigT + delta - pw/2, 1/prf);
% the radars: realistic edges, instant edges, realistic edges seen only through -35 dB
% sidelobes of a rotating antenna (main beam 5 s away)
radarArgs = {{'RiseTime', 0.1e-6}, {'RiseTime', 0}, ...
             {'RiseTime', 0.1e-6, 'ScanPeriod', 10, 'BeamTime', 5, 'BeamWidth', 39e-3, 'SidelobeDB', -35}};
radarName = {'0.1 us edges', 'instant edges', 'sidelobes -35'};
% cases: [radar, excision, FreqGuard, periodic mask]
caseList = [1 0 0 0; 1 1 0 0; 1 1 3 0; 1 1 5 0; 1 1 7 0; 2 0 0 0; 2 1 0 0; 2 1 5 0; ...
            1 1 7 1; 3 0 0 0; 3 1 7 0; 3 1 7 1];
runCases = 1:size(caseList, 1);                      % e.g. 9:12 for the periodic-mask cases

fprintf(['runLockedRadar: radar %.0f MHz, %.0f us, PRF %g Hz (= %d f0), +20 dB, one pulse per turn at ' ...
         '+%.0f us from the pulse peak; %.1f dB pulsar, L = %.3g s\n'], fR/1e6, pw*1e6, prf, ...
    round(prf / ephem.f0), sigT*1e6, snrDB, L);
tAll = tic;

% ---- 1. Receiver runs ---------------------------------------------------------------------
chanRef = receiverChannels(info_disp, fullfile(lkDir, "ref"), struct([]), rxArgs, iqArgs, chArgs, fs, fLO, fLow, fHigh);
chanTot = cell(1, numel(radarArgs)); chanR = chanTot;
for e = unique(caseList(runCases, 1)).'
    radar = rfiSource('pulsed', 'Freq', fR, 'PulseWidth', pw, 'PRF', prf, 'ChirpBW', 1e6, ...
        'StartTime', tStart, 'INRdB', 20, radarArgs{e}{:}, 'Label', 'locked radar');
    chanTot{e} = receiverChannels(info_disp, fullfile(lkDir, "tot" + e), radar, rxArgs, iqArgs, ...
        chArgs, fs, fLO, fLow, fHigh);
    chanR{e} = differenceChannels(chanTot{e}, chanRef, fullfile(lkDir, "radar" + e));
end
fprintf('  receiver runs: %.0f s\n', toc(tAll));

% ---- 2.-4. Cases ------------------------------------------------------------------------------
rhoDB = [-5 -10 -20 -30 -40 -50 -54 -60];
res = struct('edges', {}, 'excision', {}, 'guard', {}, 'periodic', {}, 'blanked', {}, 'Lpeak', {}, ...
    'Lflank', {}, 'pulse', {}, 'dMeas', {}, 'dPred', {}, 'dSyn', {}, 'rho1us', {}, 'L', {});
for c = runCases
    e = caseList(c, 1); exc = caseList(c, 2) == 1; g = caseList(c, 3); per = caseList(c, 4) == 1;
    tCase = tic;
    if exc
        [mask, info_rfi] = detectRFI(chanTot{e}, excisionArgs{:}, 'FreqGuard', g, 'Verbose', false);
        if per                                       % + every predicted pulse (B6)
            [perRows, info_per] = periodicRFI(chanTot{e}, info_rfi, periodicArgs{:}, 'Verbose', false);
            mask = [mask; perRows];
            fprintf('  periodicRFI: %d emitter(s)%s\n', numel(info_per.emitters), ...
                sprintf(', PRF %.4f Hz from %d pulses', [[info_per.emitters.prf]; [info_per.emitters.nInliers]]));
        end
    else
        mask = zeros(0, 3);
    end
    % reference first (it also makes the data weights of the mask), then the radar
    % alone and the full data with the same mask and weights
    [fRef, ifRef, det, info_w] = processSet(chanRef, mask, [], lkDir, ephem, refFreq, f_out, nBin, subintPeriods);
    fR_  = processSet(chanR{e},   mask, info_w, lkDir, ephem, refFreq, f_out, nBin, subintPeriods);
    fTot = processSet(chanTot{e}, mask, info_w, lkDir, ephem, refFreq, f_out, nBin, subintPeriods);

    % leftover radar profile in units of each channel's baseline, channels averaged
    aC = median(fRef.profTotal, 1, 'omitnan');                 % noise baseline per channel
    Lc = fR_.profTotal ./ aC;
    Lphi = mean(Lc, 2, 'omitnan');                             % [nBin x 1]
    P5 = mean((fRef.profTotal - aC) ./ aC, 2, 'omitnan');      % the -5 dB pulsar, same units
    flank = abs(mod(fRef.phase(:) - sigT * ephem.f0 + 0.5, 1) - 0.5) < 0.02;

    % TOAs at -5 dB: full data vs reference, and reference + leftover
    toaArgs = {'Bnoise', det.noise.Bnoise, 'Weighting', weighting, 'Verbose', false};
    tRef  = estimateTOA(fRef, ifRef, template, toaArgs{:});
    tTot  = estimateTOA(fTot, ifRef, template, toaArgs{:});
    tPred = estimateTOA(addFold(fRef, fR_, 1), ifRef, template, toaArgs{:});
    v = tRef.valid & tTot.valid & tPred.valid;
    dMeas = median(tTot.toa(v) - tRef.toa(v));
    dPred = median(tPred.toa(v) - tRef.toa(v));

    % weak SNR: noise-free rho * template + leftover, the real estimator
    dSyn = zeros(size(rhoDB));
    for k = 1:numel(rhoDB)
        rho = 10^(rhoDB(k)/10);
        f0  = syntheticFold(fRef, aC, rho, template);
        t0s = estimateTOA(f0, ifRef, template, toaArgs{:});
        t1s = estimateTOA(addFold(f0, fR_, 1), ifRef, template, toaArgs{:});
        vv = t0s.valid & t1s.valid;
        dSyn(k) = median(t1s.toa(vv) - t0s.toa(vv));
    end
    ab = abs(dSyn);
    k1 = find(ab >= 1e-6, 1);                                  % first rho (from bright) with |bias| >= 1 us
    if isempty(k1), rho1 = NaN;
    elseif k1 == 1, rho1 = rhoDB(1);
    else
        rho1 = interp1(log10(ab(k1-1:k1)), rhoDB(k1-1:k1), -6);
    end
    blanked = 0;
    if exc, blanked = mean(info_w.blankedFraction); end
    res(end+1) = struct('edges', radarName{e}, 'excision', exc, 'guard', g, 'periodic', per, ...
        'blanked', blanked, ...
        'Lpeak', max(abs(Lphi)), 'Lflank', max(abs(Lphi(flank))), 'pulse', max(P5), ...
        'dMeas', dMeas, 'dPred', dPred, 'dSyn', dSyn, 'rho1us', rho1, 'L', Lphi); %#ok<SAGROW>
    fprintf('  %s, excision %d, frequency guard %d, periodic mask %d: %.0f s\n', radarName{e}, exc, ...
        g, per, toc(tCase));
end
fprintf('runLockedRadar: done in %.0f s\n\n', toc(tAll));

% ---- Table --------------------------------------------------------------------------------------
fprintf(['leftover L = folded radar power / noise baseline (max over phase; on the flank); the pulsar ' ...
         'peak at %.0f dB is %.3g in these units (rho = %.3g)\n'], snrDB, res(1).pulse, 10^(snrDB/10));
fprintf('%-14s %4s %6s %4s %9s %10s %10s %11s %11s   %s   %s\n', 'radar', 'exc', 'guard', 'per', ...
    'blanked', 'L max', 'L flank', 'dTOA meas', 'dTOA pred', ['bias [us] at ' mat2str(rhoDB) ' dB'], '1 us at');
for k = 1:numel(res)
    r = res(k);
    fprintf('%-14s %4d %6d %4d %8.4f%% %10.2e %10.2e %9.3f us %9.3f us   %s   %5.1f dB\n', r.edges, ...
        r.excision, r.guard, r.periodic, 100*r.blanked, r.Lpeak, r.Lflank, r.dMeas*1e6, r.dPred*1e6, ...
        mat2str(round(r.dSyn*1e6, 3)), r.rho1us);
end
save(fullfile(lkDir, "runLockedRadar_results.mat"), 'res', 'rhoDB', 'tStart', 'prf', 'caseList', 'runCases');

if doPlot
    figure('Name', 'Locked radar: leftover in the fold');
    ph = fRef.phase(:); sel = abs(mod(ph + 0.5, 1) - 0.5) < 0.1;
    phs = mod(ph(sel) + 0.5, 1) - 0.5;
    [phs, o] = sort(phs);
    tiledlayout(2, 1);
    nexttile; hold on
    for k = 1:numel(res)
        Lk = res(k).L(sel); plot(phs, Lk(o), 'DisplayName', ...
            sprintf('%s, exc %d, guard %d, per %d', res(k).edges, res(k).excision, res(k).guard, ...
            res(k).periodic));
    end
    tpl = template(:); tpl = tpl(sel); tpl = tpl(o);
    plot(phs, 10^(snrDB/10) * tpl, 'k--', 'DisplayName', sprintf('pulsar %.0f dB', snrDB));
    ylabel('power / baseline'); legend('Location', 'best'); title('leftover radar profile vs pulsar');
    set(gca, 'YScale', 'linear'); grid on
    nexttile; hold on
    for k = 1:numel(res)
        Lk = abs(res(k).L(sel)); semilogy(phs, Lk(o) + eps, 'DisplayName', ...
            sprintf('%s, exc %d, guard %d, per %d', res(k).edges, res(k).excision, res(k).guard, ...
            res(k).periodic));
    end
    semilogy(phs, 10^(-54/10) * tpl, 'k:', 'LineWidth', 1.5);
    set(gca, 'YScale', 'log'); ylabel('|power| / baseline'); xlabel('pulse phase [turns]');
    title('log scale; dotted: pulsar at -54 dB'); grid on
end


% ======================================================================
%  Local functions
% ======================================================================
function info_chan = receiverChannels(info_disp, base, rfi, rxArgs, iqArgs, chArgs, fs, fLO, fLow, fHigh)
% Receiver noise + rfi at RF -> IQ -> filterbank; the RF and IQ files are deleted.
base = char(base);
info_rx = addNoiseAndRFI(info_disp.file, [base '_rx.dat'], info_disp.actualFsOut, rxArgs{:}, 'RFI', rfi);
info_IQ = applyIQmodulation(info_rx.file, [base '_iq.dat'], info_rx.actualFsOut, fs, fLO, iqArgs{:});
delete(info_rx.file);
info_chan = channelizeIQ(info_IQ.file, [base '_chan'], info_IQ.actualFsOut, info_IQ.fLO, fLow, fHigh, ...
    chArgs{:}, 'T0', info_IQ.t0);
delete(info_IQ.file);
end

function info = differenceChannels(infoA, infoB, base)
% Channel files A - B (cf32), e.g. the radar alone; info = infoA with the new files.
base = char(base);
info = infoA; info.file = string(base); info.chanFiles = strings(infoA.nChan, 1);
for j = 1:infoA.nChan
    a = readRaw(infoA.chanFiles(j)); b = readRaw(infoB.chanFiles(j));
    info.chanFiles(j) = sprintf('%s_ch%03d.dat', base, j);
    fid = fopen(info.chanFiles(j), 'w', 'ieee-le'); fwrite(fid, a - b, 'single'); fclose(fid);
end
end

function [fold, info_fold, info_det, info_w] = processSet(info_chan, mask, info_w, lkDir, ephem, refFreq, f_out, nBin, subintPeriods)
% Main's channel path on one set of channel files: blank (if mask), dedisperse,
% detect, fold with DataWeights (given, or made here from the mask) or, without a
% mask, NoiseCoeffs. Products in lkDir (overwritten by the next set).
if ~isempty(mask)
    info_chan = blankChannels(info_chan, mask, fullfile(lkDir, "blanked"), 'Verbose', false);
end
info_dc  = dedisperseChannels(info_chan, fullfile(lkDir, "dedisp"), ephem.DM, 'RefFreq', refFreq, 'Verbose', false);
info_det = detectChannels(info_dc, fullfile(lkDir, "power.dat"), info_dc.fs / round(info_dc.fs / f_out), 'Verbose', false);
if ~isempty(mask) && isempty(info_w)
    info_w = blankingWeights(info_chan, info_dc, info_det, mask, fullfile(lkDir, "weights.dat"), 'Verbose', false);
end
if isempty(info_w)
    noiseArgs = {'NoiseCoeffs', [info_det.noise.V, info_det.noise.X]};
else
    noiseArgs = {'DataWeights', info_w};
end
[info_fold, fold] = foldProfile(info_det, fullfile(lkDir, "fold.mat"), ephem.f0, 'F1', ephem.F1, ...
    'TRef', ephem.TRef, 'NBin', nBin, 'SubintPeriods', subintPeriods, noiseArgs{:}, ...
    'SaveFile', false, 'Verbose', false);
end

function f = addFold(f, g, scale)
% Fold f plus scale x the power of fold g (same weights): sum, prof, profTotal.
f.sum = f.sum + scale * g.sum;
f.prof = f.prof + scale * g.prof;
f.profTotal = f.profTotal + scale * g.profTotal;
end

function f = syntheticFold(f, aC, rho, template)
% Noise-free fold with the weights of f: channel c = aC(c) * (1 + rho * template).
nSub = size(f.prof, 2); nChan = size(f.prof, 3);
p = reshape(aC, 1, 1, nChan) .* (1 + rho * template(:));        % [nBin x 1 x nChan]
p = repmat(p, 1, nSub, 1);
w = f.weight; if size(w, 3) == 1, w = repmat(w, 1, 1, nChan); end
p(w == 0) = NaN;
f.prof = p;
f.sum = p .* w; f.sum(w == 0) = 0;
f.profTotal = reshape(sum(f.sum, 2), [], nChan) ./ reshape(sum(w, 2), [], nChan);
end

function x = readRaw(file)
fid = fopen(file, 'r', 'ieee-le');
c = onCleanup(@() fclose(fid));
x = fread(fid, [2 Inf], 'single=>single');
end
