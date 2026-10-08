function testPeriodicRFI()
%TESTPERIODICRFI  Unit tests for periodicRFI (B6): find periodic RFI, blank its predicted pulses.
%{
Run from the PulsarSimMatlab folder: run('tests/testPeriodicRFI.m'). ~1 min.
Makes its own channel files (0.05 s through addNoiseAndRFI -> applyIQmodulation
-> channelizeIQ on a zero sky signal, the same noise seed in every run, so a
run minus the noise-only run is the RFI alone; ~3 GB in tempdir, deleted on
the way). Scored with the sources' parameters (ground truth, test only).

  1. Noise only: no emitter, nothing blanked by prediction.
  2. Broadband impulses (1000/s, +20 dB, Poisson): no periodic emitter.
  3. The scenario radar (1300 MHz, 2 us, PRF 373 Hz, 0.1 us edges): one
     emitter; PRF within 4 sigma of its own error bar and to 1e-5; every
     true pulse centre within TimeTol of a predicted one; all pulses blanked.
  4. A rotating radar seen only through -35 dB sidelobes (beam 5 s away;
     +3 dB per channel, most pulses missed by detectRFI): the emitter is
     found from the detected ones; with the predicted pulses all pulses are
     blanked, and the RFI energy left (radar windows, all channels) drops
     below 1e-4 of the radar's energy. The detectRFI-only numbers shown.
  Shown: the data blanked by prediction.
Errors at the end if any check fails.
%}

root = fullfile(fileparts(mfilename('fullpath')), '..');
addpath(fullfile(root, 'functions'));
oldDir = cd(root); restoreDir = onCleanup(@() cd(oldDir));
pipelineParams;
tmp = fullfile(tempdir, 'testPeriodicRFI');
if ~isfolder(tmp), mkdir(tmp); end
cleanTmp = onCleanup(@() cleanFolder(tmp));
fails = {};

% ---------------------------------------------------------------------------------
% Data (0.05 s, receiver chain on a zero sky signal)
% ---------------------------------------------------------------------------------
tD = tic;
zeroFile = fullfile(tmp, 'zero.dat');
writeZeros(zeroFile, round(0.05 * f_in));
rx = struct('f_in', f_in, 'fs', fs, 'fLO', fLO, 'band', [fLow fHigh], 'order', filterOrder);
radar = rfiSource('pulsed', 'Freq', 1300e6, 'PulseWidth', 2e-6, 'PRF', 373, 'ChirpBW', 1e6, ...
    'INRdB', 20, 'Label', 'radar');
side  = rfiSource('pulsed', 'Freq', 1300e6, 'PulseWidth', 2e-6, 'PRF', 373, 'ChirpBW', 1e6, ...
    'INRdB', 20, 'ScanPeriod', 10, 'BeamTime', 5, 'BeamWidth', 39e-3, 'SidelobeDB', -35, ...
    'Label', 'radar sidelobes -35 dB');
imp   = rfiSource('impulse', 'Rate', 1000, 'Duration', 200e-9, 'INRdB', 20, 'Label', 'impulses');
chan0 = makeChannels(zeroFile, fullfile(tmp, 'n'), struct([]), rx);
chanI = makeChannels(zeroFile, fullfile(tmp, 'i'), imp, rx);
chanR = makeChannels(zeroFile, fullfile(tmp, 'r'), radar, rx);
chanS = makeChannels(zeroFile, fullfile(tmp, 's'), side, rx);
delete(zeroFile);
fprintf('data: 4 sets of %d channels x %d samples (%.0f s)\n', chan0.nChan, chan0.N, toc(tD));
dArgs = {'Verbose', false};

% ---------------------------------------------------------------------------------
% 1. Noise only
% ---------------------------------------------------------------------------------
[~, d0] = detectRFI(chan0, dArgs{:});
[r0, p0] = periodicRFI(chan0, d0, dArgs{:});
pass = isempty(p0.emitters) && isempty(r0);
fprintf('1. noise only: %d detected pulses, %d emitters, %d mask rows: %s\n', p0.nPulses, ...
    numel(p0.emitters), size(r0, 1), passStr(pass));
if ~pass, fails{end+1} = 'noise'; end

% ---------------------------------------------------------------------------------
% 2. Impulses
% ---------------------------------------------------------------------------------
[~, dI] = detectRFI(chanI, dArgs{:});
[rI, pI] = periodicRFI(chanI, dI, dArgs{:});
pass = isempty(pI.emitters);
fprintf(['2. impulses (1000/s): %d detected pulses, %d emitters, %d mask rows; groups without a ' ...
         'period %d: %s\n'], pI.nPulses, numel(pI.emitters), size(rI, 1), size(pI.unfitted, 1), passStr(pass));
if ~pass, fails{end+1} = 'impulses'; end

% ---------------------------------------------------------------------------------
% 3. Scenario radar
% ---------------------------------------------------------------------------------
[mR, dR] = detectRFI(chanR, dArgs{:});
[rR, pR] = periodicRFI(chanR, dR, dArgs{:});
ok = numel(pR.emitters) == 1;
if ok
    E = pR.emitters(1);
    sPRF = E.periodErr / E.period * E.prf;
    okP = abs(E.prf - 373) < 4 * sPRF && abs(E.prf / 373 - 1) < 1e-5;
    tc = truthCentres(radar, chanR);                       % [samples], continuous
    pc = (E.t0 - chanR.t0) * chanR.fs + 1 + (-5:ceil(0.05 * E.prf) + 5).' * E.period * chanR.fs;
    dT = min(abs(tc - pc.'), [], 2);
    [fOwn, ~] = pulseStats(radar, mR, chanR);
    [fAll, ~] = pulseStats(radar, [mR; rR], chanR);
    ok = okP && all(dT <= pR.timeTol) && fAll == 1;
    fprintf(['3. radar: 1 emitter, PRF %.5f Hz (+- %.2g; true 373), %d of %d pulses fit; true pulse ' ...
             'centres vs predicted: max %.1f samples (tol %g); pulses blanked: detectRFI %.0f %%, ' ...
             '+ prediction %.0f %%; blanked by prediction %.4f %%: %s\n'], E.prf, sPRF, ...
        E.nInliers, E.nGroup, max(dT), pR.timeTol, 100*fOwn, 100*fAll, ...
        100*mean(pR.blankedFraction), passStr(ok));
else
    fprintf('3. radar: %d emitters found (expected 1): FAIL\n', numel(pR.emitters));
end
if ~ok, fails{end+1} = 'radar'; end

% ---------------------------------------------------------------------------------
% 4. Rotating radar, -35 dB sidelobes only
% ---------------------------------------------------------------------------------
[mS, dS] = detectRFI(chanS, dArgs{:});
[rS, pS] = periodicRFI(chanS, dS, dArgs{:});
ok = numel(pS.emitters) == 1 && abs(pS.emitters(1).prf / 373 - 1) < 1e-5;
[fOwn, eOwn] = pulseStats(side, mS, chanS);
[fAll, eAll] = pulseStats(side, [mS; rS], chanS);
[lOwn, lAll] = leftover(side, chanS, chan0, mS, [mS; rS]);
ok = ok && fAll == 1 && lAll < 1e-4;
nE = numel(pS.emitters); prfS = NaN; if nE, prfS = pS.emitters(1).prf; end
fprintf(['4. -35 dB sidelobes: %d emitter(s), PRF %.5f Hz; pulses blanked detectRFI %.1f %% (energy ' ...
         'missed %.2f), + prediction %.1f %% (%.2f); radar energy left, all channels: %.1e -> %.1e; ' ...
         'blanked by prediction %.4f %%: %s\n'], nE, prfS, 100*fOwn, eOwn, 100*fAll, eAll, lOwn, lAll, ...
    100*mean(pS.blankedFraction), passStr(ok));
if ~ok, fails{end+1} = 'sidelobes'; end

% ---------------------------------------------------------------------------------
if isempty(fails)
    fprintf('testPeriodicRFI: ALL PASSED\n');
else
    error('testPeriodicRFI:failed', 'FAILED: %s', strjoin(fails, ', '));
end
end


% =================================================================================
function info_chan = makeChannels(zeroFile, base, rfi, rx)
% Receiver chain on a zero sky signal: noise (seed 11, std 1) + rfi -> IQ -> channels.
info_rx = addNoiseAndRFI(zeroFile, [base '_rx.dat'], rx.f_in, 'Band', rx.band, 'NoiseStd', 1, ...
    'RFI', rfi, 'Seed', 11, 'Verbose', false);
info_IQ = applyIQmodulation(info_rx.file, [base '_iq.dat'], rx.f_in, rx.fs, rx.fLO, ...
    'FilterOrder', rx.order, 'Band', rx.band, 'Verbose', false);
delete(info_rx.file);
info_chan = channelizeIQ(info_IQ.file, [base '_chan'], info_IQ.actualFsOut, info_IQ.fLO, ...
    rx.band(1), rx.band(2), 'T0', info_IQ.t0, 'Verbose', false);
delete(info_IQ.file);
end

function tc = truthCentres(d, info_chan)
% True pulse centres of a pulsed source in channel samples (1-based, continuous).
fsC = info_chan.fs; Nc = info_chan.N;
tk = d.startTime + (0 : floor(Nc / fsC * d.prf)).' / d.prf;
tc = (tk + d.pulseWidth/2 - info_chan.t0) * fsC + 1;
tc = tc(tc >= 1 & tc <= Nc);
end

function [frac, missedE] = pulseStats(d, mask, info_chan)
% Fraction of the pulses blanked in all the source's own channels (within ChanWidth/2 +
% ChirpBW/2 of its frequency), and the fraction of its energy (antenna gain) in the others.
fsC = info_chan.fs; Nc = info_chan.N;
ch = find(abs(info_chan.chanFreqs - d.freq) < info_chan.chanWidth/2 + d.chirpBW/2).';
sOf = @(t) round((t - info_chan.t0) * fsC) + 1;
tk = d.startTime + (0 : floor(Nc / fsC * d.prf)).' / d.prf;
tk = tk(sOf(tk) >= 1 & sOf(tk + d.pulseWidth) <= Nc);
g = gainOf(d, tk + d.pulseWidth/2);
hit = true(size(tk));
for j = ch
    cv = coveredOf(mask, j, Nc);
    for i = 1:numel(tk)
        hit(i) = hit(i) && all(cv(sOf(tk(i)) : sOf(tk(i) + d.pulseWidth)));
    end
end
frac = mean(hit);
missedE = sum(g(~hit)) / sum(g);
end

function [lA, lB] = leftover(d, chanX, chan0, maskA, maskB)
% Fraction of the source's energy (RFI alone = chanX - chan0) in +-20 us windows around
% its pulses, all channels, that masks A and B leave.
fsC = chanX.fs; Nc = chanX.N;
tc = truthCentres(d, chanX);
hw = round(20e-6 * fsC);
inW = false(1, Nc);
for i = 1:numel(tc), inW(max(1, round(tc(i)) - hw) : min(Nc, round(tc(i)) + hw)) = true; end
eT = 0; eA = 0; eB = 0;
for j = 1:chanX.nChan
    Pr = abs(readIQ(chanX.chanFiles(j)) - readIQ(chan0.chanFiles(j))).^2;
    eT = eT + sum(Pr(inW));
    eA = eA + sum(Pr(inW & ~coveredOf(maskA, j, Nc)));
    eB = eB + sum(Pr(inW & ~coveredOf(maskB, j, Nc)));
end
lA = eA / eT; lB = eB / eT;
end

function g = gainOf(d, t)
% Antenna gain of a source at times t (1 without rotation).
g = ones(size(t));
if isfinite(d.scanPeriod)
    T = d.scanPeriod;
    dt = mod(t - d.beamTime + T/2, T) - T/2;
    g = max(exp(-4*log(2) * (dt / d.beamWidth).^2), 10^(d.sidelobeDB/10));
end
end

function cv = coveredOf(mask, j, N)
r = mask(mask(:, 1) == j, 2:3);
r = [max(r(:, 1), 1), min(r(:, 2), N)];
r = r(r(:, 1) <= r(:, 2), :);
d = accumarray([r(:, 1); r(:, 2) + 1], [ones(size(r, 1), 1); -ones(size(r, 1), 1)], [N + 1, 1]);
cv = cumsum(d(1:N)).' > 0;
end

function writeZeros(file, N)
fid = fopen(file, 'w', 'ieee-le');
c = onCleanup(@() fclose(fid));
blk = 4e6;
for n0 = 0:blk:N - 1
    fwrite(fid, zeros(1, min(blk, N - n0), 'single'), 'single');
end
end

function z = readIQ(file)
fid = fopen(file, 'r', 'ieee-le');
c = onCleanup(@() fclose(fid));
x = fread(fid, [2 Inf], 'single=>double');
z = complex(x(1, :), x(2, :));
end

function s = passStr(p)
if p, s = 'PASS'; else, s = 'FAIL'; end
end

function cleanFolder(d)
f = dir(fullfile(d, '*'));
for i = 1:numel(f)
    if ~f(i).isdir, delete(fullfile(d, f(i).name)); end
end
end
