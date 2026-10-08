function testRfiGating()
%TESTRFIGATING  Unit tests for B5c-1: rotating antenna (rfiSource 'ScanPeriod' etc.).
%{
Run from the PulsarSimMatlab folder: run('tests/testRfiGating.m'). ~10 s.
Needs data/mc/rfiRef_preB5a.mat (tests/makeRfiReference.m).

The source power is multiplied by g(t) = max(exp(-4 ln2 (dt/BeamWidth)^2),
10^(SidelobeDB/10)), dt = time from the nearest beam passage (BeamTime +
k ScanPeriod). A fast toy scan is used so that two passages fit in 2.5 ms
of 4 GHz data: ScanPeriod 1 ms, BeamTime 0.6 ms, BeamWidth 100 us, sidelobes
-25 dB.

  1. Regression: no rotation (default) gives the pre-B5a reference output
     (all five RFI types, radar 'RiseTime' 0) bit for bit.
  2. Defaults and checks: ScanPeriod Inf by default; a rotation without
     BeamWidth, BeamWidth >= ScanPeriod, SidelobeDB > 0 are refused.
  3. Radar (PRF 20 kHz, 0.5 us pulses, 50 pulses): the energy of every pulse
     with rotation / without = g averaged over that pulse, to 1e-4 — in
     both main beams, on the slopes and on the -25 dB floor.
  4. Carrier: the local mean power (100 ns windows, ~135 carrier cycles) with
     rotation / without = g averaged over the window, to 2e-3.
Errors at the end if any check fails.
%}

root = fullfile(fileparts(mfilename('fullpath')), '..');
addpath(fullfile(root, 'functions'));
oldDir = cd(root); restoreDir = onCleanup(@() cd(oldDir));
pipelineParams;
tmp = fullfile(tempdir, 'testRfiGating');
if ~isfolder(tmp), mkdir(tmp); end
cleanTmp = onCleanup(@() cleanFolder(tmp));
fails = {};

% ---------------------------------------------------------------------------------
% 1. Regression: no rotation = the code before B5a / B5c
% ---------------------------------------------------------------------------------
R = load(fullfile(dataDir, 'mc', 'rfiRef_preB5a.mat')); ref = R.ref;
src = [ ...                                          % as tests/makeRfiReference.m
    rfiSource('bpsk',    'Freq', 1575.42e6, 'ChipRate', 1.023e6, 'INRdB', 0), ...
    rfiSource('bpsk',    'Freq', 1227.60e6, 'ChipRate', 10.23e6, 'INRdB', 0), ...
    rfiSource('pulsed',  'Freq', 1300e6, 'PulseWidth', 2e-6, 'PRF', 373, 'ChirpBW', 1e6, ...
              'StartTime', 5e-6, 'INRdB', 20, 'RiseTime', 0), ...
    rfiSource('cw',      'Freq', 1351.3e6, 'Drift', 1e9, 'INRdB', 0), ...
    rfiSource('impulse', 'Rate', 2e5, 'Duration', 200e-9, 'INRdB', 20)];
y = runRx(tmp, ref.N, ref.fs, ref.band, src, ref.seed);
pass = isequal(y(:), ref.y(:));
fprintf('1. no rotation (default), all five RFI types: output = pre-B5a reference bit for bit: %s\n', ...
    passStr(pass));
if ~pass, fails{end+1} = 'regression'; end

% ---------------------------------------------------------------------------------
% 2. Defaults and checks
% ---------------------------------------------------------------------------------
sD = rfiSource('cw', 'Freq', 1351.3e6);
bad = { {'ScanPeriod', 1}, {'ScanPeriod', 1, 'BeamWidth', 1}, ...
        {'ScanPeriod', 1, 'BeamWidth', 0.1, 'SidelobeDB', 3} };
refused = 0;
for i = 1:numel(bad)
    try
        rfiSource('cw', 'Freq', 1351.3e6, bad{i}{:});
    catch
        refused = refused + 1;
    end
end
pass = isinf(sD.scanPeriod) && refused == numel(bad);
fprintf('2. default ScanPeriod %g; refused %d of %d bad rotations: %s\n', sD.scanPeriod, refused, ...
    numel(bad), passStr(pass));
if ~pass, fails{end+1} = 'defaults'; end

% ---------------------------------------------------------------------------------
% 3. Radar pulses through two beam passages
% ---------------------------------------------------------------------------------
fsR = ref.fs; N = round(2.5e-3 * fsR); seed = 7;
scan = {'ScanPeriod', 1e-3, 'BeamTime', 0.6e-3, 'BeamWidth', 100e-6, 'SidelobeDB', -25};
gOf  = @(t) max(exp(-4*log(2) * ((mod(t - 0.6e-3 + 0.5e-3, 1e-3) - 0.5e-3) / 100e-6).^2), 10^(-25/10));
pw = 0.5e-6; prf = 20e3; t0 = 10e-6;
radar = @(rot) rfiSource('pulsed', 'Freq', 1300e6, 'PulseWidth', pw, 'PRF', prf, 'StartTime', t0, ...
    'INRdB', 20, 'RiseTime', 0, rot{:});
y0   = double(runRx(tmp, N, fsR, ref.band, struct([]), seed));
xFix = double(runRx(tmp, N, fsR, ref.band, radar({}), seed)) - y0;
xRot = double(runRx(tmp, N, fsR, ref.band, radar(scan), seed)) - y0;
tS = (0:N-1) / fsR;
tk = t0 + (0:floor((tS(end) - t0 - pw) * prf)) / prf;     % pulse starts
ratio = zeros(size(tk)); gk = ratio;
for k = 1:numel(tk)
    in = tS >= tk(k) & tS < tk(k) + pw;
    ratio(k) = sum(xRot(in).^2) / sum(xFix(in).^2);
    gk(k) = mean(gOf(tS(in)));
end
dR = max(abs(ratio ./ gk - 1));
[~, kHalf] = min(abs(abs(mod(tk + pw/2 - 0.6e-3 + 0.5e-3, 1e-3) - 0.5e-3) - 50e-6));
pass = dR < 1e-4;
dtHalf = mod(tk(kHalf) + pw/2 - 0.6e-3 + 0.5e-3, 1e-3) - 0.5e-3;
fprintf(['3. radar, %d pulses: energy with / without rotation = mean g over the pulse to %.1e; ' ...
         'peak %.4f (two passages); pulse nearest the -3 dB point, %+.1f us from the beam centre: ' ...
         '%.3f (g %.3f); floor %.5f (10^-2.5 = %.5f): %s\n'], numel(tk), dR, max(ratio), ...
    dtHalf*1e6, ratio(kHalf), gk(kHalf), min(ratio), 10^-2.5, passStr(pass));
if ~pass, fails{end+1} = 'radar envelope'; end

% ---------------------------------------------------------------------------------
% 4. Carrier
% ---------------------------------------------------------------------------------
cw = @(rot) rfiSource('cw', 'Freq', 1351.3e6, 'INRdB', 0, rot{:});
cFix = double(runRx(tmp, N, fsR, ref.band, cw({}), seed)) - y0;
cRot = double(runRx(tmp, N, fsR, ref.band, cw(scan), seed)) - y0;
w  = 400;                                            % 100 ns windows
nw = floor(N / w);
pF = mean(reshape(cFix(1:nw*w).^2, w, nw), 1);
pR = mean(reshape(cRot(1:nw*w).^2, w, nw), 1);
gw = mean(reshape(gOf(tS(1:nw*w)), w, nw), 1);
dC = max(abs((pR ./ pF) ./ gw - 1));
pass = dC < 2e-3;
fprintf('4. carrier: local power with / without rotation = mean g over 100 ns to %.1e: %s\n', ...
    dC, passStr(pass));
if ~pass, fails{end+1} = 'carrier envelope'; end

% ---------------------------------------------------------------------------------
if isempty(fails)
    fprintf('testRfiGating: ALL PASSED\n');
else
    error('testRfiGating:failed', 'FAILED: %s', strjoin(fails, ', '));
end
end


% =================================================================================
function y = runRx(tmp, N, fs, band, src, seed)
% addNoiseAndRFI on N zero samples, receiver noise std 1: the output samples (row).
zeroFile = fullfile(tmp, 'zero.dat'); outFile = fullfile(tmp, 'rx.dat');
fid = fopen(zeroFile, 'w', 'ieee-le'); fwrite(fid, zeros(1, N, 'single'), 'single'); fclose(fid);
addNoiseAndRFI(zeroFile, outFile, fs, 'Band', band, 'NoiseStd', 1, 'RFI', src, 'Seed', seed, ...
    'SaveInfo', false, 'Verbose', false);
fid = fopen(outFile, 'r', 'ieee-le'); y = fread(fid, Inf, 'single=>single'); fclose(fid);
y = y.';
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
