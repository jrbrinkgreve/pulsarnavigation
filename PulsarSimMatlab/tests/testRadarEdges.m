function testRadarEdges()
%TESTRADAREDGES  Unit tests for B5a: radar pulses with raised-cosine edges ('RiseTime').
%{
Run from the PulsarSimMatlab folder: run('tests/testRadarEdges.m'). ~5 s.
Needs data/mc/rfiRef_preB5a.mat (tests/makeRfiReference.m, made with the
code of f4181f5, before B5a).

  1. Regression: all five RFI types with 'RiseTime' 0 for the radar give the
     reference output (noise + RFI, 20 us at 4 GHz) bit for bit.
  2. Defaults: rfiSource 'pulsed' gets RiseTime 0.1e-6, the other types 0;
     RiseTime > PulseWidth is refused.
  3. Envelope (RiseTime 0.1 us): the radar alone (noise + radar minus noise
     only, same seed) equals an independent formula: amplitude x raised-
     cosine edges centred on the nominal edges (50 % there) x the chirped
     carrier; pulse energy = (amp^2/2) (PulseWidth - RiseTime/4).
  4. Spectrum of one 2 us pulse, fraction of its energy farther than F from
     the carrier, vs the edge model (two step edges, each smoothed by the
     half-sine kernel H = cos(pi f T)/(1 - 4 f^2 T^2) of the ramp): instant
     edges 1/(pi^2 F PulseWidth) (power ~ 1/f^2); raised-cosine edges
     int_F |H|^2/(pi^2 f^2) df / PulseWidth (~ 1/(160 pi^2 F^5 T^4 PulseWidth),
     power ~ 1/f^6, for F >> 1/(2 T) = 5 MHz); within 20 % at F = 10, 25,
     50 MHz (the chirp moves the edges by +-0.5 MHz: ~7 % at 10 MHz).
Errors at the end if any check fails.
%}

root = fullfile(fileparts(mfilename('fullpath')), '..');
addpath(fullfile(root, 'functions'));
oldDir = cd(root); restoreDir = onCleanup(@() cd(oldDir));
pipelineParams;
tmp = fullfile(tempdir, 'testRadarEdges');
if ~isfolder(tmp), mkdir(tmp); end
cleanTmp = onCleanup(@() cleanFolder(tmp));
fails = {};

% ---------------------------------------------------------------------------------
% 1. Regression: RiseTime 0 = the code before B5a
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
fprintf('1. RiseTime 0, all five RFI types: output = pre-B5a reference bit for bit: %s\n', passStr(pass));
if ~pass, fails{end+1} = 'regression'; end

% ---------------------------------------------------------------------------------
% 2. Defaults
% ---------------------------------------------------------------------------------
sP = rfiSource('pulsed', 'Freq', 1300e6, 'PulseWidth', 2e-6, 'PRF', 373);
sC = rfiSource('cw', 'Freq', 1351.3e6);
refused = false;
try
    rfiSource('pulsed', 'Freq', 1300e6, 'PulseWidth', 2e-6, 'PRF', 373, 'RiseTime', 3e-6);
catch
    refused = true;
end
pass = sP.riseTime == 0.1e-6 && sC.riseTime == 0 && refused;
fprintf('2. defaults: pulsed RiseTime %.2g s, cw %g; RiseTime > PulseWidth refused %d: %s\n', ...
    sP.riseTime, sC.riseTime, refused, passStr(pass));
if ~pass, fails{end+1} = 'defaults'; end

% ---------------------------------------------------------------------------------
% 3. Envelope (radar alone = noise + radar minus noise only)
% ---------------------------------------------------------------------------------
fsR = ref.fs; N = 40000; pw = 2e-6; tr = 0.1e-6; t0 = 2e-6; fR = 1300e6; bw = 1e6;
rad = rfiSource('pulsed', 'Freq', fR, 'PulseWidth', pw, 'PRF', 373, 'ChirpBW', bw, ...
    'StartTime', t0, 'INRdB', 20, 'RiseTime', tr);
x = double(runRx(tmp, N, fsR, ref.band, rad, 5)) - double(runRx(tmp, N, fsR, ref.band, struct([]), 5));
amp = sqrt(2 * 10^(20/10) * 1^2 * diff(ref.band) / (fsR/2));   % as addNoiseAndRFI, NoiseStd 1
nn  = 0:N-1; u = nn / fsR - t0;                      % time since the nominal start
env = zeros(1, N);
env(u >= tr/2 & u <= pw - tr/2) = 1;
a = u > -tr/2 & u < tr/2;        env(a) = 0.5 * (1 - cos(pi * (u(a) + tr/2) / tr));
b = u > pw - tr/2 & u < pw + tr/2; env(b) = 0.5 * (1 + cos(pi * (u(b) - pw + tr/2) / tr));
xe = amp * env .* cos(2*pi * (nn * fR / fsR + 0.5 * (bw/pw) * (u - pw/2).^2));
dMax = max(abs(x - xe)) / amp;
E = sum(x.^2) / fsR; Ep = amp^2/2 * (pw - tr/4);
half = interp1(u(u > -tr & u < tr), env(u > -tr & u < tr), 0);
pass = dMax < 1e-5 && abs(E / Ep - 1) < 2e-3;
fprintf(['3. RiseTime 0.1 us: radar = formula to %.1e of the amplitude (float32 rounding); envelope ' ...
         'at the nominal edge %.3f; energy / (amp^2/2 (PW - RiseTime/4)) %.5f: %s\n'], dMax, half, ...
    E / Ep, passStr(pass));
if ~pass, fails{end+1} = 'envelope'; end

% ---------------------------------------------------------------------------------
% 4. Spectrum: energy far from the carrier
% ---------------------------------------------------------------------------------
F = [10 25 50] * 1e6;
frac = zeros(2, numel(F)); pred = frac;
trs = [0 tr];
for i = 1:2
    s = rfiSource('pulsed', 'Freq', fR, 'PulseWidth', pw, 'PRF', 373, 'ChirpBW', bw, ...
        'StartTime', t0, 'INRdB', 20, 'RiseTime', trs(i));
    xr = double(runRx(tmp, N, fsR, ref.band, s, 5)) - double(runRx(tmp, N, fsR, ref.band, struct([]), 5));
    X  = abs(fft(xr)).^2;
    f  = (0:N-1) * fsR / N;
    pos = f > 0 & f < fsR/2;                         % one side (real signal)
    for k = 1:numel(F)
        frac(i, k) = sum(X(pos & abs(f - fR) > F(k))) / sum(X(pos));
    end
end
% two edges, each a step (1/(2 pi f)) smoothed by the half-sine kernel of the
% raised-cosine ramp, H = cos(pi f T) / (1 - 4 f^2 T^2) (H = 1 for T = 0); the
% interference of the two edges averages to 1/2 over the 1/PW ripple
H2 = @(f, T) (cos(pi * f * T) ./ (1 - 4 * f.^2 * T^2)).^2;
% (on a 10 kHz grid to 1 GHz: integral() to Inf fails on the oscillating cos^2)
for k = 1:numel(F)
    pred(1, k) = 1 / (pi^2 * F(k) * pw);
    fg = F(k):1e4:1e9;
    pred(2, k) = trapz(fg, H2(fg, tr) ./ (pi^2 * fg.^2)) / pw;
end
r = frac ./ pred;
pass = all(abs(r(:) - 1) < 0.2);
fprintf(['4. energy fraction beyond %s MHz of the carrier: instant edges %s (predicted %s); ' ...
         'RiseTime 0.1 us %s (predicted %s): %s\n'], mat2str(F/1e6), mat2str(frac(1, :), 2), ...
    mat2str(pred(1, :), 2), mat2str(frac(2, :), 2), mat2str(pred(2, :), 2), passStr(pass));
if ~pass, fails{end+1} = 'spectrum'; end

% ---------------------------------------------------------------------------------
if isempty(fails)
    fprintf('testRadarEdges: ALL PASSED\n');
else
    error('testRadarEdges:failed', 'FAILED: %s', strjoin(fails, ', '));
end
end


% =================================================================================
function y = runRx(tmp, N, fs, band, src, seed)
% addNoiseAndRFI on N zero samples, receiver noise std 1: the output samples.
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
