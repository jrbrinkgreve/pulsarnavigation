function testRfiNoise()
%TESTRFINOISE  Unit tests for B5d: the 'noise' RFI type (band-limited Gaussian, LTE-like).
%{
Run from the PulsarSimMatlab folder: run('tests/testRfiNoise.m'). ~20 s.
Needs data/mc/rfiRef_preB5a.mat (tests/makeRfiReference.m).
The RFI alone = (noise + RFI) - (noise only), same seed. Source: 20 MHz at
1472 MHz (inside LTE band 32), INR 0 dB, 2.5 ms of 4 GHz data.

  1. Regression: the other RFI types (no 'noise' source) give the pre-B5a
     reference output bit for bit.
  2. Block size: BlockSize 1e6+3 (odd, not a multiple of the 65,536-sample
     random substreams) = one block, bit for bit.
  3. Power = INR x the receiver-noise power in the band, within 4 sigma of
     1/sqrt(Bandwidth x T) (the number of independent samples).
  4. Spectrum (averaged periodogram, 244 kHz bins, Hann): flat in the band
     (|f - Freq| < 0.45 Bandwidth: level = power / Bandwidth to 2 %, bin
     scatter as expected for the averaging); out of band (> Bandwidth/2 +
     2 MHz) below 1e-5 of the in-band level.
  5. Gaussian: kurtosis 3 within 4 sigma (sqrt(24 / (2 Bandwidth T))).
Errors at the end if any check fails.
%}

root = fullfile(fileparts(mfilename('fullpath')), '..');
addpath(fullfile(root, 'functions'));
oldDir = cd(root); restoreDir = onCleanup(@() cd(oldDir));
pipelineParams;
tmp = fullfile(tempdir, 'testRfiNoise');
if ~isfolder(tmp), mkdir(tmp); end
cleanTmp = onCleanup(@() cleanFolder(tmp));
fails = {};

% ---------------------------------------------------------------------------------
% 1. Regression
% ---------------------------------------------------------------------------------
R = load(fullfile(dataDir, 'mc', 'rfiRef_preB5a.mat')); ref = R.ref;
src = [ ...                                          % as tests/makeRfiReference.m
    rfiSource('bpsk',    'Freq', 1575.42e6, 'ChipRate', 1.023e6, 'INRdB', 0), ...
    rfiSource('bpsk',    'Freq', 1227.60e6, 'ChipRate', 10.23e6, 'INRdB', 0), ...
    rfiSource('pulsed',  'Freq', 1300e6, 'PulseWidth', 2e-6, 'PRF', 373, 'ChirpBW', 1e6, ...
              'StartTime', 5e-6, 'INRdB', 20, 'RiseTime', 0), ...
    rfiSource('cw',      'Freq', 1351.3e6, 'Drift', 1e9, 'INRdB', 0), ...
    rfiSource('impulse', 'Rate', 2e5, 'Duration', 200e-9, 'INRdB', 20)];
y = runRx(tmp, ref.N, ref.fs, ref.band, src, ref.seed, 4e6);
pass = isequal(y(:), ref.y(:));
fprintf('1. other RFI types: output = pre-B5a reference bit for bit: %s\n', passStr(pass));
if ~pass, fails{end+1} = 'regression'; end

% ---------------------------------------------------------------------------------
% 2.-5. The noise source
% ---------------------------------------------------------------------------------
fsR = ref.fs; T = 2.5e-3; N = round(T * fsR); seed = 9;
fc = 1472e6; bw = 20e6;
lte = rfiSource('noise', 'Freq', fc, 'Bandwidth', bw, 'INRdB', 0, 'Label', 'LTE-like');
yA = runRx(tmp, N, fsR, ref.band, lte, seed, 4e7);            % one block
yB = runRx(tmp, N, fsR, ref.band, lte, seed, 1e6 + 3);        % odd blocks
pass = isequal(yA, yB);
fprintf('2. BlockSize 1e6+3 vs one block: identical %d: %s\n', pass, passStr(pass));
if ~pass, fails{end+1} = 'block size'; end

y0 = double(runRx(tmp, N, fsR, ref.band, struct([]), seed, 4e7));
x  = double(yA) - y0;                                         % the RFI alone
PnBand = 1^2 * diff(ref.band) / (fsR/2);                      % NoiseStd 1
Pm = mean(x.^2); Pp = 10^(0/10) * PnBand;
sP = 1 / sqrt(bw * T);
pass = abs(Pm / Pp - 1) < 4 * sP;
fprintf('3. power / (INR x noise power in the band) %.4f (sigma %.4f): %s\n', Pm / Pp, sP, passStr(pass));
if ~pass, fails{end+1} = 'power'; end

nS = 2^14; nSeg = floor(N / nS);
win = 0.5 - 0.5 * cos(2*pi*(0:nS-1).' / nS);
Xs = reshape(x(1:nS*nSeg), nS, nSeg) .* win;
S  = mean(abs(fft(Xs)).^2, 2) / sum(win.^2) / fsR;           % two-sided PSD [power/Hz]
f  = (0:nS-1).' * fsR / nS;
inB  = abs(f - fc) < 0.45 * bw;
outB = abs(f - fc) > bw/2 + 2e6 & f < fsR/2 & f > 0;
lvl  = 2 * mean(S(inB));                                      % one-sided level
flat = abs(lvl / (Pm / bw) - 1);
scat = std(S(inB)) / mean(S(inB));
out  = max(S(outB)) / mean(S(inB));
pass = flat < 0.02 && scat < 3 / sqrt(nSeg) && out < 1e-5;
fprintf(['4. spectrum: in-band level / (power / Bandwidth) %.4f, bin scatter %.3f (averaging of %d ' ...
         'segments: %.3f), max out of band %.1e of in-band: %s\n'], lvl / (Pm / bw), scat, nSeg, ...
    1 / sqrt(nSeg), out, passStr(pass));
if ~pass, fails{end+1} = 'spectrum'; end

kurt = mean(x.^4) / mean(x.^2)^2;
sK = sqrt(24 / (2 * bw * T));
pass = abs(kurt - 3) < 4 * sK;
fprintf('5. kurtosis %.4f (Gaussian 3, sigma %.4f): %s\n', kurt, sK, passStr(pass));
if ~pass, fails{end+1} = 'Gaussian'; end

% ---------------------------------------------------------------------------------
if isempty(fails)
    fprintf('testRfiNoise: ALL PASSED\n');
else
    error('testRfiNoise:failed', 'FAILED: %s', strjoin(fails, ', '));
end
end


% =================================================================================
function y = runRx(tmp, N, fs, band, src, seed, blk)
% addNoiseAndRFI on N zero samples, receiver noise std 1, block size blk: output (row).
zeroFile = fullfile(tmp, 'zero.dat'); outFile = fullfile(tmp, 'rx.dat');
fid = fopen(zeroFile, 'w', 'ieee-le'); fwrite(fid, zeros(1, N, 'single'), 'single'); fclose(fid);
addNoiseAndRFI(zeroFile, outFile, fs, 'Band', band, 'NoiseStd', 1, 'RFI', src, 'Seed', seed, ...
    'BlockSize', blk, 'SaveInfo', false, 'Verbose', false);
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
