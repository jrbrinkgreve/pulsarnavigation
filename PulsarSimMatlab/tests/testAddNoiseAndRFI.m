function testAddNoiseAndRFI()
%TESTADDNOISEANDRFI  Unit tests for addNoiseAndRFI (stage 3: receiver noise and RFI at RF).
%{
Run from the PulsarSimMatlab folder: run('tests/testAddNoiseAndRFI.m'). ~4 s.
Writes its own files to tempdir (~0.1 GB, deleted at the end); needs no data files.
Ground truth is the help texts' definitions, computed here:
  SNR = A^2 / sigma_n^2 (pulse-peak PSD over noise PSD, white signal and noise at f_in),
  P_band = sigma_n^2 * B / (fs/2) (the band's share of white noise),
  INR = RFI power while on / P_band: carrier amplitude sqrt(2*INR*P_band); a white burst
  amp sqrt(INR)*sigma_n, of which the band holds the same share B/(fs/2) as of the noise.
The RFI alone = (noise + RFI) - (noise only), same seed. The band 1.2-1.6 GHz at 4 GHz,
1 ms files (4,000,000 samples) unless noted.

  1. Format and info: N, 4*N bytes, format fields, band, snrDB, noiseStd = A/sqrt(rho),
     noisePSD = sigma_n^2/fs, noisePowerBand = P_band, seed, no RFI; the saved _info.mat =
     the returned info. 'NoiseStd' gives snrDB = 10 log10(A^2/sigma^2). SNRdB Inf, or no
     noise option at all: the output is the input, bit for bit.
  2. The SNR definition (input: white Gaussian of std A = 2, a pulsar held at its peak;
     SNRdB -5): the noise's variance = sigma_n^2 within 4 sigma (sqrt(2/N)); the in-band
     power ratio input / noise = rho within 4 sigma (sqrt(2/K), K bins in the band); the
     noise is white (8 sub-bands of 250 MHz, each within 4 sigma); the noise does not depend
     on the input (zero input and signal input give the same noise to float32 rounding).
  3. RFI power against INR (receiver noise std 1), each type alone:
     cw (+10 dB, 1351.3 MHz, phase 0.7, drift 1 GHz/s) = amp*cos(2*pi*(f*n/fs + drift*t^2/2)
     + phase) sample for sample (1e-6 of amp), mean power INR*P_band within 1e-4;
     bpsk (-5 dB, 1575.42 MHz, 1.023 Mchip/s): RFI / (amp*carrier) = +-1 (1e-4) and constant
     within each chip, chips balanced within 4 sigma, power INR*P_band within 1e-3;
     impulse (+10 dB, 1e5 bursts/s, 200 ns): nonzero exactly in the reported bursts (800
     samples each), burst count = rate * T within 4 sigma (Poisson), energy = INR*sigma_n^2
     x the burst coverage within 4 sigma, in-band share = the share of bins in the band
     (0.2) within 4 sigma. Not repeated: 'pulsed' (testRadarEdges 3, amplitude from INR) and
     'noise' (testRfiNoise 3).
  4. Block size (2.5 ms, noise + all five types, one of them with a rotating antenna, a
     30 Mchip/s bpsk that crosses a 65,536-chip substream): BlockSize 4e6 (default, 3
     blocks), 1,000,003 (10 blocks) and one block give the same file, bit for bit.
  5. Prediction (info.prediction, reference pulse: 0.5 ms power FWHM, 10 pulses): against
     adaptive quadrature of the help text's radiometer integrals within 1e-6, at -54 dB and
     0 dB (self-noise), and in 'amplitude' mode (profile std sigma/sqrt(2)); at -54 dB also
     the weak-pulse forms SNR = rho*sqrt(B*sigma*sqrt(pi)) (= 388*rho) and sigma_TOA =
     sqrt(2*sigma/(B*sqrt(pi)))/rho within 1e-5 (self-noise moves them by ~3e-6 there);
     folded values = per pulse x / sqrt(10).
  6. Input checks: inFile = outFile, a band upside down or above fs/2, SNRdB without an
     amplitude, RFI without receiver noise, an RFI frequency outside (0, fs/2) and a
     missing file are refused, each with its own id.
Errors at the end if any check fails.
%}

root = fullfile(fileparts(mfilename('fullpath')), '..');
addpath(fullfile(root, 'functions'));
oldDir = cd(root); restoreDir = onCleanup(@() cd(oldDir));
tmp = tempname(tempdir);                         % own folder per run
mkdir(tmp);
cleanTmp = onCleanup(@() rmdir(tmp, 's'));
fails = {};

fs = 4e9; band = [1.2e9 1.6e9]; B = diff(band);  % the receiver band at the RF rate (testParams)
share = B / (fs/2);                              % the band's share of white noise: 0.2
N = 4e6;                                         % 1 ms
fSig = fullfile(tmp, 'sig.dat'); fZero = fullfile(tmp, 'zero.dat'); fOut = fullfile(tmp, 'out.dat');
rx = @(fin, varargin) addNoiseAndRFI(fin, fOut, fs, 'Band', band, 'Verbose', false, ...
    'SaveInfo', false, varargin{:});
A = 2;
rs = RandStream('mt19937ar', 'Seed', 5);
x0 = single(A * randn(rs, N, 1));                % white, std A: a pulsar held at its peak
writeF32(fSig, x0);
writeF32(fZero, zeros(N, 1, 'single'));

% ---------------------------------------------------------------------------------
% 1. Format and info
% ---------------------------------------------------------------------------------
rho = 10^(-5/10);
sig = A / sqrt(rho);                             % sigma_n from the SNR definition
i1 = addNoiseAndRFI(fSig, fOut, fs, 'Band', band, 'SNRdB', -5, 'SignalAmp', A, 'Seed', 3, ...
    'Verbose', false);
y1 = readF32(fOut);
f = dir(fOut);
S = load(fullfile(tmp, 'out_info.mat'));
ok = [i1.N == N, f.bytes == 4*N, strcmp(i1.precision, 'single'), strcmp(i1.byteOrder, 'ieee-le'), ...
      ~i1.isComplex, i1.fs == fs, i1.actualFsOut == fs, isequal(i1.band, band), i1.snrDB == -5, ...
      abs(i1.noiseStd / sig - 1) < 1e-14, abs(i1.noisePSD / (sig^2 / fs) - 1) < 1e-14, ...
      abs(i1.noisePowerBand / (sig^2 * share) - 1) < 1e-14, i1.seed == 3, isempty(i1.rfi), ...
      isequal(rmfield(i1, 'infoFile'), S.info)];
iS = rx(fSig, 'NoiseStd', 3, 'SignalAmp', A);
ok(end+1) = iS.noiseStd == 3 && abs(iS.snrDB - 10*log10(A^2 / 9)) < 1e-12;
iInf = rx(fSig, 'SNRdB', Inf, 'SignalAmp', A); yInf = readF32(fOut);
iNone = rx(fSig); yNone = readF32(fOut);
ok(end+1) = isequal(yInf, x0) && isequal(yNone, x0) && iInf.noiseStd == 0 && iNone.noiseStd == 0 && ...
            iNone.snrDB == Inf;
pass = all(ok);
fprintf(['1. format and info: N %d, %d bytes, sigma_n %.6f (A/sqrt(rho) %.6f), P_band %.6f; NoiseStd 3 -> ' ...
         'snrDB %.4f; SNRdB Inf and no noise give the input bit for bit: %d of %d checks: %s\n'], ...
    i1.N, f.bytes, i1.noiseStd, sig, i1.noisePowerBand, iS.snrDB, nnz(ok), numel(ok), passStr(pass));
if ~pass, fails{end+1} = 'format / info'; end

% ---------------------------------------------------------------------------------
% 2. The SNR definition
% ---------------------------------------------------------------------------------
rx(fZero, 'SNRdB', -5, 'SignalAmp', A, 'Seed', 3);
n0 = readF32(fOut);                              % the noise alone (zero input)
n1 = double(y1) - double(x0);                    % the noise of the signal run
okIndep = all(abs(n1 - double(n0)) <= double(eps(y1)) / 2);   % y1 = x0 + noise, rounded once
zV = (mean(double(n0).^2) / sig^2 - 1) / sqrt(2/N);
X  = fft(double(x0));
Nn = fft(double(n0));
fk = (0:N/2).' * (fs/N);
kb = find(fk >= band(1) & fk <= band(2));        % the bins in the band (1-based)
K  = numel(kb);
ratio = sum(abs(X(kb)).^2) / sum(abs(Nn(kb)).^2);
zR = (ratio / rho - 1) / sqrt(2/K);
e  = round(linspace(1, N/2, 9));                 % 8 sub-bands of 0..fs/2, without DC and Nyquist
P2 = abs(Nn).^2;
zW = zeros(1, 8);
for j = 1:8
    idx = (e(j):e(j+1) - 1) + 1;
    zW(j) = (sum(P2(idx)) / (numel(idx) * N * sig^2) - 1) * sqrt(numel(idx));   % E|N_k|^2 = N sigma^2
end
clear X Nn P2 n1 y1
pass = abs(zV) < 4 && abs(zR) < 4 && all(abs(zW) < 4) && okIndep;
fprintf(['2. SNR definition: noise variance / sigma_n^2 %.5f (%+.2f sigma); in-band power ratio ' ...
         'input / noise %.5f vs rho %.5f (%+.2f sigma, %d bins); 8 sub-bands within %.2f sigma; ' ...
         'noise independent of the input %d: %s\n'], mean(double(n0).^2) / sig^2, zV, ratio, rho, zR, ...
    K, max(abs(zW)), okIndep, passStr(pass));
if ~pass, fails{end+1} = 'SNR definition'; end

% ---------------------------------------------------------------------------------
% 3. RFI power against INR
% ---------------------------------------------------------------------------------
P = 1^2 * share;                                 % receiver noise std 1: P_band
base = {'NoiseStd', 1, 'Seed', 7};
rx(fZero, base{:});
yN = double(readF32(fOut));                      % noise only
n = (0:N-1).'; t = n / fs;

INRc = 10; fc = 1351.3e6; ph = 0.7; dr = 1e9;    % cw
rx(fZero, base{:}, 'RFI', rfiSource('cw', 'Freq', fc, 'INRdB', 10, 'Phase', ph, 'Drift', dr));
r = double(readF32(fOut)) - yN;
ampC = sqrt(2 * INRc * P);
dC = max(abs(r - ampC * cos(2*pi * (fc * n / fs + 0.5 * dr * t.^2) + ph))) / ampC;
pC = mean(r.^2) / (INRc * P);
okC = dC < 1e-6 && abs(pC - 1) < 1e-4;

INRb = 10^(-5/10); fb = 1575.42e6; cr = 1.023e6; % bpsk
rx(fZero, base{:}, 'RFI', rfiSource('bpsk', 'Freq', fb, 'ChipRate', cr, 'INRdB', -5));
r = double(readF32(fOut)) - yN;
ampB = sqrt(2 * INRb * P);
c = cos(2*pi * fb * n / fs);
sel = abs(c) > 0.1;                              % away from the carrier's zero crossings
q = r(sel) ./ (ampB * c(sel));                   % the chips
m = floor(n(sel) * (cr / fs));                   % chip index of each sample
sMin = accumarray(m + 1, sign(q), [], @min);
sMax = accumarray(m + 1, sign(q), [], @max);
nChip = numel(sMin);
dB = max(abs(abs(q) - 1));
zB = (mean(sMax > 0) - 0.5) / (0.5 / sqrt(nChip));
pB = mean(r.^2) / (INRb * P);
okB = dB < 1e-4 && isequal(sMin, sMax) && abs(zB) < 4 && abs(pB - 1) < 1e-3;

INRi = 10; rate = 1e5; dur = 200e-9;             % impulse
iI = rx(fZero, base{:}, 'RFI', rfiSource('impulse', 'Rate', rate, 'Duration', dur, 'INRdB', 10));
r = double(readF32(fOut)) - yN;
ev0 = iI.rfi.evStart(:); evL = iI.rfi.evLen(:);
cover = zeros(N, 1);                             % number of bursts covering each sample
for k = 1:numel(ev0)
    a = ev0(k) + 1; b = min(N, ev0(k) + evL(k));
    cover(a:b) = cover(a:b) + 1;
end
on = cover > 0;
nEv = numel(ev0);
ampI = sqrt(INRi) * 1;
zE = (sum(r.^2) / (ampI^2 * sum(cover)) - 1) / (sqrt(2 * sum(cover.^2)) / sum(cover));
R  = fft(r);
sh = 2 * sum(abs(R(kb)).^2) / N / sum(r.^2);     % in-band share of the burst power (Parseval)
s0 = 2 * K / N;                                  % white: the share of the bins in the band
zS = (sh - s0) / sqrt(2 * s0 * (1 - s0) / nnz(on));
clear R
okI = all(r(~on) == 0) && nnz(r(on)) >= 0.9999 * nnz(on) && all(evL == round(dur * fs)) && ...
      abs(nEv - rate * N / fs) < 4 * sqrt(rate * N / fs) && abs(zE) < 4 && abs(zS) < 4;
pass = okC && okB && okI;
fprintf(['3. RFI against INR: cw = formula to %.1e of amp, power / (INR P_band) %.6f: %s; bpsk chips ' ...
         '+-1 to %.1e, constant in all %d chips %d, +1 fraction %+.2f sigma, power %.5f: %s; impulse ' ...
         '%d bursts (rate x T = %g), zero outside %d, energy %+.2f sigma, in-band share %.4f (%.4f, ' ...
         '%+.2f sigma), in-band power while on / (INR P_band) %.4f: %s\n'], dC, pC, passStr(okC), dB, ...
    nChip, isequal(sMin, sMax), zB, pB, passStr(okB), nEv, rate * N / fs, all(r(~on) == 0), zE, sh, ...
    s0, zS, sh * sum(r.^2) / nnz(on) / (INRi * P), passStr(okI));
if ~pass, fails{end+1} = 'RFI against INR'; end
clear yN r c q m sel cover on t n

% ---------------------------------------------------------------------------------
% 4. Block size
% ---------------------------------------------------------------------------------
N4 = 1e7;                                        % 2.5 ms
f4 = fullfile(tmp, 'sig4.dat');
writeF32(f4, single(randn(rs, N4, 1)));
rfi4 = [rfiSource('cw', 'Freq', 1351.3e6, 'INRdB', 0, 'Drift', 1e9, 'ScanPeriod', 1e-3, ...
                  'BeamTime', 0.6e-3, 'BeamWidth', 1e-4, 'SidelobeDB', -25), ...
        rfiSource('bpsk', 'Freq', 1575.42e6, 'ChipRate', 1.023e6, 'INRdB', 0), ...
        rfiSource('bpsk', 'Freq', 1.3e9, 'ChipRate', 30e6, 'INRdB', 0), ...   % 75,000 chips
        rfiSource('pulsed', 'Freq', 1.25e9, 'PulseWidth', 2e-6, 'PRF', 2e4, 'ChirpBW', 1e6, 'INRdB', 20), ...
        rfiSource('impulse', 'Rate', 2e5, 'Duration', 200e-9, 'INRdB', 20), ...
        rfiSource('noise', 'Freq', 1472e6, 'Bandwidth', 20e6, 'INRdB', 0)];
blk = {{}, {'BlockSize', 1000003}, {'BlockSize', N4}};
y4 = cell(1, numel(blk));
for b = 1:numel(blk)
    rx(f4, 'NoiseStd', 1, 'Seed', 11, 'RFI', rfi4, blk{b}{:});
    y4{b} = readF32(fOut);
end
same = [isequal(y4{1}, y4{2}), isequal(y4{1}, y4{3})];
pass = all(same);
fprintf(['4. block size (noise + %d sources, 2.5 ms): BlockSize 1,000,003 = default %d, one block = ' ...
         'default %d: %s\n'], numel(rfi4), same, passStr(pass));
if ~pass, fails{end+1} = 'block size'; end
clear y4

% ---------------------------------------------------------------------------------
% 5. Prediction
% ---------------------------------------------------------------------------------
sg = 0.5e-3 / (2*sqrt(2*log(2)));                % 0.5 ms power FWHM (10 ms, 5 % duty)
g  = struct('A', A, 'sigma', sg, 'envelopeMode', 'power', 'pulseCenters', ((0:9) + 0.5) * 10e-3);
gA = g; gA.envelopeMode = 'amplitude';
fTiny = fullfile(tmp, 'tiny.dat');
writeF32(fTiny, zeros(1000, 1, 'single'));
cases5 = {'-54 dB', -54, g, sg; '0 dB (self-noise)', 0, g, sg; '-54 dB, amplitude mode', -54, gA, sg / sqrt(2)};
pass = true;
for c5 = 1:size(cases5, 1)
    [label, snr, gi, sp] = cases5{c5, :};
    i5 = rx(fTiny, 'SNRdB', snr, 'SignalInfo', gi);
    pr = i5.prediction;
    rh = 10^(snr/10);
    [sQ, tQ] = radiometer(rh, sp, B);
    d = [pr.snrPulse / sQ - 1, pr.toaErrPulse / tQ - 1];
    okP = all(abs(d) < 1e-6) && abs(pr.snrFolded / (pr.snrPulse * sqrt(10)) - 1) < 1e-14 && ...
          abs(pr.toaErrFolded * sqrt(10) / pr.toaErrPulse - 1) < 1e-14 && pr.nPulses == 10 && ...
          pr.bandwidth == B && pr.rho == rh;
    msg = '';
    if snr == -54
        w = [pr.snrPulse / (rh * sqrt(B * sp * sqrt(pi))) - 1, pr.toaErrPulse / (sqrt(2*sp / (B*sqrt(pi))) / rh) - 1];
        okP = okP && all(abs(w) < 1e-5);
        msg = sprintf(', weak-pulse forms %+.1e / %+.1e, SNR per pulse = %.2f rho', w, pr.snrPulse / rh);
    end
    pass = pass && okP;
    fprintf(['5. prediction, %s: SNR per pulse %.5g, TOA error per pulse %.5g us; against quadrature ' ...
             '%+.1e / %+.1e%s: %s\n'], label, pr.snrPulse, pr.toaErrPulse * 1e6, d, msg, passStr(okP));
end
if ~pass, fails{end+1} = 'prediction'; end

% ---------------------------------------------------------------------------------
% 6. Input checks
% ---------------------------------------------------------------------------------
cw = rfiSource('cw', 'Freq', 1.35e9);
bad = @(fin, fo, varargin) errId(@() addNoiseAndRFI(fin, fo, fs, 'Verbose', false, 'SaveInfo', false, varargin{:}));
ids = {bad(fZero, fZero, 'NoiseStd', 1), ...
       bad(fZero, fOut, 'NoiseStd', 1, 'Band', [1.6e9 1.2e9]), ...
       bad(fZero, fOut, 'NoiseStd', 1, 'Band', [1.2e9 2.1e9]), ...
       bad(fZero, fOut, 'SNRdB', -5), ...
       bad(fZero, fOut, 'RFI', cw), ...
       bad(fZero, fOut, 'NoiseStd', 1, 'RFI', rfiSource('cw', 'Freq', 2.5e9)), ...
       bad(fullfile(tmp, 'none.dat'), fOut, 'NoiseStd', 1)};
want = strcat('addNoiseAndRFI:', {'sameFile', 'band', 'band', 'amp', 'ref', 'freq', 'notFound'});
pass = isequal(ids, want);
fprintf('6. input checks: %d of %d bad inputs refused with the right id: %s\n', ...
    nnz(strcmp(ids, want)), numel(want), passStr(pass));
if ~pass, fails{end+1} = 'input checks'; end

% ---------------------------------------------------------------------------------
if isempty(fails)
    fprintf('testAddNoiseAndRFI: ALL PASSED\n');
else
    error('testAddNoiseAndRFI:failed', 'FAILED: %s', strjoin(fails, ', '));
end
end


% =================================================================================
function [snr, toa] = radiometer(rho, sp, B)
% The help text's flat-band radiometer integrals for a unit-peak Gaussian power profile of
% std sp, by adaptive quadrature (independent of the code's trapz grid):
% SNR^2 = B int (rho p)^2/(1 + rho p)^2 dt, 1/sigma_TOA^2 = B int (rho p')^2/(1 + rho p)^2 dt.
p   = @(t) exp(-0.5 * (t / sp).^2);
q   = @(f) integral(f, -20*sp, 20*sp, 'RelTol', 1e-10, 'AbsTol', 0);
snr = sqrt(B * q(@(t) (rho * p(t)).^2 ./ (1 + rho * p(t)).^2));
toa = 1 / sqrt(B * q(@(t) (rho * t / sp^2 .* p(t)).^2 ./ (1 + rho * p(t)).^2));
end

function writeF32(f, x)
fid = fopen(f, 'w', 'ieee-le'); fwrite(fid, x, 'single'); fclose(fid);
end

function y = readF32(f)
fid = fopen(f, 'r', 'ieee-le'); y = fread(fid, Inf, 'single=>single'); fclose(fid);
end

function id = errId(f)
% The identifier of the error f() raises ('' if none).
id = '';
try
    f();
catch e
    id = e.identifier;
end
end

function s = passStr(p)
if p, s = 'PASS'; else, s = 'FAIL'; end
end
