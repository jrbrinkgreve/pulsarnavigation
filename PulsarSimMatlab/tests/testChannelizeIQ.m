function testChannelizeIQ()
%TESTCHANNELIZEIQ  Unit tests for channelizeIQ.
%{
Run from the PulsarSimMatlab folder: run('tests/testChannelizeIQ.m'). ~10 s.
Synthetic tests write small files to tempdir; test 3 uses the receiver IQ
file (data/test/test_rx_IQ.dat, tests/makeTestData.m) if it exists.
Tests 4-7 (27 Oct 2026, review R-C1) cover grids and options tests 1-3 never ran.

  1. Tones: every channel sample must equal the exact prediction
     A*exp(i*phi) * Hc(df) * exp(i*2*pi*df*t_m), df = f - fc_j, t_m = m*D/fs.
     Checks frequency mapping, LO phase reference, timing (group delay
     removed) and gain at once. Also prints the passband/stopband response.
  2. White noise: channel variance = sum(h^2) (unit-variance input), all
     channels equal; channel PSD flat within +-ChanWidth/2.
  3. Real data: sum over channels of the power inside +-ChanWidth/2 equals
     the IQ power inside [fLow, fHigh] over the same time span (Parseval).
  4. LO position: channel centres sit at (k_j + beta)*fs/K from the LO; tests 1-3
     have beta = 1/2 (fLO 1.3 GHz). The tones of test 1 plus two at the band edges,
     with fLO 1.3015625 GHz (on a channel centre: beta = 0) and 1.301 GHz (off the
     grid: beta = 0.18): channels 1, 2, 6-8, 39-42, 99-101, 127, 128 = the formula of
     test 1 within 1e-4; the same channel centres, K, D and supported range; and the
     same channel files as with fLO 1.3 GHz to float32 rounding (relative rms 1e-5):
     a channel depends only on the RF signal and its centre, not on the LO.
  5. Oversampling 2 (D = 128, 6.25 MHz): prototype length from the Kaiser formula
     (1285 taps), flat to ChanWidth/2 and stopband beyond fs_out - ChanWidth/2 as in
     test 1; N and the supported range for D = 128; tones = the formula within 1e-4,
     one of them at +4 MHz from channel 50 (between fs_out/2 and the stopband: it
     aliases to -2.25 MHz, outside +-ChanWidth/2).
  6. BlockOut and T0 (complex noise, 2^16 samples, 342 outputs per channel): BlockOut
     1, 7, 100 and 341 (two blocks) give all 128 channel files byte for byte as the
     default (one block); T0 0.123 s changes info.t0 only.
  7. Input checks: fs/ChanWidth or K/Oversampling not an integer, fLow >= fHigh, a band
     that is not a whole number of channels, a band (with its transitions) outside
     fLO +- fs/2, a missing file and an empty file are refused, each with its own id.
Errors at the end if any check fails.
%}

addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'functions'));
fs = 800e6; fLO = 1.3e9; fLow = 1.2e9; fHigh = 1.6e9;
tmp = fullfile(tempdir, 'testChannelizeIQ');
if ~isfolder(tmp), mkdir(tmp); end
fails = {};

% ---------------------------------------------------------------------------------
% 1. Tones
% ---------------------------------------------------------------------------------
dF = 3.125e6;
fc = @(j) fLow + (j - 0.5) * dF;
tones = [ ... % RF frequency, amplitude, phase
    fc(40),            0.7,  0.3
    fc(41) + 1.2e6,    1.1, -1.0
    fc(100) - 1.5e6,   0.5,  2.0
    fc(7) + 0.37e6,    0.9,  1.234];
N = 2^20;
n = (0:N-1).';
x = zeros(N, 1);
for t = 1:size(tones, 1)
    x = x + tones(t,2) * exp(1i*(2*pi*(tones(t,1) - fLO)/fs * n + tones(t,3)));
end
fileIn = fullfile(tmp, 'tones.dat');
writeCF32(fileIn, x);
clear x n
info = channelizeIQ(fileIn, fullfile(tmp, 'tones_chan'), fs, fLO, fLow, fHigh, 'Verbose', false);
D = info.decimation; h = info.prototype; G = (numel(h) - 1)/2;
Hc = @(df) sum(h .* exp(-2i*pi*df*((0:numel(h)-1) - G)/fs));   % zero-phase response

ok = info.fullySupported(1):info.fullySupported(2);
tm = (ok - 1).' * D / fs;
worst = 0;
for j = [6 7 8 39 40 41 42 99 100 101 10]
    y = readCF32(info.chanFiles(j));
    y = y(ok);
    yExp = zeros(size(y));
    for t = 1:size(tones, 1)
        df = tones(t,1) - info.chanFreqs(j);
        yExp = yExp + tones(t,2)*exp(1i*tones(t,3)) * Hc(df) * exp(2i*pi*df*tm);
    end
    err = max(abs(y - yExp));
    worst = max(worst, err);
    fprintf('  tone test ch %3d: max |y - pred| = %.2e (max |pred| %.3g)\n', j, err, max(abs(yExp)));
end
pass = worst < 1e-4;
fprintf('1. tones: worst error %.2e (limit 1e-4): %s\n', worst, passStr(pass));
if ~pass, fails{end+1} = 'tones'; end

% Prototype response
dfP = [0 0.5 1.0 1.2 1.5 dF/2/1e6] * 1e6;
dfS = linspace(info.stopbandEdge, info.fs/2 + info.stopbandEdge, 2000);
gP  = arrayfun(Hc, dfP);
gS  = max(abs(arrayfun(Hc, dfS)));
fprintf('   passband |Hc| at %s MHz: %s\n', mat2str(dfP/1e6, 4), mat2str(abs(gP), 7));
fprintf('   stopband max |Hc| beyond %.4g MHz: %.2e (%.1f dB)\n', info.stopbandEdge/1e6, gS, 20*log10(gS));
pass = all(abs(abs(gP) - 1) < 2e-4) && gS < 2e-4;
fprintf('   prototype flat to dF/2 and >= ~74 dB stopband: %s\n', passStr(pass));
if ~pass, fails{end+1} = 'prototype'; end

% ---------------------------------------------------------------------------------
% 2. White noise
% ---------------------------------------------------------------------------------
N = 2^23;
rs = RandStream('mt19937ar', 'Seed', 7);
x = (randn(rs, N, 1) + 1i*randn(rs, N, 1)) / sqrt(2);       % E|x|^2 = 1
fileIn = fullfile(tmp, 'noise.dat');
writeCF32(fileIn, x);
clear x
info = channelizeIQ(fileIn, fullfile(tmp, 'noise_chan'), fs, fLO, fLow, fHigh, 'Verbose', false);
ok = info.fullySupported(1):info.fullySupported(2);
expVar = sum(info.prototype.^2);
pw = zeros(info.nChan, 1);
nSeg = 512; psd = zeros(nSeg, 1); nAvg = 0;
win = 0.5 - 0.5*cos(2*pi*(0:nSeg-1).'/nSeg);                 % Hann
for j = 1:info.nChan
    y = readCF32(info.chanFiles(j));
    y = y(ok);
    pw(j) = mean(abs(y).^2);
    nS = floor(numel(y)/nSeg);
    Y = fft(reshape(y(1:nS*nSeg), nSeg, nS) .* win);
    psd = psd + sum(abs(Y).^2, 2); nAvg = nAvg + nS;
end
% Periodogram normalized by sum(w^2): E = PSD [per Hz] * fs_out. The input
% (unit variance over fs) has PSD 1/fs; unit DC gain keeps that in the
% passband, so the expected level there is fs_out/fs = 1/D.
psd = fftshift(psd / nAvg / sum(win.^2));
fAx = ((0:nSeg-1).' - nSeg/2) * info.fs / nSeg;
inP = abs(fAx) <= dF/2 - info.fs/nSeg;                       % stay off the edge bin
lvl  = mean(psd(inP)) * info.decimation;                     % should be 1
flat = psd(inP) / mean(psd(inP));
rel  = pw / expVar;
nEff = numel(ok) * dF / info.fs;                              % ~independent samples per channel
fprintf(['2. noise: channel power / sum(h^2): mean %.5f, min %.4f, max %.4f ' ...
         '(1-sigma per channel ~%.4f)\n'], mean(rel), min(rel), max(rel), 1/sqrt(nEff));
fprintf('   passband PSD level x D = %.4f (expect 1); within +-dF/2 min/max vs mean %.4f / %.4f\n', ...
    lvl, min(flat), max(flat));
pass = abs(mean(rel) - 1) < 0.01 && all(abs(rel - 1) < 6/sqrt(nEff)) && ...
       abs(lvl - 1) < 0.01 && max(abs(flat - 1)) < 0.03;
fprintf('   noise checks: %s\n', passStr(pass));
if ~pass, fails{end+1} = 'noise'; end

% ---------------------------------------------------------------------------------
% 3. Real data (the frozen receiver IQ, data/test/), Parseval over [fLow, fHigh]
% ---------------------------------------------------------------------------------
dataDir = fullfile(fileparts(mfilename('fullpath')), '..', 'data', 'test');   % tests/makeTestData.m
fileIQ  = fullfile(dataDir, 'test_rx_IQ.dat');
if isfile(fileIQ)
    infoIQ = loadInfo(fileIQ);
    info = channelizeIQ(fileIQ, fullfile(dataDir, 'chan', 'test_rx_IQ_chan'), ...
        infoIQ.actualFsOut, infoIQ.fLO, fLow, fHigh, 'T0', infoIQ.t0);
    D = info.decimation;
    % Common span: channel samples m = mA..mB (fully supported), IQ samples m*D ...
    nCh = 2^16;                                               % channel samples used
    mA  = info.fullySupported(1) - 1;                         % 0-based
    nIQ = nCh * D;
    xIQ = readCF32(fileIQ, mA*D, nIQ);
    X   = fft(double(xIQ)); clear xIQ
    fIQ = infoIQ.fLO + [0:nIQ/2-1, -nIQ/2:-1].' * infoIQ.actualFsOut / nIQ;
    pIQ = sum(abs(X(fIQ >= fLow & fIQ < fHigh)).^2) / nIQ^2;  % band power (mean |x|^2 units)
    clear X fIQ
    % Band part of mean|x|^2 = sum over in-band bins |X|^2/n^2 (Parseval). A unit-gain
    % channel keeps the PSD, so its part inside +-dF/2 is that slice of mean|x|^2:
    % the channel slices must add up to the IQ band power (no factor D).
    pj  = zeros(info.nChan, 1);
    fCh = [0:nCh/2-1, -nCh/2:-1].' * info.fs / nCh;
    inB = fCh >= -dF/2 & fCh < dF/2;
    for j = 1:info.nChan
        y = readCF32(info.chanFiles(j), mA, nCh);
        Y = fft(double(y));
        pj(j) = sum(abs(Y(inB)).^2) / nCh^2;
    end
    pCh = sum(pj);
    ratio = pCh / pIQ;
    fprintf(['3. real data: sum of channel band powers / IQ band power = %.5f; ' ...
             'channel powers min/mean/max relative %.3f / 1 / %.3f\n'], ratio, ...
        min(pj)/mean(pj), max(pj)/mean(pj));
    pass = abs(ratio - 1) < 2e-3;
    fprintf('   Parseval check (limit 2e-3): %s\n', passStr(pass));
    if ~pass, fails{end+1} = 'real-data Parseval'; end
else
    fprintf('3. real data: %s not found, skipped\n', fileIQ);
end

% ---------------------------------------------------------------------------------
% 4. LO position: the grid offset beta (review R-C1)
% ---------------------------------------------------------------------------------
HcOf = @(h, df) sum(h .* exp(-2i*pi*df*((0:numel(h)-1) - (numel(h)-1)/2)/fs));   % as Hc of test 1
tones4 = [tones; fc(1) + 0.4e6, 0.8, -0.6; fc(128) - 0.9e6, 0.6, 0.9];
chk4 = [1 2 6 7 8 39 40 41 42 99 100 101 127 128];
fLOs = [1.3e9, 1.3015625e9, 1.301e9];
N = 2^20;
fileIn = fullfile(tmp, 'lo.dat');
Y4 = cell(1, numel(fLOs)); I4 = cell(1, numel(fLOs));
pass = true;
for k = 1:numel(fLOs)
    writeCF32(fileIn, toneInput(tones4, fLOs(k), fs, N));
    I4{k} = channelizeIQ(fileIn, fullfile(tmp, sprintf('lo%d_chan', k)), fs, fLOs(k), fLow, fHigh, ...
        'Verbose', false, 'SaveInfo', false);
    [err4, Y4{k}] = toneError(I4{k}, tones4, chk4, HcOf, fs);
    a = (fLow - fLOs(k)) / dF + 0.5;               % first channel centre from the LO, in channels
    beta = a - floor(a + 1e-9);
    okG = I4{k}.fLOIn == fLOs(k) && isequal(I4{k}.chanFreqs, I4{1}.chanFreqs) && I4{k}.K == I4{1}.K && ...
          I4{k}.decimation == I4{1}.decimation && I4{k}.N == I4{1}.N && ...
          isequal(I4{k}.fullySupported, I4{1}.fullySupported);
    okK = err4 < 1e-4 && okG;
    dLO = 0;
    if k > 1
        dLO = sqrt(sum(abs(double(Y4{k}(:)) - double(Y4{1}(:))).^2) / sum(abs(double(Y4{1}(:))).^2));
        okK = okK && dLO < 1e-5;
    end
    pass = pass && okK;
    fprintf(['4. LO %.7f GHz (beta %.2f): worst |y - pred| over %d channels %.2e (< 1e-4); same grid %d; ' ...
             'relative rms from the fLO 1.3 GHz files %.1e (< 1e-5): %s\n'], fLOs(k) / 1e9, beta, numel(chk4), ...
        err4, okG, dLO, passStr(okK));
end
if ~pass, fails{end+1} = 'LO position'; end
clear Y4

% ---------------------------------------------------------------------------------
% 5. Oversampling 2
% ---------------------------------------------------------------------------------
tones5 = [tones; fc(50) + 4.0e6, 0.9, 0.4];      % channel 50: between fs_out/2 and its stopband
writeCF32(fileIn, toneInput(tones5, fLO, fs, N));
info = channelizeIQ(fileIn, fullfile(tmp, 'ov2_chan'), fs, fLO, fLow, fHigh, 'Oversampling', 2, ...
    'Verbose', false, 'SaveInfo', false);
K5 = fs / dF; D5 = K5 / 2; fsO = fs / D5;
M5 = ceil((80 - 8) / (2.285 * 2*pi * (fsO - dF) / fs)) + 1;   % the help text's Kaiser length, odd
M5 = M5 + (mod(M5, 2) == 0);
G5 = (M5 - 1) / 2;
okI = info.K == K5 && info.decimation == D5 && info.fs == fsO && info.oversampling == 2 && ...
      info.filterLen == M5 && info.stopbandEdge == fsO - dF/2 && info.N == ceil(N / D5) && ...
      isequal(info.fullySupported, [ceil(G5/D5) + 1, floor((N - 1 - G5)/D5) + 1]);
gP = arrayfun(@(d) HcOf(info.prototype, d), [0 0.5 1.0 1.2 1.5 dF/2/1e6] * 1e6);
gS = max(abs(arrayfun(@(d) HcOf(info.prototype, d), linspace(info.stopbandEdge, info.fs/2 + info.stopbandEdge, 2000))));
okP = all(abs(abs(gP) - 1) < 2e-4) && gS < 2e-4;
err5 = toneError(info, tones5, [1 7 40 41 49 50 51 52 100 128], HcOf, fs);
pass = okI && okP && err5 < 1e-4;
fprintf(['5. oversampling 2: K %d, D %d, %.4g MHz, %d taps (formula %d), flat to dF/2 within %.1e, stopband ' ...
         'max %.1f dB beyond %.4g MHz; tones worst |y - pred| %.2e (< 1e-4); the +4 MHz tone lands in ' ...
         'channel 50 at %.2f MHz with |Hc| %.3f: %s\n'], info.K, info.decimation, info.fs / 1e6, info.filterLen, ...
    M5, max(abs(abs(gP) - 1)), 20*log10(gS), info.stopbandEdge / 1e6, err5, (4.0e6 - fsO) / 1e6, ...
    abs(HcOf(info.prototype, 4.0e6)), passStr(pass));
if ~pass, fails{end+1} = 'oversampling 2'; end

% ---------------------------------------------------------------------------------
% 6. BlockOut and T0
% ---------------------------------------------------------------------------------
rs6 = RandStream('mt19937ar', 'Seed', 8);
fileN6 = fullfile(tmp, 'noise6.dat');
writeCF32(fileN6, (randn(rs6, 2^16, 1) + 1i*randn(rs6, 2^16, 1)) / sqrt(2));
ch = @(name, varargin) channelizeIQ(fileN6, fullfile(tmp, ['bo_' name]), fs, fLO, fLow, fHigh, ...
    'Verbose', false, 'SaveInfo', false, varargin{:});
iRef = ch('ref');
ref  = fileBytes(iRef.chanFiles);
bos  = [1 7 100 iRef.N - 1];
same = false(size(bos)); nb = zeros(size(bos));
for k = 1:numel(bos)
    ib = ch(sprintf('%d', bos(k)), 'BlockOut', bos(k));
    same(k) = isequal(fileBytes(ib.chanFiles), ref);
    nb(k) = ib.nBlocks;
end
iT = ch('t0', 'T0', 0.123);
okT = isequal(fileBytes(iT.chanFiles), ref) && iT.t0 == 0.123;
pass = all(same) && okT && iRef.nBlocks == 1;
fprintf(['6. BlockOut %s (%s blocks) vs the default (%d block, %d samples per channel): all %d channels byte ' ...
         'for byte %s; T0 0.123 s: files unchanged and info.t0 set %d: %s\n'], mat2str(bos), mat2str(nb), ...
    iRef.nBlocks, iRef.N, iRef.nChan, mat2str(same), okT, passStr(pass));
if ~pass, fails{end+1} = 'BlockOut / T0'; end

% ---------------------------------------------------------------------------------
% 7. Input checks
% ---------------------------------------------------------------------------------
fEmpty = fullfile(tmp, 'empty.dat');
fclose(fopen(fEmpty, 'w'));
cz = @(fi, varargin) errId(@() channelizeIQ(fi, fullfile(tmp, 'bad_chan'), fs, varargin{:}, ...
    'Verbose', false, 'SaveInfo', false));
ids = {cz(fileN6, fLO, fLow, fHigh, 'ChanWidth', 3e6), cz(fileN6, fLO, fLow, fHigh, 'Oversampling', 1.5), ...
       cz(fileN6, fLO, fHigh, fLow), cz(fileN6, fLO, fLow, 1.6005e9), cz(fileN6, 1.0e9, fLow, fHigh), ...
       cz(fullfile(tmp, 'none.dat'), fLO, fLow, fHigh), cz(fEmpty, fLO, fLow, fHigh)};
want = strcat('channelizeIQ:', {'grid', 'grid', 'band', 'band', 'band', 'notFound', 'empty'});
pass = isequal(ids, want);
fprintf('7. input checks: %d of %d bad inputs refused with the right id: %s\n', ...
    nnz(strcmp(ids, want)), numel(want), passStr(pass));
if ~pass, fails{end+1} = 'input checks'; end

% ---------------------------------------------------------------------------------
if isempty(fails)
    fprintf('testChannelizeIQ: ALL PASSED\n');
else
    error('testChannelizeIQ:failed', 'FAILED: %s', strjoin(fails, ', '));
end
end


% =================================================================================
function s = passStr(p)
if p, s = 'PASS'; else, s = 'FAIL'; end
end

function writeCF32(file, x)
fid = fopen(file, 'w', 'ieee-le');
fwrite(fid, [real(x(:)).'; imag(x(:)).'], 'single');
fclose(fid);
end

function x = toneInput(tones, fLO, fs, N)
%TONEINPUT  Complex baseband (LO fLO) of the RF tones [f, A, phi] (rows): sum A e^(i(2 pi (f - fLO) n/fs + phi)).
n = (0:N-1).';
x = zeros(N, 1);
for t = 1:size(tones, 1)
    x = x + tones(t,2) * exp(1i*(2*pi*(tones(t,1) - fLO)/fs * n + tones(t,3)));
end
end

function [worst, Y] = toneError(info, tones, chk, HcOf, fs)
%TONEERROR  All channels of a channelizeIQ run (one column each) and the largest error, over the
% channels chk and their fully supported samples, against the formula of test 1:
% A e^(i phi) Hc(df) e^(i 2 pi df m D/fs), df = f - fc_j, Hc from the reported prototype.
Y = complex(zeros(info.N, info.nChan, 'single'));
for j = 1:info.nChan
    Y(:, j) = readCF32(info.chanFiles(j));
end
ok = info.fullySupported(1):info.fullySupported(2);
tm = (ok - 1).' * info.decimation / fs;
worst = 0;
for j = chk
    yExp = zeros(numel(ok), 1);
    for t = 1:size(tones, 1)
        df = tones(t,1) - info.chanFreqs(j);
        yExp = yExp + tones(t,2) * exp(1i*tones(t,3)) * HcOf(info.prototype, df) * exp(2i*pi*df*tm);
    end
    worst = max(worst, max(abs(double(Y(ok, j)) - yExp)));
end
end

function b = fileBytes(files)
%FILEBYTES  The raw bytes of each file (cell, one per file), for byte-for-byte comparisons.
b = cell(numel(files), 1);
for k = 1:numel(files)
    fid = fopen(files(k), 'r'); b{k} = fread(fid, Inf, 'uint8=>uint8'); fclose(fid);
end
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

function y = readCF32(file, first, count)
%READCF32  Complex samples (column); optional 0-based first sample and count.
fid = fopen(file, 'r', 'ieee-le');
c = onCleanup(@() fclose(fid));
if nargin > 1
    fseek(fid, first * 8, 'bof');
    raw = fread(fid, [2 count], 'single=>single');
else
    raw = fread(fid, [2 Inf], 'single=>single');
end
y = complex(raw(1, :), raw(2, :)).';
end
