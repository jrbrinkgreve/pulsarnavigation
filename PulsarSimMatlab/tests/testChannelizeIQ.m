function testChannelizeIQ()
%TESTCHANNELIZEIQ  Unit tests for channelizeIQ.
%{
Run from the PulsarSimMatlab folder: run('tests/testChannelizeIQ.m'). ~10 s.
Synthetic tests write small files to tempdir; test 3 uses the receiver IQ
file of main.m (data/test_rx_IQ.dat) if it exists.

  1. Tones: every channel sample must equal the exact prediction
     A*exp(i*phi) * Hc(df) * exp(i*2*pi*df*t_m), df = f - fc_j, t_m = m*D/fs.
     Checks frequency mapping, LO phase reference, timing (group delay
     removed) and gain at once. Also prints the passband/stopband response.
  2. White noise: channel variance = sum(h^2) (unit-variance input), all
     channels equal; channel PSD flat within +-ChanWidth/2.
  3. Real data: sum over channels of the power inside +-ChanWidth/2 equals
     the IQ power inside [fLow, fHigh] over the same time span (Parseval).
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
% 3. Real data (receiver IQ of main.m), Parseval over [fLow, fHigh]
% ---------------------------------------------------------------------------------
dataDir = fullfile(fileparts(mfilename('fullpath')), '..', 'data');
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
