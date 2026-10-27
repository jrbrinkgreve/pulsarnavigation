function testApplyDispersionStream()
%TESTAPPLYDISPERSIONSTREAM  Unit tests for applyDispersionStream (stage 2: interstellar dispersion).
%{
Run from the PulsarSimMatlab folder: run('tests/testApplyDispersionStream.m'). ~3 s.
Writes its own files to tempdir (~0.2 GB, deleted at the end); needs no data files.
Ground truth is the transfer function of the help text, computed here in double precision,
independent of the kernel code: from the input file to the output file (the stream removes
the kernel's s0 samples again)
  H_eff(f) = W(f) * exp(i*phi(f)),   phi(f) = 2*pi*K*DM*(f - fRef)^2 / (f*fRef^2),
whose group delay is tau(f) - tau(fRef), tau(f) = K*DM/f^2: low frequencies arrive later.
W is 1 in the band, with sin^2 edges 8 MHz wide inside fLow and fHigh, and 0 outside.
A forward-plus-inverse round trip cannot catch a wrong sign (the first kernel had one, and
its round trip still worked), so every check compares with this formula.
The receiver band 1.2-1.6 GHz at 4 GHz. DM 0.05 (63 us smear, kernel ~272k samples) where
the size does not matter; the reference DM 5 (6.3 ms, 25.2M samples) for the delays.

  1. Format and info (one impulse, DM 0.05): N samples, 4*N bytes, format fields, refFreq
     default fHigh, dispersionConst, smearTime = tau(fLow) - tau(fHigh), bulkDelayRef =
     tau(fRef), edgeWidth, a block size that fits (Nfft >= blockLen + M - 1, nBlocks covers
     N), the saved _info.mat = the returned info.
  2. Transfer function (one impulse, DM 0.05): the FFT of the output, re-referenced to the
     impulse's sample, against H_eff on the whole 0..fs/2 grid: relative error energy
     sum|dH|^2 / sum|H|^2 below 1e-9 (kernel truncation ~1e-12 plus float32 rounding; one
     sample off in s0, the wrong sign or a wrong constant phase each give an error of
     order 1). RefFreq = fHigh (default), 1.4 GHz (part of the band arrives before the
     impulse), fLow; DM = 0 (H = W: band-pass only); EdgeFrac 0.01 and 0.05 (the edge width,
     and the guard 20/edgeW that follows it). Printed: the largest |H| outside the band, and
     the measured group delay at 1.25 ... 1.55 GHz against tau(f) - tau(fRef).
  3. Delays at DM 5: five Gaussian tone bursts (sigma 0.5 us) at 1.25, 1.3, 1.4, 1.5 and
     1.55 GHz, all at one instant, leave 5.17 ... 0.53 ms later, each in its own window.
     Energy-centroid shift = the power-weighted mean group delay over the burst's own
     spectrum (exact where |H| = 1) within 0.01 sample (2.5 ps); burst energy unchanged
     (|H| = 1) within 1e-5. Printed: the curvature term, shift - (tau(f0) - tau(fRef)) (a
     few samples), which is why the plain tau(f0) is not the prediction.
  4. Overlap-add = linear convolution (white noise, 4,000,001 samples, DM 0.05). Reference:
     one FFT in double, zero-padded to >= N + M, with H_eff (zeros before sample 1, anything
     delayed past the end dropped). Automatic block size, BlockLen = M - 1 (the minimum, 15
     blocks), 300,007 (odd), 400,000 (a last block of 1 sample) and N (one block): each
     within relative rms 1e-5 and max 1e-4 of the reference, the first and last M samples
     (carry start, flush) within rms 1e-5 too. Not bit for bit between block lengths: each
     uses another FFT size, so the rounding differs. And a file of 5001 samples, shorter
     than the kernel's zero lag s0 = 10,000 (all its output comes from the flush): within
     rms 1e-5 of its own reference.
  5. Power: info.kernelEnergy = 2*(B - 1.25*edgeW)/fs (each sin^4 edge averages 3/8) within
     1e-6; the output power of check 4's noise away from the file ends = kernelEnergy
     within 4 sigma (sigma from the exact sums of |H|^2 and |H|^4: ~1.7e-3).
  6. Leakage: every run above stays below 1e-10 (MaxLeakage warns at 1e-8); GuardTime
     2/edgeW (a tenth of the default) warns (leakage ~4e-8), and check 2's error energy then
     equals info.kernelLeakage within 1 % (the cut tails are the error; the error floor of
     the default runs, ~1e-12, is 3e-5 of it).
  7. Input checks: fLow >= fHigh, fHigh >= fs/2, RefFreq outside the band (both sides),
     EdgeFrac > 0.5, inFile = outFile, a missing file and BlockLen < M - 1 are refused,
     each with its own id.
Errors at the end if any check fails.
%}

root = fullfile(fileparts(mfilename('fullpath')), '..');
addpath(fullfile(root, 'functions'));
oldDir = cd(root); restoreDir = onCleanup(@() cd(oldDir));
tmp = tempname(tempdir);                         % own folder per run
mkdir(tmp);
cleanTmp = onCleanup(@() rmdir(tmp, 's'));
fails = {};

fs = 4e9; fLow = 1.2e9; fHigh = 1.6e9;           % the receiver band at the RF rate (testParams)
K = 4.148808e3;                                  % s MHz^2 pc^-1 cm^3 (help text)
edgeW = 0.02 * (fHigh - fLow);                   % default EdgeFrac 0.02: 8 MHz
tau = @(f, DM) K * 1e12 * DM ./ f.^2;            % s, f in Hz
fin = fullfile(tmp, 'in.dat'); fout = fullfile(tmp, 'out.dat');
ads = @(DM, varargin) applyDispersionStream(fin, fout, DM, fs, fLow, fHigh, ...
    'Verbose', false, 'SaveInfo', false, varargin{:});
leaks = [];                                      % kernelLeakage of every default-guard run (check 6)
warned = false;                                  % any leakage warning in those runs

% ---------------------------------------------------------------------------------
% 1. Format and info
% ---------------------------------------------------------------------------------
DM2 = 0.05; N2 = 2^20; k0 = 2^19;
x2 = zeros(N2, 1, 'single'); x2(k0 + 1) = 1;     % one impulse at 0-based sample k0
writeF32(fin, x2);
i1 = applyDispersionStream(fin, fout, DM2, fs, fLow, fHigh, 'Verbose', false);
f = dir(fout);
S = load(fullfile(tmp, 'out_info.mat'));
M1 = i1.kernelLen;
leaks(end+1) = i1.kernelLeakage;
ok = [i1.N == N2, f.bytes == 4*N2, i1.samplesWritten == N2, strcmp(i1.precision, 'single'), ...
      strcmp(i1.byteOrder, 'ieee-le'), ~i1.isComplex, i1.fs == fs, i1.actualFsOut == fs, ...
      i1.DM == DM2, i1.fLow == fLow, i1.fHigh == fHigh, i1.refFreq == fHigh, i1.dispersionConst == K, ...
      abs(i1.smearTime / (tau(fLow, DM2) - tau(fHigh, DM2)) - 1) < 1e-12, ...
      abs(i1.bulkDelayRef / tau(fHigh, DM2) - 1) < 1e-12, i1.edgeWidth == edgeW, i1.autoBlockSize, ...
      i1.Nfft >= i1.blockLen + M1 - 1, (i1.nBlocks - 1) * i1.blockLen < N2, N2 <= i1.nBlocks * i1.blockLen, ...
      isequal(rmfield(i1, 'infoFile'), S.info)];
pass = all(ok);
fprintf(['1. format and info: N %d, %d bytes, refFreq %.4g GHz (default fHigh), smear %.6f us, ' ...
         'bulk delay %.6f us, kernel %d samples, Nfft 2^%d, %d block(s), info file: %d of %d checks: %s\n'], ...
    i1.N, f.bytes, i1.refFreq / 1e9, i1.smearTime * 1e6, i1.bulkDelayRef * 1e6, M1, log2(i1.Nfft), ...
    i1.nBlocks, nnz(ok), numel(ok), passStr(pass));
if ~pass, fails{end+1} = 'format / info'; end

% ---------------------------------------------------------------------------------
% 2. Transfer function from one impulse
% ---------------------------------------------------------------------------------
cases2 = {'RefFreq = fHigh (default)', DM2, fHigh, {},                  0.02
          'RefFreq = 1.4 GHz',          DM2, 1.4e9, {'RefFreq', 1.4e9},   0.02
          'RefFreq = fLow',             DM2, fLow,  {'RefFreq', fLow},    0.02
          'DM = 0 (band-pass only)',    0,   fHigh, {},                  0.02
          'EdgeFrac 0.01',              DM2, fHigh, {'EdgeFrac', 0.01},  0.01
          'EdgeFrac 0.05',              DM2, fHigh, {'EdgeFrac', 0.05},  0.05};
f0 = [1.25 1.3 1.4 1.5 1.55] * 1e9;              % where the group delay is printed
df = fs / N2;
kb = round(f0 / df);                             % 0-based bins
fk = (0:N2/2).' * df;
for c = 1:size(cases2, 1)
    [label, dm, fRef, opts2, eF] = cases2{c, :};
    lastwarn('');
    info = ads(dm, opts2{:});
    [~, wId] = lastwarn;
    warned = warned || strcmp(wId, 'applyDispersionStream:leakage');
    leaks(end+1) = info.kernelLeakage; %#ok<AGROW>
    [errE, Hm] = transferError(readF32(fout), k0, fs, dm, fRef, fLow, fHigh, eF * (fHigh - fLow));
    pos = Hm(1:N2/2 + 1);
    outB = max(abs(pos(fk < fLow | fk > fHigh)));
    tauM = -angle(Hm(kb + 2) .* conj(Hm(kb + 1))) / (2*pi*df);    % phase step across one bin
    tauP = tau((kb + 0.5) * df, dm) - tau(fRef, dm);              % at the bin's midpoint
    dTau = max(abs(tauM(:) - tauP(:)));
    pass = errE < 1e-9;
    fprintf(['2. transfer function, %s: error energy %.2e (< 1e-9), largest |H| outside the band ' ...
             '%.1e, group delay at 1.25..1.55 GHz %s us, largest difference from tau(f) - tau(fRef) ' ...
             '%.2g ps: %s\n'], label, errE, outB, sprintf('%.4f ', tauM * 1e6), dTau * 1e12, passStr(pass));
    if ~pass, fails{end+1} = ['transfer function, ' label]; end %#ok<AGROW>
end

% ---------------------------------------------------------------------------------
% 3. Delays at DM 5: tone bursts
% ---------------------------------------------------------------------------------
DM3 = 5; N3 = 23e6; c3 = 1e6; sigN = 2000;       % burst centre (0-based sample), sigma 0.5 us
f3 = [1.25 1.3 1.4 1.5 1.55] * 1e9;
m = (-10*sigN:10*sigN).';
xb = exp(-0.5 * (m / sigN).^2) .* cos(2*pi * (f3 / fs) .* abs(m));   % one burst per column, symmetric
x3 = zeros(N3, 1, 'single');
x3(c3 + 1 + m) = single(sum(xb, 2));
writeF32(fin, x3);
clear x3
lastwarn('');
i3 = ads(DM3);
[~, wId] = lastwarn;
warned = warned || strcmp(wId, 'applyDispersionStream:leakage');
leaks(end+1) = i3.kernelLeakage;
y3 = readF32(fout);
Np = 2^20;                                       % spectrum grid of one burst (>= 2 x its length)
fp = (0:Np/2).' * (fs / Np);
band = fp >= fLow & fp <= fHigh;                 % the bursts have no power outside (Gaussian tails)
[shiftMs, curv, err3, eRatio] = deal(zeros(1, numel(f3)));
for j = 1:numel(f3)
    X = fft(xb(:, j), Np);
    P = abs(X(1:Np/2 + 1)).^2;
    shift = sum(P(band) .* (tau(fp(band), DM3) - tau(fHigh, DM3))) / sum(P(band));   % s
    cP = c3 + shift * fs;                        % predicted output centroid (0-based sample)
    n = round(cP) + (-250000:250000).';          % +-62.5 us (the dispersed burst: sigma <= 5 us)
    e = double(y3(n + 1)).^2;
    err3(j) = sum(n .* e) / sum(e) - cP;
    eRatio(j) = sum(e) / sum(xb(:, j).^2);
    shiftMs(j) = shift * 1e3;
    curv(j) = (shift - (tau(f3(j), DM3) - tau(fHigh, DM3))) * fs;
end
clear y3
pass = max(abs(err3)) < 0.01 && max(abs(eRatio - 1)) < 1e-5;
fprintf(['3. delays at DM 5 (kernel %d samples): bursts at 1.25..1.55 GHz shifted %s ms; ' ...
         'curvature term %s samples; centroid minus prediction %s samples (< 0.01); ' ...
         'energy ratio - 1 at most %.1e (< 1e-5): %s\n'], i3.kernelLen, sprintf('%.4f ', shiftMs), ...
    sprintf('%.2f ', curv), sprintf('%+.1e ', err3), max(abs(eRatio - 1)), passStr(pass));
if ~pass, fails{end+1} = 'delays at DM 5'; end

% ---------------------------------------------------------------------------------
% 4. Overlap-add = linear convolution
% ---------------------------------------------------------------------------------
N4 = 4000001;
rs = RandStream('mt19937ar', 'Seed', 4);
x4 = single(randn(rs, N4, 1));
writeF32(fin, x4);
lastwarn('');
i4 = ads(DM2);
[~, wId] = lastwarn;
warned = warned || strcmp(wId, 'applyDispersionStream:leakage');
leaks(end+1) = i4.kernelLeakage;
yAuto = readF32(fout);
M4 = i4.kernelLen;
Np4 = 2^nextpow2(N4 + M4);
H4 = hEffGrid(Np4, fs, DM2, fHigh, fLow, fHigh, edgeW);
yRef = ifft(fft(double(x4), Np4) .* H4, 'symmetric');
yRef = yRef(1:N4);
rmsRef = sqrt(mean(yRef.^2));
bls = [NaN, M4 - 1, 300007, 400000, N4];         % NaN: the automatic block size
pass = true;
for b = 1:numel(bls)
    if isnan(bls(b))
        y = yAuto; info = i4; label = 'automatic';
    else
        info = ads(DM2, 'BlockLen', bls(b)); y = readF32(fout); label = sprintf('%d', bls(b));
    end
    d = double(y) - yRef;
    r = [sqrt(mean(d.^2)), max(abs(d)), sqrt(mean(d(1:M4).^2)), sqrt(mean(d(end-M4+1:end).^2))] / rmsRef;
    okB = r(1) < 1e-5 && r(2) < 1e-4 && r(3) < 1e-5 && r(4) < 1e-5;
    pass = pass && okB;
    fprintf(['4. overlap-add, BlockLen %s: %d blocks (last %d samples), relative rms %.2e, max %.2e, ' ...
             'first M %.2e, last M %.2e: %s\n'], label, info.nBlocks, ...
        N4 - (info.nBlocks - 1) * info.blockLen, r, passStr(okB));
end
N5 = 5001;                                       % shorter than s0: all output from the flush
x5 = single(randn(rs, N5, 1));
writeF32(fin, x5);
i5 = ads(DM2);
y5 = readF32(fout);
Np5 = 2^nextpow2(N5 + M4);
y5Ref = ifft(fft(double(x5), Np5) .* hEffGrid(Np5, fs, DM2, fHigh, fLow, fHigh, edgeW), 'symmetric');
y5Ref = y5Ref(1:N5);
r5 = sqrt(mean((double(y5) - y5Ref).^2)) / sqrt(mean(y5Ref.^2));
okS = numel(y5) == N5 && i5.kernelZeroLag > N5 && r5 < 1e-5;
pass = pass && okS;
leaks(end+1) = i5.kernelLeakage;
fprintf('4. overlap-add, a file of %d samples (s0 = %d, all output from the flush): relative rms %.2e: %s\n', ...
    N5, i5.kernelZeroLag, r5, passStr(okS));
if ~pass, fails{end+1} = 'overlap-add'; end

% ---------------------------------------------------------------------------------
% 5. Power
% ---------------------------------------------------------------------------------
E0 = 2 * (fHigh - fLow - 1.25*edgeW) / fs;       % sum of h^2 = output/input power, white input
C0 = sum(abs(H4).^2) / Np4;                      % the same from the formula's grid
S4 = sum(abs(H4).^4) / Np4;                      % sum over lags of the autocovariance squared
inner = (M4 + 1 : N4 - M4).';                    % full kernel support, away from the file ends
pw = mean(double(yAuto(inner)).^2);
sRel = sqrt(2 * S4 / (numel(inner) * C0^2));
z5 = (pw / i4.kernelEnergy - 1) / sRel;
pass = abs(i4.kernelEnergy / E0 - 1) < 1e-6 && abs(z5) < 4;
fprintf(['5. power: kernelEnergy %.8f, 2(B - 1.25 edgeW)/fs %.8f (formula grid %.8f); output power ' ...
         '%.6f = kernelEnergy x (1 %+.2e), %.2f sigma (sigma %.2e): %s\n'], i4.kernelEnergy, E0, C0, ...
    pw, pw / i4.kernelEnergy - 1, z5, sRel, passStr(pass));
if ~pass, fails{end+1} = 'power'; end
clear x4 yAuto yRef y d

% ---------------------------------------------------------------------------------
% 6. Leakage
% ---------------------------------------------------------------------------------
writeF32(fin, x2);                               % the impulse of checks 1-2
lastwarn('');
i6 = ads(DM2, 'GuardTime', 2 / edgeW);
[~, wId6] = lastwarn;
err6 = transferError(readF32(fout), k0, fs, DM2, fHigh, fLow, fHigh, edgeW);
ratio = err6 / i6.kernelLeakage;
pass = all(leaks < 1e-10) && ~warned && strcmp(wId6, 'applyDispersionStream:leakage') && ...
       abs(ratio - 1) < 0.01;
fprintf(['6. leakage: default guard at most %.1e over %d runs (< 1e-10), warnings %d; GuardTime 2/edgeW: ' ...
         'leakage %.2e, warns %d (the warning above), error energy %.2e = %.5f x leakage: %s\n'], ...
    max(leaks), numel(leaks), ...
    warned, i6.kernelLeakage, strcmp(wId6, 'applyDispersionStream:leakage'), err6, ratio, passStr(pass));
if ~pass, fails{end+1} = 'leakage'; end

% ---------------------------------------------------------------------------------
% 7. Input checks
% ---------------------------------------------------------------------------------
bad = @(varargin) errId(@() applyDispersionStream(varargin{:}, 'Verbose', false, 'SaveInfo', false));
ids = {bad(fin, fout, DM2, fs, fHigh, fLow), ...
       bad(fin, fout, DM2, 3e9, fLow, fHigh), ...
       bad(fin, fout, DM2, fs, fLow, fHigh, 'RefFreq', 1.7e9), ...
       bad(fin, fout, DM2, fs, fLow, fHigh, 'RefFreq', 1.1e9), ...
       bad(fin, fout, DM2, fs, fLow, fHigh, 'EdgeFrac', 0.6), ...
       bad(fin, fin, DM2, fs, fLow, fHigh), ...
       bad(fullfile(tmp, 'none.dat'), fout, DM2, fs, fLow, fHigh), ...
       bad(fin, fout, DM2, fs, fLow, fHigh, 'BlockLen', M1 - 2)};
want = strcat('applyDispersionStream:', {'band', 'nyquist', 'ref', 'ref', 'edge', 'sameFile', ...
       'notFound', 'blockLen'});
pass = isequal(ids, want);
fprintf('7. input checks: %d of %d bad inputs refused with the right id: %s\n', ...
    nnz(strcmp(ids, want)), numel(want), passStr(pass));
if ~pass, fails{end+1} = 'input checks'; end

% ---------------------------------------------------------------------------------
if isempty(fails)
    fprintf('testApplyDispersionStream: ALL PASSED\n');
else
    error('testApplyDispersionStream:failed', 'FAILED: %s', strjoin(fails, ', '));
end
end


% =================================================================================
function H = hEff(f, DM, fRef, fLow, fHigh, edgeW)
% The help text's transfer function W(f)*exp(i*phi(f)) at frequencies f >= 0 (Hz), from the
% formula alone (not the kernel code): 0 outside (fLow, fHigh), sin^2 edges inside.
KDM = 4.148808e3 * 1e12 * DM;
in = f > fLow & f < fHigh;
W = double(in);
lo = in & f < fLow + edgeW;  W(lo) = sin(pi/2 * (f(lo) - fLow) / edgeW).^2;
hi = in & f > fHigh - edgeW; W(hi) = sin(pi/2 * (fHigh - f(hi)) / edgeW).^2;
H = zeros(size(f));
H(in) = W(in) .* exp(1i * 2*pi * KDM * (f(in) - fRef).^2 ./ (f(in) * fRef^2));
end

function Hf = hEffGrid(Np, fs, varargin)
% hEff on the full Np-point FFT grid (column), conjugate-symmetric (real impulse response).
Hp = hEff((0:Np/2).' * (fs / Np), varargin{:});
Hf = [Hp; conj(Hp(end-1:-1:2))];
end

function [errE, Hm] = transferError(y, k0, fs, varargin)
% The response to an impulse at 0-based sample k0: its FFT re-referenced to k0 (the measured
% H_eff) and the relative error energy against hEff on the same grid.
N  = numel(y);
k  = (0:N-1).';
Hm = fft(double(y)) .* exp(1i * 2*pi * mod(k * k0, N) / N);
Hf = hEffGrid(N, fs, varargin{:});
errE = sum(abs(Hm - Hf).^2) / sum(abs(Hf).^2);
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
