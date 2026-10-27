function testApplyInverseDispersion()
%TESTAPPLYINVERSEDISPERSION  Unit tests for applyInverseDispersion (stage 5: coherent dedispersion).
%{
Run from the PulsarSimMatlab folder: run('tests/testApplyInverseDispersion.m'). ~2 s.
Writes its own files to tempdir (~0.1 GB, deleted at the end); needs no data files.
Ground truth is the help text's filter, computed here in double precision: per baseband
bin, RF f = fLO + f_bb inside [fLow, fHigh],
  H(f_bb) = W(f) * exp(-i*2*pi*K*DM*(f - fRef)^2 / (f*fRef^2)),
the conjugate of the forward chirp (same sin^2 edges W), group delay -(tau(f) - tau(fRef)),
0 outside the band; overlap nFuture = ceil((tau(fLow) - tau(fRef)) fs) + guard, nPast =
ceil((tau(fRef) - tau(fHigh)) fs) + guard (each >= guard), fully supported nPast+1 ..
N-nFuture. main's settings: 800 MHz IQ, fLO 1.3 GHz, band 1.2-1.6 GHz. DM 0.05 on the full
band (63 us smear, 50,420 samples), DM 25 on the channel path.

  1. Format and info (one impulse): N, 8*N bytes, format fields, refFreq default fHigh,
     bulkDelayRef = tau(fRef), edgeWidth, nPast / nFuture / fullySupported as above,
     refOutsideBand, the saved _info.mat = the returned info.
  2. Transfer function (one impulse, 2^19 samples): the FFT of the output re-referenced to
     the impulse's sample against H on the whole -fs/2..fs/2 grid, relative error energy
     < 1e-9; RefFreq fHigh (default), 1.4 GHz, fLow, 1.8 GHz (above the band: a pure advance,
     nPast = guard), 1.1 GHz (below: a pure delay, nFuture = guard; both with
     AllowRefOutsideBand) and DM 0 (H = W); nPast / nFuture as above in every case. Printed:
     the measured group delay at 1.25 ... 1.55 GHz against -(tau(f) - tau(fRef)).
  3. Round trip at RF (white noise at 4 GHz, 1 ms, zero near the file ends so that no part
     of a sweep is cut off): applyDispersionStream -> applyIQmodulation (main's settings) ->
     applyInverseDispersion = the IQ of the input filtered by W^2 (the chirps cancel), from
     the input's FFT in double (band kept, factor 2, shifted by fLO, x W^2): relative rms
     < 1e-4 over the whole file (the IQ low-pass ripple is <= 9e-6 in the band).
  4. Channel path: dedisperseChannels on four fake 3.125 MHz channels (1, 64, 65, 128 of the
     grid, at 4/3 x 3.125 MHz), one impulse each, DM 25, RefFreq 1.4 GHz inside the band (the
     edge between channels 64 and 65): each channel = H with its own band (error energy
     < 1e-9); channels above 1.4 GHz need the past (nPast > guard), those below the future
     (nFuture > guard), as the formulas; fullySupported = [max(s1 + nPast), min(s2 - nFuture)].
     DM 0: each channel = its taper W, nPast = nFuture = guard.
  5. Overlap-save = linear filtering (complex white noise, 2^20 samples): the automatic Nfft,
     the minimum 2^nextpow2(2*overlap), twice that and one block, each against one long FFT
     in double with H (zeros before sample 1 and after the end): relative rms < 1e-5
     (float32 rounding plus the filter's leakage, ~1e-12 in energy).
  6. Leakage (the filter energy outside the overlap lags): every default run below 1e-10
     (MaxLeakage warns at 1e-6) without a warning. info.leakage = the same quantity from the
     formula's impulse response on the run's own Nfft grid, in double, within 1 % + 1e-13:
     at the default guard (~1.2e-12), GuardTime 2/edgeW (~3.6e-8) and 0.1/edgeW (~1.2e-6,
     which warns :leakage; silenced, lastwarn still records it). Check 5's noise error stays
     <= 2 x the leakage in each case: a leaked lag can at most pick one unrelated sample (most
     do not: the buffer still holds their true sample, except at the start of each block).
     Until 27 Oct 2026 the code summed in single and reported 0 or multiples of 6e-8.
  7. Input checks: inFile = outFile, fLow >= fHigh, a band outside fLO +- fs/2, RefFreq
     outside the band without AllowRefOutsideBand, RefFreq <= 0, EdgeFrac > 0.5, a missing
     file, an empty file, Nfft <= the overlap and a too small MaxMemoryGB are refused, each
     with its own id.
Errors at the end if any check fails.
%}

root = fullfile(fileparts(mfilename('fullpath')), '..');
addpath(fullfile(root, 'functions'));
oldDir = cd(root); restoreDir = onCleanup(@() cd(oldDir));
tmp = tempname(tempdir);                         % own folder per run
mkdir(tmp);
cleanTmp = onCleanup(@() rmdir(tmp, 's'));
fails = {};

fs = 800e6; fLO = 1.3e9; fLow = 1.2e9; fHigh = 1.6e9;   % main.m / testParams
K = 4.148808e3;                                  % s MHz^2 pc^-1 cm^3 (help text)
tau = @(f, DM) K * 1e12 * DM ./ f.^2;            % s, f in Hz
edgeW = 0.02 * (fHigh - fLow);                   % default EdgeFrac: 8 MHz
g = ceil(20 / edgeW * fs);                       % default guard: 2.5 us = 2000 samples
nOv = @(DM, fRef, lo, hi, fsx, gx) deal(max(ceil((tau(fRef, DM) - tau(hi, DM)) * fsx), 0) + gx, ...
                                        max(ceil((tau(lo, DM) - tau(fRef, DM)) * fsx), 0) + gx);   % [nPast, nFuture]
fIn = fullfile(tmp, 'in.dat'); fOut = fullfile(tmp, 'out.dat');
ded = @(DM, varargin) applyInverseDispersion(fIn, fOut, fs, fLO, DM, fLow, fHigh, 'Verbose', false, ...
    'SaveInfo', false, varargin{:});
rs = RandStream('mt19937ar', 'Seed', 13);
DM2 = 0.05;
leaks = [];                                      % leakage of every default-guard run (check 6)
warned = false;

% ---------------------------------------------------------------------------------
% 1. Format and info
% ---------------------------------------------------------------------------------
N2 = 2^19; k0 = 2^18;
z2 = complex(zeros(N2, 1, 'single')); z2(k0 + 1) = 1;   % one impulse at 0-based sample k0
writeCF32(fIn, z2);
lastwarn('');
i1 = applyInverseDispersion(fIn, fOut, fs, fLO, DM2, fLow, fHigh, 'Verbose', false);
[~, wId] = lastwarn; warned = warned || strcmp(wId, 'applyInverseDispersion:leakage');
leaks(end+1) = i1.leakage;
f = dir(fOut);
S = load(fullfile(tmp, 'out_info.mat'));
[nP, nF] = nOv(DM2, fHigh, fLow, fHigh, fs, g);
ok = [i1.N == N2, f.bytes == 8*N2, i1.samplesRead == N2, i1.samplesWritten == N2, ...
      strcmp(i1.precision, 'single'), strcmp(i1.byteOrder, 'ieee-le'), i1.isComplex, i1.fs == fs, ...
      i1.actualFsOut == fs, i1.fLO == fLO, i1.DM == DM2, i1.fLow == fLow, i1.fHigh == fHigh, ...
      i1.refFreq == fHigh, ~i1.refOutsideBand, i1.dispersionConst == K, ...
      abs(i1.bulkDelayRef / tau(fHigh, DM2) - 1) < 1e-12, i1.edgeWidth == edgeW, i1.nPast == nP, ...
      i1.nFuture == nF, isequal(i1.fullySupported, [nP + 1, N2 - nF]), i1.Nfft > nP + nF, ...
      isequal(rmfield(i1, 'infoFile'), S.info)];
pass = all(ok);
fprintf(['1. format and info: N %d, %d bytes, refFreq %.4g GHz (default), bulk delay %.6f us, nPast %d ' ...
         '(guard %d), nFuture %d, fully supported %d..%d, Nfft 2^%d: %d of %d checks: %s\n'], i1.N, f.bytes, ...
    i1.refFreq / 1e9, i1.bulkDelayRef * 1e6, i1.nPast, g, i1.nFuture, i1.fullySupported, log2(i1.Nfft), ...
    nnz(ok), numel(ok), passStr(pass));
if ~pass, fails{end+1} = 'format / info'; end

% ---------------------------------------------------------------------------------
% 2. Transfer function from one impulse
% ---------------------------------------------------------------------------------
cases2 = {'RefFreq = fHigh (default)',     DM2, fHigh, {}
          'RefFreq = 1.4 GHz',             DM2, 1.4e9, {'RefFreq', 1.4e9}
          'RefFreq = fLow',                DM2, fLow,  {'RefFreq', fLow}
          'RefFreq = 1.8 GHz (above)',     DM2, 1.8e9, {'RefFreq', 1.8e9, 'AllowRefOutsideBand', true}
          'RefFreq = 1.1 GHz (below)',     DM2, 1.1e9, {'RefFreq', 1.1e9, 'AllowRefOutsideBand', true}
          'DM = 0 (band-pass only)',       0,   fHigh, {}};
f0 = [1.25 1.3 1.4 1.5 1.55] * 1e9;              % where the group delay is printed
df = fs / N2;
kbS = round((f0 - fLO) / df);                    % signed baseband bins
for c = 1:size(cases2, 1)
    [label, dm, fRef, opts2] = cases2{c, :};
    lastwarn('');
    info = ded(dm, opts2{:});
    [~, wId] = lastwarn; warned = warned || strcmp(wId, 'applyInverseDispersion:leakage');
    leaks(end+1) = info.leakage; %#ok<AGROW>
    [errE, Hm] = transferError(readCF32(fOut), k0, fs, dm, fRef, fLO, fLow, fHigh, edgeW);
    tauM = -angle(Hm(mod(kbS + 1, N2) + 1) .* conj(Hm(mod(kbS, N2) + 1))) / (2*pi*df);
    tauP = -(tau(fLO + (kbS + 0.5) * df, dm) - tau(fRef, dm));
    [nP, nF] = nOv(dm, fRef, fLow, fHigh, fs, g);
    okO = info.nPast == nP && info.nFuture == nF && info.refOutsideBand == (fRef < fLow || fRef > fHigh);
    pass = errE < 1e-9 && okO;
    fprintf(['2. transfer function, %s: error energy %.2e (< 1e-9), nPast %d, nFuture %d (formula %d), group ' ...
             'delay at 1.25..1.55 GHz %s us, largest difference from -(tau(f) - tau(fRef)) %.2g ps: %s\n'], ...
        label, errE, info.nPast, info.nFuture, okO, sprintf('%.4f ', tauM * 1e6), ...
        max(abs(tauM(:) - tauP(:))) * 1e12, passStr(pass));
    if ~pass, fails{end+1} = ['transfer function, ' label]; end %#ok<AGROW>
end

% ---------------------------------------------------------------------------------
% 3. Round trip at RF
% ---------------------------------------------------------------------------------
f_in = 4e9; Nx = 4e6; NB = Nx / 5;               % 1 ms at 4 GHz -> 800,000 IQ samples
x = randn(rs, Nx, 1) .* endTaper(Nx, 3e5, 1e5);  % zero for 0.3e6 samples at each end (75 us)
xs = single(x); clear x
fRF = fullfile(tmp, 'rf.dat'); fDisp = fullfile(tmp, 'disp.dat'); fIQ = fullfile(tmp, 'iq.dat');
writeF32(fRF, xs);
applyDispersionStream(fRF, fDisp, DM2, f_in, fLow, fHigh, 'Verbose', false, 'SaveInfo', false);
applyIQmodulation(fDisp, fIQ, f_in, fs, fLO, 'Verbose', false, 'SaveInfo', false);
lastwarn('');
iR = applyInverseDispersion(fIQ, fOut, fs, fLO, DM2, fLow, fHigh, 'Verbose', false, 'SaveInfo', false);
[~, wId] = lastwarn; warned = warned || strcmp(wId, 'applyInverseDispersion:leakage');
leaks(end+1) = iR.leakage;
zR = readCF32(fOut);
X = fft(double(xs)); clear xs
m = (ceil(fLow * Nx / f_in) : floor(fHigh * Nx / f_in)).';   % positive RF bins in the band (0-based)
Zb = zeros(NB, 1);
Zb(mod(m - fLO * Nx / f_in, NB) + 1) = X(m + 1) .* wTaper(m * (f_in / Nx), fLow, fHigh, edgeW).^2;
clear X
zRef = (2 * NB / Nx) * ifft(Zb);                 % 2 x the positive-frequency part, at fLO, x W^2
e3 = sqrt(mean(abs(zR - zRef).^2) / mean(abs(zRef).^2));
pass = numel(zR) == NB && e3 < 1e-4;
fprintf(['3. round trip RF -> dispersion -> IQ -> dedispersion (DM %.2f): = IQ of the input x W^2 to ' ...
         'relative rms %.2e (< 1e-4) over all %d samples: %s\n'], DM2, e3, numel(zR), passStr(pass));
if ~pass, fails{end+1} = 'round trip'; end
clear zR zRef Zb

% ---------------------------------------------------------------------------------
% 4. Channel path (dedisperseChannels, RefFreq inside the band)
% ---------------------------------------------------------------------------------
dF = 3.125e6; fsC = dF * 4/3; N4 = 2^18; k4 = 2^17; DM4 = 25; fR4 = 1.4e9;
jj = [1 64 65 128];
fc = fLow + (jj(:) - 0.5) * dF;                  % channel centres of the 128-channel grid
chanFiles = strings(numel(jj), 1);
z4 = complex(zeros(N4, 1, 'single')); z4(k4 + 1) = 1;
for c = 1:numel(jj)
    chanFiles(c) = fullfile(tmp, sprintf('chan%03d.dat', jj(c)));
    writeCF32(chanFiles(c), z4);
end
s4 = [101, N4 - 50];                             % the channelizer's supported range (made up)
info_chan = struct('nChan', numel(jj), 'chanWidth', dF, 'chanFreqs', fc, 'chanFiles', chanFiles, ...
    'fs', fsC, 'N', N4, 'fullySupported', s4, 'file', "fake", 'fLow', fLow, 'fHigh', fHigh, 't0', 0);
gC = ceil(20 / (0.02 * dF) * fsC);               % per-channel guard: 320 us
pass = true;
for dm = [DM4, 0]
    lastwarn('');
    d4 = dedisperseChannels(info_chan, fullfile(tmp, sprintf('dd%g', dm)), dm, 'RefFreq', fR4, ...
        'SaveInfo', false, 'Verbose', false);
    [~, wId] = lastwarn; warned = warned || strcmp(wId, 'applyInverseDispersion:leakage');
    leaks = [leaks, d4.leakage(:).']; %#ok<AGROW>
    [nP, nF, eC] = deal(zeros(numel(jj), 1));
    for c = 1:numel(jj)
        lo = fc(c) - dF/2; hi = fc(c) + dF/2;
        [nP(c), nF(c)] = nOv(dm, fR4, lo, hi, fsC, ceil(20 / (0.02 * (hi - lo)) * fsC));
        eC(c) = transferError(readCF32(d4.chanFiles(c)), k4, fsC, dm, fR4, fc(c), lo, hi, 0.02 * (hi - lo));
    end
    above = fc > fR4;
    okC = all(eC < 1e-9) && isequal(d4.nPast, nP) && isequal(d4.nFuture, nF) && ...
          isequal(d4.fullySupported, [max(s4(1) + nP), min(s4(2) - nF)]);
    if dm > 0
        okC = okC && all(nP(above) > gC) && all(nF(above) == gC) && all(nF(~above) > gC) && all(nP(~above) == gC);
    else
        okC = okC && all(nP == gC) && all(nF == gC);
    end
    pass = pass && okC;
    fprintf(['4. channel path, DM %g, RefFreq 1.4 GHz, channels %s: error energy at most %.2e (< 1e-9), ' ...
             'nPast %s, nFuture %s (guard %d), fully supported %d..%d: %s\n'], dm, sprintf('%d ', jj), max(eC), ...
        sprintf('%d ', nP), sprintf('%d ', nF), gC, d4.fullySupported, passStr(okC));
end
if ~pass, fails{end+1} = 'channel path'; end

% ---------------------------------------------------------------------------------
% 5. Overlap-save = linear filtering
% ---------------------------------------------------------------------------------
N5 = 2^20;
z5 = double(single(complex(randn(rs, N5, 1), randn(rs, N5, 1)) / sqrt(2)));
writeCF32(fIn, z5);
[nP5, nF5] = nOv(DM2, fHigh, fLow, fHigh, fs, g);
nO = nP5 + nF5;
Np = 2^nextpow2(N5 + nO + 1);                    % zero padding >= the filter's reach on both sides
y5 = ifft(fft(z5, Np) .* hInv(gridF(Np, fs), DM2, fHigh, fLO, fLow, fHigh, edgeW));
y5 = y5(1:N5);                                   % linear filtering, zeros before and after the file
P5 = mean(abs(y5).^2);
nMin = 2^nextpow2(2 * nO);
cases5 = {'automatic', {}; sprintf('2^%d (minimum)', log2(nMin)), {'Nfft', nMin}; ...
          sprintf('2^%d', log2(2*nMin)), {'Nfft', 2*nMin}; sprintf('2^%d (one block)', nextpow2(N5 + nO)), ...
          {'Nfft', 2^nextpow2(N5 + nO)}};
pass = true;
for c = 1:size(cases5, 1)
    lastwarn('');
    info = ded(DM2, cases5{c, 2}{:});
    [~, wId] = lastwarn; warned = warned || strcmp(wId, 'applyInverseDispersion:leakage');
    leaks(end+1) = info.leakage; %#ok<AGROW>
    e5 = sqrt(mean(abs(readCF32(fOut) - y5).^2) / P5);
    if c == 1, iAuto = info; eAuto = e5^2; end   % for check 6
    okF = e5 < 1e-5;
    pass = pass && okF;
    fprintf('5. overlap-save, Nfft %s: %d blocks, relative rms from one double FFT %.2e (< 1e-5): %s\n', ...
        cases5{c, 1}, info.nBlocks, e5, passStr(okF));
end
if ~pass, fails{end+1} = 'overlap-save'; end

% ---------------------------------------------------------------------------------
% 6. Leakage
% ---------------------------------------------------------------------------------
L0 = leakFormula(iAuto, fs, DM2, fHigh, fLO, fLow, fHigh, edgeW);
okL = abs(iAuto.leakage - L0) <= 0.01 * L0 + 1e-13 && eAuto <= 2 * iAuto.leakage;
pass = all(leaks < 1e-10) && ~warned && okL;
fprintf(['6. leakage, default guard: at most %.1e over %d runs (< 1e-10), warnings %d; automatic run %.3e, ' ...
         'formula %.3e; noise error energy %.2e = %.4f x leakage (<= 2): %s\n'], max(leaks), numel(leaks), warned, ...
    iAuto.leakage, L0, eAuto, eAuto / iAuto.leakage, passStr(pass));
for gt = [2 0.1]
    s = warning('off', 'applyInverseDispersion:leakage');
    lastwarn('');
    i6 = ded(DM2, 'GuardTime', gt / edgeW);
    [~, wId6] = lastwarn;
    warning(s);
    [nP6, nF6] = nOv(DM2, fHigh, fLow, fHigh, fs, ceil(gt / edgeW * fs));
    L6 = leakFormula(i6, fs, DM2, fHigh, fLO, fLow, fHigh, edgeW);
    e6 = mean(abs(readCF32(fOut) - y5).^2) / P5;   % y5 is the full filter: it does not depend on the guard
    warns = strcmp(wId6, 'applyInverseDispersion:leakage');
    okG = i6.nPast == nP6 && i6.nFuture == nF6 && abs(i6.leakage - L6) <= 0.01 * L6 + 1e-13 && ...
          e6 <= 2 * i6.leakage && warns == (L6 > 1e-6);
    pass = pass && okG;
    fprintf(['6. leakage, GuardTime %g/edgeW (%d samples): %.3e, formula %.3e, warns %d (expected %d); noise ' ...
             'error energy %.2e = %.4f x leakage (<= 2): %s\n'], gt, ceil(gt / edgeW * fs), i6.leakage, L6, ...
        warns, L6 > 1e-6, e6, e6 / i6.leakage, passStr(okG));
end
if ~pass, fails{end+1} = 'leakage'; end
clear z5 y5

% ---------------------------------------------------------------------------------
% 7. Input checks
% ---------------------------------------------------------------------------------
writeCF32(fIn, z2);
fEmpty = fullfile(tmp, 'empty.dat');
fclose(fopen(fEmpty, 'w'));
bad = @(fi, fo, lo, hi, varargin) errId(@() applyInverseDispersion(fi, fo, fs, fLO, DM2, lo, hi, ...
    'Verbose', false, 'SaveInfo', false, varargin{:}));
ids = {bad(fIn, fIn, fLow, fHigh), bad(fIn, fOut, fHigh, fLow), bad(fIn, fOut, fLow, 1.75e9), ...
       bad(fIn, fOut, fLow, fHigh, 'RefFreq', 1.8e9), bad(fIn, fOut, fLow, fHigh, 'RefFreq', -1, ...
       'AllowRefOutsideBand', true), bad(fIn, fOut, fLow, fHigh, 'EdgeFrac', 0.6), ...
       bad(fullfile(tmp, 'none.dat'), fOut, fLow, fHigh), bad(fEmpty, fOut, fLow, fHigh), ...
       bad(fIn, fOut, fLow, fHigh, 'Nfft', i1.nPast + i1.nFuture), bad(fIn, fOut, fLow, fHigh, 'MaxMemoryGB', 1e-4)};
want = strcat('applyInverseDispersion:', {'sameFile', 'band', 'band', 'ref', 'ref', 'edge', 'notFound', ...
       'empty', 'nfft', 'memory'});
pass = isequal(ids, want);
fprintf('7. input checks: %d of %d bad inputs refused with the right id: %s\n', ...
    nnz(strcmp(ids, want)), numel(want), passStr(pass));
if ~pass, fails{end+1} = 'input checks'; end

% ---------------------------------------------------------------------------------
if isempty(fails)
    fprintf('testApplyInverseDispersion: ALL PASSED\n');
else
    error('testApplyInverseDispersion:failed', 'FAILED: %s', strjoin(fails, ', '));
end
end


% =================================================================================
function W = wTaper(f, fLow, fHigh, edgeW)
% The help text's band taper at RF frequencies f: 1 inside, sin^2 edges of width edgeW
% inside fLow and fHigh, 0 outside.
in = f > fLow & f < fHigh;
W = double(in);
lo = in & f < fLow + edgeW;  W(lo) = sin(pi/2 * (f(lo) - fLow) / edgeW).^2;
hi = in & f > fHigh - edgeW; W(hi) = sin(pi/2 * (fHigh - f(hi)) / edgeW).^2;
end

function H = hInv(fbb, DM, fRef, fLO, fLow, fHigh, edgeW)
% The help text's inverse filter at baseband frequencies fbb (RF f = fLO + fbb), from the
% formula alone (not the code): W(f) * exp(-i 2 pi K DM (f - fRef)^2 / (f fRef^2)).
f = fLO + fbb;
KDM = 4.148808e3 * 1e12 * DM;
H = wTaper(f, fLow, fHigh, edgeW) .* exp(-1i * 2*pi * KDM * (f - fRef).^2 ./ (f * fRef^2));
end

function f = gridF(N, fs)
% Signed frequencies of the N FFT bins (column), -fs/2 <= f < fs/2.
b = (0:N-1).';
b(b >= N/2) = b(b >= N/2) - N;
f = b * (fs / N);
end

function [errE, Hm] = transferError(y, k0, fs, DM, fRef, fLO, fLow, fHigh, edgeW)
% The response to an impulse at 0-based sample k0: its FFT re-referenced to k0 (the measured
% H) and the relative error energy against hInv on the same grid.
N  = numel(y);
k  = (0:N-1).';
Hm = fft(y) .* exp(1i * 2*pi * mod(k * k0, N) / N);
Hf = hInv(gridF(N, fs), DM, fRef, fLO, fLow, fHigh, edgeW);
errE = sum(abs(Hm - Hf).^2) / sum(abs(Hf).^2);
end

function L = leakFormula(info, fs, DM, fRef, fLO, fLow, fHigh, edgeW)
% The filter energy outside the overlap lags, from the formula on the run's own Nfft grid,
% in double: the impulse response h = ifft(H), lags 0..nPast and -nFuture..-1 kept.
N = info.Nfft;
h = ifft(hInv(gridF(N, fs), DM, fRef, fLO, fLow, fHigh, edgeW));
inside = sum(abs(h(1:info.nPast+1)).^2) + sum(abs(h(N-info.nFuture+1:N)).^2);
L = max(0, 1 - inside / sum(abs(h).^2));
end

function w = endTaper(N, z, r)
% 0 for z samples at each end, raised-cosine ramps of r samples, 1 in between (column).
w = zeros(N, 1);
u = (0:r-1).' / r;
w(z+1 : z+r) = 0.5 * (1 - cos(pi * u));
w(z+r+1 : N-z-r) = 1;
w(N-z-r+1 : N-z) = flipud(w(z+1 : z+r));
end

function writeF32(f, x)
fid = fopen(f, 'w', 'ieee-le'); fwrite(fid, x, 'single'); fclose(fid);
end

function writeCF32(f, z)
fid = fopen(f, 'w', 'ieee-le'); fwrite(fid, [real(z(:)).'; imag(z(:)).'], 'single'); fclose(fid);
end

function z = readCF32(f)
% Interleaved I,Q float32 -> complex double column.
fid = fopen(f, 'r', 'ieee-le'); v = fread(fid, [2 Inf], 'single=>double'); fclose(fid);
z = complex(v(1, :), v(2, :)).';
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
