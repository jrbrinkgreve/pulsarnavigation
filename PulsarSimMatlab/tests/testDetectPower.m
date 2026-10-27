function testDetectPower()
%TESTDETECTPOWER  Unit tests for detectPower (stage 6: square-law detection).
%{
Run from the PulsarSimMatlab folder: run('tests/testDetectPower.m'). ~2 s.
Writes its own files to tempdir (~0.1 GB, deleted at the end); needs no data files.
Ground truth is the help text, computed here: bin k (1-based) = the mean of |z|^2 = I^2 + Q^2
over input samples (k-1)*binLen+1 .. k*binLen, binLen = round(fs/f_out), at the centroid time
t_k = T0 + ((k-1)*binLen + (binLen-1)/2)/fs; a trailing partial bin is dropped. For complex
Gaussian noise with autocovariance C(L) the variance of a bin is exactly
(1/binLen) sum_L (1 - |L|/binLen) |C(L)|^2. main's settings: 800 MHz IQ, fLO 1.3 GHz,
f_out 1 MHz (binLen 800).

  1. Format and info (800,123 random samples, T0 1.234 us, FullySupported [2002, 790000]):
     N = floor(N_in/binLen) = 1000, 123 samples dropped, binLen, rates, binDt, binTime0 =
     T0 + (binLen-1)/(2 fs), the single-channel fields, chanMeanPower = the mean of the bins,
     fullySupportedBins = the bins wholly inside the sample range (found by brute force), the
     saved _info.mat = the returned info.
  2. Exactness: every bin = the mean of I^2 + Q^2 over its samples, computed here in double,
     to float32 rounding (1e-6); BlockSamples 1001 (one bin per block) and 300,007 give the
     same file as the default (one block here), byte for byte.
  3. Timing: a Gaussian power burst (sigma 20 us, at a fractional sample) with T0 1.234 us:
     the power-weighted centroid of the bins at their reported times = the burst's own
     centroid within 1e-3 bin (the symmetric boxcar does not bias it).
  4. Statistics (2^23 samples, 10,485 bins): white complex noise (variance 2): bin mean =
     sigma^2 and bin variance = sigma^4/binLen, each within 4 sigma, neighbouring bins
     uncorrelated; noise shaped by the dedispersion taper W (band 1.2-1.6 GHz, 8 MHz sin^2
     edges, made by FFT here): mean = mean(W^2), variance = the exact sum above from C(L) =
     ifft(W^2), lag-1 correlation as predicted, each within 4 sigma. Printed: the radiometer
     approximation P^2/(B_noise dt), B_noise = (int W^2)^2 / int W^4.
  5. Sub-band mode ('NChan' > 1; unused by main and the tests, review list D1-26): 8 channels
     over the whole baseband add up to the total power (Parseval) within 1e-5; a tone on an
     FFT bin of channel 3 of 8 (1.2-1.6 GHz) puts its power A^2 there and nothing elsewhere;
     no Band and empty channels are refused (the resolution warning silenced).
  6. Warnings and input checks: f_out 0.7 MHz (not fs/integer) warns :rate (silenced;
     lastwarn still records it); inFile = outFile, a missing file, a file shorter than one bin
     and a Band outside fLO +- fs/2 are refused, each with its own id.
Errors at the end if any check fails.
%}

root = fullfile(fileparts(mfilename('fullpath')), '..');
addpath(fullfile(root, 'functions'));
oldDir = cd(root); restoreDir = onCleanup(@() cd(oldDir));
tmp = tempname(tempdir);                         % own folder per run
mkdir(tmp);
cleanTmp = onCleanup(@() rmdir(tmp, 's'));
fails = {};

fs = 800e6; fLO = 1.3e9; f_out = 1e6;            % main.m / testParams
binLen = round(fs / f_out);                      % 800
fIn = fullfile(tmp, 'in.dat'); fOut = fullfile(tmp, 'out.dat');
dp = @(varargin) detectPower(fIn, fOut, fs, fLO, f_out, 'Verbose', false, 'SaveInfo', false, varargin{:});
rs = RandStream('mt19937ar', 'Seed', 17);

% ---------------------------------------------------------------------------------
% 1. Format and info
% ---------------------------------------------------------------------------------
Nx1 = 1000 * binLen + 123; T0 = 1.234e-6; s1 = [2002, 790000];
z1 = single(complex(randn(rs, Nx1, 1), randn(rs, Nx1, 1)));
writeCF32(fIn, z1);
lastwarn('');
i1 = detectPower(fIn, fOut, fs, fLO, f_out, 'T0', T0, 'FullySupported', s1, 'Verbose', false);
[~, wId1] = lastwarn;
p1 = readF32(fOut);
f = dir(fOut);
S = load(fullfile(tmp, 'out_info.mat'));
k = (1:1000).';
inside = find((k - 1) * binLen + 1 >= s1(1) & k * binLen <= s1(2));   % bins wholly inside [s1(1), s1(2)]
ok = [i1.N == 1000, numel(p1) == 1000, f.bytes == 4000, i1.samplesDropped == 123, ...
      i1.samplesUsed == 1000 * binLen, i1.binLen == binLen, i1.actualFout == f_out, i1.fs == f_out, ...
      i1.fsIn == fs, i1.binDt == binLen / fs, abs(i1.binTime0 - (T0 + (binLen - 1) / (2*fs))) < 1e-18, ...
      i1.T0 == T0, i1.nChan == 1, i1.fLO == fLO, i1.chanFreqs == fLO, ...
      isequal(i1.chanEdges, [fLO - fs/2, fLO + fs/2]), i1.frameLen == binLen, ...
      abs(i1.chanMeanPower / mean(double(p1)) - 1) < 1e-6, ...
      isequal(i1.fullySupportedBins, [inside(1), inside(end)]), isempty(wId1), ...
      isequal(rmfield(i1, 'infoFile'), S.info)];
pass = all(ok);
fprintf(['1. format and info: %d bins of %d samples (%d dropped), %.6g Hz, bin 1 at %.6f us, fully supported ' ...
         'bins %d..%d (brute force %d..%d): %d of %d checks: %s\n'], i1.N, i1.binLen, i1.samplesDropped, ...
    i1.actualFout, i1.binTime0 * 1e6, i1.fullySupportedBins, inside(1), inside(end), nnz(ok), numel(ok), passStr(pass));
if ~pass, fails{end+1} = 'format / info'; end

% ---------------------------------------------------------------------------------
% 2. Exactness
% ---------------------------------------------------------------------------------
q = double(real(z1(1:1000*binLen))).^2 + double(imag(z1(1:1000*binLen))).^2;
pRef = mean(reshape(q, binLen, []), 1).';        % the mean of I^2 + Q^2 per bin, in double
e2 = max(abs(double(p1) - pRef) ./ pRef);
ref = fileBytes(fOut);
bsz = [1001, 300007];
same = false(size(bsz)); nBlk = zeros(size(bsz));
for b = 1:numel(bsz)
    ib = detectPower(fIn, fOut, fs, fLO, f_out, 'T0', T0, 'FullySupported', s1, 'BlockSamples', bsz(b), ...
        'Verbose', false, 'SaveInfo', false);
    same(b) = isequal(fileBytes(fOut), ref);
    nBlk(b) = ceil(ib.N / max(1, round(bsz(b) / binLen)));
end
pass = e2 < 1e-6 && all(same);
fprintf(['2. exactness: largest relative difference from the mean of I^2 + Q^2 per bin %.1e (< 1e-6); ' ...
         'BlockSamples %s (%s blocks) = the default byte for byte %s: %s\n'], e2, mat2str(bsz), mat2str(nBlk), ...
    mat2str(same), passStr(pass));
if ~pass, fails{end+1} = 'exactness'; end
clear z1 q pRef

% ---------------------------------------------------------------------------------
% 3. Timing
% ---------------------------------------------------------------------------------
Nx3 = 200 * binLen;
n = (0:Nx3-1).';
nc = 81234.37; sgN = 20e-6 * fs;                 % centre (fractional sample) and sigma (16,000 samples)
z3 = single(sqrt(exp(-0.5 * ((n - nc) / sgN).^2)) .* exp(2i*pi * 0.1234 * n));   % |z|^2 = the burst
writeCF32(fIn, z3);
i3 = dp('T0', T0);
P3 = double(readF32(fOut));
tk = i3.binTime0 + (0:i3.N-1).' * i3.binDt;
pw = double(real(z3)).^2 + double(imag(z3)).^2;  % the burst as written
tTrue = T0 + sum(n .* pw) / sum(pw) / fs;
dT = (sum(tk .* P3) / sum(P3) - tTrue) / i3.binDt;
pass = abs(dT) < 1e-3;
fprintf('3. timing: burst centroid %.6f us (T0 %.3f us); bin centroid minus it %+.1e bin (< 1e-3): %s\n', ...
    tTrue * 1e6, T0 * 1e6, dT, passStr(pass));
if ~pass, fails{end+1} = 'timing'; end
clear z3 n pw

% ---------------------------------------------------------------------------------
% 4. Statistics
% ---------------------------------------------------------------------------------
N4 = 2^23;
sig2 = 2;
writeCF32(fIn, single(sqrt(sig2/2) * complex(randn(rs, N4, 1), randn(rs, N4, 1))));
dp();
P = double(readF32(fOut));
nB = numel(P);
zM = (mean(P) / sig2 - 1) * sqrt(nB * binLen);   % the mean of nB*binLen exponential samples
zV = (var(P) / (sig2^2 / binLen) - 1) / sqrt(2 / (nB - 1));
r1 = corr1(P);
okW = abs(zM) < 4 && abs(zV) < 4 && abs(r1) * sqrt(nB) < 4;
fprintf(['4. statistics, white noise (%d bins): mean / sigma^2 %.5f (%+.2f sigma), variance / (sigma^4/binLen) ' ...
         '%.4f (%+.2f sigma), lag-1 correlation %+.4f (+- %.4f): %s\n'], nB, mean(P) / sig2, zM, ...
    var(P) / (sig2^2 / binLen), zV, r1, 1 / sqrt(nB), passStr(okW));
% noise shaped by the dedispersion taper W
fb = gridF(N4, fs);
W2 = wTaper(fLO + fb, 1.2e9, 1.6e9, 8e6).^2;     % |W|^2 per FFT bin
z4 = ifft(fft(complex(randn(rs, N4, 1), randn(rs, N4, 1)) / sqrt(2)) .* sqrt(W2));
writeCF32(fIn, single(z4)); clear z4
dp();
P = double(readF32(fOut));
nB = numel(P);
C = ifft(W2);                                    % autocovariance C(L), L = 0..N4-1 (circular)
L = (0:binLen-1).';
vPred = (abs(C(1))^2 + 2 * sum((1 - L(2:end) / binLen) .* abs(C(L(2:end) + 1)).^2)) / binLen;
d = (1:2*binLen-1).';                            % lags between the samples of neighbouring bins
cPred = sum(min(d, 2*binLen - d) .* abs(C(d + 1)).^2) / binLen^2;
mPred = mean(W2);
Bn = fs * mean(W2)^2 / mean(W2.^2);              % noise bandwidth (int W^2)^2 / int W^4
zM = (mean(P) / mPred - 1) / (sqrt(vPred / nB) / mPred);
zV = (var(P) / vPred - 1) / sqrt(2 / (nB - 1));
r1 = corr1(P);
okS = abs(zM) < 4 && abs(zV) < 4 && abs(r1 - cPred / vPred) * sqrt(nB) < 4;
fprintf(['4. statistics, taper-shaped noise: mean / mean(W^2) %.5f (%+.2f sigma), variance / exact %.4f ' ...
         '(%+.2f sigma; radiometer P^2/(B_noise dt) with B_noise %.2f MHz is %.4f of exact), lag-1 correlation ' ...
         '%+.4f (predicted %+.4f, +- %.4f): %s\n'], mean(P) / mPred, zM, var(P) / vPred, zV, Bn / 1e6, ...
    mPred^2 / (Bn * binLen / fs) / vPred, r1, cPred / vPred, 1 / sqrt(nB), passStr(okS));
pass = okW && okS;
if ~pass, fails{end+1} = 'statistics'; end
clear P C W2 fb

% ---------------------------------------------------------------------------------
% 5. Sub-band mode (NChan > 1)
% ---------------------------------------------------------------------------------
Nx5 = 500 * binLen;
writeCF32(fIn, single(complex(randn(rs, Nx5, 1), randn(rs, Nx5, 1))));
dp();
pTot = double(readF32(fOut));
i5 = dp('NChan', 8, 'Band', [fLO - fs/2, fLO + fs/2]);
p5 = reshape(double(readF32(fOut)), 8, []);
e5 = max(abs(sum(p5, 1).' - pTot) ./ pTot);
A5 = 1.5; fT = 1.325e9;                          % the centre of channel 3 of 8 over 1.2-1.6 GHz
writeCF32(fIn, single(A5 * exp(2i*pi * (fT - fLO) / fs * (0:Nx5-1).')));
iT = dp('NChan', 8, 'Band', [1.2e9 1.6e9]);
pT = reshape(double(readF32(fOut)), 8, []);
okT = abs(mean(pT(3, :)) / A5^2 - 1) < 1e-5 && max(max(pT([1:2 4:8], :))) < 1e-8 * A5^2 && ...
      mod((fT - fLO) / (fs / iT.frameLen), 1) == 0;
s = warning('off', 'detectPower:resolution');
ids5 = {errId(@() dp('NChan', 4)), errId(@() dp('NChan', 800, 'Band', [1.2e9 1.6e9]))};
warning(s);
ok5 = isequal(ids5, {'detectPower:band', 'detectPower:emptyChannel'});
pass = e5 < 1e-5 && okT && ok5;
fprintf(['5. sub-band mode: 8 channels over the whole baseband (frame %d) add up to the total power within %.1e ' ...
         '(< 1e-5); a tone at 1.325 GHz (frame %d, on an FFT bin): channel 3 holds %.6f of A^2, the others at most ' ...
         '%.1e of A^2; no Band and empty channels refused %d: %s\n'], i5.frameLen, e5, iT.frameLen, ...
    mean(pT(3, :)) / A5^2, max(max(pT([1:2 4:8], :))) / A5^2, ok5, passStr(pass));
if ~pass, fails{end+1} = 'sub-band mode'; end

% ---------------------------------------------------------------------------------
% 6. Warnings and input checks
% ---------------------------------------------------------------------------------
wR = warnId(@() detectPower(fIn, fOut, fs, fLO, 0.7e6, 'Verbose', false, 'SaveInfo', false), 'detectPower:rate');
fShort = fullfile(tmp, 'short.dat');
writeCF32(fShort, complex(ones(binLen - 1, 1, 'single')));
bad = @(fi, fo, varargin) errId(@() detectPower(fi, fo, fs, fLO, f_out, 'Verbose', false, 'SaveInfo', false, varargin{:}));
ids = {bad(fIn, fIn), bad(fullfile(tmp, 'none.dat'), fOut), bad(fShort, fOut), ...
       bad(fIn, fOut, 'NChan', 2, 'Band', [0.8e9 1.2e9])};
want = strcat('detectPower:', {'sameFile', 'notFound', 'short', 'band'});
pass = strcmp(wR, 'detectPower:rate') && isequal(ids, want);
fprintf(['6. warnings and input checks: f_out 0.7 MHz warns :rate %d; %d of %d bad inputs refused with the ' ...
         'right id: %s\n'], strcmp(wR, 'detectPower:rate'), nnz(strcmp(ids, want)), numel(want), passStr(pass));
if ~pass, fails{end+1} = 'warnings / input checks'; end

% ---------------------------------------------------------------------------------
if isempty(fails)
    fprintf('testDetectPower: ALL PASSED\n');
else
    error('testDetectPower:failed', 'FAILED: %s', strjoin(fails, ', '));
end
end


% =================================================================================
function r = corr1(x)
% Lag-1 correlation coefficient of a series.
x = x - mean(x);
r = sum(x(1:end-1) .* x(2:end)) / sum(x.^2);
end

function W = wTaper(f, fLow, fHigh, edgeW)
% The dedispersion stage's band taper at RF frequencies f: 1 inside, sin^2 edges of width
% edgeW inside fLow and fHigh, 0 outside (applyInverseDispersion's help text).
in = f > fLow & f < fHigh;
W = double(in);
lo = in & f < fLow + edgeW;  W(lo) = sin(pi/2 * (f(lo) - fLow) / edgeW).^2;
hi = in & f > fHigh - edgeW; W(hi) = sin(pi/2 * (fHigh - f(hi)) / edgeW).^2;
end

function f = gridF(N, fs)
% Signed frequencies of the N FFT bins (column), -fs/2 <= f < fs/2.
b = (0:N-1).';
b(b >= N/2) = b(b >= N/2) - N;
f = b * (fs / N);
end

function writeCF32(f, z)
fid = fopen(f, 'w', 'ieee-le'); fwrite(fid, [real(z(:)).'; imag(z(:)).'], 'single'); fclose(fid);
end

function y = readF32(f)
fid = fopen(f, 'r', 'ieee-le'); y = fread(fid, Inf, 'single=>single'); fclose(fid);
end

function b = fileBytes(f)
fid = fopen(f, 'r'); b = fread(fid, Inf, 'uint8=>uint8'); fclose(fid);
end

function id = warnId(f, expected)
% The identifier of the last warning f() raises ('' if none), with the expected one silenced
% (lastwarn still records a disabled warning).
s = warning('off', expected);
restore = onCleanup(@() warning(s));
lastwarn('');
f();
[~, id] = lastwarn;
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
