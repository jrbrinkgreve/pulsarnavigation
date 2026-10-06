%TESTDEDISPERSECHANNELS  Unit tests for dedisperseChannels (+ applyInverseDispersion option).
%{
Run from the PulsarSimMatlab folder: run('tests/testDedisperseChannels.m').
Needs main.m's receiver files (data/test_rx_IQ.dat, data/test_IQ_dedispersed.dat,
seed 43 reference). Writes channel files to data/chan/ (~1.3 GB) and
temporary files to tempdir. Takes ~1-2 min.

  1. Regression: applyInverseDispersion with default options (the full-band
     path of main.m) reproduces data/test_IQ_dedispersed.dat bit for bit
     after the AllowRefOutsideBand change.
  2. Commutation: channelize(full-band dedispersed) and
     dedisperseChannels(channelize(IQ)) apply the same linear filters in a
     different order, so inside each channel's flat band (|f| <= dF/2 - edge)
     their spectra must agree: complex correlation |rho| ~ 1, phase ~ 0, equal
     power. Checks chirp segments, inter-channel delays (RefFreq outside the
     channel band) and phase reference at once. Outer channels see the
     full-band path's 8 MHz edge taper, so only interior channels are tested.
  3. End to end: per-channel detectPower, channel powers summed, fold, TOAs
     and detection; compared with the full-band chain on the same data
     (same noise): TOA differences; SNR against the expected loss from the
     channel-edge tapers (simulation only: the simulated pulsar also passed
     the forward-dispersion taper, which overlaps the full-band dedispersion
     taper; on real data both paths collect the same int W^2 = 390.0 MHz);
     fold noise ratio within its statistical scatter.
  4. Detected-power noise per time bin: in a 3 MHz channel the power has a
     correlation time ~1/B_ch = 0.33 us, not negligible against 0.96 us bins.
     Exact prediction from the channel spectrum |W_ch|^2: variance below the
     radiometer value m^2/(B*dt) and a positive lag-1 correlation (they
     nearly cancel within a phase bin of ~5 time bins). Checked on the
     off-pulse bins of the summed channel power.
Errors at the end if any check fails.
%}

root = fullfile(fileparts(mfilename('fullpath')), '..');
addpath(fullfile(root, 'functions'));
oldDir = cd(root); restoreDir = onCleanup(@() cd(oldDir));
pipelineParams;                                   % fLow, fHigh, refFreq, ephem, f_out, nBin, ...
tmp = fullfile(tempdir, 'testDedisperseChannels');
if ~isfolder(tmp), mkdir(tmp); end
chanDir = fullfile(dataDir, 'chan');
fails = {};

info_IQ  = loadInfo(fileIQ);
info_ref = loadInfo(fileDedisp);

% ---------------------------------------------------------------------------------
% 1. Regression of the full-band path
% ---------------------------------------------------------------------------------
fileRegr = fullfile(tmp, 'regr_dedisp.dat');
d0 = applyInverseDispersion(info_IQ.file, fileRegr, info_IQ.actualFsOut, info_IQ.fLO, ...
    ephem.DM, fLow, fHigh, 'RefFreq', refFreq, 'Verbose', false, 'SaveInfo', false);
same = filesEqual(fileRegr, fileDedisp);
sameSizes = isequal([d0.nPast d0.nFuture d0.Nfft], [info_ref.nPast info_ref.nFuture info_ref.Nfft]);
delete(fileRegr);
pass = same && sameSizes;
fprintf('1. regression: full-band dedispersion bit-identical %d, nPast/nFuture/Nfft equal %d: %s\n', ...
    same, sameSizes, passStr(pass));
if ~pass, fails{end+1} = 'regression'; end

% ---------------------------------------------------------------------------------
% 2. Commutation: channelize(dedisperse) vs dedisperse-per-channel(channelize)
% ---------------------------------------------------------------------------------
info_chan = channelizeIQ(info_IQ.file, fullfile(chanDir, 'test_rx_IQ_chan'), ...
    info_IQ.actualFsOut, info_IQ.fLO, fLow, fHigh, 'T0', info_IQ.t0);
info_dc = dedisperseChannels(info_chan, fullfile(chanDir, 'test_dedisp_chan'), ephem.DM, ...
    'RefFreq', refFreq);
info_fbc = channelizeIQ(fileDedisp, fullfile(chanDir, 'test_dedispFB_chan'), ...
    info_ref.fs, info_ref.fLO, fLow, fHigh, 'T0', 0, 'Verbose', false);

nSeg = 2^17; m0 = 50000;                          % inside both supported ranges
assert(m0 + 1 >= info_dc.fullySupported(1) && m0 + nSeg <= info_dc.fullySupported(2));
D = info_chan.decimation;
assert(m0*D >= info_ref.fullySupported(1) + info_chan.filterLen && ...
       (m0 + nSeg)*D <= info_ref.fullySupported(2) - info_chan.filterLen);
dF  = info_chan.chanWidth;
fAx = [0:nSeg/2-1, -nSeg/2:-1].' * info_chan.fs / nSeg;
inF = abs(fAx) <= dF/2 - info_dc.edgeWidth;
rho = zeros(info_chan.nChan, 1); pr = zeros(info_chan.nChan, 1);
for j = 1:info_chan.nChan
    A = fft(double(readCF32(info_fbc.chanFiles(j), m0, nSeg)));   % dedisperse, then channelize
    B = fft(double(readCF32(info_dc.chanFiles(j), m0, nSeg)));    % channelize, then dedisperse
    A = A(inF); B = B(inF);
    rho(j) = sum(B .* conj(A)) / sqrt(sum(abs(A).^2) * sum(abs(B).^2));
    pr(j)  = sum(abs(B).^2) / sum(abs(A).^2);
end
nEdge = ceil(info_ref.edgeWidth / dF);             % channels under the 8 MHz full-band taper
inner = (nEdge + 1 : info_chan.nChan - nEdge).';
fprintf(['2. commutation (channels %d..%d): |rho| min %.6f, |phase| max %.2e rad, ' ...
         'power ratio %.5f..%.5f\n'], inner(1), inner(end), min(abs(rho(inner))), ...
    max(abs(angle(rho(inner)))), min(pr(inner)), max(pr(inner)));
fprintf('   outer channels (full-band taper): |rho| %s, phase max %.2e rad\n', ...
    mat2str(abs(rho([1:nEdge, end-nEdge+1:end])).', 4), max(abs(angle(rho))));
pass = min(abs(rho(inner))) > 0.9999 && max(abs(angle(rho))) < 1e-3 && ...
       all(abs(pr(inner) - 1) < 2e-3);
fprintf('   commutation checks: %s\n', passStr(pass));
if ~pass, fails{end+1} = 'commutation'; end

% ---------------------------------------------------------------------------------
% 3. End to end: TOAs from the summed channel power vs the full-band chain
% ---------------------------------------------------------------------------------
template = gaussianTemplate(nBin, ephem.profileFWHM);
foldArgs = {'F1', ephem.F1, 'TRef', ephem.TRef, 'NBin', nBin, 'SubintPeriods', subintPeriods, ...
            'SaveFile', false, 'Verbose', false};

% Full-band reference (same data, same settings as main.m)
detR = detectPower(fileDedisp, fullfile(tmp, 'ref_env.dat'), info_ref.fs, info_ref.fLO, f_out, ...
    'FullySupported', info_ref.fullySupported, 'SaveInfo', false, 'Verbose', false);
[ifR, foldR] = foldProfile(detR, fullfile(tmp, 'ref_fold.mat'), ephem.f0, foldArgs{:});
BnR = noiseBandwidth(fLow, fHigh, info_ref.edgeWidth);
toaR = estimateTOA(foldR, ifR, template, 'Bnoise', BnR, 'Verbose', false);
detnR = detectPulsar(foldR, ifR, template, 'Bnoise', BnR, 'Verbose', false);

% Channelized: power per channel, summed (channel powers add up to the band power).
% Bin rate: the channel rate (4.1667 MHz) allows fs/4 = 1.0417 MHz, not 1 MHz.
fOutC = info_dc.fs / round(info_dc.fs / f_out);
pSum = [];
for j = 1:info_dc.nChan
    dj = detectPower(info_dc.chanFiles(j), fullfile(tmp, 'chan_env.dat'), info_dc.fs, ...
        info_dc.chanFreqs(j), fOutC, 'FullySupported', info_dc.fullySupported, ...
        'T0', info_dc.t0, 'SaveInfo', false, 'Verbose', false);
    p = readF32(dj.file);
    if isempty(pSum), pSum = zeros(size(p)); detC = dj; end
    pSum = pSum + double(p);
end
detC.file = fullfile(tmp, 'chan_sum_env.dat');
writeF32(detC.file, pSum);
[ifC, foldC] = foldProfile(detC, fullfile(tmp, 'chan_fold.mat'), ephem.f0, foldArgs{:});
toaC = estimateTOA(foldC, ifC, template, 'Bnoise', info_dc.BnoiseTotal, 'Verbose', false);
detnC = detectPulsar(foldC, ifC, template, 'Bnoise', info_dc.BnoiseTotal, 'Verbose', false);

[tr, iR, iC] = intersect(toaR.turnRef(toaR.valid), toaC.turnRef(toaC.valid));
tR = toaR.toa(toaR.valid); tC = toaC.toa(toaC.valid);
sR = toaR.toaErr(toaR.valid);
dT = tC(iC) - tR(iR);
snrRatio = median(toaC.snr(toaC.valid)) / median(toaR.snr(toaR.valid));
% Expected SNR ratio (simulation): signal sees W_f^2*W_ch^2 vs W_f^2*W_i^2, noise
% int W_ch^2 = int W_i^2; W_f = forward-dispersion taper (ground truth, test only)
info_disp = loadInfo(fileDispersed);
snrExp = taperedPower(fLow, fHigh, info_disp.edgeWidth, dF, info_dc.edgeWidth) / ...
         taperedPower(fLow, fHigh, info_disp.edgeWidth, fHigh - fLow, info_ref.edgeWidth);
nOff = 0.9 * nBin;                                   % ~off-pulse bins in the noise check
sigNR = sqrt(2 / nOff);                              % scatter of a variance estimate
fprintf(['3. end to end (%d common TOAs; bins %.4g us vs %.4g us): median SNR %.2f vs %.2f ' ...
         '(ratio %.4f, expected from tapers %.4f); median sigma %.3f vs %.3f us\n'], ...
    numel(tr), detC.binDt*1e6, detR.binDt*1e6, median(toaC.snr(toaC.valid)), ...
    median(toaR.snr(toaR.valid)), snrRatio, snrExp, median(toaC.toaErr(toaC.valid))*1e6, ...
    median(sR)*1e6);
fprintf(['   TOA difference channelized - full band: mean %+.3f us, rms %.3f us ' ...
         '(%.3f of sigma); total fold offset %+.3f vs %+.3f us\n'], mean(dT)*1e6, ...
    rms(dT)*1e6, rms(dT)/median(sR), toaC.total.timeOffset*1e6, toaR.total.timeOffset*1e6);
fprintf(['   fold noise ratio (off-pulse var / radiometer): %.4f vs %.4f (1-sigma ~%.3f each); ' ...
         'Bnoise %.4g vs %.4g MHz\n'], detnC.total.noiseRatio, detnR.total.noiseRatio, sigNR, ...
    info_dc.BnoiseTotal/1e6, BnR/1e6);
pass = numel(tr) >= 5 && abs(snrRatio/snrExp - 1) < 0.015 && rms(dT) < 0.2*median(sR) && ...
       abs(mean(dT)) < 0.3e-6 && abs(detnC.total.noiseRatio - 1) < 3*sigNR;
fprintf('   end-to-end checks: %s\n', passStr(pass));
if ~pass, fails{end+1} = 'end-to-end'; end

% ---------------------------------------------------------------------------------
% 4. Detected-power noise per time bin (summed channels, off-pulse)
% ---------------------------------------------------------------------------------
nPerBin = round(info_dc.fs / fOutC);
[vrExp, r1Exp] = binNoisePrediction(dF, info_dc.edgeWidth, info_dc.fs, nPerBin);
k  = (1:numel(pSum)).';
tk = detC.binTime0 + (k - 1) * detC.binDt;
ph = mod((tk - ephem.TRef) * ephem.f0 + 0.5, 1) - 0.5;          % pulse at phase 0
sb = detC.fullySupportedBins;
ok = abs(ph) > 0.15 & k >= sb(1) & k <= sb(2);
q  = pSum(ok); m = mean(q);
vr = var(q) / (m^2 / (info_dc.BnoiseTotal * detC.binDt));
x  = pSum - m; pair = ok(1:end-1) & ok(2:end);
r1 = mean(x([pair; false]) .* x([false; pair])) / var(q);
sigV = sqrt(2 / nnz(ok)); sigR1 = 1 / sqrt(nnz(pair));
fprintf(['4. time-bin noise (%d off-pulse bins of %d samples): var / radiometer %.4f ' ...
         '(pred %.4f +- %.4f), lag-1 corr %.4f (pred %.4f +- %.4f)\n'], nnz(ok), nPerBin, ...
    vr, vrExp, sigV*vrExp, r1, r1Exp, sigR1);
pass = abs(vr - vrExp) < 4*sigV*vrExp && abs(r1 - r1Exp) < 4*sigR1;
fprintf('   time-bin noise checks: %s\n', passStr(pass));
if ~pass, fails{end+1} = 'time-bin noise'; end

% ---------------------------------------------------------------------------------
if isempty(fails)
    fprintf('testDedisperseChannels: ALL PASSED\n');
else
    error('testDedisperseChannels:failed', 'FAILED: %s', strjoin(fails, ', '));
end


% =================================================================================
function s = passStr(p)
if p, s = 'PASS'; else, s = 'FAIL'; end
end

function P = taperedPower(fLow, fHigh, edgeF, chanW, edgeC)
%TAPEREDPOWER  int W_f^2 * W_c^2 df: W_f sin^2 edges (width edgeF) at the band edges,
% W_c sin^2 edges (width edgeC) at the edges of every chanW-wide channel.
f  = linspace(fLow, fHigh, 4000001);
W  = sin2Taper(f, fLow, fHigh, edgeF);
fc = fLow + floor((f - fLow) / chanW) * chanW;              % channel start of each f
fc(f >= fHigh) = fHigh - chanW;
Wc = sin2Taper(f - fc, 0, chanW, edgeC);
P  = trapz(f, W.^2 .* Wc.^2);
end

function W = sin2Taper(f, a, b, e)
W  = ones(size(f));
lo = f < a + e;  W(lo) = sin(pi/2 * (f(lo) - a) / e).^2;
hi = f > b - e;  W(hi) = sin(pi/2 * (b - f(hi)) / e).^2;
end

function [vr, r1] = binNoisePrediction(dF, edgeW, fs, n)
%BINNOISEPREDICTION  Variance (relative to the radiometer value) and lag-1 correlation
% of time-bin averages of |y|^2 over n samples, for Gaussian y with spectrum
% S = W^2 (sin^2 edges of width edgeW inside +-dF/2), sampled at fs.
% R(l) = normalized autocorrelation of y; cov of bin powers at bin lag L
% = (1/n^2) sum_{a,b} |R(a - b + L*n)|^2; radiometer value = 1/(B*dt) with
% B = (int S)^2 / int S^2 = fs / sum_l |R(l)|^2 and dt = n/fs.
Nf = 2^18;
f  = ((0:Nf-1).' - Nf/2) * fs / Nf;
S  = zeros(Nf, 1); in = abs(f) <= dF/2;
S(in) = sin2Taper(f(in), -dF/2, dF/2, edgeW).^2;
R  = fft(ifftshift(S));                       % R(l) up to a constant, l = 0..Nf-1 (circular)
R  = conj(R) / R(1);                          % sign convention irrelevant for |R|^2
R2 = abs(R).^2;
lagR2 = @(l) R2(mod(l, Nf) + 1);
C  = @(L) sum(arrayfun(@(a) sum(lagR2(a - (0:n-1) + L*n)), 0:n-1)) / n^2;
radiometer = sum(R2) / n;                     % = fs/(B*n) in units of m^2
vr = C(0) / radiometer;
r1 = C(1) / C(0);
end

function eq = filesEqual(f1, f2)
%FILESEQUAL  Byte-wise comparison in chunks.
a = dir(f1); b = dir(f2);
eq = ~isempty(a) && ~isempty(b) && a.bytes == b.bytes;
if ~eq, return; end
fa = fopen(f1, 'r'); fb = fopen(f2, 'r');
c = onCleanup(@() cellfun(@fclose, {fa, fb}));
while eq
    xa = fread(fa, 2^24, 'uint8=>uint8');
    xb = fread(fb, 2^24, 'uint8=>uint8');
    eq = isequal(xa, xb);
    if isempty(xa), break; end
end
end

function y = readCF32(file, first, count)
fid = fopen(file, 'r', 'ieee-le');
c = onCleanup(@() fclose(fid));
fseek(fid, first * 8, 'bof');
raw = fread(fid, [2 count], 'single=>single');
y = complex(raw(1, :), raw(2, :)).';
end

function p = readF32(file)
fid = fopen(file, 'r', 'ieee-le');
c = onCleanup(@() fclose(fid));
p = fread(fid, Inf, 'single=>single');
end

function writeF32(file, p)
fid = fopen(file, 'w', 'ieee-le');
fwrite(fid, single(p), 'single');
fclose(fid);
end
