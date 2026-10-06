function testDetectChannels()
%TESTDETECTCHANNELS  Unit tests for detectChannels and powerCovariance (unit 3b).
%{
Run from the PulsarSimMatlab folder: run('tests/testDetectChannels.m'). ~30 s.
Uses the dedispersed channels of tests/testDedisperseChannels.m
(data/chan/test_dedisp_chan_*); makes them first if they are missing.

  1. powerCovariance (from the channel spectrum) against the same statistics
     computed from the actual impulse response h of a channel's dedispersion
     filter: C(L) = (1/n^2) sum_{a,b} |c(a - b + L*n)|^2, c(tau) = sum_l h(l)
     conj(h(l+tau)) (the chirp drops out). V, X(1..Lmax) must agree.
  2. Layout: row j of the [nChan x nBins] file is bit-identical to detectPower
     run on channel j alone (channels 1, 64, 128).
  3. foldProfile reads the 128-channel file; the sum over channels of its
     profiles equals the fold of the summed channel power (float rounding).
  4. Measured noise of the detected bins (off-pulse, pooled over 128 channels):
     variance and lag-1..3 covariances vs V and X(L); variance of sums of 5
     consecutive bins (~ one 4.9 us phase bin) vs V + 2*sum_L (1 - L/5) X(L).
Errors at the end if any check fails.
%}

root = fullfile(fileparts(mfilename('fullpath')), '..');
addpath(fullfile(root, 'functions'));
oldDir = cd(root); restoreDir = onCleanup(@() cd(oldDir));
pipelineParams;
chanDir = fullfile(dataDir, 'chan');
tmp = fullfile(tempdir, 'testDetectChannels');
if ~isfolder(tmp), mkdir(tmp); end
fails = {};

baseDC = char(fullfile(chanDir, 'test_dedisp_chan'));     % char: [baseDC '_info.mat'] must join text
if isfile([baseDC '_info.mat'])
    info_dc = loadInfo(baseDC);
else
    info_IQ = loadInfo(fileIQ);
    info_chan = channelizeIQ(info_IQ.file, fullfile(chanDir, 'test_rx_IQ_chan'), ...
        info_IQ.actualFsOut, info_IQ.fLO, fLow, fHigh, 'T0', info_IQ.t0);
    info_dc = dedisperseChannels(info_chan, baseDC, ephem.DM, 'RefFreq', refFreq);
end
fOutC = info_dc.fs / round(info_dc.fs / f_out);
n = round(info_dc.fs / fOutC);

% ---------------------------------------------------------------------------------
% 1. powerCovariance vs the real dedispersion filter
% ---------------------------------------------------------------------------------
nc = powerCovariance(info_dc.chanWidth, info_dc.edgeWidth, info_dc.fs, n);
fc = info_dc.chanFreqs(1); dF = info_dc.chanWidth;
Nimp = 2^16; n0 = Nimp/2;
imp = zeros(Nimp, 1); imp(n0 + 1) = 1;
writeCF32(fullfile(tmp, 'imp.dat'), imp);
d = applyInverseDispersion(fullfile(tmp, 'imp.dat'), fullfile(tmp, 'imp_out.dat'), info_dc.fs, ...
    fc, ephem.DM, fc - dF/2, fc + dF/2, 'RefFreq', refFreq, 'AllowRefOutsideBand', true, ...
    'EdgeFrac', info_dc.edgeWidth/dF, 'SaveInfo', false, 'Verbose', false);
yi = double(readCF32(fullfile(tmp, 'imp_out.dat')));
hv = yi(n0 + 1 + (-d.nFuture:d.nPast));
c  = conv(hv, conj(flipud(hv)));                   % c(tau) = sum_l h(l+tau) conj(h(l)), centred
c  = c / c(numel(hv));                              % normalize: c(0) = 1 at index numel(hv)
c2 = abs(c).^2; mid = numel(hv);
radH = sum(c2) / n;
dd = (-(n-1):(n-1)).'; cnt = n - abs(dd);
CL = @(L) sum(cnt .* c2(mid + dd + L*n)) / n^2;
VH = CL(0) / radH;
XH = arrayfun(@(L) CL(L) / radH, 1:nc.Lmax);
dV = abs(VH - nc.V); dX = max(abs(XH - nc.X));
fprintf(['1. noise statistics, spectrum vs filter (channel 1, %d samples/bin): V %.5f vs %.5f, ' ...
         'X(1..3) %s vs %s, max |dX| %.1e; Lmax %d, captured %.4f\n'], n, nc.V, VH, ...
    mat2str(nc.X(1:3), 4), mat2str(XH(1:3), 4), dX, nc.Lmax, nc.captured);
pass = dV < 1e-3 && dX < 1e-3 && nc.captured >= 0.998;
fprintf('   %s\n', passStr(pass));
if ~pass, fails{end+1} = 'powerCovariance vs filter'; end

% ---------------------------------------------------------------------------------
% 2. Layout of the interleaved power file
% ---------------------------------------------------------------------------------
fileP = fullfile(chanDir, 'test_dedisp_chan_power.dat');
info_det = detectChannels(info_dc, fileP, fOutC);
Pall = readF32(fileP);
Pall = reshape(Pall, info_det.nChan, []);
same = true;
for j = [1 64 128]
    dj = detectPower(info_dc.chanFiles(j), fullfile(tmp, 'one.dat'), info_dc.fs, ...
        info_dc.chanFreqs(j), fOutC, 'FullySupported', info_dc.fullySupported, ...
        'T0', info_dc.t0, 'SaveInfo', false, 'Verbose', false);
    same = same && isequal(readF32(dj.file).', Pall(j, :)) && dj.N == info_det.N && ...
        isequal(dj.fullySupportedBins, info_det.fullySupportedBins);
end
noTmp = isempty(dir(fullfile(chanDir, 'test_dedisp_chan_power_tmpch*.dat')));
fprintf('2. layout: rows bit-identical to detectPower per channel %d; temporary files removed %d: %s\n', ...
    same, noTmp, passStr(same && noTmp));
if ~(same && noTmp), fails{end+1} = 'layout'; end

% ---------------------------------------------------------------------------------
% 3. foldProfile on the 128-channel file vs the fold of the summed power
% ---------------------------------------------------------------------------------
foldArgs = {'F1', ephem.F1, 'TRef', ephem.TRef, 'NBin', nBin, 'SubintPeriods', subintPeriods, ...
            'SaveFile', false, 'Verbose', false};
[ifA, foldA] = foldProfile(info_det, fullfile(tmp, 'fA.mat'), ephem.f0, foldArgs{:});
detS = info_det; detS.nChan = 1; detS.file = fullfile(tmp, 'sum.dat');
writeF32(detS.file, sum(double(Pall), 1));
[~, foldS] = foldProfile(detS, fullfile(tmp, 'fS.mat'), ephem.f0, foldArgs{:});
dFold = max(abs(sum(foldA.profTotal, 2) - foldS.profTotal)) / mean(foldS.profTotal);
fprintf('3. fold of 128 channels (%d x %d x %d): channel sum vs fold of summed power, max rel diff %.1e: %s\n', ...
    ifA.NBin, ifA.nSub, ifA.nChan, dFold, passStr(dFold < 1e-5));
if dFold >= 1e-5, fails{end+1} = 'fold'; end

% ---------------------------------------------------------------------------------
% 4. Measured noise of detected bins vs V and X(L)
% ---------------------------------------------------------------------------------
k  = (1:info_det.N).';
tk = info_det.binTime0 + (k - 1) * info_det.binDt;
ph = mod((tk - ephem.TRef) * ephem.f0 + 0.5, 1) - 0.5;
sb = info_det.fullySupportedBins;
off = abs(ph) > 0.15 & k >= sb(1) & k <= sb(2);
rad = nc.radiometer;
P   = double(Pall).';                                % nBins x nChan
mu  = mean(P(off, :), 1);
r   = P ./ mu - 1;                                   % relative fluctuation per channel
r(~off, :) = NaN;
nOff = nnz(off) * info_det.nChan;
Vm = mean(r(off, :).^2, 'all') / rad;
Xm = zeros(1, 3);
for L = 1:3
    pr = r(1:end-L, :) .* r(1+L:end, :);
    Xm(L) = mean(pr(~isnan(pr))) / rad;
end
K = 5;                                               % sums of 5 bins (~ one phase bin)
nBlk = floor(info_det.N / K);
okB = all(reshape(off(1:nBlk*K), K, nBlk), 1).';
Bs  = squeeze(sum(reshape(P(1:nBlk*K, :) ./ mu - 1, K, nBlk, []), 1));   % nBlk x nChan
VK  = mean(Bs(okB, :).^2, 'all') / (K * rad);
VKp = nc.V + 2*sum((1 - (1:min(K-1, nc.Lmax))/K) .* nc.X(1:min(K-1, nc.Lmax)));
sV = 1.4*sqrt(2/nOff); sX = 1.2/sqrt(nOff); sK = 1.4*sqrt(2/(nnz(okB)*info_det.nChan));
fprintf(['4. measured (%d off-pulse bins x %d channels): V %.4f (pred %.4f +- %.4f), ' ...
         'X(1..3) %s (pred %s +- %.4f)\n'], nnz(off), info_det.nChan, Vm, nc.V, sV*nc.V, ...
    mat2str(Xm, 4), mat2str(nc.X(1:3), 4), sX*nc.V);
fprintf(['   variance of %d-bin sums / (%d x radiometer): %.4f (pred %.4f +- %.4f; ' ...
         'independent-bin model 1)\n'], K, K, VK, VKp, sK*VKp);
pass = abs(Vm - nc.V) < 4*sV*nc.V && all(abs(Xm - nc.X(1:3)) < 4*sX*nc.V) && ...
       abs(VK - VKp) < 4*sK*VKp;
fprintf('   %s\n', passStr(pass));
if ~pass, fails{end+1} = 'measured noise'; end

% ---------------------------------------------------------------------------------
if isempty(fails)
    fprintf('testDetectChannels: ALL PASSED\n');
else
    error('testDetectChannels:failed', 'FAILED: %s', strjoin(fails, ', '));
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

function y = readCF32(file)
fid = fopen(file, 'r', 'ieee-le');
c = onCleanup(@() fclose(fid));
raw = fread(fid, [2 Inf], 'single=>single');
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
