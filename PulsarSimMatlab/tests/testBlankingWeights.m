function testBlankingWeights()
%TESTBLANKINGWEIGHTS  Unit tests for blankingWeights (A3a): exact data weights from a mask.
%{
Run from the PulsarSimMatlab folder: run('tests/testBlankingWeights.m').
Uses the seed-43 channel files of tests/testDetectChannels.m (data/chan/
test_rx_IQ_chan_*, test_dedisp_chan_*, test_dedisp_chan_power.dat). Writes
~1.5 GB of temporary files to tempdir (deleted on the way and at the end).

  1. Exactness: W, V, X(1..7) of the weight file vs a direct calculation per
     bin with the full sample covariance matrix C = H*K*H' (H from the
     channel's own dedispersion filter, an impulse through
     applyInverseDispersion; K = diag(keep)) for ~150 bins each in
     channels 1, 64, 128 (bottom, middle, top: different delays to 1.6 GHz)
     with blanks of 1, 9, 50 (overlapping), 300 and 3001 samples; W with the
     channel autocorrelation as u' R u (independent matrix form).
  2. Consistency: (a) unblanked V0, X0 from each channel's filter = the
     spectrum constants of detectChannels (powerCovariance); (b) empty mask:
     the fold with these data weights = the fold with NoiseCoeffs (A1a) in
     every channel; (c) far from blanks a blanked channel has its V0, X0;
     (d) a messy mask (overlapping, unsorted, outside the file) gives the
     same file as the equivalent clean mask; bins without data have W = 0;
     (e) the channel noise autocorrelation (prototype spectrum aliased at the
     channel rate, used for the mean) = measured on the raw channel samples.
  3. End to end on the seed-43 data: blank copies of the channel IQ files
     BEFORE dedispersion, then dedisperseChannels, detectChannels,
     blankingWeights, foldProfile('DataWeights'), estimateTOA, detectPulsar.
     Mask: random 20 us blanks (5 %, every 4th channel), broadband 1 us
     impulses (all channels), radar 2 us at 373 Hz (channels 31-34), radar
     2 us at 400 Hz (locked to the 100 Hz pulsar; channels 50-53), 30 us
     blanks on the pulse itself in every turn (channels 90-100), 100 us
     blanks every 1.1 ms (channels 110-115), a 30 ms gap (channel 64),
     channel 10 blanked completely.
     (a) means: detected power E[P_b] = W * E[P_u] per channel and per W
         group down to W = 1e-6 (white-input W shown: ~4 % off there); folded
         profile unbiased (inverse-variance weighted mean) off-pulse (all
         channels, and per blanking type and per valid-fraction group) and
         ON the pulse in the pulse-blanked channels (the naive fold,
         blanking ignored, shown);
     (b) noise: measured variance and lag-1/2 covariances of off-pulse
         profile bins vs weight2 / weightX, in groups of the valid fraction
         of the phase bin (1 .. heavily blanked), within 4 sigma, as a
         fraction of the variance; the simple w^2 and w models shown;
     (c) channel use: channel 10 never, channel 64 left out where its gap
         lies, all as the exclusion rule says;
     (d) TOAs vs unblanked consistent; noise ratio of detectPulsar = 1, and
         equal to that of the unblanked fold with the same channels per
         sub-int (paired, much tighter).
Errors at the end if any check fails.
%}

root = fullfile(fileparts(mfilename('fullpath')), '..');
addpath(fullfile(root, 'functions'));
oldDir = cd(root); restoreDir = onCleanup(@() cd(oldDir));
pipelineParams;
tmp = fullfile(tempdir, 'testBlankingWeights');
if ~isfolder(tmp), mkdir(tmp); end
cleanTmp = onCleanup(@() cleanFolder(tmp));
fails = {};

chanDir   = fullfile(dataDir, 'chan');
info_chan = loadInfo(fullfile(chanDir, 'test_rx_IQ_chan'));
info_dc   = loadInfo(fullfile(chanDir, 'test_dedisp_chan'));
detD      = loadInfo(fullfile(chanDir, 'test_dedisp_chan_power.dat'));
nc  = detD.noise; nCh = info_dc.nChan; Nc = info_dc.N; n = detD.binLen; nb = detD.N;
Lm  = nc.Lmax; rad = nc.radiometer; sb = detD.fullySupportedBins; fs = info_dc.fs;

% ---------------------------------------------------------------------------------
% 1. Exactness vs the direct covariance matrix
% ---------------------------------------------------------------------------------
chans = [1 64 128];
mask1 = zeros(0, 3);
for c = chans
    s0 = 150000 + 1000*c;
    mask1 = [mask1; c s0 s0; c s0+500 s0+508; c s0+2000 s0+2049; c s0+2030 s0+2079; ...
             c s0+5000 s0+5299; c s0+20000 s0+23000]; %#ok<AGROW>
end
bw1 = blankingWeights(info_chan, info_dc, detD, mask1, fullfile(tmp, 'w1.dat'), 'Verbose', false);
Y1 = readWeights(bw1);
errW = 0; errV = 0; errX = 0; nChk = 0; zeroOk = true;
for c = chans
    [h, nF, nP] = channelKernel(info_dc, c, tmp);
    m0 = sum(abs(h).^2);
    keep = keepOf(mask1, c, Nc);
    W = squeeze(Y1(c, 1, :));
    kk = (sb(1):sb(2) - Lm).';
    aff = kk(W(kk) < 0.99999);
    pick = aff(unique(round(linspace(1, numel(aff), min(130, numel(aff))))));
    edgeZ = find(diff(W == 0) ~= 0);                 % bins at the border of W = 0
    edgeZ = edgeZ(edgeZ >= sb(1) & edgeZ <= sb(2) - Lm - 1);
    pick = unique([pick; edgeZ; edgeZ + 1; kk(round(linspace(1, numel(kk), 10)))]);
    for k = pick.'
        [We, Ve, Xe] = directStats(h, nF, nP, keep, k, n, Lm, m0, rad, bw1.Rx);
        Yk = squeeze(Y1(c, :, k)).';
        if Yk(1) == 0                                 % set to 0 by MinWeight
            zeroOk = zeroOk && We < bw1.minWeight && all(Yk == 0);
            continue
        end
        errW = max(errW, abs(We - Yk(1)));
        errV = max(errV, abs(Ve - Yk(2)));
        errX = max(errX, max(abs(Xe(:) - Yk(3:end))));
        nChk = nChk + 1;
    end
end
pass = errW < 1e-6 && errV < 1e-6 && errX < 1e-6 && zeroOk;
fprintf(['1. exact vs direct covariance matrix (%d bins, channels %s): max |dW| %.1e, |dV| %.1e, ' ...
         '|dX| %.1e (float32 file); W = 0 bins consistent %d: %s\n'], nChk, mat2str(chans), errW, ...
    errV, errX, zeroOk, passStr(pass));
if ~pass, fails{end+1} = 'exact vs direct'; end

% ---------------------------------------------------------------------------------
% 2. Consistency
% ---------------------------------------------------------------------------------
dV0 = max(abs(bw1.V0 - nc.V)); dX0 = max(abs(bw1.X0 - nc.X), [], 'all');
fprintf('2a. unblanked V0, X0 from each filter vs spectrum constants: max |dV| %.1e, |dX| %.1e: %s\n', ...
    dV0, dX0, passStr(dV0 < 1e-4 && dX0 < 1e-4));
if ~(dV0 < 1e-4 && dX0 < 1e-4), fails{end+1} = 'V0/X0 vs spectrum'; end

base = {'F1', ephem.F1, 'TRef', ephem.TRef, 'NBin', nBin, 'SubintPeriods', subintPeriods, ...
        'SaveFile', false, 'Verbose', false};
bw0 = blankingWeights(info_chan, info_dc, detD, [], fullfile(tmp, 'w0.dat'), 'Verbose', false);
[~, f0w] = foldProfile(detD, 'x', ephem.f0, base{:}, 'DataWeights', bw0, 'ChunkBins', 20000);
[~, f0a] = foldProfile(detD, 'x', ephem.f0, base{:}, 'ChunkBins', 20000, 'NoiseCoeffs', ...
    double(single([bw0.V0(1), bw0.X0(1, :)])));     % same chunks: same rounding
r2 = max(abs(f0w.weight2 ./ f0a.weight2 - 1), [], 'all');
rX = max(abs(f0w.weightX(:, :, :, 1) ./ f0a.weightX(:, :, :, 1) - 1), [], 'all');
sameW = isequal(f0w.weight, repmat(f0a.weight, 1, 1, nCh)) && isequaln(f0w.prof, f0a.prof);
pass = sameW && r2 < 1e-6 && rX < 1e-6;
fprintf(['2b. empty mask: fold with data weights = NoiseCoeffs fold: weight, prof bit-identical %d; ' ...
         'weight2 / weightX rel diff %.1e / %.1e: %s\n'], sameW, r2, rX, passStr(pass));
if ~pass, fails{end+1} = 'empty mask'; end

far = (sb(1):sb(1) + 2000).';                        % blanks start at sample 150064
dFar = 0;
for c = chans
    dFar = max([dFar, max(abs(squeeze(Y1(c, 2, far)) - bw1.V0(c))), ...
        max(abs(squeeze(Y1(c, 3:end, far)).' - bw1.X0(c, :)), [], 'all'), ...
        max(abs(squeeze(Y1(c, 1, far)) - 1))]);
end
fprintf('2c. far from blanks: W = 1, V = V0, X = X0 to %.1e: %s\n', dFar, passStr(dFar < 1e-6));
if dFar >= 1e-6, fails{end+1} = 'far bins'; end

% 2e. the channel noise spectrum (prototype aliased at the channel rate) vs data
rho = bw1.Rx / bw1.Rx(1); nR = numel(rho) - 1;
chk = round(linspace(8, nCh - 8, 8));
rhoM = zeros(numel(chk), nR);
for i = 1:numel(chk)
    fid = fopen(info_chan.chanFiles(chk(i)), 'r', 'ieee-le');
    raw = fread(fid, [2 Inf], 'single=>double'); fclose(fid);
    xc = complex(raw(1, :), raw(2, :)).';
    xc = xc(info_chan.fullySupported(1):info_chan.fullySupported(2));
    P0 = mean(abs(xc).^2);
    for d = 1:nR
        rhoM(i, d) = real(mean(xc(1+d:end) .* conj(xc(1:end-d)))) / P0;
    end
end
dRho = mean(rhoM, 1) - rho(2:end).';
sRho = std(rhoM, 0, 1) / sqrt(numel(chk));
okRho = all(abs(dRho) < 4 * sRho);
fprintf(['2e. channel autocorrelation, %d lags, measured (%d channels) vs prototype spectrum: max ' ...
         '|diff| %.1e (%.1f sigma max); white input would be off by up to %.3f: %s\n'], nR, ...
    numel(chk), max(abs(dRho)), max(abs(dRho) ./ sRho), max(abs(rho(2:end))), passStr(okRho));
if ~okRho, fails{end+1} = 'channel spectrum'; end

messy = [5 3000 3100; 5 3050 3200; 5 -40 10; 5 Nc-5 Nc+100; 5 9000 9000; 5 2000 1990];
clean = [5 1 10; 5 3000 3200; 5 9000 9000; 5 Nc-5 Nc];
bwM = blankingWeights(info_chan, info_dc, detD, messy(randperm(size(messy, 1)), :), fullfile(tmp, 'wm.dat'), ...
    'Verbose', false);
bwC = blankingWeights(info_chan, info_dc, detD, clean, fullfile(tmp, 'wc.dat'), 'Verbose', false);
sameMask = isequal(readWeights(bwM), readWeights(bwC));
fprintf('2d. messy mask (overlaps, unsorted, outside, empty row) = clean mask, identical file %d: %s\n', ...
    sameMask, passStr(sameMask));
if ~sameMask, fails{end+1} = 'mask normalization'; end
deleteFiles({bw0.file, bwM.file, bwC.file, bw1.file});
clear Y1

% ---------------------------------------------------------------------------------
% 3. End to end: blank the channel voltages before dedispersion
% ---------------------------------------------------------------------------------
t3 = tic;
mask3 = makeMask3(info_chan, ephem, refFreq, nCh, Nc, fs);
info_chanB = info_chan;
info_chanB.chanFiles = blankChannelFiles(info_chan, mask3, fullfile(tmp, 'chanB'));
info_dcB = dedisperseChannels(info_chanB, fullfile(tmp, 'dcB'), ephem.DM, 'RefFreq', refFreq, ...
    'Verbose', false);
deleteFiles(cellstr(info_chanB.chanFiles));
fOutC = info_dcB.fs / round(info_dcB.fs / f_out);
detB = detectChannels(info_dcB, fullfile(tmp, 'pB.dat'), fOutC, 'Verbose', false);
deleteFiles(cellstr(info_dcB.chanFiles));
gridOk = detB.N == nb && detB.binTime0 == detD.binTime0 && isequal(detB.fullySupportedBins, sb) && ...
    isequal(info_dcB.Nfft, info_dc.Nfft) && isequal(info_dcB.nFuture, info_dc.nFuture);
bw3 = blankingWeights(info_chan, info_dcB, detB, mask3, fullfile(tmp, 'w3.dat'), 'Verbose', false);
fprintf(['3. blanked %d mask rows: %.2f %% of channel samples, valid fraction %.4f (fully supported); ' ...
         'same grid as unblanked %d; chain %.0f s\n'], size(mask3, 1), 100*mean(bw3.blankedFraction), ...
    mean(bw3.validFraction), gridOk, toc(t3));
if ~gridOk, fails{end+1} = 'grid'; end

% (a0) detected power: E[P_b] = W * E[P_u], channel by channel
Pb = reshape(readF32(detB.file), nCh, nb);
Pu = reshape(readF32(detD.file), nCh, nb);
Y3 = readWeights(bw3, 'single');
W3 = reshape(Y3(:, 1, :), nCh, nb);
ks = sb(1):sb(2);
num = sum(double(Pb(:, ks)) - double(W3(:, ks)) .* double(Pu(:, ks)), 2);
den = sum(double(W3(:, ks)) .* double(Pu(:, ks)), 2);
hasW = den > 0;
dP = num(hasW) ./ den(hasW);
okP = abs(mean(dP)) < 4 * std(dP) / sqrt(numel(dP));
fprintf(['   (a0) detected power, sum(P_b - W P_u) / sum(W P_u) per channel: mean %+.2e +- %.2e ' ...
         '(max |.| %.1e, %d channels): %s\n'], mean(dP), std(dP)/sqrt(numel(dP)), max(abs(dP)), ...
    numel(dP), passStr(okP));
if ~okP, fails{end+1} = 'detected power mean'; end
% the same per W group (nearly empty bins: where the channel spectrum matters), vs white input
bw3w = blankingWeights(info_chan, info_dcB, detB, mask3, fullfile(tmp, 'w3w.dat'), ...
    'InputSpectrum', 'white', 'Verbose', false);
Y3w = readWeights(bw3w, 'single'); W3w = reshape(Y3w(:, 1, :), nCh, nb); clear Y3w
deleteFiles({bw3w.file});
wEdges = [1e-6 1e-4 1e-3 1e-2 0.1 0.5 0.9 0.99999];
Pbd = double(Pb); Pud = double(Pu);
for g = 1:numel(wEdges) - 1
    sel = false(nCh, nb); sel(:, ks) = W3(:, ks) >= wEdges(g) & W3(:, ks) < wEdges(g + 1);
    Dc  = sum((Pbd - double(W3) .* Pud) .* sel, 2);       % per channel
    Dw  = sum((Pbd - double(W3w) .* Pud) .* sel, 2);
    Sc  = sum(double(W3) .* Pud .* sel, 2);
    nPc = sum(sel, 2); use = nPc >= 3;
    if nnz(use) < 5, continue; end
    tot = sum(Sc(use)); sg = std(Dc(use)) * sqrt(nnz(use)) / tot;
    dev = sum(Dc(use)) / tot; devW = sum(Dw(use)) / tot;
    okg = abs(dev) < 4 * sg;
    fprintf(['        W %-7.1g..%-7.1g %7d bins, %3d chan: sum(P_b)/sum(W P_u) - 1 = %+.4f +- %.4f ' ...
             '(%+.1f sigma): %s;  white-input W: %+.4f (%+.1f sigma)\n'], wEdges(g), wEdges(g+1), ...
        sum(nPc(use)), nnz(use), dev, sg, dev/sg, passStr(okg), devW, devW/sg);
    if ~okg, fails{end+1} = sprintf('detected mean W group %.1g', wEdges(g)); end %#ok<AGROW>
end
clear Pbd Pud W3w
clear Pb Pu

% folds: unblanked (A1a), blanked with exact weights, blanked naive
[ifU, fU] = foldProfile(detD, 'x', ephem.f0, base{:}, 'NoiseCoeffs', [nc.V, nc.X]);
[ifB, fB] = foldProfile(detB, 'x', ephem.f0, base{:}, 'DataWeights', bw3);
[~, fN] = foldProfile(detB, 'x', ephem.f0, base{:}, 'NoiseCoeffs', [nc.V, nc.X]);
ph   = fU.phase.';
offB = ph >= 0.15 & ph <= 0.85;
onB  = ph <= 0.02 | ph >= 0.98;
okU  = fU.weight > 0;                                % [NBin x nSub]
muC  = zeros(1, 1, nCh);                             % true level per channel (unblanked)
for c = 1:nCh
    pc = fU.prof(:, :, c);
    muC(c) = mean(pc(offB & okU));
end

% (a) folded means, inverse-variance weighted (bins with tiny W are very noisy)
wt = fB.weight.^2 ./ fB.weight2;                     % 1/var of a profile bin (relative)
selO = offB & fB.weight > 0 & okU;
wO = wt; wO(~selO) = 0;
dd = (fB.prof - fU.prof) ./ muC; dd(~selO) = 0;
dOff = squeeze(sum(dd .* wO, [1 2]) ./ sum(wO, [1 2]));
dOff = dOff(~isnan(dOff));
okOff = abs(mean(dOff)) < 4 * std(dOff) / sqrt(numel(dOff));
plc = 90:100;                                        % pulse-blanked channels
amp = zeros(numel(plc), 1); dOn = amp; dOnN = amp;
for i = 1:numel(plc)
    c = plc(i);
    amp(i) = max(fU.profTotal(:, c)) - muC(c);
    sel = onB & fB.weight(:, :, c) > 0 & okU;
    wc = wt(:, :, c);
    dB = fB.prof(:, :, c) - fU.prof(:, :, c);
    dN = fN.prof(:, :, c) - fU.prof(:, :, c);
    dOn(i)  = sum(dB(sel) .* wc(sel)) / sum(wc(sel)) / amp(i);
    dOnN(i) = mean(dN(sel)) / amp(i);
end
okOn = abs(mean(dOn)) < 4 * std(dOn) / sqrt(numel(plc));
selOn = onB & okU;
wFrac = fB.weight(:, :, plc) ./ fU.weight;
fprintf(['   (a) folded mean vs unblanked: off-pulse (all channels) %+.2e +- %.2e: %s; ON the pulse ' ...
         'in the pulse-blanked channels 90-100 (valid fraction there %.2f): %+.4f +- %.4f of the ' ...
         'pulse height (naive fold %+.3f): %s\n'], mean(dOff), std(dOff)/sqrt(numel(dOff)), ...
    passStr(okOff), mean(wFrac(repmat(selOn, 1, 1, numel(plc)))), mean(dOn), ...
    std(dOn)/sqrt(numel(plc)), mean(dOnN), passStr(okOn));
if ~okOff, fails{end+1} = 'folded mean off-pulse'; end
if ~okOn, fails{end+1} = 'folded mean on the pulse'; end
% per blanking type (an effect in a few channels must not hide in the average)
dAll = squeeze(sum(dd .* wO, [1 2]) ./ sum(wO, [1 2]));     % per channel, NaN for channel 10
tName = {'random 20 us', 'radar 373 Hz', 'radar 400 Hz locked', 'on-pulse blanks (off-pulse)', ...
         '100 us / 1.1 ms', 'impulses only'};
tChan = {1:4:nCh, 31:34, 50:53, 90:100, 110:115, []};
tChan{end} = setdiff(1:nCh, [tChan{1:end-1}, 10, 64]);
for g = 1:numel(tName)
    x = dAll(tChan{g}); x = x(~isnan(x));
    sg = std(x) / sqrt(numel(x)); okg = abs(mean(x)) < 4*sg;
    fprintf('       mean, %-27s %3d chan: %+.2e +- %.2e (%+.1f sigma): %s\n', tName{g}, numel(x), ...
        mean(x), sg, mean(x)/sg, passStr(okg));
    if ~okg, fails{end+1} = ['folded mean ' tName{g}]; end %#ok<AGROW>
end
% per valid fraction of the phase bin (all channels pooled per channel)
fr = fB.weight ./ fU.weight;
fEdges = [0 0.01 0.1 0.5 0.9 1.01];
for g = 1:numel(fEdges) - 1
    wg = wO .* (fr >= fEdges(g) & fr < fEdges(g + 1));
    x = squeeze(sum(dd .* wg, [1 2]) ./ sum(wg, [1 2]));
    nbg = squeeze(sum(wg > 0, [1 2]));
    x = x(nbg >= 5);
    if numel(x) < 3, continue; end
    sg = std(x) / sqrt(numel(x)); okg = abs(mean(x)) < 4*sg;
    fprintf('       mean, valid fraction %-5.2g..%-5.2g %3d chan, %6d bins: %+.2e +- %.2e (%+.1f sigma): %s\n', ...
        fEdges(g), min(fEdges(g+1), 1), numel(x), sum(nbg(nbg >= 5)), mean(x), sg, mean(x)/sg, passStr(okg));
    if ~okg, fails{end+1} = sprintf('folded mean fraction %.2g', fEdges(g)); end %#ok<AGROW>
end

% (b) noise in groups of the valid fraction of the phase bin, with the w^2 / w models
Y3a = Y3; Y3b = Y3;
for c = 1:nCh
    Wc = double(reshape(Y3(c, 1, :), [], 1));
    Y3a(c, 2, :) = bw3.V0(c) * Wc.^2;                % "w^2 model"
    Y3b(c, 2, :) = bw3.V0(c) * Wc;                   % "w model"
    for L = 1:Lm
        WW = Wc .* [Wc(1+L:end); zeros(L, 1)];
        Y3a(c, 2 + L, :) = bw3.X0(c, L) * WW;
        Y3b(c, 2 + L, :) = bw3.X0(c, L) * sqrt(WW);
    end
end
bwA = bw3; bwA.file = fullfile(tmp, 'w3a.dat'); writeArray(bwA.file, Y3a); clear Y3a
[~, fA] = foldProfile(detB, 'x', ephem.f0, base{:}, 'DataWeights', bwA); deleteFiles({bwA.file});
bwW = bw3; bwW.file = fullfile(tmp, 'w3b.dat'); writeArray(bwW.file, Y3b); clear Y3b Y3
[~, fW] = foldProfile(detB, 'x', ephem.f0, base{:}, 'DataWeights', bwW); deleteFiles({bwW.file});

u = (fB.prof ./ muC - 1) .* fB.weight / sqrt(rad);
u(~(offB & fB.weight > 0)) = NaN;
frac = fB.weight ./ fU.weight;                       % valid fraction of the phase bin
gEdges = [0 0.5 0.9 0.999 1.01]; gName = {'0-0.5', '0.5-0.9', '0.9-0.999', '>=0.999'};
for d = 0:2
    pr = u .* circshift(u, -d, 1);
    if d == 0
        qE = fB.weight2; qA = fA.weight2; qW = fW.weight2; lb = 'variance';
    else
        qE = fB.weightX(:, :, :, d); qA = fA.weightX(:, :, :, d); qW = fW.weightX(:, :, :, d);
        lb = sprintf('lag %d', d);
    end
    for g = 1:numel(gName)
        v = ~isnan(pr) & frac >= gEdges(g) & frac < gEdges(g + 1);
        nPerC = squeeze(sum(v, [1 2]));
        use = nPerC >= 20;
        if nnz(use) < 3, continue; end
        Sm = squeeze(sum(pr .* v, [1 2], 'omitnan'));
        SE = squeeze(sum(qE .* v, [1 2]));
        SA = squeeze(sum(qA .* v, [1 2]));
        SW = squeeze(sum(qW .* v, [1 2]));
        S0 = squeeze(sum(fB.weight2 .* v, [1 2]));     % variance scale of these bins
        Dc = Sm(use) - SE(use);
        sD = std(Dc) * sqrt(nnz(use));                  % sigma of sum(Dc)
        tot = sum(S0(use));
        dev = sum(Dc) / tot; sg = sD / tot;
        devA = (sum(Sm(use)) - sum(SA(use))) / tot;
        devW = (sum(Sm(use)) - sum(SW(use))) / tot;
        okg = abs(sum(Dc)) < 4 * sD;
        fprintf(['   (b) %-8s valid fraction %-9s %7d bins, %3d chan: (measured - exact) / variance = ' ...
                 '%+.4f +- %.4f (%+.1f sigma): %s;  w^2 model %+.3f, w model %+.3f\n'], lb, ...
            gName{g}, sum(nPerC(use)), nnz(use), dev, sg, dev/sg, passStr(okg), devA, devW);
        if ~okg, fails{end+1} = sprintf('noise %s group %s', lb, gName{g}); end %#ok<AGROW>
    end
end

% (c) channel use: estimateTOA vs the exclusion rule; channels 10 and 64
q = {'Bnoise', nc.Bnoise, 'Verbose', false};
toaB = estimateTOA(fB, ifB, pulseTemplate(nBin, ephem), q{:});
usedRule = false(nCh, ifB.nSub);
for s = 1:ifB.nSub
    Ws = reshape(fB.weight(:, s, :), nBin, nCh);
    have = any(Ws > 0, 2);
    usedRule(:, s) = all(Ws(have, :) > 0, 1).';
end
vs = toaB.valid;
ruleOk = isequal(toaB.nChanUsed(vs), sum(usedRule(:, vs), 1).');
ex64 = find(~usedRule(64, :) & vs.');
okC = ruleOk && ~any(usedRule(10, :)) && numel(ex64) >= 2 && usedRule(64, find(vs, 1, 'first')) ...
    && usedRule(64, find(vs, 1, 'last'));
fprintf(['   (c) channels used per sub-int %s (= rule %d); channel 10 never; channel 64 left out in ' ...
         'sub-ints %s (gap 30-60 ms); left out in every sub-int: %s: %s\n'], ...
    mat2str(toaB.nChanUsed(vs).'), ruleOk, mat2str(ex64), mat2str(find(~any(usedRule(:, vs), 2)).'), ...
    passStr(okC));
if ~okC, fails{end+1} = 'channel use'; end

% (d) TOAs and detection
toaU = estimateTOA(fU, ifU, pulseTemplate(nBin, ephem), q{:});
v = toaU.valid & toaB.valid;
dT = (toaB.toa(v) - toaU.toa(v)) * 1e6;
okT = abs(mean(dT)) < 4 * std(dT) / sqrt(nnz(v));
detN = detectPulsar(fB, ifB, pulseTemplate(nBin, ephem), q{:});
nr = detN.noiseRatio(detN.tested);
okN = abs(mean(nr) - 1) < 4 * std(nr) / sqrt(numel(nr));
% paired: the unblanked fold with the same channels per sub-int as the blanked one
fR = fU;
keepC = reshape(usedRule, 1, ifB.nSub, nCh);
fR.weight  = repmat(fU.weight, 1, 1, nCh) .* keepC;
fR.weight2 = repmat(fU.weight2, 1, 1, nCh) .* keepC;
fR.weightX = repmat(fU.weightX, 1, 1, nCh, 1) .* keepC;
detR = detectPulsar(fR, ifU, pulseTemplate(nBin, ephem), q{:});
tt = detN.tested & detR.tested;
dNR = detN.noiseRatio(tt) - detR.noiseRatio(tt);
okNR = abs(mean(dNR)) < 4 * std(dNR) / sqrt(nnz(tt));
fprintf(['   (d) TOAs blanked - unblanked: %d, mean %+.3f +- %.3f us (rms %.3f us; errors %.3f vs %.3f ' ...
         'us): %s; red. chi2 median %.3f; noise ratio %.4f +- %.4f: %s\n'], nnz(v), mean(dT), ...
    std(dT)/sqrt(nnz(v)), rms(dT), median(toaB.toaErr(v))*1e6, median(toaU.toaErr(v))*1e6, ...
    passStr(okT), median(toaB.redChi2(v)), mean(nr), std(nr)/sqrt(numel(nr)), passStr(okN));
fprintf(['       noise ratio blanked - unblanked with the same channels, per sub-int: %+.4f +- %.4f ' ...
         '(unblanked same channels %.4f): %s\n'], mean(dNR), std(dNR)/sqrt(nnz(tt)), ...
    mean(detR.noiseRatio(tt)), passStr(okNR));
if ~okT, fails{end+1} = 'TOAs'; end
if ~okN, fails{end+1} = 'noise ratio'; end
if ~okNR, fails{end+1} = 'noise ratio paired'; end

% ---------------------------------------------------------------------------------
if isempty(fails)
    fprintf('testBlankingWeights: ALL PASSED\n');
else
    error('testBlankingWeights:failed', 'FAILED: %s', strjoin(fails, ', '));
end
end


% =================================================================================
function t = pulseTemplate(nBin, ephem)
t = gaussianTemplate(nBin, ephem.profileFWHM);
end

function mask = makeMask3(info_chan, ephem, refFreq, nCh, Nc, fs)
% The end-to-end mask of test 3 (channel samples, 1-based).
rng(3);
sOf = @(t) round((t - info_chan.t0) * fs) + 1;
rows = {};
len = round(20e-6 * fs);                             % (i) random 20 us blanks, ~5 %
for c = 1:4:nCh
    s = 1; r = zeros(0, 3);
    while true
        s = s + round(-log(rand) * len / 0.05 * 0.95);
        if s + len > Nc, break; end
        r(end+1, :) = [c, s, s + len - 1]; %#ok<AGROW>
        s = s + len;
    end
    rows{end+1} = r; %#ok<AGROW>
end
imp = sort(randi(Nc - 10, 10, 1));                   % (ii) broadband 1 us impulses
for c = 1:nCh
    rows{end+1} = [c*ones(10, 1), imp, imp + 4]; %#ok<AGROW>
end
for c = 31:34                                        % (iii) radar 2 us at 373 Hz
    s = sOf(1e-3 : 1/373 : info_chan.t0 + Nc/fs).';
    rows{end+1} = [c*ones(numel(s), 1), s, s + 8]; %#ok<AGROW>
end
for c = 50:53                                        % (iv) radar 2 us at 400 Hz (phase-locked)
    s = sOf(0.37e-3 : 1/400 : info_chan.t0 + Nc/fs).';
    rows{end+1} = [c*ones(numel(s), 1), s, s + 8]; %#ok<AGROW>
end
Kdm = 4.148808e3 * 1e12 * ephem.DM;                  % (v) 30 us on the pulse, every turn
fc  = info_chan.chanFreqs;
half = round(15e-6 * fs);
for c = 90:100
    tP = ephem.TRef + (-2:12).' / ephem.f0 + Kdm * (1/fc(c)^2 - 1/refFreq^2);
    s = sOf(tP);
    rows{end+1} = [c*ones(numel(s), 1), s - half, s + half]; %#ok<AGROW>
end
for c = 110:115                                      % (vi) 100 us blanks every 1.1 ms
    s = sOf(0.2e-3 : 1.1e-3 : info_chan.t0 + Nc/fs).';
    rows{end+1} = [c*ones(numel(s), 1), s, s + round(100e-6 * fs) - 1]; %#ok<AGROW>
end
rows{end+1} = [64, sOf(30e-3), sOf(60e-3)];          % (vii) a 30 ms gap
rows{end+1} = [10, 1, Nc];                           % (viii) a dead channel
mask = vertcat(rows{:});
end

function files = blankChannelFiles(info_chan, mask, outBase)
% Copies of the channel IQ files with the mask applied (samples set to 0).
nCh = info_chan.nChan;
files = strings(nCh, 1);
for c = 1:nCh
    fid = fopen(info_chan.chanFiles(c), 'r', 'ieee-le');
    x = fread(fid, [2 Inf], 'single=>single');
    fclose(fid);
    keep = keepOf(mask, c, size(x, 2));
    x(:, ~keep) = 0;
    files(c) = sprintf('%s_ch%03d.dat', outBase, c);
    fid = fopen(files(c), 'w', 'ieee-le');
    fwrite(fid, x, 'single');
    fclose(fid);
end
end

function keep = keepOf(mask, c, Nc)
keep = true(Nc, 1);
r = mask(mask(:, 1) == c, 2:3);
for i = 1:size(r, 1)
    a = max(1, r(i, 1)); b = min(Nc, r(i, 2));
    if a <= b, keep(a:b) = false; end
end
end

function [h, nF, nP] = channelKernel(info_dc, c, tmp)
% The channel's dedispersion filter h(l), l = -nF..nP: an impulse through it.
nF = info_dc.nFuture(c); nP = info_dc.nPast(c);
imp = zeros(nF + nP + 1, 1); imp(nF + 1) = 1;
fi = fullfile(tmp, 'kimp.dat'); fo = fullfile(tmp, 'kout.dat');
fid = fopen(fi, 'w', 'ieee-le'); fwrite(fid, [imp.'; zeros(1, numel(imp))], 'single'); fclose(fid);
dF = info_dc.chanWidth; fcn = info_dc.chanFreqs(c);
applyInverseDispersion(fi, fo, info_dc.fs, fcn, info_dc.DM, fcn - dF/2, fcn + dF/2, ...
    'RefFreq', info_dc.refFreq, 'AllowRefOutsideBand', true, 'EdgeFrac', info_dc.edgeWidth / dF, ...
    'Nfft', info_dc.Nfft(c), 'SaveInfo', false, 'Verbose', false);
fid = fopen(fo, 'r', 'ieee-le'); raw = fread(fid, [2 Inf], 'single=>double'); fclose(fid);
h = complex(raw(1, :), raw(2, :)).';
end

function [We, Ve, Xe] = directStats(h, nF, nP, keep, k, n, Lm, m0, rad, Rx)
% W, V, X(1..Lm) of bin k from the sample covariance matrix C = H*diag(keep)*H'
% (white input, for V and X) and, for W, with the channel autocorrelation Rx:
% E|y(a)|^2 = u' * R * u with u(s) = h(a-s) keep(s), R Toeplitz from Rx.
Nc = numel(keep);
a = (k - 1)*n + (1:(Lm + 1)*n).';                  % samples of bins k .. k+Lm
j = a(1) - nP : a(end) + nF;                       % inputs that reach them
lag = a - j;
in = lag >= -nF & lag <= nP;
Hm = zeros(size(lag));
Hm(in) = h(lag(in) + nF + 1);
kj = false(1, numel(j));
ok = j >= 1 & j <= Nc; kj(ok) = keep(j(ok));
G = (Hm .* kj) * Hm' / m0;                         % C(a, b) / m0
A2 = abs(G).^2;
Rk = [flipud(Rx(2:end)); Rx(:)].';                 % lags -nR..nR
U  = Hm(1:n, :) .* kj;                             % u for the n samples of bin k
Wa = real(sum(U .* conv2(conj(U), Rk, 'same'), 2));
hr = h(:).';
We = mean(Wa) / real(sum(hr .* conv(conj(hr), Rk, 'same')));
Ve = sum(A2(1:n, 1:n), 'all') / n^2 / rad;
Xe = zeros(1, Lm);
for L = 1:Lm
    Xe(L) = sum(A2(1:n, L*n + 1:L*n + n), 'all') / n^2 / rad;
end
end

function Y = readWeights(bw, prec)
if nargin < 2, prec = 'double'; end
fid = fopen(bw.file, 'r', 'ieee-le');
c = onCleanup(@() fclose(fid));
Y = fread(fid, Inf, ['single=>' prec]);
Y = reshape(Y, bw.nChan, 2 + bw.Lmax, bw.N);
end

function writeArray(file, Y)
fid = fopen(file, 'w', 'ieee-le');
fwrite(fid, Y, 'single');
fclose(fid);
end

function p = readF32(file)
fid = fopen(file, 'r', 'ieee-le');
c = onCleanup(@() fclose(fid));
p = fread(fid, Inf, 'single=>single');
end

function s = passStr(p)
if p, s = 'PASS'; else, s = 'FAIL'; end
end

function deleteFiles(files)
for i = 1:numel(files)
    if isfile(files{i}), delete(files{i}); end
end
end

function cleanFolder(d)
f = dir(fullfile(d, '*'));
for i = 1:numel(f)
    if ~f(i).isdir, delete(fullfile(d, f(i).name)); end
end
end
