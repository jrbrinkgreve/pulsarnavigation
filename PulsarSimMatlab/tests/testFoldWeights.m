function testFoldWeights()
%TESTFOLDWEIGHTS  Unit tests for the noise covariance (A1a) and data weights (A1b) in foldProfile.
%{
Run from the PulsarSimMatlab folder: run('tests/testFoldWeights.m'). ~7 s (measured; tests 1-3 ~1.5 s).
Needs data/mc/foldRef_pre3c.mat (tests/makeFoldReference.m), the full-band
power data/test_envelope.dat (main.m) and the 128-channel power file
data/chan/test_dedisp_chan_power.dat (tests/testDetectChannels.m).

  1. Regression: without NoiseCoeffs, foldProfile reproduces the four
     reference folds made before A1 bit for bit (all fold fields and the
     info fields that existed then). The inputs must be the files the
     references were made from (same seeds -> same files).
  2. Algebra: weight2 and weightX against the explicit covariance of the
     phase-bin sums, Q = A' * C * A per sub-integration (A: assignment of
     time bins to phase bins; C: banded time-bin covariance, V on the
     diagonal and X(L) at +-L). Synthetic 0.96 us bins, P ~ 10 ms, 4096 phase
     bins (2.4 us: correlated time bins reach 4 phase bins further), F1 ~= 0,
     chunks of 1000 bins (pairs across chunk boundaries), 3 turns per sub-int
     (phase wraps inside a sub-int); linear and nearest. Plus: zero lag
     coefficients leave the fold bit-identical to the default.
  3. Physics: the 128-channel seed-43 fold with the channel coefficients
     (info.noise of detectChannels). Off-pulse profile bins, normalized as
     u_j = (prof_j/mu_c - 1) * weight_j / sqrt(radiometer), have
     E[u_j^2] = weight2_j and E[u_j*u_j+d] = weightX_j,d. Measured means vs
     the new model (pass within 4 sigma; sigma from the scatter between the
     128 independent channels) and vs the old model (independent time bins).
  4. Data weights, constant stream (W = 1, V and X of test 3 as float32):
     per-channel weights bit-identical to the NoiseCoeffs fold in every
     channel (same chunks).
  5. Data weights, random streams (3 channels, W in [0, 1] with 10 % zeros,
     random V and X(1..7)) on the geometry of test 2: weight = A' * W,
     sum = A' * (P where W > 0) (empty bins add no power) and
     weight2 / weightX = A' * C * A per sub-int and channel (C from the
     per-bin V and X); prof = sum ./ weight per channel.
  6. Fake blanking of the 128-channel data: whole detected time bins zeroed
     (15 % random per channel + every other bin in the off-pulse window at
     phase 0.30-0.36), exact stream W = keep, V*keep, X(L)*keep_k*keep_k+L.
     Mean level (vs the unblanked fold) unbiased in and outside the window
     (the naive fold, blanking ignored, shown for contrast); measured
     variance and lag covariances vs weight2 / weightX per channel, for all
     off-pulse bins and for the window alone (pass within 4 sigma).
     (Real blanking acts on voltages before dedispersion; its streams come
     from mask convolutions, A3/B.)
Errors at the end if any check fails. Writes ~0.5 GB of temporary files to
tempdir (deleted at the end).
%}

root = fullfile(fileparts(mfilename('fullpath')), '..');
addpath(fullfile(root, 'functions'));
oldDir = cd(root); restoreDir = onCleanup(@() cd(oldDir));
pipelineParams;
tmp = fullfile(tempdir, 'testFoldWeights');
if ~isfolder(tmp), mkdir(tmp); end
cleanTmp = onCleanup(@() delete(fullfile(tmp, '*.dat')));   % ~0.5 GB of weight files
fails = {};

% ---------------------------------------------------------------------------------
% 1. Regression against the folds made before A1
% ---------------------------------------------------------------------------------
refFile = fullfile(dataDir, 'mc', 'foldRef_pre3c.mat');
R = load(refFile); ref = R.ref;
detA = loadInfo(fileEnvelope);
detD = loadInfo(fullfile(dataDir, 'chan', 'test_dedisp_chan_power.dat'));
base = {'F1', ephem.F1, 'TRef', ephem.TRef, 'NBin', nBin, 'SaveFile', false, 'Verbose', false};
cs = {'A', detA, {'SubintPeriods', subintPeriods},                     'full band, linear'
      'B', detA, {'SubintPeriods', subintPeriods, 'Assign', 'nearest'}, 'full band, nearest'
      'C', detA, {'SubintPeriods', 3},                                  'full band, 3 turns per sub-int'
      'D', detD, {'SubintPeriods', subintPeriods},                      '128 channels, linear'};
for i = 1:size(cs, 1)
    [ifo, fo] = foldProfile(cs{i, 2}, 'x', ephem.f0, base{:}, cs{i, 3}{:});
    r  = ref.(cs{i, 1});
    fn = setdiff(fieldnames(r.info), {'elapsed'});
    sameFold = isequaln(fo, r.fold);
    sameInfo = all(cellfun(@(f) isequaln(ifo.(f), r.info.(f)), fn));
    fprintf('1%s. regression, %s: fold bit-identical %d, info %d: %s\n', lower(cs{i, 1}), ...
        cs{i, 4}, sameFold, sameInfo, passStr(sameFold && sameInfo));
    if ~(sameFold && sameInfo)
        fails{end+1} = ['regression ' cs{i, 1}]; %#ok<AGROW>
        dIn = dir(cs{i, 2}.file); dRef = dir(refFile);
        if dIn.datenum > dRef.datenum
            fprintf('    note: %s is newer than the reference (regenerated since?)\n', cs{i, 2}.file);
        end
    end
    if cs{i, 1} == 'D', foldOld = fo; tOld = ifo.elapsed; end
end

% ---------------------------------------------------------------------------------
% 2. weight2 / weightX against the explicit covariance A' * C * A
% ---------------------------------------------------------------------------------
dt2 = 0.96e-6; nb2 = 4096; f02 = 100.3; F12 = -0.5; TRef2 = 1.7e-3; Phi02 = 0.2; per2 = 3;
nu  = [0.9, 0.3, 0.2, 0.1, 0.06, 0.04, 0.02, 0.01];    % V, X(1..7), exaggerated
N2  = 70000;
det2 = struct('file', fullfile(tmp, 'p.dat'), 'nChan', 1, 'N', N2, 'binTime0', 0.3e-3, ...
    'binDt', dt2, 'fullySupportedBins', [11, N2 - 7], 'byteOrder', 'ieee-le', 'chanFreqs', 0);
writeF32(det2.file, rand(N2, 1));
a2 = {'F1', F12, 'TRef', TRef2, 'Phi0', Phi02, 'NBin', nb2, 'SubintPeriods', per2, ...
      'ChunkBins', 1000, 'SaveFile', false, 'Verbose', false};
k  = (det2.fullySupportedBins(1):det2.fullySupportedBins(2)).';
t  = det2.binTime0 + (k - 1) * dt2;
ph = Phi02 + f02*(t - TRef2) + 0.5*F12*(t - TRef2).^2;   % as in foldProfile
turn = round(ph);
x  = (ph - turn) * nb2;
sb = floor((turn - turn(1)) / per2) + 1;
Cd = [fliplr(nu(2:end)), nu];                          % diagonals -Lmax..Lmax
Lm = numel(nu) - 1;
j  = (1:nb2).';
for asg = ["linear", "nearest"]
    [ifo2, f2] = foldProfile(det2, 'x', f02, a2{:}, 'Assign', asg, 'NoiseCoeffs', nu);
    [~, f2d]   = foldProfile(det2, 'x', f02, a2{:}, 'Assign', asg);
    [~, f2z]   = foldProfile(det2, 'x', f02, a2{:}, 'Assign', asg, 'NoiseCoeffs', [1 0 0 0]);
    if asg == "linear"
        j0 = floor(x); a = x - j0;
        cols = [mod(j0, nb2), mod(j0 + 1, nb2)] + 1;  vals = [1 - a, a];
    else
        cols = mod(round(x), nb2) + 1;  vals = ones(size(x));
    end
    D = ifo2.covLags;
    errQ = 0; errW = 0; outside = 0;
    for s = 1:ifo2.nSub
        rr = find(sb == s); ns = numel(rr);
        A  = sparse(repmat((1:ns).', 1, size(cols, 2)), cols(rr, :), vals(rr, :), ns, nb2);
        C  = spdiags(repmat(Cd, ns, 1), -Lm:Lm, ns, ns);
        Q  = A.' * C * A;
        q0 = full(Q(sub2ind([nb2 nb2], j, j)));
        sc = max(q0);
        errQ = max(errQ, max(abs(q0 - f2.weight2(:, s))) / sc);
        for d = 1:D
            qd = full(Q(sub2ind([nb2 nb2], j, mod(j - 1 + d, nb2) + 1)));
            errQ = max(errQ, max(abs(qd - f2.weightX(:, s, 1, d))) / sc);
        end
        [ii, jj, vv] = find(Q);
        dd = mod(jj - ii, nb2);
        outside = max([outside; abs(vv(dd > D & dd < nb2 - D))]);
        wE = full(sum(A, 1)).';
        errW = max(errW, max(abs(wE - f2.weight(:, s))) / max(wE));
    end
    zeroSame = isequal(f2z.weight2, f2d.weight2) && isequal(f2z.weightX(:, :, 1, 1), f2d.weightX) && ...
        ~any(f2z.weightX(:, :, 1, 2:end), 'all') && isequaln(f2z.prof, f2d.prof) && ...
        isequal(f2.sum, f2d.sum);
    pass = errQ < 1e-12 && errW < 1e-12 && outside == 0 && zeroSame;
    fprintf(['2. %s, %d sub-ints, D = %d lags: weight2/weightX vs A''CA max rel err %.1e, ' ...
             'weight %.1e, beyond lag D %g; zero lags = default %d: %s\n'], asg, ifo2.nSub, D, ...
        errQ, errW, outside, zeroSame, passStr(pass));
    if ~pass, fails{end+1} = char("algebra " + asg); end %#ok<AGROW>
end

% ---------------------------------------------------------------------------------
% 3. Measured noise of the 128-channel fold vs new and old model
% ---------------------------------------------------------------------------------
nc = detD.noise;
[ifoN, fN] = foldProfile(detD, 'x', ephem.f0, base{:}, 'SubintPeriods', subintPeriods, ...
    'NoiseCoeffs', [nc.V, nc.X]);
nC = ifoN.nChan;
sameData = isequal(fN.sum, foldOld.sum) && isequal(fN.weight, foldOld.weight);
ok  = repmat((fN.phase.' >= 0.15 & fN.phase.' <= 0.85) & fN.weight > 0, 1, 1, nC);  % off-pulse
Pz  = fN.prof;  Pz(~ok) = 0;
mu  = sum(Pz, [1 2]) / nnz(ok(:, :, 1));                    % 1 x 1 x nC
u   = (fN.prof ./ mu - 1) .* fN.weight / sqrt(nc.radiometer);
u(~ok) = NaN;
fprintf(['3. 128-channel fold, %d off-pulse bins x %d channels, Lmax %d -> D = %d lags ' ...
         '(fold %.2f s, default %.2f s); profiles unchanged %d\n'], nnz(ok(:, :, 1)), nC, ...
    nc.Lmax, ifoN.covLags, ifoN.elapsed, tOld, sameData);
pass = sameData;
lbl = {'variance (d = 0)', 'lag 1 covariance', 'lag 2 covariance'};
for d = 0:2
    pr = u .* circshift(u, -d, 1);                         % u_j * u_j+d (circular per profile)
    v  = ~isnan(pr(:, :, 1));
    mc = squeeze(sum(pr .* v, [1 2], 'omitnan')) / nnz(v); % mean per channel
    meas = mean(mc); sig = std(mc) / sqrt(nC);
    if d == 0
        qN = fN.weight2; qO = foldOld.weight2;
    else
        qN = fN.weightX(:, :, 1, d);
        if d == 1, qO = foldOld.weightX; else, qO = zeros(size(qN)); end
    end
    pN = mean(qN(v)); pO = mean(qO(v));
    okd = abs(meas - pN) < 4*sig;
    fprintf(['   %-17s measured %.5f +- %.5f; new %.5f (%+.1f sigma); old %.5f (%+.1f sigma): %s\n'], ...
        lbl{d+1}, meas, sig, pN, (meas - pN)/sig, pO, (meas - pO)/sig, passStr(okd));
    pass = pass && okd;
end
if ~pass, fails{end+1} = 'measured fold noise'; end

% ---------------------------------------------------------------------------------
% 4. Constant data weights (W = 1, V and X as NoiseCoeffs) = the A1a fold
% ---------------------------------------------------------------------------------
nCh = detD.nChan; Nd = detD.N; Lmx = nc.Lmax;
coefS = double(single([nc.V, nc.X]));               % the values the float32 file holds
infoW4 = struct('file', fullfile(tmp, 'w_const.dat'), 'nChan', nCh, 'N', Nd, 'Lmax', Lmx, ...
    'byteOrder', 'ieee-le');
writeWeightBlocks(infoW4.file, Nd, @(kk) repmat(single([1, coefS]), nCh, 1, numel(kk)));
cb = 20000;                                          % same chunks for both folds
[i4a, f4a] = foldProfile(detD, 'x', ephem.f0, base{:}, 'SubintPeriods', subintPeriods, ...
    'NoiseCoeffs', coefS, 'ChunkBins', cb);
[i4w, f4w] = foldProfile(detD, 'x', ephem.f0, base{:}, 'SubintPeriods', subintPeriods, ...
    'DataWeights', infoW4, 'ChunkBins', cb);
one = ones(1, 1, nCh);
same = isequal(f4w.weight, f4a.weight .* one) && isequal(f4w.weight2, f4a.weight2 .* one) && ...
    isequal(f4w.weightX, repmat(f4a.weightX, 1, 1, nCh, 1)) && isequal(f4w.sum, f4a.sum) && ...
    isequaln(f4w.prof, f4a.prof) && isequaln(f4w.profTotal, f4a.profTotal);
dmax = max(abs(f4w.weight2 - f4a.weight2 .* one), [], 'all') / max(f4a.weight2, [], 'all');
fprintf(['4. constant data weights (%d chan): per-channel weights bit-identical to the A1a fold %d ' ...
         '(max rel diff weight2 %.1e); fold %.2f s vs %.2f s: %s\n'], nCh, same, dmax, ...
    i4w.elapsed, i4a.elapsed, passStr(same));
if ~same, fails{end+1} = 'constant data weights'; end

% ---------------------------------------------------------------------------------
% 5. Random data weights: weight = A'W, weight2/weightX = A'C A per channel
% ---------------------------------------------------------------------------------
nC5 = 3; Lm5 = 7;
rng(11);
Wst = rand(nC5, N2); Wst(rand(nC5, N2) < 0.1) = 0;  % some bins fully blanked
Vst = 0.5 + 0.5*rand(nC5, N2);
Xst = rand(nC5, Lm5, N2) .* (0.3 ./ (1:Lm5));
Y5  = single(cat(2, reshape(Wst, nC5, 1, N2), reshape(Vst, nC5, 1, N2), Xst));
Wst = double(reshape(Y5(:, 1, :), nC5, N2));        % what the fold reads
Vst = double(reshape(Y5(:, 2, :), nC5, N2));
Xst = double(Y5(:, 3:end, :));
det5 = det2; det5.nChan = nC5; det5.chanFreqs = zeros(nC5, 1); det5.file = fullfile(tmp, 'p5.dat');
P5 = double(single(rand(nC5, N2)));                % power, also in bins with W = 0
writeF32(det5.file, P5);
infoW5 = struct('file', fullfile(tmp, 'w5.dat'), 'nChan', nC5, 'N', N2, 'Lmax', Lm5, ...
    'byteOrder', 'ieee-le');
writeWeightBlocks(infoW5.file, N2, @(kk) Y5(:, :, kk));
for asg = ["linear", "nearest"]
    [ifo5, f5] = foldProfile(det5, 'x', f02, a2{:}, 'Assign', asg, 'DataWeights', infoW5);
    if asg == "linear"
        j0 = floor(x); a = x - j0;
        cols = [mod(j0, nb2), mod(j0 + 1, nb2)] + 1;  vals = [1 - a, a];
    else
        cols = mod(round(x), nb2) + 1;  vals = ones(size(x));
    end
    D = ifo5.covLags;
    errQ = 0; errW = 0; errS = 0; outside = 0;
    for s = 1:ifo5.nSub
        rr = find(sb == s); ns = numel(rr); kk = k(rr);
        A  = sparse(repmat((1:ns).', 1, size(cols, 2)), cols(rr, :), vals(rr, :), ns, nb2);
        for c = 1:nC5
            sE = A.' * (P5(c, kk).' .* (Wst(c, kk).' > 0));   % empty bins add no power
            errS = max(errS, max(abs(sE - f5.sum(:, s, c))) / max(sE));
            I = (1:ns).'; Jc = I; Cv = Vst(c, kk).';
            for L = 1:Lm5
                i1 = (1:ns-L).';
                xv = reshape(Xst(c, L, kk(i1)), [], 1);   % covariance of bin k with k+L
                I = [I; i1; i1 + L]; Jc = [Jc; i1 + L; i1]; Cv = [Cv; xv; xv]; %#ok<AGROW>
            end
            Q  = A.' * sparse(I, Jc, Cv, ns, ns) * A;
            q0 = full(Q(sub2ind([nb2 nb2], j, j)));
            sc = max(q0);
            errQ = max(errQ, max(abs(q0 - f5.weight2(:, s, c))) / sc);
            for d = 1:D
                qd = full(Q(sub2ind([nb2 nb2], j, mod(j - 1 + d, nb2) + 1)));
                errQ = max(errQ, max(abs(qd - f5.weightX(:, s, c, d))) / sc);
            end
            [ii, jj, vv] = find(Q);
            dd = mod(jj - ii, nb2);
            outside = max([outside; abs(vv(dd > D & dd < nb2 - D))]);
            wE = A.' * Wst(c, kk).';
            errW = max(errW, max(abs(wE - f5.weight(:, s, c))) / max(wE));
        end
    end
    profOk = isequaln(f5.prof, f5.sum ./ f5.weight .* nanWhere(f5.weight == 0));
    pass = errQ < 1e-12 && errW < 1e-12 && errS < 1e-12 && outside == 0 && profOk && ...
        isequal(size(f5.weight), [nb2, ifo5.nSub, nC5]);
    fprintf(['5. random data weights, %s, %d chan x %d sub-ints, D = %d: weight vs A''W %.1e, ' ...
             'sum vs A''(P where W > 0) %.1e, weight2/weightX vs A''CA %.1e, beyond lag D %g, ' ...
             'prof per channel %d: %s\n'], asg, nC5, ifo5.nSub, D, errW, errS, errQ, outside, ...
        profOk, passStr(pass));
    if ~pass, fails{end+1} = char("random data weights " + asg); end %#ok<AGROW>
end

% ---------------------------------------------------------------------------------
% 6. Fake blanking of whole time bins in the 128-channel data
% ---------------------------------------------------------------------------------
Pall = reshape(readF32(detD.file), nCh, Nd);
kAll = 1:Nd;
phk  = mod((detD.binTime0 + (kAll - 1)*detD.binDt - ephem.TRef) * ephem.f0 + 0.5, 1) - 0.5;
rng(7);
keep = rand(nCh, Nd) >= 0.15;                        % 15 % random, per channel
win  = phk >= 0.30 & phk < 0.36;                     % off-pulse phase window ...
keep(:, win & mod(kAll, 2) == 1) = false;            % ... every other time bin blanked
detB = detD; detB.file = fullfile(tmp, 'p_blank.dat');
writeF32(detB.file, Pall .* single(keep));
coef = [nc.V, nc.X];
infoW6 = struct('file', fullfile(tmp, 'w_blank.dat'), 'nChan', nCh, 'N', Nd, 'Lmax', Lmx, ...
    'byteOrder', 'ieee-le');
writeWeightBlocks(infoW6.file, Nd, @(kk) blankStream(kk, keep, coef));
[i6, f6] = foldProfile(detB, 'x', ephem.f0, base{:}, 'SubintPeriods', subintPeriods, ...
    'DataWeights', infoW6);
[~, f6n] = foldProfile(detB, 'x', ephem.f0, base{:}, 'SubintPeriods', subintPeriods, ...
    'NoiseCoeffs', coef);                            % naive: blanking ignored
muT = mu;                                            % true level: unblanked fold (test 3)
offB = fN.phase.' >= 0.15 & fN.phase.' <= 0.85;
winB = fN.phase.' >= 0.30 & fN.phase.' < 0.36;
fprintf('6. fake blanking (15 %% random + every other bin at phase 0.30-0.36), valid fraction %.3f\n', ...
    i6.validFraction);
% mean level: weighted fold unbiased, naive fold low by the blanked fraction
pass = true;
regions = {winB, offB & ~winB}; rName = {'window', 'rest off-pulse'};
for ir = 1:2
    sel = repmat(regions{ir}, 1, i6.nSub, nCh) & f6.weight > 0;
    rw  = squeeze(mean((f6.prof ./ muT - 1) .* nanWhere(~sel), [1 2], 'omitnan'));
    rn  = squeeze(mean((f6n.prof ./ muT - 1) .* nanWhere(~sel), [1 2], 'omitnan'));
    okm = abs(mean(rw)) < 4*std(rw)/sqrt(nCh);
    fprintf('   mean level, %-15s weighted %+.5f +- %.5f, naive %+.4f: %s\n', rName{ir}, ...
        mean(rw), std(rw)/sqrt(nCh), mean(rn), passStr(okm));
    pass = pass && okm;
end
% noise: E[u_j u_j+d] = weightX(j, d) per channel, u = (prof/mu - 1) * weight / sqrt(rad)
okB = offB & f6.weight > 0;                          % NBin x nSub x nCh
u6  = (f6.prof ./ muT - 1) .* f6.weight / sqrt(nc.radiometer);
u6(~okB) = NaN;
nName = {'all off-pulse', 'window'}; nRegion = {offB, winB}; nLags = [2, 1];
for ir = 1:2
    for d = 0:nLags(ir)
        pr = u6 .* circshift(u6, -d, 1);
        if d == 0, q = f6.weight2; lb = 'variance'; else, q = f6.weightX(:, :, :, d); lb = sprintf('lag %d', d); end
        v  = ~isnan(pr) & nRegion{ir};
        df = squeeze(sum((pr - q) .* v, [1 2], 'omitnan') ./ sum(v, [1 2]));   % per channel
        ms = squeeze(sum(pr .* v, [1 2], 'omitnan') ./ sum(v, [1 2]));
        pq = squeeze(sum(q .* v, [1 2]) ./ sum(v, [1 2]));
        sg = std(df) / sqrt(nCh);
        okd = abs(mean(df)) < 4*sg;
        fprintf('   noise, %-13s %-8s measured %.5f, predicted %.5f, diff %+.5f +- %.5f (%+.1f sigma): %s\n', ...
            nName{ir}, lb, mean(ms), mean(pq), mean(df), sg, mean(df)/sg, passStr(okd));
        pass = pass && okd;
    end
end
if ~pass, fails{end+1} = 'fake blanking'; end

% ---------------------------------------------------------------------------------
if isempty(fails)
    fprintf('testFoldWeights: ALL PASSED\n');
else
    error('testFoldWeights:failed', 'FAILED: %s', strjoin(fails, ', '));
end
end


% =================================================================================
function s = passStr(p)
if p, s = 'PASS'; else, s = 'FAIL'; end
end

function writeF32(file, p)
fid = fopen(file, 'w', 'ieee-le');
fwrite(fid, single(p), 'single');
fclose(fid);
end

function p = readF32(file)
fid = fopen(file, 'r', 'ieee-le');
c = onCleanup(@() fclose(fid));
p = fread(fid, Inf, 'single=>single');
end

function z = nanWhere(m)
%NANWHERE  NaN where m is true, 1 elsewhere (to mask values in sums with 'omitnan').
z = ones(size(m));
z(m) = NaN;
end

function writeWeightBlocks(file, N, blockFun)
%WRITEWEIGHTBLOCKS  Write a data-weight file [nChan x (2+Lmax) x N] (float32) in blocks;
% blockFun(kk) returns the single array for the time bins kk.
fid = fopen(file, 'w', 'ieee-le');
c = onCleanup(@() fclose(fid));
blk = 8192;
for k0 = 1:blk:N
    fwrite(fid, blockFun(k0:min(N, k0 + blk - 1)), 'single');
end
end

function Y = blankStream(kk, keep, coef)
%BLANKSTREAM  Exact W, V, X(L) of time bins kk when whole detected bins are zeroed:
% W = keep, V = V0*keep_k, X(L) = X0(L)*keep_k*keep_k+L.
[nCh, N] = size(keep);
Lmax = numel(coef) - 1;
kp = keep(:, kk);
Y = zeros(nCh, 2 + Lmax, numel(kk), 'single');
Y(:, 1, :) = reshape(kp, nCh, 1, []);
Y(:, 2, :) = reshape(coef(1) * kp, nCh, 1, []);
for L = 1:Lmax
    kn = false(nCh, numel(kk));
    inside = kk + L <= N;
    kn(:, inside) = keep(:, kk(inside) + L);
    Y(:, 2 + L, :) = reshape(coef(1 + L) * (kp & kn), nCh, 1, []);
end
end
