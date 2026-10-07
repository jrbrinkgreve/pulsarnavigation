function testFoldWeights()
%TESTFOLDWEIGHTS  Unit tests for the noise covariance in foldProfile (A1a).
%{
Run from the PulsarSimMatlab folder: run('tests/testFoldWeights.m'). ~1.5 s (measured).
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
Errors at the end if any check fails.
%}

root = fullfile(fileparts(mfilename('fullpath')), '..');
addpath(fullfile(root, 'functions'));
oldDir = cd(root); restoreDir = onCleanup(@() cd(oldDir));
pipelineParams;
tmp = fullfile(tempdir, 'testFoldWeights');
if ~isfolder(tmp), mkdir(tmp); end
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
