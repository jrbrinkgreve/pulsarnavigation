function [info, fold] = foldProfile(info_det, outFile, f0, opts)
%FOLDPROFILE  Fold a detected power time series into pulse profiles.
%{
Folds the output of detectPower with a pulsar phase model into NBin phase
bins, per sub-integration and per channel. Sums and weights are kept
separately, so sub-integrations and channels can be combined exactly
later (TOA estimation, noise normalization).

  [info, fold] = foldProfile(info_det, outFile, f0)
  [info, fold] = foldProfile(info_det, outFile, f0, 'TRef', t0, 'NBin', 1024, ...
                             'SubintPeriods', 1)

Inputs:
  info_det  info struct from detectPower (file, nChan, N, binTime0, binDt,
            fullySupportedBins, byteOrder).
  outFile   .mat file to save fold and info to.
  f0        [Hz] spin frequency at TRef (1/period).

Name-value options:
  'F1'             [Hz/s] spin-frequency derivative (default 0).
  'TRef'           [s] reference epoch of the phase model (default 0).
  'Phi0'           [turns] phase at TRef (default 0).
                   Phase model: phi(t) = Phi0 + f0*dt + F1*dt^2/2, dt = t - TRef.
                   Phase 0 is the centre of phase bin 1. For the synthetic
                   data, TRef = info_gen.pulseCenters(1) puts every true
                   pulse peak at phase 0, so a perfect TOA has offset 0.
  'NBin'           phase bins (default 1024). Every phase bin needs data in
                   every sub-integration: time bins of at most one phase bin
                   (binDt <= 1/(f0*NBin)), or enough turns per sub-integration
                   to fill the gaps. Otherwise warning 'foldProfile:coverage'
                   (checked on the time grid before folding): estimateTOA and
                   detectPulsar (MinCoverage 1) would reject every
                   sub-integration.
  'SubintPeriods'  turns per sub-integration (default 1).
  'SubintTime'     [s] alternative: sub-integration length, rounded to
                   whole turns (overrides SubintPeriods).
  'Assign'         'linear' (default): each time bin is split between the
                   two nearest phase bins (unbiased at sub-bin level);
                   'nearest': whole time bin to the nearest phase bin
                   (half the work; bias up to ~binDt/2 depending on how the
                   time grid falls on the phase grid).
  'NoiseCoeffs'    [V X(1) ... X(Lmax)]: noise of the detected time bins
                   relative to sigma^2 = m^2/(Bnoise*binDt) (m = mean power):
                   variance V, covariance of time bins L apart X(L)
                   (powerCovariance; detectChannels stores it in info.noise).
                   Default 1: independent time bins, right for the full band
                   (Bnoise*binDt >> 1). In a 3.125 MHz channel with 0.96 us
                   bins V = 0.888, X(1) = 0.046, Lmax = 7.
  'DataWeights'    info struct of a data-weight file for the power file, e.g.
                   after RFI blanking: fields file, nChan, N, Lmax, byteOrder.
                   float32 [nChan x (2+Lmax) x N], same time grid as the power:
                   per time bin the valid fraction W of every channel, then
                   the variance V, then the covariances X(1) ... X(Lmax) with
                   the following time bins (V, X relative to sigma^2 of that
                   channel without blanking, as NoiseCoeffs). The weights are
                   then kept per channel. A time bin with W = 0 holds no data:
                   its power is not added either (whatever is left in it,
                   e.g. filter tails below blankingWeights' MinWeight, would
                   otherwise bias phase bins that have almost no weight).
                   Replaces NoiseCoeffs. Default: none
                   (W = 1 and the constants of NoiseCoeffs, shared by all
                   channels).
  'UseSupportedOnly' fold only info_det.fullySupportedBins (default true).
  'ChunkBins'      time bins read per chunk (default 2^22; with DataWeights
                   at most 2^25 / (nChan*(2+Lmax)) to limit memory).
  'MaxMemoryGB'    limit for the fold arrays (default 2).
  'SaveFile'       save fold and info to outFile (default true).
  'Verbose'        print a summary (default true).

Outputs:
  fold  struct:
    .sum     [NBin x nSub x nChan] sum of w * power (w = assignment weight)
    .weight  [NBin x nSub x nW]    data weight: sum of w * W (W = valid fraction
                                   of the time bin; 1 without DataWeights). nW =
                                   nChan with DataWeights, else 1 (the same for
                                   every channel).
    .weight2 [NBin x nSub x nW]    noise of a phase-bin sum: var(sum_j) =
                                   sigma^2 * weight2_j, sigma^2 = m^2/(Bnoise*binDt)
                                   (see NoiseCoeffs); the variance of a phase bin
                                   of prof is sigma^2 * weight2 / weight^2. For
                                   NoiseCoeffs = 1 this is the sum of w^2.
    .weightX [NBin x nSub x nW x D] covariance of phase bins d = 1..D apart
                                   (circular): cov(sum_j, sum_j+d) = sigma^2 *
                                   weightX(j, :, :, d), so for prof bins
                                   sigma_j*sigma_j+d * weightX / (weight_j*weight_j+d).
                                   For NoiseCoeffs = 1, D = 1 and weightX is
                                   [NBin x nSub]: the sum of w_j*w_(j+1) over time
                                   bins split between bin j and bin j+1 (zero for
                                   'nearest'). D > 1 when correlated time bins
                                   reach beyond the next phase bin.
                                   Needed for correct uncertainties of anything
                                   that sums over bins. Not stored: covariance
                                   between consecutive sub-integrations (time bins
                                   close to a sub-int boundary, at phase 0.5).
    .prof    [NBin x nSub x nChan] sum ./ weight (NaN where weight = 0)
    .profTotal [NBin x nChan]      all sub-integrations combined (per channel)
    .phase   [1 x NBin]            phase of each bin centre, turns in [0,1)
    .subint  struct with per-sub-integration vectors:
               turnFirst, turnLast  integer turns covered
               turnRef   reference turn (middle); its phase 0 is at
               tRef      [s] time where phi = turnRef
               fRef      [Hz] spin frequency at tRef
               tMean     [s] mean time of the folded time bins
               nTimeBins number of time bins folded
    A TOA for sub-integration s is then tRef(s) + dphi / fRef(s), with dphi
    the fitted phase offset (turns) of its profile.
  info  metadata (model, sizes, conventions, source files).

Method:
  Per chunk: t_k = binTime0 + (k-1)*binDt (double), phi(t_k), turn =
  round(phi), frac = phi - turn in [-0.5, 0.5). Sub-integration boundaries
  therefore fall at phase 0.5, away from a pulse centred at phase 0.
  Bin index from frac*NBin; accumulation with accumarray into
  (phase bin, sub-integration). Everything is vectorized; one pass over the
  (small) detected-power file.
  Noise: var(sum_j) and cov(sum_j, sum_j+d) are sums over pairs of time
  bins (k, k') of a_kj * a_k',j+d * c(k' - k), with a the assignment weights
  and c(0) = V, c(+-L) = X(L). L = 0 gives the w^2 and w_j*w_j+1 terms
  (times V). For L = 1..Lmax every pair of time bins L apart in the same
  sub-integration adds X(L)*a*a' to its pair of phase bins; the last Lmax
  time bins of a chunk are kept for the pairs with the next chunk.
  Data weights: a blanked time bin has mean power W_k * m, so prof =
  sum(a*P) / sum(a*W) is unbiased even when blanking depends on the pulse
  phase. With DataWeights the noise sums use V_k and X_k(L) of each time
  bin (X of the earlier bin of a pair), per channel.
%}

arguments
    info_det                    struct
    outFile                     {mustBeTextScalar}
    f0                    (1,1) double {mustBePositive, mustBeFinite}
    opts.F1               (1,1) double {mustBeFinite} = 0
    opts.TRef             (1,1) double {mustBeFinite} = 0
    opts.Phi0             (1,1) double {mustBeFinite} = 0
    opts.NBin             (1,1) double {mustBeInteger, mustBePositive} = 1024
    opts.SubintPeriods    (1,1) double {mustBeInteger, mustBePositive} = 1
    opts.SubintTime             double = []
    opts.Assign                 {mustBeTextScalar} = 'linear'
    opts.NoiseCoeffs      (1,:) double {mustBeNonempty, mustBeFinite} = 1
    opts.DataWeights            = []
    opts.UseSupportedOnly (1,1) logical = true
    opts.ChunkBins        (1,1) double {mustBePositive} = 2^22
    opts.MaxMemoryGB      (1,1) double {mustBePositive} = 2
    opts.SaveFile         (1,1) logical = true
    opts.Verbose          (1,1) logical = true
end

tStart = tic;
outFile = char(outFile);
assign = lower(char(opts.Assign));
if ~any(strcmp(assign, {'linear', 'nearest'}))
    error('foldProfile:assign', 'Assign must be ''linear'' or ''nearest''.');
end
if opts.NoiseCoeffs(1) <= 0
    error('foldProfile:noise', 'NoiseCoeffs(1) (the time-bin variance V) must be positive.');
end
V     = opts.NoiseCoeffs(1);
Xc    = opts.NoiseCoeffs(2:end);                    % X(1..Lmax)
Lmax  = numel(Xc);
Nbin  = opts.NBin;
nChan = info_det.nChan;
DW    = opts.DataWeights;
haveW = ~isempty(DW);
nW    = 1;                                          % channels with their own weights
if haveW
    need = {'file', 'nChan', 'N', 'Lmax', 'byteOrder'};
    if ~isstruct(DW) || ~all(isfield(DW, need))
        error('foldProfile:weights', 'DataWeights must be a struct with fields %s.', ...
            strjoin(need, ', '));
    end
    if DW.nChan ~= nChan || DW.N ~= info_det.N
        error('foldProfile:weights', ['DataWeights (%d chan x %d bins) does not match ' ...
            'the power file (%d x %d).'], DW.nChan, DW.N, nChan, info_det.N);
    end
    if ~isequal(opts.NoiseCoeffs, 1)
        error('foldProfile:weights', 'Give either NoiseCoeffs or DataWeights, not both.');
    end
    Lmax = DW.Lmax;
    nW   = nChan;
end
nQ    = 2 + Lmax;                                   % values per bin and channel in DataWeights
dt    = info_det.binDt;
t1    = info_det.binTime0;
F1    = opts.F1;
TRef  = opts.TRef;
Phi0  = opts.Phi0;
phiOf = @(t) Phi0 + f0*(t - TRef) + 0.5*F1*(t - TRef).^2;

% ---- Time-bin range --------------------------------------------------------------
if opts.UseSupportedOnly
    kRange = info_det.fullySupportedBins;
else
    kRange = [1, info_det.N];
end
k1 = kRange(1); k2 = kRange(2);
if k2 < k1
    error('foldProfile:empty', 'No time bins to fold.');
end
tFirst = t1 + (k1 - 1)*dt;
tLast  = t1 + (k2 - 1)*dt;
if f0 + F1*(tFirst - TRef) <= 0 || f0 + F1*(tLast - TRef) <= 0
    error('foldProfile:model', 'Spin frequency is not positive over the data span.');
end

% ---- Sub-integrations ---------------------------------------------------------------
if ~isempty(opts.SubintTime)
    nPer = max(1, round(opts.SubintTime * f0));
else
    nPer = opts.SubintPeriods;
end
turnFirst = round(phiOf(tFirst));
turnLast  = round(phiOf(tLast));
nSub = floor((turnLast - turnFirst) / nPer) + 1;

% ---- Phase-bin lags of the noise covariance ---------------------------------------------
% Time bins Lmax apart are at most Lmax*dt*fMax*NBin phase bins apart; the
% linear split adds one more. D = 1 without correlated time bins.
fMax = max(f0 + F1*(tFirst - TRef), f0 + F1*(tLast - TRef));
D = ceil(Lmax * dt * fMax * Nbin * (1 + 1e-9)) + 1;
if Lmax > 0 && 2*D >= Nbin
    error('foldProfile:lags', ['Correlated time bins span %d of %d phase bins; ' ...
        'use fewer lags or fewer phase bins.'], D, Nbin);
end

% ---- Phase coverage of a complete sub-integration ---------------------------------------
% Consecutive time bins are f*dt*NBin phase bins apart. Above 1, phase bins can stay
% empty in every sub-integration, and estimateTOA / detectPulsar (MinCoverage 1)
% then reject them all (f_out sweep 9 Oct: 100 kHz, NBin 2048, 1 turn -> no TOAs).
% Checked on the time grid alone (no data): the time bins of sub-integration 2, the
% first complete one (the first and last are usually partial), assigned as in the
% fold. Skipped without a complete sub-integration or above 2^22 time bins.
if fMax * dt * Nbin > 1 && nSub >= 3
    fMin = min(f0 + F1*(tFirst - TRef), f0 + F1*(tLast - TRef));
    nCov = min(k2 - k1 + 1, ceil((2*nPer + 2) / (fMin * dt)));   % reaches past sub-int 2
    if nCov <= 2^22
        phCov   = phiOf(t1 + ((k1 : k1 + nCov - 1).' - 1) * dt);
        turnCov = round(phCov);
        xCov    = (phCov - turnCov) * Nbin;
        xCov    = xCov(floor((turnCov - turnFirst) / nPer) == 1);   % sub-integration 2
        if strcmp(assign, 'nearest')
            hit = mod(round(xCov), Nbin);
        else                                        % both neighbours with weight > 0
            j0  = floor(xCov);
            hit = [mod(j0, Nbin); mod(j0(xCov > j0) + 1, Nbin)];
        end
        nHit = numel(unique(hit));
        if nHit < Nbin
            warning('foldProfile:coverage', ...
                ['NBin = %d is too fine for %.4g us time bins: a complete sub-integration ' ...
                 '(%d turn(s)) puts data in only %d of %d phase bins, so estimateTOA and ' ...
                 'detectPulsar (MinCoverage 1) reject every sub-integration. Use NBin <= %d ' ...
                 '(time bins per turn), time bins <= %.4g us (f_out >= %.4g Hz), or more ' ...
                 'turns per sub-integration.'], Nbin, dt*1e6, nPer, nHit, Nbin, ...
                floor(1/(fMax*dt)), 1e6/(fMax*Nbin), fMax*Nbin);
        end
    end
end

memGB = Nbin * nSub * (2*nChan + nW*(2 + D)) * 8 / 1e9;
if memGB > opts.MaxMemoryGB
    error('foldProfile:memory', ...
        ['Fold arrays need ~%.2f GB (%d bins x %d sub-integrations x %d chan). ' ...
         'Use longer sub-integrations, fewer bins, or raise MaxMemoryGB.'], ...
        memGB, Nbin, nSub, nChan);
end

if opts.Verbose
    fprintf(['foldProfile: bins %d..%d (%.4g s), f0 = %.10g Hz, F1 = %.3g Hz/s, ' ...
             'NBin = %d (%.3g us), %d sub-int(s) of %d turn(s), %s assignment\n'], ...
        k1, k2, tLast - tFirst + dt, f0, F1, Nbin, 1e6/(f0*Nbin), nSub, nPer, assign);
    if haveW
        fprintf(['foldProfile: data weights per channel from %s (Lmax %d) -> ' ...
                 'phase-bin covariance over %d lag(s)\n'], DW.file, Lmax, D);
    elseif Lmax > 0
        fprintf(['foldProfile: correlated time bins (V %.4f, X(1) %.4f, Lmax %d) -> ' ...
                 'phase-bin covariance over %d lag(s)\n'], V, Xc(1), Lmax, D);
    end
end

% ---- Accumulate ------------------------------------------------------------------------
nCell   = Nbin * nSub;
S       = zeros(nCell, nChan);
Wt      = zeros(nCell, nW);
W2      = zeros(nCell, nW);
WX      = zeros(nCell, nW, D);
off     = nCell * (0:nW-1);                         % channel offsets into Wt, W2, WX
nA      = 1 + strcmp(assign, 'linear');             % phase bins per time bin
tail    = struct('jA', zeros(0, nA), 'wA', zeros(0, nA), 'sub', zeros(0, 1), ...
                 'X', zeros(0, nW, Lmax));
tSum    = zeros(nSub, 1);
nTB     = zeros(nSub, 1);
sumW    = 0;                                        % for the valid fraction

[fid, msg] = fopen(info_det.file, 'r', info_det.byteOrder);
if fid == -1
    error('foldProfile:open', 'Could not open "%s": %s', info_det.file, msg);
end
cleanupIn = onCleanup(@() fclose(fid)); %#ok<NASGU>
if fseek(fid, (k1 - 1) * nChan * 4, 'bof') ~= 0
    error('foldProfile:seek', 'Could not seek in "%s".', info_det.file);
end
chunkMax = opts.ChunkBins;
if haveW
    [fidW, msg] = fopen(DW.file, 'r', DW.byteOrder);
    if fidW == -1
        error('foldProfile:open', 'Could not open "%s": %s', DW.file, msg);
    end
    cleanupW = onCleanup(@() fclose(fidW)); %#ok<NASGU>
    if fseek(fidW, (k1 - 1) * nChan * nQ * 4, 'bof') ~= 0
        error('foldProfile:seek', 'Could not seek in "%s".', DW.file);
    end
    chunkMax = min(chunkMax, max(1, floor(2^25 / (nChan * nQ))));
end

kPos = k1;
while kPos <= k2
    n = min(chunkMax, k2 - kPos + 1);
    X = fread(fid, [nChan n], 'single=>double');
    if size(X, 2) < n
        error('foldProfile:read', 'Unexpected end of "%s".', info_det.file);
    end
    if haveW                                         % per bin and channel: W, V, X(1..Lmax)
        Y = fread(fidW, [nChan*nQ n], 'single=>double');
        if size(Y, 2) < n
            error('foldProfile:read', 'Unexpected end of "%s".', DW.file);
        end
        Y  = permute(reshape(Y, nChan, nQ, n), [3 1 2]);   % n x nChan x nQ
        Wk = Y(:, :, 1);
        Vk = Y(:, :, 2);
        Xk = Y(:, :, 3:end);                         % n x nChan x Lmax
        sumW = sumW + sum(Wk, 'all');
        X(Wk.' == 0) = 0;                            % empty bins add no power
    else                                             % the same constants for every bin
        Wk = 1;
        Vk = V;
        Xk = Xc;
        sumW = sumW + n;
    end
    k  = (kPos : kPos + n - 1).';
    t  = t1 + (k - 1) * dt;
    ph = phiOf(t);
    turn = round(ph);
    x  = (ph - turn) * Nbin;                         % in [-Nbin/2, Nbin/2]
    sub = floor((turn - turnFirst) / nPer) + 1;
    base = Nbin * (sub - 1);

    tSum = tSum + accumarray(sub, t, [nSub 1]);
    nTB  = nTB  + accumarray(sub, 1, [nSub 1]);

    if strcmp(assign, 'nearest')
        L = mod(round(x), Nbin) + 1 + base;
        Wt = addTo(Wt, L + off, Wk);
        W2 = addTo(W2, L + off, Vk);
        for c = 1:nChan
            S(:, c) = S(:, c) + accumarray(L, X(c, :).', [nCell 1]);
        end
    else
        j0 = floor(x);
        a  = x - j0;
        L  = [mod(j0, Nbin) + 1 + base; mod(j0 + 1, Nbin) + 1 + base];
        w  = [1 - a; a];
        if haveW, Wk2 = [Wk; Wk]; Vk2 = [Vk; Vk]; else, Wk2 = Wk; Vk2 = Vk; end
        Wt = addTo(Wt, L + off, w .* Wk2);
        W2 = addTo(W2, L + off, Vk2 .* w.^2);
        WX = addTo(WX, mod(j0, Nbin) + 1 + base + off, Vk .* (1 - a) .* a);
        for c = 1:nChan
            v = X(c, :).';
            S(:, c) = S(:, c) + accumarray(L, w .* [v; v], [nCell 1]);
        end
    end
    if Lmax > 0                                      % correlated time bins: lag pairs
        if nA == 1
            jA = mod(round(x), Nbin);  wA = ones(n, 1);
        else
            jA = [mod(j0, Nbin), mod(j0 + 1, Nbin)];  wA = [1 - a, a];
        end
        [W2, WX, tail] = addLagPairs(W2, WX, tail, jA, wA, sub, Xk, haveW, Nbin);
    end
    kPos = kPos + n;
end

% ---- Assemble fold struct ----------------------------------------------------------------
fold = struct();
fold.sum     = reshape(S,  Nbin, nSub, nChan);
fold.weight  = reshape(Wt, Nbin, nSub, nW);
fold.weight2 = reshape(W2, Nbin, nSub, nW);
fold.weightX = reshape(WX, Nbin, nSub, nW, D);
wFull = repmat(fold.weight, 1, 1, nChan / nW);
fold.prof = fold.sum ./ wFull;
fold.prof(wFull == 0) = NaN;
wTot = reshape(sum(fold.weight, 2), Nbin, nW);
fold.profTotal = reshape(sum(fold.sum, 2), Nbin, nChan) ./ wTot;
fold.profTotal(repmat(wTot == 0, 1, nChan / nW)) = NaN;
fold.phase = (0:Nbin-1) / Nbin;

sI = (1:nSub).';
tf = turnFirst + (sI - 1) * nPer;
tl = min(tf + nPer - 1, turnLast);
tr = round((tf + tl) / 2);
d  = tr - Phi0;                                     % phase to advance from TRef
dtRef = 2*d ./ (f0 + sqrt(f0^2 + 2*F1*d));          % stable root of phi(t) = tr
fold.subint = struct( ...
    'turnFirst', tf, 'turnLast', tl, 'turnRef', tr, ...
    'tRef', TRef + dtRef, 'fRef', f0 + F1*dtRef, ...
    'tMean', tSum ./ max(nTB, 1), 'nTimeBins', nTB);

% ---- Info ------------------------------------------------------------------------------------
info = struct();
info.file         = outFile;
info.sourceFile   = info_det.file;
info.f0           = f0;
info.F1           = F1;
info.TRef         = TRef;
info.Phi0         = Phi0;
info.phaseModel   = 'phi(t) = Phi0 + f0*(t-TRef) + F1*(t-TRef)^2/2';
info.phaseConvention = 'bin j (1-based) centred at phase (j-1)/NBin; turn = round(phi)';
info.NBin         = Nbin;
info.nSub         = nSub;
info.subintPeriods = nPer;
info.nChan        = nChan;
info.chanFreqs    = info_det.chanFreqs;
info.assign       = assign;
info.noiseCoeffs  = opts.NoiseCoeffs;               % [V X(1..Lmax)] of the time bins
info.covLags      = D;                              % phase-bin lags in weightX
info.noiseModel   = ['var(sum_j) = sigma^2*weight2, cov(sum_j, sum_j+d) = ' ...
                     'sigma^2*weightX(:,:,:,d), sigma^2 = m^2/(Bnoise*binDt)'];
info.dataWeights  = '';                             % data-weight file (DataWeights)
if haveW, info.dataWeights = DW.file; end
info.nWeightChan  = nW;                             % size of dimension 3 of the weights
info.validFraction = sumW / ((k2 - k1 + 1) * nW);   % mean W over the folded bins
info.binDt        = dt;                             % time-bin length of the input
info.binsFolded   = [k1, k2];
info.units        = 'power as in detectPower (mean per input sample)';
info.elapsed      = toc(tStart);

if opts.SaveFile
    outDir = fileparts(outFile);
    if ~isempty(outDir) && ~isfolder(outDir)
        mkdir(outDir);
    end
    save(outFile, 'fold', 'info');
end
if opts.Verbose
    fprintf('foldProfile: folded %d time bins into %d x %d x %d in %.2f s\n', ...
        k2 - k1 + 1, Nbin, nSub, nChan, info.elapsed);
    if haveW
        fprintf('foldProfile: valid fraction %.4f (mean W)\n', info.validFraction);
    end
end
end


% =========================================================================================
function [W2, WX, tail] = addLagPairs(W2, WX, tail, jA, wA, sub, Xk, perBin, Nbin)
%ADDLAGPAIRS  Add the covariance of time bins 1..Lmax apart to weight2 and weightX.
% jA [n x nA] phase bins (0-based) and wA [n x nA] assignment weights of the
% chunk's time bins, sub [n x 1] their sub-integration. Xk: the covariances,
% either constants [1 x Lmax] (perBin false) or per time bin and channel
% [n x nW x Lmax] (X of bin k with bin k+L). tail holds the last Lmax time
% bins of the previous chunk, so every pair (k, k+L) is counted once, in the
% chunk that contains k+L. Pairs in different sub-integrations are skipped
% (covariance between sub-int profiles is not stored).
[nCell, nW, D] = size(WX);
off = nCell * (0:nW-1);
nT = size(tail.jA, 1);
J  = [tail.jA; jA];
Wa = [tail.wA; wA];
Sb = [tail.sub; sub];
if perBin
    Xs = [tail.X; Xk];
    Lmax = size(Xk, 3);
else
    Lmax = numel(Xk);
end
m  = size(J, 1);
nA = size(J, 2);
for lag = 1:Lmax
    i1 = (max(1, nT - lag + 1) : m - lag).';
    i1 = i1(Sb(i1) == Sb(i1 + lag));                 % same sub-integration only
    i2 = i1 + lag;
    base = Nbin * (Sb(i1) - 1);
    if perBin, Xl = Xs(i1, :, lag); else, Xl = Xk(lag); end
    for p = 1:nA
        for q = 1:nA
            jp = J(i1, p);  jq = J(i2, q);
            c  = Xl .* Wa(i1, p) .* Wa(i2, q);       % [pairs x nW] (or x 1)
            d  = mod(jq - jp, Nbin);                 % phase-bin lag, circular
            s0 = d == 0;                             % same phase bin: (k, k+L) and (k+L, k)
            fw = d >= 1 & d <= D;                    % k+L in a later phase bin
            bw = d >= Nbin - D;                      % k+L in an earlier one (phase wrap / split)
            if ~all(s0 | fw | bw)
                error('foldProfile:lags', 'Correlated time bins reach beyond %d phase bins.', D);
            end
            W2 = addTo(W2, jp(s0) + 1 + base(s0) + off, 2*c(s0, :));
            cellX = [jp(fw) + 1 + base(fw) + nCell*nW*(d(fw) - 1); ...
                     jq(bw) + 1 + base(bw) + nCell*nW*(Nbin - d(bw) - 1)];
            WX = addTo(WX, cellX + off, [c(fw, :); c(bw, :)]);
        end
    end
end
keep = max(1, m - Lmax + 1) : m;
tail = struct('jA', J(keep, :), 'wA', Wa(keep, :), 'sub', Sb(keep), ...
              'X', zeros(0, nW, Lmax));
if perBin, tail.X = Xs(keep, :, :); end
end


% =========================================================================================
function A = addTo(A, idx, vals)
%ADDTO  A(idx) += vals, repeated indices summed (idx: linear indices into A;
% vals: same size as idx, or a scalar).
A(:) = A(:) + accumarray(idx(:), vals(:), [numel(A) 1]);
end