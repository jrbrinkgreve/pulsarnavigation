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
  'NBin'           phase bins (default 1024).
  'SubintPeriods'  turns per sub-integration (default 1).
  'SubintTime'     [s] alternative: sub-integration length, rounded to
                   whole turns (overrides SubintPeriods).
  'Assign'         'linear' (default): each time bin is split between the
                   two nearest phase bins (unbiased at sub-bin level);
                   'nearest': whole time bin to the nearest phase bin
                   (half the work; bias up to ~binDt/2 depending on how the
                   time grid falls on the phase grid).
  'UseSupportedOnly' fold only info_det.fullySupportedBins (default true).
  'ChunkBins'      time bins read per chunk (default 2^22).
  'MaxMemoryGB'    limit for the fold arrays (default 2).
  'SaveFile'       save fold and info to outFile (default true).
  'Verbose'        print a summary (default true).

Outputs:
  fold  struct:
    .sum     [NBin x nSub x nChan] sum of w * power
    .weight  [NBin x nSub]         sum of w (same for every channel)
    .weight2 [NBin x nSub]         sum of w^2 (for noise: the variance of a
                                   phase bin is sigma_timebin^2 * weight2 / weight^2)
    .weightX [NBin x nSub]         sum of w_j*w_(j+1) over time bins split between
                                   bin j and bin j+1 (circular); covariance of
                                   neighbouring phase bins is
                                   sigma_j*sigma_j+1 * weightX_j / (weight_j*weight_j+1).
                                   Zero for 'nearest'. Needed for correct
                                   uncertainties of anything that sums over bins.
    .prof    [NBin x nSub x nChan] sum ./ weight (NaN where weight = 0)
    .profTotal [NBin x nChan]      all sub-integrations combined
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
Nbin  = opts.NBin;
nChan = info_det.nChan;
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

memGB = Nbin * nSub * (2*nChan + 2) * 8 / 1e9;
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
end

% ---- Accumulate ------------------------------------------------------------------------
nCell   = Nbin * nSub;
S       = zeros(nCell, nChan);
Wt      = zeros(nCell, 1);
W2      = zeros(nCell, 1);
WX      = zeros(nCell, 1);
tSum    = zeros(nSub, 1);
nTB     = zeros(nSub, 1);

[fid, msg] = fopen(info_det.file, 'r', info_det.byteOrder);
if fid == -1
    error('foldProfile:open', 'Could not open "%s": %s', info_det.file, msg);
end
cleanupIn = onCleanup(@() fclose(fid)); %#ok<NASGU>
if fseek(fid, (k1 - 1) * nChan * 4, 'bof') ~= 0
    error('foldProfile:seek', 'Could not seek in "%s".', info_det.file);
end

kPos = k1;
while kPos <= k2
    n = min(opts.ChunkBins, k2 - kPos + 1);
    X = fread(fid, [nChan n], 'single=>double');
    if size(X, 2) < n
        error('foldProfile:read', 'Unexpected end of "%s".', info_det.file);
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
        Wt = Wt + accumarray(L, 1, [nCell 1]);
        W2 = W2 + accumarray(L, 1, [nCell 1]);
        for c = 1:nChan
            S(:, c) = S(:, c) + accumarray(L, X(c, :).', [nCell 1]);
        end
    else
        j0 = floor(x);
        a  = x - j0;
        L  = [mod(j0, Nbin) + 1 + base; mod(j0 + 1, Nbin) + 1 + base];
        w  = [1 - a; a];
        Wt = Wt + accumarray(L, w,    [nCell 1]);
        W2 = W2 + accumarray(L, w.^2, [nCell 1]);
        WX = WX + accumarray(mod(j0, Nbin) + 1 + base, (1 - a) .* a, [nCell 1]);
        for c = 1:nChan
            v = X(c, :).';
            S(:, c) = S(:, c) + accumarray(L, w .* [v; v], [nCell 1]);
        end
    end
    kPos = kPos + n;
end

% ---- Assemble fold struct ----------------------------------------------------------------
fold = struct();
fold.sum     = reshape(S,  Nbin, nSub, nChan);
fold.weight  = reshape(Wt, Nbin, nSub);
fold.weight2 = reshape(W2, Nbin, nSub);
fold.weightX = reshape(WX, Nbin, nSub);
wFull = repmat(fold.weight, 1, 1, nChan);
fold.prof = fold.sum ./ wFull;
fold.prof(wFull == 0) = NaN;
wTot = sum(fold.weight, 2);
fold.profTotal = reshape(sum(fold.sum, 2), Nbin, nChan) ./ wTot;
fold.profTotal(wTot == 0, :) = NaN;
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
end
end