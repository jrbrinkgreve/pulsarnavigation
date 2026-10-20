function prof = rfiProfile(info_det, info_w, opts)
%RFIPROFILE  Per-channel statistics of the detected power with RFI and excision (C4b calibration).
%{
Measures, from a voltage-chain run WITHOUT pulsar (zero sky + receiver noise + RFI,
channel path with excision), what the RFI leaves in each channel's detected power after
excision, so simulateFold ('RFI') can reproduce it in folds of any length:

  level        L_c  mean power over the receiver-noise level (1 = clean; constant-envelope
                    and noise-like RFI that excision leaves in raise it)
  noiseCoeffs  V_c, X_c(1..Lmax): variance and covariances of the detected time bins,
               relative to the radiometer value mu_c^2/(BnoiseChan*binDt) with the channel's
               own mean mu_c - the same normalization as detectChannels' clean model
               (V 0.888, X(1) 0.046, ...). A channel partly covered by noise-like RFI has a
               narrower spectrum, so its V_c and X_c are larger than the model: the
               observer's error bars are too small there (B5e).
  validFraction mean data weight W (blankingWeights): 1 - the blanked fraction.

Time-domain quantities, so they transfer to any pulsar, phase-bin count or sub-int length.
Only the fully valid time bins (W >= MinValid) enter mu_c, V_c and X_c (a partly blanked bin
has mean W*mu and other noise; its weight loss is in validFraction); lag pairs need both
bins valid. The receiver-noise level is the median channel mean (most channels carry no
RFI), or 'RefLevel'. Time-variable RFI (a radar beam passage) is averaged over the run: for
a profile of the steady part keep beam passages outside the run (runRFIProfile).

  prof = rfiProfile(info_det, info_w)
  prof = rfiProfile(info_det, [], 'Label', 'noise only')      % no excision

Inputs:
  info_det  detectChannels' info of the run (power file [nChan x N] float32; noise =
            the clean model: V, X, Bnoise; binDt; fullySupportedBins)
  info_w    blankingWeights' info of the run (weights file [nChan x (2+Lmax) x N]: W, V,
            X per bin), or [] without excision (every bin valid)

Name-value options:
  'Label'      text stored in prof.label (default '')
  'RefLevel'   receiver-noise power per channel (default [] = median of the channel means)
  'MinValid'   W at or above which a time bin counts as fully valid (default 1 - 1e-6)
  'ChunkBins'  time bins read per chunk (default 2^18)
  'Verbose'    print a summary (default true)

Output prof: label, nChan, chanFreqs, chanWidth, binDt, BnoiseChan, Lmax, level [1 x nChan],
  noiseCoeffs [nChan x (1+Lmax)] (measured V_c X_c(L)), modelCoeffs [1 x (1+Lmax)] (clean
  model of the run), validFraction [1 x nChan], nValidBins [1 x nChan], refLevel, duration
  [s] (fully supported span), sourceFiles, made (date).
%}

arguments
    info_det struct
    info_w   = []
    opts.Label      {mustBeTextScalar} = ''
    opts.RefLevel   double = []
    opts.MinValid   (1,1) double {mustBeInRange(opts.MinValid, 0, 1)} = 1 - 1e-6
    opts.ChunkBins  (1,1) double {mustBeInteger, mustBePositive} = 2^18
    opts.Verbose    (1,1) logical = true
end

tStart = tic;
nChan = info_det.nChan;
Xm    = info_det.noise.X(:).';
Lmax  = numel(Xm);
Bn    = info_det.noise.Bnoise;
dt    = info_det.binDt;
haveW = ~isempty(info_w);
if haveW
    if info_w.nChan ~= nChan || info_w.N ~= info_det.N
        error('rfiProfile:weights', 'The weights file does not match the power file.');
    end
    nQ = 2 + info_w.Lmax;
end
k1 = info_det.fullySupportedBins(1);
k2 = info_det.fullySupportedBins(2);

[fid, msg] = fopen(info_det.file, 'r', info_det.byteOrder);
if fid == -1, error('rfiProfile:open', 'Could not open "%s": %s', info_det.file, msg); end
cP = onCleanup(@() fclose(fid));
fseek(fid, (k1 - 1) * nChan * 4, 'bof');
if haveW
    [fidW, msg] = fopen(info_w.file, 'r', info_w.byteOrder);
    if fidW == -1, error('rfiProfile:open', 'Could not open "%s": %s', info_w.file, msg); end
    cW = onCleanup(@() fclose(fidW));
    fseek(fidW, (k1 - 1) * nChan * nQ * 4, 'bof');
end

% ---- Pass 1: means over the valid bins, valid fraction ---------------------------------------
s1 = zeros(nChan, 1); n1 = zeros(nChan, 1); sW = zeros(nChan, 1);
forEachChunk(@pass1);
mu = s1 ./ n1;
validFraction = sW / (k2 - k1 + 1);

% ---- Pass 2: variance and lag covariances about mu (pairs of valid bins) ----------------------
sV = zeros(nChan, 1); sX = zeros(nChan, Lmax); nX = zeros(nChan, Lmax);
tailD = zeros(nChan, 0); tailV = false(nChan, 0);   % last Lmax bins of the previous chunk
fseek(fid, (k1 - 1) * nChan * 4, 'bof');
if haveW, fseek(fidW, (k1 - 1) * nChan * nQ * 4, 'bof'); end
forEachChunk(@pass2);

rad = mu.^2 / (Bn * dt);                            % radiometer variance per channel
V   = sV ./ n1 ./ rad;
X   = sX ./ max(nX, 1) ./ rad;
ref = opts.RefLevel;
if isempty(ref), ref = median(mu); end

prof = struct();
prof.label        = char(opts.Label);
prof.nChan        = nChan;
prof.chanFreqs    = reshape(info_det.chanFreqs, 1, []);
prof.chanWidth    = info_det.chanWidth;
prof.binDt        = dt;
prof.BnoiseChan   = Bn;
prof.Lmax         = Lmax;
prof.level        = (mu / ref).';
prof.noiseCoeffs  = [V, X];
prof.modelCoeffs  = [info_det.noise.V, Xm];
prof.validFraction = validFraction.';
prof.nValidBins   = n1.';
prof.refLevel     = ref;
prof.duration     = (k2 - k1 + 1) * dt;
prof.sourceFiles  = {info_det.file, ''};
if haveW, prof.sourceFiles{2} = info_w.file; end
prof.made         = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm'));

if opts.Verbose
    r = (prof.noiseCoeffs(:, 1) + 2*sum(prof.noiseCoeffs(:, 2:end), 2)) / ...
        (prof.modelCoeffs(1) + 2*sum(prof.modelCoeffs(2:end)));       % long-run noise / model
    fprintf(['rfiProfile: %d channels, %.3g s; level %.3f-%.3f (above 1.05: %d channels); ' ...
             'long-run noise / model %.3f-%.3f (above 1.02: %d); valid fraction min %.4f, mean %.5f; ' ...
             'model V %.4f, median measured V %.4f; %.1f s\n'], nChan, prof.duration, ...
        min(prof.level), max(prof.level), nnz(prof.level > 1.05), min(r), max(r), nnz(r > 1.02), ...
        min(validFraction), mean(validFraction), prof.modelCoeffs(1), median(V), toc(tStart));
end

% ===========================================================================================
    function forEachChunk(fn)
        kPos = k1;
        chunk = opts.ChunkBins;
        if haveW, chunk = min(chunk, max(1, floor(2^25 / (nChan * nQ)))); end   % memory
        while kPos <= k2
            n = min(chunk, k2 - kPos + 1);
            P = fread(fid, [nChan n], 'single=>double');
            if haveW
                Y = fread(fidW, [nChan*nQ n], 'single=>double');
                W = Y(1:nChan, :);                  % first nChan values of each bin: W
            else
                W = ones(nChan, n);
            end
            fn(P, W);
            kPos = kPos + n;
        end
    end

    function pass1(P, W)
        ok = W >= opts.MinValid;
        s1 = s1 + sum(P .* ok, 2);
        n1 = n1 + sum(ok, 2);
        sW = sW + sum(W, 2);
    end

    function pass2(P, W)
        ok = W >= opts.MinValid;
        D  = (P - mu) .* ok;
        sV = sV + sum(D.^2, 2);
        Dx = [tailD, D]; Vx = [tailV, ok];          % with the previous chunk's tail
        o  = size(tailD, 2);
        for L = 1:Lmax
            a = (1 + o - L) : (size(Dx, 2) - L);    % pairs (k, k+L) with k+L in this chunk
            a = a(a >= 1);
            sX(:, L) = sX(:, L) + sum(Dx(:, a) .* Dx(:, a + L), 2);
            nX(:, L) = nX(:, L) + sum(Vx(:, a) & Vx(:, a + L), 2);
        end
        keep  = max(1, size(Dx, 2) - Lmax + 1) : size(Dx, 2);
        tailD = Dx(:, keep); tailV = Vx(:, keep);
    end
end
