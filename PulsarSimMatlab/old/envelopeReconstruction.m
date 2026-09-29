function info = envelopeReconstruction(inFile, outFile, fs, f_out, blockLen, mode, filterOrder)
%{
Reduce complex IQ baseband data down to its real-valued pulse envelope
(or power), discarding carrier/phase content entirely, with an
optional further decimation to f_out. Fully streamed, memory bounded
by blockLen and filter length.

  info = envelopeReconstruction(inFile, outFile, fs, f_out, fLO, blockLen, mode, filterOrder)

Inputs:
  inFile      path to a complex baseband file: interleaved float32 I,Q
              pairs (I0,Q0,I1,Q1,...), e.g. from applyIQmodulation.
  outFile     path to write the real-valued envelope/power output to
              (single float32 per sample, NOT interleaved -- it's real
              now).
  fs          [Hz] sample rate of the complex data in inFile.
  f_out       [Hz] desired output sample rate. Must be <= fs (this
              function only decimates -- see note on fLO below for why
              going the other way isn't what you want here). If f_out
              doesn't divide fs by an integer, the nearest achievable
              rate is used and reported in info.actualFout.
  blockLen    (opt) complex samples read per block. Default 1e6.
  mode        (opt) 'power' (default) or 'magnitude'.
              'power'     -> y = I^2 + Q^2 (the standard "detection"
                             step in real pulsar backends -- this is
                             the physically meaningful quantity for
                             pulse search/folding, and matches how the
                             Gaussian-noise-pulse envelopes were
                             generated in the first place).
              'magnitude' -> y = sqrt(I^2 + Q^2), same shape, different
                             scaling/units.
  filterOrder (opt) anti-alias FIR filter order used only if decimating
              (f_out < fs). Default 512. Same order/history tradeoff
              as elsewhere in this pipeline (see prior discussion of
              applyIQmodulation's filterOrder).

Output:
  info        struct: .samplesRead (complex, from inFile),
              .samplesWritten (real, to outFile), .decimationFactor,
              .actualFout, .mode.

--- Why fLO is unused, and what actually replaces the template's
    "y = x .* exp{j*2*pi*fLO*t}" step ---

The pulse ENVELOPE is exactly the magnitude of the complex baseband,
and magnitude is invariant to any phase rotation: |x*exp(j*theta)| =
|x| for ANY theta, including the inverse-LO phase. So multiplying by
exp(+j*2*pi*fLO*t) changes nothing about the quantity you actually
want -- it's not a missing step, it's a no-op for this goal, and is
left out entirely.

It would also be numerically wrong if used: fs here is the DECIMATED
baseband rate (e.g. 200 MHz), while fLO sits up near the original RF
band (e.g. 1.4 GHz) -- far above fs's own Nyquist. A discrete-time
complex exponential exp(j*2*pi*fLO*n/fs) at that ratio does not
represent true fLO-Hz content anymore; it's already aliased to some
other apparent frequency once you're sampling at fs. Reconstructing
real content up near the original RF band would require upsampling
back to (at least) the original f_in first, which is a much heavier
operation and contradicts the stated goal of reducing data size.

--- What this function actually does (the corrected pipeline) ---

  read inFile (interleaved I/Q) -> cut into blocks -> for each block:
    1. DETECT: collapse each complex sample to a real power or
       magnitude value (this IS "without the frequency content" --
       carrier phase is discarded completely here).
    2. If f_out < fs (the genuinely missing step in the template):
       anti-alias FIR low-pass the real detected sequence via
       overlap-save (same carry-forward-the-tail pattern as
       applyDispersionStream/applyIQmodulation), then decimate.
    3. fwrite the real result.

  If f_out == fs, step 2 is skipped entirely -- just detect and write,
  block by block.
%}

if nargin < 6 || isempty(blockLen), blockLen = 1e6; end
if nargin < 7 || isempty(mode), mode = 'power'; end
if nargin < 8 || isempty(filterOrder), filterOrder = 512; end

if f_out > fs
    error(['envelopeReconstruction: f_out (%.6g Hz) > fs (%.6g Hz). ', ...
        'This function only decimates the envelope/power trace; it does ', ...
        'not upsample or reconstruct RF content (see fLO note in the ', ...
        'docstring for why that would require a different, heavier operation).'], ...
        f_out, fs);
end

switch lower(mode)
    case 'power'
        detectFcn = @(iq) real(iq).^2 + imag(iq).^2;
    case 'magnitude'
        detectFcn = @(iq) sqrt(real(iq).^2 + imag(iq).^2);
    otherwise
        error('envelopeReconstruction: unknown mode "%s" (use ''power'' or ''magnitude'').', mode);
end

% --- Decimation factor for the envelope/power trace ---
D = round(fs / f_out);
D = max(D, 1);
actualFout = fs / D;
if abs(actualFout - f_out) > 1e-6*f_out
    warning(['envelopeReconstruction: f_out (%.6g Hz) does not divide fs ', ...
        '(%.6g Hz) by an integer factor. Using nearest achievable rate ', ...
        '%.6g Hz (decimation factor %d) instead.'], f_out, fs, actualFout, D);
end
needsFilter = (D > 1);

% --- Determine total complex sample count from file size ---
fInfo = dir(inFile);
if isempty(fInfo)
    error('envelopeReconstruction: inFile "%s" not found.', inFile);
end
if mod(fInfo.bytes, 8) ~= 0
    warning('envelopeReconstruction: inFile size is not a multiple of 8 bytes (interleaved I/Q pairs); truncating trailing partial sample.');
end
Nx = floor(fInfo.bytes / 8); % 8 bytes per complex sample (2 x float32)

fidIn = fopen(inFile, 'r');
if fidIn == -1
    error('envelopeReconstruction: could not open inFile "%s".', inFile);
end
fidOut = fopen(outFile, 'w');
if fidOut == -1
    fclose(fidIn);
    error('envelopeReconstruction: could not open outFile "%s" for writing.', outFile);
end
cleanupIn  = onCleanup(@() fclose(fidIn));  %#ok<NASGU>
cleanupOut = onCleanup(@() fclose(fidOut)); %#ok<NASGU>

if needsFilter
    marginFactor = 0.9;
    Wn = marginFactor * (actualFout/2) / (fs/2); % normalized to fs's Nyquist
    Wn = min(max(Wn, 1e-6), 0.999);
    h_lpf = single(fir1(filterOrder, Wn));
    M = numel(h_lpf);
    if blockLen <= 4*M
        warning(['envelopeReconstruction: blockLen (%d) is not much larger ', ...
            'than the filter length (%d). Overlap-save is still correct, ', ...
            'but efficiency drops.'], blockLen, M);
    end
    Nfft = 2^nextpow2(blockLen + M - 1);
    H = fft(h_lpf, Nfft);
    history = zeros(1, M-1, 'single'); % "to store for next round" (real now)
end

startIdx = 1;
samplesRead = 0;
samplesWritten = 0;
blockNum = 0;

while startIdx <= Nx
    nWant = min(blockLen, Nx - startIdx + 1);
    raw = fread(fidIn, 2*nWant, 'single=>single').'; % interleaved I,Q
    nPairs = floor(numel(raw)/2);
    if nPairs == 0
        break;
    end
    if mod(numel(raw), 2) ~= 0
        warning('envelopeReconstruction: dropped a trailing unpaired float at end of file.');
    end
    endIdx = startIdx + nPairs - 1;

    iqBlock = complex(raw(1:2:2*nPairs-1), raw(2:2:2*nPairs));
    detected = single(detectFcn(iqBlock)); % real, length nPairs

    if ~needsFilter
        % No decimation requested -- write the detected trace directly.
        fwrite(fidOut, detected, 'single');
        samplesWritten = samplesWritten + numel(detected);
    else
        block = [history, detected]; % length M-1 + nPairs
        if numel(block) <= Nfft
            X = fft(block, Nfft);
            convOut = real(ifft(X .* H));
        else
            NfftLocal = 2^nextpow2(numel(block));
            Hlocal = fft(h_lpf, NfftLocal);
            X = fft(block, NfftLocal);
            convOut = real(ifft(X .* Hlocal));
        end
        validOut = convOut(M:M+nPairs-1); % discard first M-1 corrupted (causal, no centering needed for a smoothing/anti-alias filter)
        history = block(end-(M-2):end);   % carry last M-1 detected samples forward

        globalIdx0 = startIdx - 1; % 0-based
        keepMask = mod(globalIdx0 + (0:nPairs-1), D) == 0;
        decimated = single(validOut(keepMask));

        if ~isempty(decimated)
            fwrite(fidOut, decimated, 'single');
            samplesWritten = samplesWritten + numel(decimated);
        end
    end

    samplesRead = samplesRead + nPairs;
    startIdx = endIdx + 1;
    blockNum = blockNum + 1;
end

info = struct( ...
    'samplesRead', samplesRead, ...
    'samplesWritten', samplesWritten, ...
    'decimationFactor', D, ...
    'actualFout', actualFout, ...
    'mode', mode);

fprintf(['envelopeReconstruction: read %d complex samples @ %.4g Hz, wrote ', ...
    '%d real (%s) samples @ %.4g Hz (D=%d, %d blocks) to %s\n'], ...
    samplesRead, fs, samplesWritten, mode, actualFout, D, blockNum, outFile);

end