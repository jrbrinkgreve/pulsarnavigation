function y = applyDispersion(x, DM, f_in, fLow, fHigh, blockLen)
%{
Applying interstellar dispersion to a raw generated pulsar signal, via
time-block-based overlap-add convolution with edge handling.

  y = applyDispersion(x, DM, f_in, fLow, fHigh, blockLen)

Inputs:
  x        raw pulsar signal (real, e.g. from generatePulsarSignal),
           row or column vector.
  DM       [pc/cm^3] dispersion measure.
  f_in       [Hz] sample rate of x (must match how x was generated).
  fLow     [Hz] (opt) low edge of the signal band. Default 1.2e9,
           matching generatePulsarSignal's assumed RF band -- override
           if you generated at a different frequency.
  fHigh    [Hz] (opt) high edge of the signal band. Default 1.6e9.
  blockLen (opt) number of samples of x processed per block. Default
           1e6. Must be well above the dispersion kernel length for the
           overlap-add bookkeeping below to make sense (a warning is
           printed if it isn't).

Output:
  y        dispersed signal, SAME LENGTH as x, with the Gaussian pulses
           still centered at the same nominal times as in x (the
           kernel's own group delay is removed by construction -- see
           "centering" below).

--- Method (overlap-add, matching the intended structure) ---

The dispersion impulse response h(t) has support wider than a single
processing block (that's the whole reason edge handling is needed).
Convolving block-by-block naively drops the part of each block's
response that spills into the next block:

  ---------t0       -      t0+dt
  ----------|----------------|---------
  ---------/ \              / \
  --------/   \      H     /   \
  -------/     \          /     \
  ------|-------|--------|-------|-----
            ^    to write    ^
            |      now       |
          to sum / to store for next round

So each block's LOCAL full convolution (length nBlock+M-1, M = length
of h) is split in two: the first nBlock samples are this block's final,
writeable output (after adding in whatever tail was carried over FROM
the previous block); the last M-1 samples are not yet final -- they
overlap the next block's response and must be carried forward and added
in before that portion can be written. That carried vector is exactly
the "to store for next round" piece in the diagram.

Centering: h is built symmetric (see makeDispersionKernel, which
fftshift-centers it) and forced to odd length M here, so halfM =
(M-1)/2 is its exact group delay in samples on each side. The
overlap-add loop tracks a running "global full-convolution index" and
writes each finished sample into y at (globalIdx - halfM) -- i.e. the
block algorithm computes the same thing conv(x,h,'same') would, just
without ever materializing the full N+M-1 length array or a single
huge FFT.

Edges:
  - True start of x: initial carry is zero (no data exists before
    x(1); this is the correct physical edge, not an approximation).
  - True end of x: after the last block, the remaining carry (length
    M-1) still contains valid, not-yet-written samples -- these are
    flushed into y explicitly after the loop, since (thanks to the
    halfM centering shift) roughly the first half of that leftover tail
    still lands inside y's valid range.
%}

if nargin < 4 || isempty(fLow),  fLow  = 1.2e9; end
if nargin < 5 || isempty(fHigh), fHigh = 1.6e9; end
if nargin < 6 || isempty(blockLen), blockLen = 1e6; end

x = single(x(:).'); % row vector, single precision
Nx = numel(x);

% --- Build the dispersion kernel h(t) ---
h = makeDispersionKernel(DM, fLow, fHigh, f_in, 'hann'); % real, single, power-of-2 length
h = h(:).';
disp("Constructed dispersion kernel")
if mod(numel(h), 2) == 0
    h = h(1:end-1); % force odd length for exact symmetric centering
end
M = numel(h);
halfM = (M-1)/2;

if blockLen <= 4*M
    warning(['applyDispersion: blockLen (%d) is not much larger than the ', ...
        'kernel length (%d). Overlap-add is still correct, but efficiency ', ...
        'drops sharply and the "block-wise" edge picture stops being ', ...
        'meaningful. Consider a larger blockLen or a smaller-DM/narrower-', ...
        'band kernel.'], blockLen, M);
end

y = zeros(1, Nx, 'single');
carry = zeros(1, M-1, 'single'); % "to store for next round" (see diagram)

Nfft = 2^nextpow2(blockLen + M - 1); % fixed FFT size, reused every block
H = fft(h, Nfft);

for startIdx = 1:blockLen:Nx
    endIdx = min(startIdx + blockLen - 1, Nx);
    nBlock = endIdx - startIdx + 1;

    xBlock = x(startIdx:endIdx);

    if nBlock + M - 1 <= Nfft
        Xb = fft(xBlock, Nfft);
        convFull = ifft(Xb .* H);
        convFull = real(convFull(1:nBlock + M - 1)); % local full conv, this block only
    else
        % last block bigger than planned (shouldn't happen with fixed
        % blockLen, but guard against it) -- fall back to a
        % block-specific FFT size
        NfftLocal = 2^nextpow2(nBlock + M - 1);
        Xb = fft(xBlock, NfftLocal);
        Hlocal = fft(h, NfftLocal);
        convFull = real(ifft(Xb .* Hlocal));
        convFull = convFull(1:nBlock + M - 1);
    end

    % Add in the carry from the previous block (the "to sum" step)
    convFull(1:M-1) = convFull(1:M-1) + carry;

    outputNow = convFull(1:nBlock);          % finalized for this block
    carry = convFull(nBlock+1:end);          % "to store for next round"

    % Map local block (global full-conv index = startIdx..endIdx) to y,
    % shifted by halfM to remove the kernel's group delay (centering).
    fullIdxRange = startIdx:endIdx;
    yIdxRange = fullIdxRange - halfM;
    valid = yIdxRange >= 1 & yIdxRange <= Nx;
    y(yIdxRange(valid)) = outputNow(valid);
end

% Flush the final leftover carry (true end-of-signal edge): these
% samples are final (nothing more will ever be added to them) and,
% after the halfM shift, some of them still fall inside y.
finalFullIdx = Nx+1 : Nx+M-1;
finalYIdx = finalFullIdx - halfM;
validFinal = finalYIdx >= 1 & finalYIdx <= Nx;
y(finalYIdx(validFinal)) = carry(validFinal);
y = y';

end











%helper function
%================================================================================

function [kernel, info] = makeDispersionKernel(DM, fLow, fHigh, f_in, taper)
%MAKEDISPERSIONKERNEL Build a real FIR kernel that applies interstellar-
% medium (ISM) dispersion to a real-valued RF-sampled signal, for use
% with chirpFilterStream (from the earlier matched-filter work).
%
%   [kernel, info] = makeDispersionKernel(DM, fLow, fHigh, f_in, taper)
%
% Inputs:
%   DM      [pc/cm^3] dispersion measure
%   fLow    [Hz]      lowest frequency of the signal band (e.g. 1.2e9)
%   fHigh   [Hz]      highest frequency of the signal band (e.g. 1.6e9)
%   f_in      [Hz]      sample rate of the raw RF data (e.g. 4e9)
%   taper   (opt)     'hann' (default) or 'none' -- windows the kernel
%                      in time to suppress Gibbs ringing at its edges.
%
% Output:
%   kernel  real single-precision FIR kernel. Feed this straight into
%           chirpFilterStream(inFile, outFile, f_in, kernel, blockSize)
%           to disperse a real-valued raw RF file exactly like the
%           chirp-matched-filter case, just with a different kernel.
%   info    struct with .smearTime (s), .kernelLen, .fRefHz, etc., for
%           sanity-checking before you commit to a run.
%
% Convention: reference frequency is fHigh (top of band arrives with
% zero extra delay; lower frequencies lag behind it), matching how real
% dispersed pulses look on a filterbank. To DEdisperse instead of
% disperse, conjugate H below (flip the sign in the exponent) when
% building the kernel.
%
% IMPORTANT -- kernel size scales with DM and bandwidth:
%   Kernel length (in samples) is roughly f_in * dispersion smear time
%   across the band. For large DM and/or wide, multi-GHz bandwidth
%   sampled at multi-GHz rates, this can run into hundreds of millions
%   of samples -- a real cost, not a bug. This is exactly why real
%   pulsar back-ends channelize the band into many narrow sub-bands and
%   coherently disperse/dedisperse each channel separately (each
%   channel's smear, and hence its kernel, shrinks roughly with the
%   square of the channel bandwidth). If info.kernelLen comes back huge,
%   channelize first rather than trying to build one giant kernel.

if nargin < 5 || isempty(taper)
    taper = 'hann';
end

D = 4148.808; % s * MHz^2 * pc^-1 * cm^3

fLow_MHz  = fLow  / 1e6;
fHigh_MHz = fHigh / 1e6;
fRef_MHz  = fHigh_MHz;          % reference: top of band, zero delay there

tau = @(f_MHz) D * DM ./ (f_MHz.^2);   % delay [s] at absolute RF freq (MHz)
tauRef = tau(fRef_MHz);

smearTime = tau(fLow_MHz) - tauRef;    % total delay spread across the band [s]

marginFactor = 1.5; % headroom so the kernel isn't truncated right at the edge
Nfft = 2^nextpow2(ceil(marginFactor * smearTime * f_in));
Nfft = max(Nfft, 64); % floor, for degenerate tiny-DM cases

% Build a Hermitian-symmetric H(f) over the full [0, f_in) bin grid so
% ifft(H) comes out real.
binFreqsHz = (0:Nfft-1) * (f_in / Nfft);
H = ones(1, Nfft); % default: no phase shift (used at/near DC, and Nyquist)

half = floor(Nfft/2);
for k = 2:half              % skip DC (k=1); handle positive freqs
    f_Hz = binFreqsHz(k);
    f_MHz = f_Hz / 1e6;
    relDelay = tau(f_MHz) - tauRef;         % seconds
    H(k) = exp(-1i * 2*pi * f_Hz * relDelay);
end
% Mirror to negative-frequency bins for conjugate (Hermitian) symmetry
for k = 2:half
    mirrorIdx = Nfft - k + 2;
    if mirrorIdx <= Nfft && mirrorIdx ~= k
        H(mirrorIdx) = conj(H(k));
    end
end
% Nyquist bin (if Nfft even) must be real
if mod(Nfft,2) == 0
    H(half+1) = real(H(half+1));
end

kernelFull = real(ifft(H));      % Hermitian H -> real kernel (numerically)
kernelFull = fftshift(kernelFull); % center the delay structure in time

switch lower(taper)
    case 'hann'
        w = 0.5 - 0.5*cos(2*pi*(0:Nfft-1)/(Nfft-1));
    case 'none'
        w = ones(1, Nfft);
    otherwise
        error('Unknown taper "%s".', taper);
end

kernel = single(kernelFull .* w);
kernel = kernel / sum(abs(kernel)); % keep the filter roughly unit-gain

info = struct( ...
    'smearTime', smearTime, ...
    'kernelLen', Nfft, ...
    'fRefHz', fHigh, ...
    'fLowHz', fLow, ...
    'fHighHz', fHigh, ...
    'DM', DM);

fprintf(['makeDispersionKernel: DM=%.2f pc/cm^3, band %.3f-%.3f GHz, ', ...
    'smear=%.3f ms, kernel length=%d samples (%.1f MB, single)\n'], ...
    DM, fLow/1e9, fHigh/1e9, smearTime*1e3, Nfft, Nfft*4/1e6);

end