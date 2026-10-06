function nc = powerCovariance(chanWidth, edgeWidth, fs, nPerBin, opts)
%POWERCOVARIANCE  Exact noise statistics of detected time bins in a narrow channel.
%{
After square-law detection, a time bin is the mean of |y|^2 over nPerBin
samples. For Gaussian y (receiver noise, and the pulsar's own noise-like
signal) with power spectrum S(f), relative to the radiometer value
m^2/(B*dt) (B = (int S)^2/int S^2, dt = nPerBin/fs):

  variance of a bin          V   = C(0) / rad
  covariance with bin k + L  X_L = C(L) / rad,   L = 1..Lmax
  C(L) = (1/n^2) * sum_{a,b = 0..n-1} |R(a - b + L*n)|^2,  rad = sum_l |R(l)|^2 / n

with R(l) the normalized autocorrelation of y at lag l samples (the inverse
Fourier transform of S, R(0) = 1). Summed over all lags the bin covariances
give exactly the radiometer value: V + 2*sum_L X_L = 1. That is why the
simple model "independent bins, var = m^2/(B*dt)" is right for sums over
many bins (long phase bins) but not bin by bin when B*dt is small: in a
3.125 MHz channel with 0.96 us bins (B*dt ~ 3) V = 0.888 and X_1 = 0.046,
whereas in the 400 MHz full band (B*dt ~ 400) V = 1 and X ~ 0.

Spectrum: S = W^2 with W the sin^2 band-edge taper of width edgeWidth inside
+-chanWidth/2 (the per-channel dedispersion taper of dedisperseChannels; the
channelizer prototype is flat to 1e-4 there). The chirp of the dedispersion
filter is an all-pass phase and drops out of R; it matters only when parts
of the input are blanked (then the statistics vary from bin to bin and come
from the mask, see tests/expBlankingVariance.m).

  nc = powerCovariance(chanWidth, edgeWidth, fs, nPerBin)
  nc = powerCovariance(..., 'Tolerance', 0.002, 'MaxLag', 256)

Inputs:
  chanWidth  [Hz] useful channel width (ChanWidth of channelizeIQ).
  edgeWidth  [Hz] taper width at each channel edge (info_dc.edgeWidth).
  fs         [Hz] channel sample rate.
  nPerBin    samples per detected time bin.

Name-value options:
  'Tolerance'  Lmax = smallest L with V + 2*sum_{1..L} X >= 1 - Tolerance
               (default 0.002: the lags kept carry 99.8 % of the variance of
               a long sum).
  'MaxLag'     upper limit for Lmax (default min(256, what the frequency grid
               resolves without wrap-around)).
  'NFreq'      frequency grid points for R (default 2^18).

Output nc: V, X (1 x Lmax), Lmax, captured (= V + 2*sum(X)), Bnoise [Hz],
  binDt [s], radiometer (= 1/(Bnoise*binDt)), nPerBin, description.
%}

arguments
    chanWidth  (1,1) double {mustBePositive}
    edgeWidth  (1,1) double {mustBePositive}
    fs         (1,1) double {mustBePositive}
    nPerBin    (1,1) double {mustBeInteger, mustBePositive}
    opts.Tolerance (1,1) double {mustBePositive} = 0.002
    opts.MaxLag          double = []
    opts.NFreq     (1,1) double {mustBeInteger, mustBePositive} = 2^18
end

if chanWidth > fs
    error('powerCovariance:band', 'chanWidth must not exceed the sample rate.');
end
n  = nPerBin;
Nf = opts.NFreq;
maxLag = opts.MaxLag;
if isempty(maxLag), maxLag = min(256, floor(Nf / (2*n)) - 2); end
if maxLag < 1 || (maxLag + 1) * n >= Nf/2
    error('powerCovariance:grid', 'NFreq too small for MaxLag (lags wrap around).');
end

% Spectrum on the sampled-frequency grid and its autocorrelation
f  = ((0:Nf-1).' - Nf/2) * fs / Nf;
in = abs(f) <= chanWidth/2;
W  = zeros(Nf, 1);
W(in) = 1;
lo = in & f < -chanWidth/2 + edgeWidth;  W(lo) = sin(pi/2 * (f(lo) + chanWidth/2) / edgeWidth).^2;
hi = in & f >  chanWidth/2 - edgeWidth;  W(hi) = sin(pi/2 * (chanWidth/2 - f(hi)) / edgeWidth).^2;
S  = W.^2;
R  = fft(ifftshift(S));                  % R(l+1) ~ autocorrelation at lag l (circular)
R2 = abs(R / R(1)).^2;

rad = sum(R2) / n;                       % radiometer value in units of m^2 (= fs/(B*n))
d   = (-(n-1):(n-1)).';
cnt = n - abs(d);                        % number of sample pairs (a, b) with a - b = d
C   = zeros(maxLag + 1, 1);
for L = 0:maxLag
    C(L+1) = sum(cnt .* R2(mod(d + L*n, Nf) + 1)) / n^2;
end
V = C(1) / rad;
X = C(2:end).' / rad;
cum  = V + 2*cumsum(X);
Lmax = find(cum >= 1 - opts.Tolerance, 1);
if isempty(Lmax)
    warning('powerCovariance:lags', ...
        'Lags up to %d carry only %.4f of the variance; raise MaxLag.', maxLag, cum(end));
    Lmax = maxLag;
end

nc = struct();
nc.V          = V;
nc.X          = X(1:Lmax);
nc.Lmax       = Lmax;
nc.captured   = V + 2*sum(nc.X);
nc.Bnoise     = fs / sum(R2);
nc.binDt      = n / fs;
nc.radiometer = rad;
nc.nPerBin    = n;
nc.chanWidth  = chanWidth;
nc.edgeWidth  = edgeWidth;
nc.description = ['time-bin power noise relative to m^2/(Bnoise*binDt): ' ...
                  'var = V, cov(k, k+L) = X(L); V + 2*sum(X) = captured (-> 1)'];
end
