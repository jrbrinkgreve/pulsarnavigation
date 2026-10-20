function [info, fold, truth] = simulateFold(f0, rhoDB, nTurns, opts)
%SIMULATEFOLD  Folded profiles drawn directly from the fold's noise model (fast simulator, block C).
%{
Hours of observation cannot be simulated at voltage level (~61 GB per second).
This draws the fold that foldProfile would make (linear assignment) straight
from its statistics, so foldToTOA / estimateTOA / detectPulsar run on it
unchanged. C1 (20 Oct 2026): full-band equivalent - one channel, independent
time bins, constant flux, the ephemeris equal to the truth except for a given
phase offset, continuous observation. C4a (20 Oct): 'NChan' channels (the channel
path: per-channel folds, correlated time bins via 'NoiseCoeffs', the true pulse
fraction per channel 'ChannelGain'); NChan = 1 (default) is C1, bit for bit.
C4b (20 Oct): 'RFI' - a per-channel profile of what RFI leaves after excision
(rfiProfile, measured on a voltage-chain run without pulsar): raised level, true
noise coefficients (the observer keeps the clean model), blanked fraction.
C3a (20 Oct): an observation schedule - 'Passes' (start/end of each pass, gaps
between them) and 'SubintTime' in seconds - with absolute times. C3b (20 Oct):
interstellar scintillation, a flux factor per sub-int and channel ('FluxScale' given,
or a random field from 'Scintillation': scintle bandwidth and time). 'NPol' (20 Oct):
dual polarization (total intensity). Truth side, like the generator: the
observer gets only fold and info (which hold nothing foldProfile's would not).

  [info, fold, truth] = simulateFold(f0, rhoDB, nTurns)
  [info, fold, truth] = simulateFold(100, -20, 72000, 'SubintPeriods', 6000, 'Seed', 3)
  [info, fold, truth] = simulateFold(f0, rhoDB, [], 'Passes', [0 4; 24 28]*3600, 'SubintTime', 60)

Example, J0437-4715 (scenario page, ATNF): f0 = 173.6879 Hz (P 5.757452 ms), power-
profile FWHM 0.141 ms = 0.0245 turns, NBin 512 (11.2 us bins = FWHM/12.5, enough by
estimateTOA's resolution rule; 20 h x 128 channels in 60 s sub-ints ~1.3 GB of folds),
rhoDB -51.0 (worst case) or -47.4 (nominal) in this convention (deep analysis 13 Oct §8;
the scenario page's -53.7 / -50.1 use W_eq = 2*W50). The scenario assumes two
polarizations: 'NPol', 2, and the observer passes Bnoise x 2 to foldToTOA.

Inputs:
  f0       [Hz] spin frequency (true and ephemeris)
  rhoDB    [dB] pulse-peak signal PSD / receiver-noise PSD in the band: the
           simulator's snrDB convention (peak of the Gaussian pulse; not the
           scenario page's W_eq = 2*W50 rule, which is 2.7 dB lower for the
           same pulsar). -Inf: no pulsar (H0).
  nTurns   turns observed from t = TRef - 0.5/f0 (turn 0); nSub = floor(nTurns /
           SubintPeriods) complete sub-ints. [] with 'Passes'.

Name-value options:
  'FWHM'          [turns] FWHM of the Gaussian power profile (default 0.05)
  'Components'    [nComp x 3] alternative to FWHM (20 Oct 2026): the power profile as a
                  sum of Gaussians, rows [amplitude, centre (turns), FWHM (turns)], scaled
                  to peak 1 (so rhoDB stays the pulse peak over the noise); phase 0 is the
                  components' reference point (the pulse sits at PhaseOffset; a perfect
                  TOA fit with gaussianTemplate(NBin, components) returns it). E.g.
                  B0329+54 at 1410 MHz: five components (arXiv:1703.03582, Table 2).
                  One row [1 0 FWHM] = 'FWHM' (to 1e-12)
  'PhaseOffset'   [turns] true pulse phase against the ephemeris (the pulse
                  peak is at phase PhaseOffset; a perfect TOA fit returns it):
                  scalar or one value per sub-int (default 0)
  'NBin'          phase bins (default 2048)
  'SubintPeriods' turns per sub-int (default 1)
  'SubintTime'    [s] alternative: sub-int length, rounded to whole turns (as foldProfile)
  'Passes'        [nPass x 2] start and end [s] of each observing pass (default []: the
                  nTurns turns from turn 0). Only complete turns inside a pass, grouped
                  into complete sub-ints of SubintPeriods turns per pass; gaps between
                  passes hold no data. Sub-int times are absolute (tRef = TRef +
                  turnRef/f0); truth.pass gives each sub-int's pass.
  'FluxScale'     true flux factor F of the pulsar (truth): [nSub x 1] or [nSub x NChan],
                  multiplies the pulse in each sub-int (and channel); default none (1)
  'Scintillation' struct, alternative to FluxScale: F = |g|^2 from a complex Gaussian
                  field g(t, nu) with mean |g|^2 = 1 (exponential statistics: fully
                  modulated scintles, no refractive variation), Gaussian correlations
                  with the intensity ACF falling to 1/e at the scintle time and to 1/2 at
                  the scintle bandwidth. Fields: bandwidth [Hz] (HWHM in frequency),
                  timescale [s] (1/e half-width in time), band [fLow fHigh] (default
                  [1.2e9 1.6e9]; NChan equal channels across it; one channel = the band
                  average). Evaluated at each sub-int's mid time, so sub-ints should be
                  short against the scintle time (warning above timescale/5). Own random
                  stream (the noise draws are unchanged). J0437 at 1.4 GHz (IAR 2020,
                  rescaled from Keith et al. 2013): bandwidth ~740 MHz, timescale ~2290 s,
                  varying notably between epochs.
  'Bnoise'        [Hz] noise-equivalent bandwidth of the detected power of one
                  channel: scalar or one per channel (default 391.6e6, the full
                  band with its dedispersion taper; channel path 3.0596e6, the
                  channel's detectChannels info.noise.Bnoise); Inf: no noise
                  (the mean profile alone)
  'NPol'          polarizations added into total intensity (default 1; 20 Oct 2026):
                  the fold is their average - baseline and pulse per polarization as
                  rhoDB says (an unpolarized pulsar: rho = S_peak/SEFD in each), the
                  noise variance divided by NPol (independent receiver noise in each;
                  with 'RFI' the fluctuations RFI leaves are taken as independent too).
                  The same as Bnoise x NPol, which is what the observer passes to
                  foldToTOA / estimateTOA / detectPulsar (fold and weights unchanged);
                  unlike Bnoise it also works with 'RFI', whose Bnoise is fixed by the
                  calibration run. NPol 1 = the single stream of the voltage chain.
  'BinDt'         [s] time-bin length of the detected power (default 1e-6; channel
                  path 0.96e-6)
  'NChan'         channels (default 1: full-band equivalent)
  'ChannelGain'   [1 x NChan] true fraction of the pulse in each channel (truth:
                  e.g. the forward band-edge taper in the edge channels; default 1)
  'NoiseCoeffs'   [V X(1) ... X(Lmax)] of the detected time bins, as foldProfile's
                  (default 1: independent; channel path detectChannels' info.noise
                  [V X], e.g. 0.8877 0.0458 0.0051 ...)
  'MaxMemoryGB'   limit for the fold arrays (default 2; NBin*nSub*(2*NChan + nW*(2 +
                  lags))*8 bytes, nW = 1, or NChan with 'RFI')
  'RFI'           [] (default) or a profile from rfiProfile (runRFIProfile): per
                  channel level L_c, measured noise coefficients and valid fraction
                  v_c. NChan must match; BinDt, Bnoise and NoiseCoeffs are then the
                  profile's (binDt, BnoiseChan, modelCoeffs: the observer's clean
                  model) unless given equal. Asymptotic weights only.
  'ObserverNoise' with 'RFI' (20 Oct 2026): 'model' (default) - the fold's weight2 /
                  weightX come from the clean model coefficients, as an observer that
                  assumes flat channels has them; 'measured' - from the profile's measured
                  coefficients per channel, as an observer that calibrates each channel's
                  noise on its own data has them (per-channel noise calibration, block B
                  item; the chain's counterpart is foldProfile 'ChannelNoise'). The noise
                  itself is drawn with the measured coefficients either way.
  'SignalGain'    pulse power factor g: peak power above baseline = g*rho
                  (default 1; the voltage chain's band tapers give slightly
                  less - pass expectedPowerModel's sigScale/(Pn*rho) to compare)
  'TRef'          [s] ephemeris reference time of phase 0 (default 0.5/f0, as
                  pipelineParams: phase 0 at the first pulse centre)
  'Weights'       [] (default): asymptotic weights, below. Or a fold struct
                  (foldProfile) with shared weights (no DataWeights), at least
                  nSub sub-ints and NBin bins: its weight, weight2, weightX are
                  used (validation against the voltage chain).
  'Seed'          random seed (default 1; own RandStream, the global one is
                  not touched)
  'Verbose'       print a summary (default true)

Model (per sub-int, phase bin j at phase phi_j = (j-1)/NBin):
  mean    m_cj = 1 + g*g_c*rho * pt(phi_j - PhaseOffset) in channel c (g_c =
          ChannelGain). Baseline 1 = the receiver-noise power of the channel. pt is the unit-peak periodic Gaussian (as
          generatePulsarSignal; or the peak-normalized sum of 'Components')
          averaged over a time bin and over the linear
          assignment triangle (+-1 phase bin), exact in Fourier space:
            pt_k = p_k * sinc(k*BinDt*f0) * sinc(k/NBin)^2.
  weights for time bins spread uniformly in phase (asymptotic), with
          n = SubintPeriods/(NBin*f0*BinDt) time bins per phase bin:
            weight = n, weight2 = 2n/3, weightX = n/6
          (exact means of the linear split; neighbour correlation 0.25).
          With correlated time bins (c(0) = V, c(+-L) = X(L), consecutive time
          bins delta = f0*BinDt*NBin phase bins apart):
            weight2 = n * sum_L c(L) K(L*delta),
            weightX(d) = n * sum_L c(L) K(d - L*delta),  d = 1..D,
          K(s) the overlap of two assignment triangles s bins apart (2/3 - s^2 +
          |s|^3/2 for |s| < 1, (2 - |s|)^3/6 for 1 <= |s| < 2); D as foldProfile.
  RFI     (C4b) channel c: mean m_cj = L_c + g*g_c*rho*pt; weights x v_c (per
          channel: the fold gets nW = NChan, as foldProfile with DataWeights); the
          fold's weight2/weightX from the clean model coefficients (what the observer
          assumes), the noise drawn with the measured ones (n*v_c*sum_L c_c(L)K(..)):
          so a channel whose RFI narrows the spectrum is noisier than the observer's
          model, a strong carrier's channel quieter (B5e); with 'ObserverNoise'
          'measured' the fold's weights use the measured ones too. Blanked bins are taken as
          spread evenly over pulse phase (RFI not locked to the pulsar).
  noise   Gaussian, var(prof_j) = s_j^2 = m_j^2/(Bnoise*BinDt) * weight2/weight^2,
          cov(prof_j, prof_j+1) = m_j*m_j+1/(Bnoise*BinDt) * weightX/(weight_j*weight_j+1),
          circular, lags 1..D - the model of estimateTOA and detectPulsar (m^2,
          so the pulsar's own noise is included); channels independent.
          Asymptotic weights: a circulant correlation matrix, drawn by FFT
          (eigenvalues 1 + 2*sum_d r_d*cos(2*pi*k*d/NBin), r_d = weightX(d)/weight2;
          C1: 1 + 0.5*cos(2*pi*k/NBin)); copied weights: sparse Cholesky factor per
          sub-int.

Outputs:
  info   the fields of foldProfile's info (f0, F1 = 0, TRef, Phi0 = 0, NBin,
         nSub, subintPeriods, nChan, assign 'linear', noiseCoeffs, covLags D,
         nWeightChan 1, validFraction 1, binDt, ...); sourceFile
         'simulateFold', units 'receiver-noise power = 1', simulated = true,
         weights ('asymptotic' or 'copied')
  fold   prof, sum [NBin x nSub x NChan], weight, weight2 [NBin x nSub],
         weightX [NBin x nSub x 1 x D], profTotal [NBin x NChan], phase, subint
         (turnFirst, turnLast, turnRef, tRef, fRef, tMean, nTimeBins) -
         foldProfile's fields and conventions (shared weights, as NoiseCoeffs)
  truth  T (= 1/f0), pulseCenters [nSub x 1] (true arrival time of each sub-
         int's reference pulse: tRef + PhaseOffset/f0; with T the fields
         validateTOA reads), phaseOffset, rhoDB, rho, signalGain,
         channelGain, FWHM, components, Bnoise, nPol, seed, peak (g*rho, pulse peak above baseline
         for g_c = 1); with 'RFI' also rfiLabel, level, validFraction,
         noiseRatio (long-run true / model noise per channel); pass (each sub-int's
         pass); with FluxScale or Scintillation fluxScale [nSub x NChan]

Notes: sub-ints are independent (foldProfile does not store the covariance
across sub-int boundaries either). Coverage: with fewer than one time bin per
phase bin per sub-int a real fold leaves phase bins empty - warning
'simulateFold:coverage' (the asymptotic weights then do not apply).
%}

arguments
    f0     (1,1) double {mustBePositive, mustBeFinite}
    rhoDB  (1,1) double {mustBeReal, mustBeNonNan}
    nTurns double {mustBeInteger, mustBeNonnegative}
    opts.FWHM          (1,1) double {mustBePositive, mustBeLessThan(opts.FWHM, 0.5)} = 0.05
    opts.Components    double = []
    opts.PhaseOffset   double {mustBeFinite} = 0
    opts.NBin          (1,1) double {mustBeInteger, mustBePositive} = 2048
    opts.SubintPeriods (1,1) double {mustBeInteger, mustBePositive} = 1
    opts.SubintTime    double = []
    opts.Passes        (:,2) double {mustBeFinite} = zeros(0, 2)
    opts.Bnoise        (1,:) double {mustBePositive} = 391.6e6
    opts.NPol          (1,1) double {mustBeInteger, mustBePositive} = 1
    opts.BinDt         (1,1) double {mustBePositive, mustBeFinite} = 1e-6
    opts.NChan         (1,1) double {mustBeInteger, mustBePositive} = 1
    opts.ChannelGain   double = []
    opts.NoiseCoeffs   (1,:) double {mustBeNonempty, mustBeFinite} = 1
    opts.MaxMemoryGB   (1,1) double {mustBePositive} = 2
    opts.RFI           = []
    opts.ObserverNoise {mustBeMember(opts.ObserverNoise, {'model', 'measured'})} = 'model'
    opts.FluxScale     double {mustBeNonnegative} = []
    opts.Scintillation = []
    opts.SignalGain    (1,1) double {mustBeNonnegative, mustBeFinite} = 1
    opts.TRef          (1,1) double = NaN
    opts.Weights       = []
    opts.Seed          (1,1) double {mustBeInteger, mustBeNonnegative} = 1
    opts.Verbose       (1,1) logical = true
end

tStart = tic;
if rhoDB == Inf
    error('simulateFold:rho', 'rhoDB must be finite or -Inf (no pulsar).');
end
N    = opts.NBin;
nPer = opts.SubintPeriods;
if ~isempty(opts.SubintTime), nPer = max(1, round(opts.SubintTime * f0)); end
TRef = opts.TRef;
if isnan(TRef), TRef = 0.5 / f0; end
% ---- Schedule: complete sub-ints of nPer turns (turn k spans phase k-1/2 .. k+1/2) ---------
if isempty(opts.Passes)
    if ~isscalar(nTurns) || nTurns < 1
        error('simulateFold:turns', 'nTurns must be a positive integer (or [] with ''Passes'').');
    end
    nSub = floor(nTurns / nPer);
    if nSub < 1
        error('simulateFold:short', 'nTurns (%d) is shorter than one sub-int (%d turns).', nTurns, nPer);
    end
    tfAll = ((1:nSub).' - 1) * nPer;                % first turn of each sub-int
    passOf = ones(nSub, 1);
else
    if ~isempty(nTurns)
        error('simulateFold:turns', 'Give either nTurns or ''Passes'', not both (nTurns = []).');
    end
    P = opts.Passes;
    if any(P(:, 2) <= P(:, 1)) || any(P(2:end, 1) < P(1:end-1, 2))
        error('simulateFold:passes', 'Passes must be [start end] rows, increasing and not overlapping.');
    end
    tfAll = zeros(0, 1); passOf = zeros(0, 1);
    for ip = 1:size(P, 1)
        k1 = ceil(f0 * (P(ip, 1) - TRef) + 0.5);    % complete turns inside the pass
        k2 = floor(f0 * (P(ip, 2) - TRef) - 0.5);
        nS = floor((k2 - k1 + 1) / nPer);
        tfAll = [tfAll; k1 + (0:nS-1).' * nPer];     %#ok<AGROW>
        passOf = [passOf; ip * ones(nS, 1)];         %#ok<AGROW>
    end
    nSub = numel(tfAll);
    if nSub < 1
        error('simulateFold:short', 'No pass holds a complete sub-int (%d turns).', nPer);
    end
end
dt   = opts.BinDt;
rho  = 10^(rhoDB / 10);                             % 0 for -Inf (no pulsar)
peak = opts.SignalGain * rho;
dPhi = opts.PhaseOffset(:);
if isscalar(dPhi)
    dPhi = repmat(dPhi, nSub, 1);
elseif numel(dPhi) ~= nSub
    error('simulateFold:offset', 'PhaseOffset must be a scalar or have nSub = %d values.', nSub);
end

% ---- RFI profile: receiver/processing settings of the calibrating run -------------------------
RP = opts.RFI;
haveRFI = ~isempty(RP);
if haveRFI
    need = {'nChan', 'level', 'noiseCoeffs', 'modelCoeffs', 'validFraction', 'binDt', 'BnoiseChan'};
    if ~isstruct(RP) || ~all(isfield(RP, need))
        error('simulateFold:rfi', 'RFI must be a profile from rfiProfile (fields %s).', strjoin(need, ', '));
    end
    if opts.NChan ~= RP.nChan
        error('simulateFold:rfi', 'NChan (%d) must equal the profile''s channels (%d).', opts.NChan, RP.nChan);
    end
    if ~isempty(opts.Weights)
        error('simulateFold:rfi', 'RFI works with the asymptotic weights only.');
    end
    if isequal(opts.BinDt, 1e-6), opts.BinDt = RP.binDt; end
    if isequal(opts.Bnoise, 391.6e6), opts.Bnoise = RP.BnoiseChan; end
    if isequal(opts.NoiseCoeffs, 1), opts.NoiseCoeffs = RP.modelCoeffs; end
    if abs(opts.BinDt - RP.binDt) > 1e-9 * RP.binDt || any(abs(opts.Bnoise - RP.BnoiseChan) > 1e-6 * RP.BnoiseChan) ...
            || ~isequal(size(opts.NoiseCoeffs), size(RP.modelCoeffs)) ...
            || any(abs(opts.NoiseCoeffs - RP.modelCoeffs) > 1e-9)
        error('simulateFold:rfi', ['BinDt, Bnoise and NoiseCoeffs must be the profile''s (%.4g s, ' ...
            '%.6g Hz, its modelCoeffs): it was measured with them.'], RP.binDt, RP.BnoiseChan);
    end
    dt = opts.BinDt;
end

% ---- Channels, time-bin correlation, memory ---------------------------------------------
nChan = opts.NChan;
gC = opts.ChannelGain;
if isempty(gC), gC = ones(1, nChan); end
gC = reshape(gC, 1, []);
if numel(gC) ~= nChan || any(~isfinite(gC)) || any(gC < 0)
    error('simulateFold:gain', 'ChannelGain must have one value >= 0 per channel (%d).', nChan);
end
Bn = reshape(opts.Bnoise, 1, []);
if ~any(numel(Bn) == [1 nChan])
    error('simulateFold:bnoise', 'Bnoise must be a scalar or one value per channel (%d).', nChan);
end
nc = opts.NoiseCoeffs;
if nc(1) <= 0
    error('simulateFold:noise', 'NoiseCoeffs(1) (the time-bin variance V) must be positive.');
end
Lmax  = numel(nc) - 1;
delta = f0 * dt * N;                                % phase bins between consecutive time bins
D = 1;
if Lmax > 0, D = ceil(Lmax * delta * (1 + 1e-9)) + 1; end   % as foldProfile
if ~isempty(opts.Weights), D = size(opts.Weights.weightX, 4); end
if 2*D >= N
    error('simulateFold:lags', 'Correlated time bins span %d of %d phase bins.', D, N);
end
nW = 1 + (nChan - 1) * haveRFI;                     % channels with their own weights
memGB = N * nSub * (2*nChan + nW*(2 + D)) * 8 / 1e9;
if memGB > opts.MaxMemoryGB
    error('simulateFold:memory', ['Fold arrays need ~%.2f GB (%d bins x %d sub-ints x %d chan). ' ...
        'Use longer sub-ints, fewer bins or channels, or raise MaxMemoryGB.'], memGB, N, nSub, nChan);
end

% ---- Weights ------------------------------------------------------------------------------
nPhase = nPer / (N * f0 * dt);                      % time bins per phase bin per sub-int
if isempty(opts.Weights)
    wSource = 'asymptotic';
    if nPhase < 1
        warning('simulateFold:coverage', ['%.3g time bins per phase bin per sub-int: a real ' ...
            'fold leaves phase bins empty, so the asymptotic weights do not apply. Use fewer ' ...
            'phase bins, shorter time bins or more turns per sub-int.'], nPhase);
    end
    if isequal(nc, 1)                               % independent time bins (C1)
        w2c = 2/3; wxc = 1/6; r = 0.25;
    else                                            % correlated time bins (channel path)
        L   = -Lmax:Lmax;
        cL  = [fliplr(nc(2:end)), nc(1), nc(2:end)];
        w2c = sum(cL .* Ktri(L * delta));
        wxc = arrayfun(@(d) sum(cL .* Ktri(d - L * delta)), 1:D);
        r   = wxc / w2c;
    end
    W  = nPhase * ones(N, nSub);
    W2 = w2c * nPhase * ones(N, nSub);
    WX = repmat(reshape(wxc * nPhase, 1, 1, 1, D), N, nSub);
    if haveRFI                                      % true coefficients per channel (measured)
        Lr  = -Lmax:Lmax;
        w2t = zeros(nChan, 1); wxt = zeros(nChan, D);
        for c = 1:nChan
            cc = RP.noiseCoeffs(c, :);
            cLc = [fliplr(cc(2:end)), cc(1), cc(2:end)];
            w2t(c) = sum(cLc .* Ktri(Lr * delta));
            wxt(c, :) = arrayfun(@(d) sum(cLc .* Ktri(d - Lr * delta)), 1:D);
        end
        vC = reshape(RP.validFraction, [], 1);      % [nChan x 1]
        lvC = reshape(RP.level, 1, 1, []);
    end
else
    wSource = 'copied';
    F = opts.Weights;
    if ~isstruct(F) || ~all(isfield(F, {'weight', 'weight2', 'weightX'}))
        error('simulateFold:weights', 'Weights must be a fold struct (weight, weight2, weightX).');
    end
    if size(F.weight, 1) ~= N || size(F.weight, 2) < nSub || size(F.weight, 3) ~= 1
        error('simulateFold:weights', ['Weights must be a fold with shared weights, NBin = %d ' ...
            'bins and at least nSub = %d sub-ints.'], N, nSub);
    end
    W  = F.weight(:, 1:nSub);
    W2 = F.weight2(:, 1:nSub);
    WX = F.weightX(:, 1:nSub, 1, :);
    if any(W <= 0, 'all')
        error('simulateFold:weights', 'The copied weights have empty phase bins.');
    end
end

% ---- Mean profile: unit-peak periodic Gaussian, smoothed, in Fourier space -------------
k   = [0:ceil(N/2)-1, -floor(N/2):-1].';          % signed harmonics (FFT order)
sincF = @(x) (x == 0) + sin(pi*x) ./ (pi*x + (x == 0));
if isempty(opts.Components)
    sig = opts.FWHM / (2*sqrt(2*log(2)));           % [turns]
    kk  = -ceil(10/sig):ceil(10/sig);               % peak of the periodic sum (the
    peakNorm = sum(exp(-2*pi^2*sig^2*kk.^2));       %   generator's normalization)
    Pk  = exp(-2*pi^2*sig^2*k.^2) / peakNorm ...    % Fourier coefficients of the profile
          .* sincF(k*dt*f0) .* sincF(k/N).^2;       %   (per turn), time bin and assignment
else                                                % sum of Gaussian components
    Cp = opts.Components;
    if size(Cp, 2) ~= 3 || any(Cp(:, 3) <= 0) || ~all(isfinite(Cp), 'all')
        error('simulateFold:components', ...
            'Components must be [nComp x 3] rows [amplitude, centre, FWHM] (turns, FWHM > 0).');
    end
    sgC = Cp(:, 3) / (2*sqrt(2*log(2)));
    Ck  = zeros(size(k));                           % Fourier series of the periodic sum
    for i = 1:size(Cp, 1)
        Ck = Ck + Cp(i, 1) * sgC(i) * sqrt(2*pi) * exp(-2*pi^2*sgC(i)^2*k.^2) .* exp(-2i*pi*k*Cp(i, 2));
    end
    pg = (0:2^16-1).' / 2^16; pf = zeros(size(pg)); % its peak on a fine grid
    for i = 1:size(Cp, 1)
        for m = -2:2
            pf = pf + Cp(i, 1) * exp(-0.5*((pg - Cp(i, 2) - m) / sgC(i)).^2);
        end
    end
    Pk  = Ck / max(pf) .* sincF(k*dt*f0) .* sincF(k/N).^2;
end

% ---- Draw per block of sub-ints -------------------------------------------------------------
rs   = RandStream('mt19937ar', 'Seed', opts.Seed);
prof = zeros(N, nSub, nChan);
c0   = 1 ./ (Bn * dt * opts.NPol);                  % radiometer variance factor (0 for Inf),
                                                    %   averaged over NPol polarizations
c0r  = reshape(c0, 1, 1, []);
pkC  = reshape(peak * gC, 1, 1, nChan);             % pulse peak per channel
% ---- Flux factor per sub-int and channel (scintillation; C3b) ------------------------------
haveF = ~isempty(opts.FluxScale) || ~isempty(opts.Scintillation);
Fsc = [];
if ~isempty(opts.FluxScale) && ~isempty(opts.Scintillation)
    error('simulateFold:flux', 'Give either FluxScale or Scintillation, not both.');
elseif ~isempty(opts.FluxScale)
    Fsc = opts.FluxScale;
    if size(Fsc, 1) ~= nSub || ~any(size(Fsc, 2) == [1 nChan])
        error('simulateFold:flux', 'FluxScale must be [nSub x 1] or [nSub x NChan] (%d x %d).', nSub, nChan);
    end
    Fsc = repmat(Fsc, 1, nChan / size(Fsc, 2));
elseif ~isempty(opts.Scintillation)
    tMid = TRef + (2*tfAll + nPer - 1) / (2*f0);    % sub-int mid times
    if nPer / f0 > opts.Scintillation.timescale / 5
        warning('simulateFold:scint', ['Sub-ints (%.3g s) are not short against the scintle time ' ...
            '(%.3g s): the flux factor is taken at mid sub-int.'], nPer / f0, opts.Scintillation.timescale);
    end
    Fsc = scintField(tMid, nChan, opts.Scintillation, opts.Seed);
end
if strcmp(wSource, 'asymptotic')
    lam   = 1 + 2 * sum(r .* cos(2*pi*k*(1:D)/N), 2);   % eigenvalues of the correlation matrix
    sqLam = sqrt(max(lam, 0));
end
if haveRFI                                          % per channel: true correlation, scale
    lamC  = 1 + 2 * (cos(2*pi*k*(1:D)/N) * (wxt ./ w2t).');    % [N x nChan]
    sqLam = reshape(sqrt(max(lamC, 0)), N, 1, nChan);
    sdC   = reshape(sqrt(c0(:) .* vC .* nPhase .* w2t) ./ (vC * nPhase), 1, 1, nChan);  % s/m
end
blk = max(1, min(512, floor(2^24 / (N * nChan))));
for b0 = 1:blk:nSub
    s  = b0 : min(b0 + blk - 1, nSub);
    nb = numel(s);
    if haveF                                        % pulse x the flux factor of each sub-int
        m = 1 + (pkC .* reshape(Fsc(s, :), 1, nb, nChan)) * N .* ...
            real(ifft(Pk .* exp(-2i*pi*k*dPhi(s).'), [], 1));
    else
        m = 1 + pkC * N .* real(ifft(Pk .* exp(-2i*pi*k*dPhi(s).'), [], 1));  % [N x nb x nChan]
    end
    if haveRFI, m = m - 1 + lvC; end                % baseline L_c (RFI left after excision)
    if all(c0 == 0)
        prof(:, s, :) = m;
        continue
    end
    e = reshape(randn(rs, N, nb * nChan), N, nb, nChan);
    if strcmp(wSource, 'asymptotic')
        z = real(ifft(sqLam .* fft(e, [], 1), [], 1));
    else
        z = zeros(N, nb, nChan);
        for i = 1:nb
            Lc = corrFactor(W2(:, s(i)), WX(:, s(i), 1, :));
            z(:, i, :) = reshape(Lc * reshape(e(:, i, :), N, nChan), N, 1, nChan);
        end
    end
    if haveRFI
        sd = m .* sdC;                              % true noise (measured coefficients)
    else
        sd = m .* sqrt(c0r .* W2(:, s)) ./ W(:, s); % s_j per channel
    end
    prof(:, s, :) = m + sd .* z;
end

% ---- Assemble fold struct (foldProfile's fields) -------------------------------------------
if haveRFI                                          % per-channel weights x v_c (observer's model)
    vR = reshape(vC, 1, 1, nChan);
    W  = W  .* vR;
    if strcmp(opts.ObserverNoise, 'measured')       % calibrated observer: the measured
        W2 = nPhase * ones(N, nSub) .* reshape(w2t .* vC, 1, 1, nChan);           %   coefficients
        WX = nPhase * ones(N, nSub) .* reshape(wxt .* vC, 1, 1, nChan, D);
    else
        W2 = W2 .* vR;
        WX = permute(WX, [1 2 3 4]) .* vR;          % [N x nSub x nChan x D]
    end
    if any(vC == 0)                                 % fully blanked channels: no data
        prof(:, :, vC == 0) = NaN;
    end
end
fold = struct();
fold.sum     = prof .* W;
fold.sum(isnan(fold.sum)) = 0;
fold.weight  = W;
fold.weight2 = W2;
fold.weightX = WX;
fold.prof    = prof;
fold.profTotal = reshape(sum(fold.sum, 2), N, nChan) ./ reshape(sum(W, 2), N, []);
fold.phase   = (0:N-1) / N;
tf = tfAll;
tl = tf + nPer - 1;
tr = round((tf + tl) / 2);
fold.subint = struct('turnFirst', tf, 'turnLast', tl, 'turnRef', tr, ...
    'tRef', TRef + tr / f0, 'fRef', f0 * ones(nSub, 1), ...
    'tMean', TRef + (tf + tl) / (2*f0), 'nTimeBins', round(nPer / (f0 * dt)) * ones(nSub, 1));

% ---- Info (only what foldProfile's info holds; no truth) -----------------------------------
info = struct();
info.file         = '';
info.sourceFile   = 'simulateFold';
info.f0           = f0;
info.F1           = 0;
info.TRef         = TRef;
info.Phi0         = 0;
info.phaseModel   = 'phi(t) = Phi0 + f0*(t-TRef) + F1*(t-TRef)^2/2';
info.phaseConvention = 'bin j (1-based) centred at phase (j-1)/NBin; turn = round(phi)';
info.NBin         = N;
info.nSub         = nSub;
info.subintPeriods = nPer;
info.passes       = opts.Passes;                    % the observing schedule (observer knowledge)
info.nChan        = nChan;
info.chanFreqs    = [];
info.assign       = 'linear';
info.noiseCoeffs  = nc;
if haveRFI && strcmp(opts.ObserverNoise, 'measured')
    info.noiseCoeffs = RP.noiseCoeffs;              % per channel [nChan x (1+Lmax)]
end
info.covLags      = D;
info.noiseModel   = ['var(sum_j) = sigma^2*weight2, cov(sum_j, sum_j+d) = ' ...
                     'sigma^2*weightX(:,:,:,d), sigma^2 = m^2/(Bnoise*binDt)'];
info.dataWeights  = '';
info.nWeightChan  = nW;
info.validFraction = 1;
if haveRFI
    info.dataWeights   = 'simulated (RFI profile)';
    info.validFraction = mean(vC);
end
info.binDt        = dt;
info.binsFolded   = [];
info.units        = 'receiver-noise power = 1';
info.simulated    = true;
info.weights      = wSource;
info.elapsed      = toc(tStart);

truth = struct('T', 1/f0, 'pulseCenters', fold.subint.tRef + dPhi / f0, ...
    'phaseOffset', dPhi, 'rhoDB', rhoDB, 'rho', rho, 'signalGain', opts.SignalGain, ...
    'channelGain', gC, 'pass', passOf, 'FWHM', opts.FWHM, 'components', opts.Components, ...
    'Bnoise', opts.Bnoise, 'nPol', opts.NPol, ...
    'seed', opts.Seed, 'peak', peak);
if haveF, truth.fluxScale = Fsc; end
if haveRFI
    truth.rfiLabel = RP.label; truth.level = RP.level; truth.validFraction = RP.validFraction;
    truth.noiseRatio = ((RP.noiseCoeffs(:, 1) + 2*sum(RP.noiseCoeffs(:, 2:end), 2)) / ...
        (RP.modelCoeffs(1) + 2*sum(RP.modelCoeffs(2:end)))).';
end

if opts.Verbose
    fprintf(['simulateFold: %d sub-int(s) of %d turn(s) (%.4g s of data), NBin %d, %d channel(s), ' ...
             'rho %g dB (peak %.3g of the baseline), %s weights (%.3g time bins per phase bin, ' ...
             '%d lag(s)), %.2f s\n'], nSub, nPer, nSub * nPer / f0, N, nChan, rhoDB, peak, ...
        wSource, nPhase, D, info.elapsed);
end
end


% =====================================================================================
function L = corrFactor(w2, wx)
%CORRFACTOR  Lower Cholesky factor of the circular banded correlation matrix with
% r_jd = wx_jd / sqrt(w2_j * w2_j+d) between phase bins j and j+d, d = 1..D (copied weights).
N = numel(w2);
wx = reshape(wx, N, []);
i = (1:N).'; I = i; J = i; V = ones(N, 1);
for d = 1:size(wx, 2)
    j = mod(i - 1 + d, N) + 1;                      % bin j+d, circular
    r = wx(:, d) ./ sqrt(w2 .* circshift(w2, -d));
    I = [I; i; j]; J = [J; j; i]; V = [V; r; r];    %#ok<AGROW>
end
R = sparse(I, J, V, N, N);
L = chol(R, 'lower');
end

function K = Ktri(s)
%KTRI  Overlap of two linear-assignment triangles (half-width 1 bin) s bins apart.
s = abs(s);
K = zeros(size(s));
a = s < 1;            K(a) = 2/3 - s(a).^2 + s(a).^3 / 2;
b = s >= 1 & s < 2;   K(b) = (2 - s(b)).^3 / 6;
end

function F = scintField(t, nChan, S, seed)
%SCINTFIELD  Flux factor F = |g|^2 [numel(t) x nChan] of a complex Gaussian field g(t, nu)
% with Gaussian correlations: intensity ACF exp(-(dt/timescale)^2) in time (1/e at the
% scintle time) and exp(-ln2*(dnu/bandwidth)^2) in frequency (1/2 at the scintle
% bandwidth). g = white complex noise convolved with Gaussian kernels of std sigma/sqrt(2)
% (so the field's own correlation has std sigma) on a padded uniform grid; time by linear
% interpolation; channels = mean of F over their frequency points. Frequency grid: at least
% 64 points across the band and a step of at most a quarter of the scintle bandwidth (20 Oct
% 2026; was one point per channel at 128 channels, which over-states a channel's fading when
% the scintles are not much wider than the channel). Own random stream.
need = {'bandwidth', 'timescale'};
if ~isstruct(S) || ~all(isfield(S, need))
    error('simulateFold:scint', 'Scintillation must be a struct with fields bandwidth, timescale.');
end
band = [1.2e9 1.6e9];
if isfield(S, 'band'), band = S.band; end
sigT = S.timescale;                                 % std of the field's time correlation
sigF = S.bandwidth / sqrt(log(2));                  % ... and of its frequency correlation
chanW  = diff(band) / nChan;
nPerCh = max([1, ceil(64 / nChan), ceil(4 * chanW / S.bandwidth)]);   % frequency points per
if nChan * nPerCh > 2^16                            %   channel: step <= bandwidth/4
    warning('simulateFold:scintGrid', ['Scintle bandwidth %.3g Hz needs %d frequency points; ' ...
        'limited to %d (the channel averages over-state the fading).'], S.bandwidth, ...
        nChan * nPerCh, nChan * floor(2^16 / nChan));
    nPerCh = max(1, floor(2^16 / nChan));
end
nuS = band(1) + ((1:nChan*nPerCh) - 0.5) * diff(band) / (nChan*nPerCh);
dNu = diff(band) / (nChan*nPerCh);
dT  = sigT / 10;
padT = ceil(3 * sigT / dT); padF = ceil(3 * sigF / dNu);
tG  = (min(t) - padT*dT : dT : max(t) + padT*dT + dT).';
kT  = exp(-((-padT:padT) * dT).^2 / (2 * (sigT/sqrt(2))^2)); kT = kT / sqrt(sum(kT.^2));
kF  = exp(-((-padF:padF) * dNu).^2 / (2 * (sigF/sqrt(2))^2)); kF = kF / sqrt(sum(kF.^2));
nT  = numel(tG) + 2*padT; nF = numel(nuS) + 2*padF;
rs  = RandStream('mrg32k3a', 'Seed', seed);
e   = (randn(rs, nT, nF) + 1i * randn(rs, nT, nF)) / sqrt(2);
g   = conv2(kT(:), kF(:), e, 'valid');               % [numel(tG) x numel(nuS)]
gt  = interp1(tG, g, t(:), 'linear');               % at the sub-int times
Fp  = abs(gt).^2;
F   = squeeze(mean(reshape(Fp, numel(t), nPerCh, nChan), 2));
F   = reshape(F, numel(t), nChan);
end
