%% Synthetic Radio Pulsar Navigation Data Generator & Pipeline Testbed
% ------------------------------------------------------------------------
% Purpose
%   Generates raw "search-mode" filterbank data (power vs. frequency
%   channel vs. time) for several radio pulsars, including:
%     - the intrinsic periodic pulse profile
%     - interstellar dispersion (frequency-dependent delay)
%     - radiometer (thermal) noise, via the radiometer equation
%     - a hidden navigation-induced timing offset + receiver clock bias
%   It then runs that data through a first-draft receiver pipeline:
%     dedispersion -> folding -> sub-bin TOA estimation -> position fix
%   and reports the recovered spacecraft position vs. the injected truth.
%
%   Use this as a testbed: swap in real pulsar parameters (ATNF catalogue),
%   sweep SNR/DM/RFI, and replace pieces with higher-fidelity models as
%   your pipeline matures. Run section-by-section (Ctrl+Enter per %% cell)
%   to inspect intermediate results while you debug.
%
%   NOTE: parameter values below (P0, DM, flux) are representative,
%   rounded figures for well-known millisecond pulsars, not
%   ephemeris-grade. Replace with precise values from the ATNF Pulsar
%   Catalogue (https://www.atnf.csiro.au/research/pulsar/psrcat/) before
%   drawing any quantitative conclusions.
% ------------------------------------------------------------------------

clear; clc; close all;
rng(1);   % reproducibility

%% ---------------- 1. USER-CONFIGURABLE PARAMETERS ----------------
c = 299792458;   % speed of light [m/s]

% --- Pulsar catalog: [name, P0 (s), DM (pc/cm^3), W50/P0, S1400 (mJy), RA (deg), Dec (deg)] ---
psr(1) = struct('name','J0437-4715', 'P0',5.757451937e-3, 'DM',2.64476,  'W50_frac',0.020, 'S1400_mJy',149.0, 'raDeg',69.3156,  'decDeg',-47.2527);
psr(2) = struct('name','B1937+21',   'P0',1.557806449e-3, 'DM',71.0398,  'W50_frac',0.020, 'S1400_mJy',12.0,  'raDeg',294.9106, 'decDeg', 21.5831);
psr(3) = struct('name','B1821-24',   'P0',3.054965123e-3, 'DM',119.866,  'W50_frac',0.030, 'S1400_mJy',3.1,   'raDeg',276.1334, 'decDeg',-24.8697);
psr(4) = struct('name','J0034-0534', 'P0',1.877275895e-3, 'DM',13.7647,  'W50_frac',0.030, 'S1400_mJy',1.7,   'raDeg',8.5993,   'decDeg', -5.5745);
Npsr = numel(psr);

% --- Receiver / telescope parameters ---
fc_MHz    = 1400;    % band center frequency [MHz]
BW_MHz    = 400;     % total receiver bandwidth [MHz]
Nchan     = 512;       % number of filterbank channels
dt_native = 10e-6;    % native (post-detection) sample time [s]
Tobs      = 10.0;      % observation length per pulsar [s]
Tsys_K       = 0.1;   % system temperature [K]
Gain_KperJy  = 2.0;   % telescope gain [K/Jy]
Npol         = 2;     % summed polarizations

% --- Ground truth the pipeline must recover ---
r_true_km   = [1.2e6, -8.4e5, 3.1e5];   % spacecraft position rel. to SSB [km]
clockBias_s = 3.7e-6;                    % unknown receiver clock offset [s]

%% ---------------- 2. DERIVED QUANTITIES ----------------
chanBW_MHz  = BW_MHz/Nchan;
f_edges_MHz = linspace(fc_MHz-BW_MHz/2, fc_MHz+BW_MHz/2, Nchan+1);
f_chan_MHz  = (f_edges_MHz(1:end-1) + f_edges_MHz(2:end))/2;
Nsamp  = round(Tobs/dt_native);
t_axis = (0:Nsamp-1)*dt_native;

%% ---------------- 3. GENERATE SYNTHETIC DATA PER PULSAR ----------------
data = struct([]);
for k = 1:Npsr
    los = radecToUnitVector(psr(k).raDeg, psr(k).decDeg);
    navDelay_s  = -dot(los, r_true_km*1e3)/c;     % geometric delay from position
    toaOffset_s = navDelay_s + clockBias_s;        % hidden quantity the pipeline must find

    filterbank = simulatePulsarFilterbank(psr(k), f_chan_MHz, chanBW_MHz, ...
                    dt_native, Nsamp, toaOffset_s, Tsys_K, Gain_KperJy, Npol);

    data(k).name  = psr(k).name;
    data(k).los   = los;
    data(k).filterbank  = filterbank;      % [Nchan x Nsamp] raw power data
    data(k).f_chan_MHz  = f_chan_MHz;
    data(k).dt    = dt_native;
    data(k).truthTOAoffset_s = toaOffset_s;
    data(k).P0    = psr(k).P0;
    data(k).W50_frac = psr(k).W50_frac;
end

%% ---------------- 4. RECEIVER PIPELINE: DEDISPERSE -> FOLD -> TOA ----------------
% Here DM is assumed known; for a real search pipeline you would run a
% trial-DM grid and pick the DM that maximizes folded S/N (see notes below).
DM_trial = [psr.DM];
measTOA  = zeros(1,Npsr);

for k = 1:Npsr
    dedisp = incoherentDedisperse(data(k).filterbank, data(k).f_chan_MHz, DM_trial(k), data(k).dt);
    Nbins  = 128;
    [profile, template] = foldProfile(dedisp, data(k).dt, data(k).P0, Nbins, data(k).W50_frac);
    dtoa = fftPhaseTOA(profile, template, data(k).P0);
    measTOA(k) = dtoa;
    fprintf('%-12s  true TOA offset = %8.3f us | recovered = %8.3f us | error = %7.3f us\n', ...
        data(k).name, data(k).truthTOAoffset_s*1e6, dtoa*1e6, (dtoa-data(k).truthTOAoffset_s)*1e6);
end

%% ---------------- 5. NAVIGATION SOLUTION (position + clock bias) ----------------
% Each pulsar TOA gives one linear equation:  -los_i . r + c*clockBias = c*TOA_i
% With >=4 pulsars of independent directions this is solvable (like GPS
% pseudorange multilateration, but using only line-of-sight *direction*,
% since pulsars are effectively at infinite range).
H = zeros(Npsr,4);
y = zeros(Npsr,1);
for k = 1:Npsr
    H(k,:) = [-data(k).los, 1];
    y(k)   = measTOA(k)*c;
end
sol = H\y;                       % least squares (exact if Npsr==4)
r_est_km        = sol(1:3)'/1e3;
clockBias_est_s = sol(4)/c;

fprintf('\nTrue position      [km]: %s\n', mat2str(r_true_km,6));
fprintf('Estimated position [km]: %s\n', mat2str(r_est_km,6));
fprintf('Position error      [km]: %.4f\n', norm(r_est_km-r_true_km));
fprintf('Clock bias true/est [us]: %.3f / %.3f\n', clockBias_s*1e6, clockBias_est_s*1e6);

%% ---------------- 6. DIAGNOSTIC PLOTS ----------------
figure;
imagesc(t_axis*1e3, data(1).f_chan_MHz, data(1).filterbank); axis xy;
xlabel('Time [ms]'); ylabel('Frequency [MHz]');
title(sprintf('%s: dispersed filterbank (dynamic spectrum)', data(1).name));
colorbar;

dedisp1 = incoherentDedisperse(data(1).filterbank, data(1).f_chan_MHz, DM_trial(1), data(1).dt);
figure;
plot(t_axis*1e3, dedisp1);
xlabel('Time [ms]'); ylabel('Power [Jy, approx.]');
title(sprintf('%s: dedispersed time series', data(1).name));

%% ================= LOCAL FUNCTIONS =================

function u = radecToUnitVector(raDeg, decDeg)
    % Unit vector toward a source given RA/Dec (deg), in an equatorial frame.
    ra = deg2rad(raDeg); dec = deg2rad(decDeg);
    u = [cos(dec)*cos(ra), cos(dec)*sin(ra), sin(dec)];
end

function fb = simulatePulsarFilterbank(psrp, f_chan_MHz, chanBW_MHz, dt, Nsamp, toaOffset_s, Tsys_K, Gain_KperJy, Npol)
    % Builds an [Nchan x Nsamp] search-mode power filterbank: a dispersed,
    % noisy periodic pulse train for one pulsar.
    Nchan = numel(f_chan_MHz);
    t = (0:Nsamp-1)*dt;
    sigma_phase = psrp.W50_frac/(2*sqrt(2*log(2)));   % Gaussian sigma, in phase units (fraction of P0)

    k_DM     = 4.148808e3;                 % ms . MHz^2 / (pc/cm^3), dispersion constant
    SEFD_Jy  = Tsys_K/Gain_KperJy;          % system equivalent flux density [Jy]
    Smean_Jy = psrp.S1400_mJy/1000;         % mean flux density [Jy]
    avgProf  = sqrt(2*pi)*sigma_phase;      % mean of a unit-height wrapped-Gaussian profile over 1 period

    fb = zeros(Nchan, Nsamp);
    for ci = 1:Nchan
        delay_s = (k_DM*psrp.DM/f_chan_MHz(ci)^2)/1000;         % ISM dispersion delay at this channel [s]
        phase   = mod((t - toaOffset_s - delay_s)/psrp.P0, 1);
        prof    = wrappedGaussianProfile(phase, sigma_phase);    % peak-normalized pulse shape
        signal_Jy = (Smean_Jy/avgProf) * prof;                   % scaled so time-average == Smean_Jy

        noise_sigma_Jy = SEFD_Jy/sqrt(Npol*chanBW_MHz*1e6*dt);   % radiometer noise per raw sample
        fb(ci,:) = signal_Jy + noise_sigma_Jy*randn(1,Nsamp);
    end
end

function y = wrappedGaussianProfile(phase, sigma)
    % Periodic (wrapped) Gaussian pulse shape, peak-normalized to 1.
    y = zeros(size(phase));
    for n = -3:3
        y = y + exp(-0.5*((phase-n)/sigma).^2);
    end
end

function dedisp = incoherentDedisperse(fb, f_chan_MHz, DM, dt)
    % Shift-and-add incoherent dedispersion: aligns each channel to the
    % top-of-band reference frequency using integer-sample shifts.
    % NOTE: circshift wraps at the data edges; in a real pipeline, pad or
    % trim the first/last max-delay-worth of samples to avoid edge artifacts.
    Nchan = size(fb,1); Nsamp = size(fb,2);
    k_DM = 4.148808e3;
    fref = max(f_chan_MHz);
    dedisp = zeros(1,Nsamp);
    for ci = 1:Nchan
        delay_s = (k_DM*DM/f_chan_MHz(ci)^2 - k_DM*DM/fref^2)/1000;
        shiftSamples = round(delay_s/dt);
        dedisp = dedisp + circshift(fb(ci,:), -shiftSamples);
    end
end

function [profile, template] = foldProfile(x, dt, P0, Nbins, W50_frac)
    % Synchronous folding at the known topocentric period, plus a
    % noiseless template of the same pulse shape for TOA matching.
    Nsamp = numel(x);
    t = (0:Nsamp-1)*dt;
    phaseBin = floor(mod(t/P0,1)*Nbins) + 1;
    profile = accumarray(phaseBin', x', [Nbins,1]) ./ accumarray(phaseBin', 1, [Nbins,1]);
    profile = profile(:)';

    phi = ((0:Nbins-1)+0.5)/Nbins;
    sigma_phase = W50_frac/(2*sqrt(2*log(2)));
    template = wrappedGaussianProfile(phi, sigma_phase);

    template = template - mean(template);
    profile  = profile  - mean(profile);
end

function dtoa = fftPhaseTOA(profile, template, P0)
    % Sub-bin TOA via FFT cross-correlation + parabolic interpolation of
    % the peak. (For higher precision, replace with the standard
    % Fourier phase-gradient fit of Taylor 1992.)
    Nbins = numel(profile);
    cc = real(ifft(fft(profile) .* conj(fft(template))));
    [~, idx] = max(cc);

    if idx == 1, ip1 = Nbins; im1 = 2;
    elseif idx == Nbins, ip1 = 1; im1 = Nbins-1;
    else, ip1 = idx+1; im1 = idx-1;
    end
    y0 = cc(im1); y1 = cc(idx); y2 = cc(ip1);
    delta = 0.5*(y0-y2)/(y0-2*y1+y2+eps);

    peakBin = idx - 1 + delta;              % zero-based bin offset
    if peakBin > Nbins/2, peakBin = peakBin - Nbins; end
    dtoa = peakBin/Nbins * P0;
end

%% ---------------- EXTENSION IDEAS ----------------
% - DM search: replace DM_trial with a grid search that maximizes folded S/N,
%   to test your pipeline's dedispersion-search stage, not just its tracking stage.
% - RFI: add narrowband spikes / broadband bursts to specific channels/times
%   to test flagging and excision logic.
% - Scintillation: multiply each channel's signal by a slowly-varying gain
%   (e.g. correlated log-normal process) to mimic ISM scintillation.
% - Precision sweep: vary Tsys/Gain/BW/Tobs and check that recovered TOA
%   error tracks the expected ~ W50/SNR scaling -- a good pipeline sanity check.
% - Real geometry: replace the flat RA/Dec unit-vector model with actual
%   barycentric ephemerides (e.g. via a SPICE toolkit) for navigation-grade fidelity.
% - Dynamics: turn r_true into a moving trajectory and feed the TOA
%   residuals into an EKF/UKF for continuous orbit determination, rather
%   than a single-epoch least-squares fix.
% - Coherent dedispersion: for the narrowest millisecond-pulsar profiles,
%   incoherent (post-detection) dedispersion smears the pulse within a
%   channel; coherent (baseband, phase-preserving) dedispersion removes
%   this and is standard for precision timing.
% - Validation: cross-check TOA/whitening behavior against established
%   pulsar timing software (PINT, TEMPO2) and, once available, real
%   archival filterbank/PSRFITS data.