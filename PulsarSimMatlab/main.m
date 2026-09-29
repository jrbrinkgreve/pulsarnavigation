%{
The main file for pulsar navigation data cleaning with synthetic data
Main thought: Pulsar signals arrive at a known location at very accurate
timings! Therefore, when we arrive such a signal earlier, it means we are
closer to the pulsar in that direction. Combining multiple pulsars means we
can locate ourselves therefore.
General structure
Data generation:
- 'Clean' signal from pulsar on receiver:
- Interstellar dispersion
- Pulsar-Earth relative movement
- Receiver noise
- Earth-RF interference
- Clock jitter addition (note: look up how to set this up)
- Polarisation?
SNR-improving processing:
- Channel estimation / Interstellar channel deconvolution
- Temporally folding observations
- 3x3-array spatial steering
- Conversion of sampling rate to pulsar "template" sampling rate
- Matched filtering / FFT-correlation based shift detection
- ... potential other ways of improving
    - Doppler-based processing
Detection
- NP-type detector? to study still how to implement this exactly
- Phase detection - -
Pulsar and phase detection --> localisation
- Barycentric corrections with orbits
Validation
- 'Checks' against ground truth? 
%}

% Parameters
% --- Pulsar / signal generation
T          = 10e-3;     % [s]  pulsar period
f_in       = 4e9;       % [Hz] generation sampling rate
A          = 1;         % [V]  noise std dev at pulse peak
L          = 1000e-3;    % [s]  total signal length
dutycycle  = 5;         % [%]  pulse FWHM as percent of T
genEnvMode = 'power';   % 'power': power profile has FWHM = dutycycle% of T
seed       = 42;        % RNG seed for reproducible runs

% --- Interstellar dispersion
DM    = 5;              % [pc cm^-3] dispersion measure
fLow  = 1.2e9;          % [Hz] lower edge of observed RF band
fHigh = 1.6e9;          % [Hz] upper edge of observed RF band

% --- Receiver front end (IQ downconversion)
fLO         = 1.4e9;    % [Hz] local oscillator
fs          = 1e9;      % [Hz] requested sample rate out of IQ stage
filterOrder = 2048;

% --- Envelope reconstruction
detMode = 'power';      % envelope detection mode ('mode' would shadow MATLAB's mode())
f_out   = 1e6;          % [Hz] envelope sampling rate

% --- Files
dataDir       = "data";
fileRaw       = fullfile(dataDir, "test.dat");
fileDispersed = fullfile(dataDir, "test_dispersed.dat");
fileIQ        = fullfile(dataDir, "test_dispersed_IQ_modulated.dat");
fileDedisp    = fullfile(dataDir, "test_IQ_dedispersed.dat");
fileEnvelope  = fullfile(dataDir, "test_envelope.dat");
fileFold      = fullfile(dataDir, "test_fold.mat");








%Synthetic data generation

% 1. Pulsar signal generation
info_gen = generatePulsarSignal(fileRaw, T, f_in, A, L, dutycycle, 'Seed', seed, 'EnvelopeMode', genEnvMode, 'Verbose', true);

% 2. Interstellar dispersion
info_disp = applyDispersionStream(info_gen.file, fileDispersed, DM, info_gen.actualFsOut, fLow, fHigh);

% 2.5: plotting for checking
%plotDispersionCheck(info_gen, info_disp, 0, 20e-3);   % first 20 ms

% 3. IQ downconversion
info_IQ = applyIQmodulation(info_disp.file, fileIQ, info_disp.actualFsOut, fs, fLO, 'FilterOrder', filterOrder, 'Band', [fLow fHigh]);








%Processing pipeline

% 4. Coherent dedispersion
info_dedisp = applyInverseDispersion(info_IQ.file, fileDedisp, info_IQ.actualFsOut, info_IQ.fLO, info_disp.DM, info_disp.fLow, info_disp.fHigh, 'RefFreq', info_disp.refFreq);

%4.5 visual de-dispersion check
%[~, check] = plotIQCheck(info_gen, info_disp, info_IQ, info_dedisp, 0, 20e-3);

% 5. Square-law detection
info_det = detectPower(info_dedisp.file, fileEnvelope, info_dedisp.fs, info_dedisp.fLO, f_out, 'FullySupported', info_dedisp.fullySupported);

%5.5: plotting
[~, check_det] = plotDetectedPower(info_det, info_gen, 0, 20e-3, 'InfoDisp', info_disp, 'InfoIQ', info_IQ, 'InfoDedisp', info_dedisp);


%6. Folding
[info_fold, fold] = foldProfile(info_det, fileFold, 1/T,'TRef', info_gen.pulseCenters(1), 'NBin', 1024, 'SubintPeriods', 1);

%6.5 Visual check: fold
[~, check_fold] = plotFoldCheck(fold, info_fold, info_gen, 'InfoDisp', info_disp, 'InfoIQ', info_IQ, 'InfoDedisp', info_dedisp);

