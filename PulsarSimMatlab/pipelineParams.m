%PIPELINEPARAMS  Parameters of the synthetic pipeline (script).
%{
Shared by main.m and runMonteCarlo.m, so every script uses the same
values. Run as a script (`pipelineParams;`): the variables land directly
in the caller's workspace. Needs functions/ on the path (rfiSource) and
MATLAB's current folder = PulsarSimMatlab (dataDir is relative).

Later, for real data: split into simulation parameters (pulsar truth,
f_in, L, seed, SNR, RFI) and observation parameters (receiver, processing).
%}

% --- Pulsar (ground truth for generation)
T          = 10e-3;     % [s]  pulsar period
A          = 1;         % [V]  pulsar noise std at pulse peak (sets the signal scale)
dutycycle  = 5;         % [%]  pulse FWHM as percent of T
genEnvMode = 'power';   % 'power': power profile has FWHM = dutycycle% of T
DM         = 5;         % [pc cm^-3] dispersion measure

% --- Simulation
f_in = 4e9;             % [Hz] RF sampling rate
L    = 100e-3;               % [s]  total signal length
seed = 42;              % RNG seed of the pulsar signal

% --- Receiver (known to the observer)
fLow        = 1.2e9;    % [Hz] lower edge of observed RF band
fHigh       = 1.6e9;    % [Hz] upper edge of observed RF band
fLO         = 1.3e9;    % [Hz] local oscillator
fs          = 800e6;    % [Hz] IQ sample rate
filterOrder = 2048;     % IQ low-pass order

% --- Receiver noise and interference (added at RF, after dispersion)
% snrDB: pulse-peak signal PSD / receiver-noise PSD in the band (S_peak/SEFD).
% -20 dB -> ~3.8 SNR per pulse (0.5 ms pulses, 400 MHz); Inf -> no noise.
% addNoiseAndRFI prints what a value means (SNR per pulse / folded, TOA error).
snrDB     = -5;
noiseSeed = seed + 1;   % different seed -> new noise, same pulsar realization
rfiOn     = false;      % validate noise-only first, then switch RFI on
rfiSelect = [];         % sources of the list below to switch on, e.g. 4 = carrier only
                        % (1 L1, 2 L2, 3 radar, 4 carrier, 5 impulses); [] = all

% Illustrative L-band RFI scenario; INR = power relative to all receiver
% noise in the 400 MHz band (while the source is on).
rfi = [ ...
    rfiSource('bpsk',    'Freq', 1575.42e6, 'ChipRate', 1.023e6, 'INRdB', -10, 'Label', 'GNSS L1 C/A-like'), ...
    rfiSource('bpsk',    'Freq', 1227.60e6, 'ChipRate', 10.23e6, 'INRdB', -15, 'Label', 'GNSS L2 P(Y)-like'), ...
    rfiSource('pulsed',  'Freq', 1300e6, 'PulseWidth', 2e-6, 'PRF', 373, 'ChirpBW', 1e6, ...
              'INRdB', 20, 'Label', 'L-band radar'), ...
    rfiSource('cw',      'Freq', 1351.3e6, 'INRdB', -5, 'Label', 'spurious carrier'), ... % not 1350 (ch 48/49 edge)
    rfiSource('impulse', 'Rate', 50, 'Duration', 200e-9, 'INRdB', 20, 'Label', 'broadband impulses')];
rfiScenario = rfi;      % the full list (runRFITest picks from it)
if ~isempty(rfiSelect), rfi = rfi(rfiSelect); end
if ~rfiOn || isinf(snrDB), rfi = struct([]); end

% --- Processing
refFreq       = fHigh;  % [Hz] dedispersion reference frequency (same as generation)
f_out         = 1e6;    % [Hz] detected-power bin rate
nBin          = 2048;   % phase bins in the fold
subintPeriods = 1;     % turns per sub-integration; aim for SNR >= 6-7 per sub-int
                       % (threshold sweep 6 Oct: 92 % / 100 % usable TOAs at SNR 6 / 7.5;
                       % turns needed N = (7 / SNR per turn)^2, SNR per turn ~ 373*rho here)
frontEnd      = 'channels';  % 'channels' (default since 8 Oct): the channelized front end
                             % (channelizeIQ -> dedisperseChannels -> detectChannels; where RFI
                             % excision will go, block B), or 'fullband': the reference path
weighting     = 'optimal';   % channel path: how estimateTOA / detectPulsar combine channels:
                             % 'optimal' (weights 1/variance, all channels) or 'equal'
chanWidth     = 3.125e6;     % [Hz] channel width of the channelized front end
excision      = true;        % channel path: RFI excision before dedispersion (block B):
                             % detectRFI -> blankChannels -> ... -> blankingWeights -> fold
                             % with DataWeights (default since 8 Oct, Jasper); false = none
excisionArgs  = {};          % name-value options for detectRFI, e.g. {'PFA', 1e-6}
                             % (defaults: windows 1-16 samples, PFA 1e-6, guard 10)

% --- Observer knowledge (ephemeris) - in real use from a pulsar catalogue
% Synthetic case: spin frequency from T; phase 0 at the first pulse centre
% (T/2), referenced to refFreq.
ephem.f0   = 1/T;
ephem.F1   = 0;
ephem.TRef = 0.5 * T;
ephem.DM   = DM;
ephem.profileFWHM = dutycycle/100;   % [turns] power-profile FWHM, for the template
                                     % (genEnvMode 'amplitude' would need /sqrt(2))

% --- Files
dataDir       = "data";
fileRaw       = fullfile(dataDir, "test.dat");
fileDispersed = fullfile(dataDir, "test_dispersed.dat");
fileRx        = fullfile(dataDir, "test_rx.dat");
fileIQ        = fullfile(dataDir, "test_rx_IQ.dat");
fileDedisp    = fullfile(dataDir, "test_IQ_dedispersed.dat");
fileEnvelope  = fullfile(dataDir, "test_envelope.dat");
fileFold      = fullfile(dataDir, "test_fold.mat");
% channelized front end (frontEnd = 'channels'); the unit tests in tests/ use these files too
chanDir        = fullfile(dataDir, "chan");
fileChanIQ     = fullfile(chanDir, "test_rx_IQ_chan");             % channelizeIQ (base name)
% with excision: own files, so the unblanked ones (the tests' reference) stay
fileChanBlank   = fullfile(chanDir, "test_rx_IQ_chan_blanked");          % blankChannels (base name)
fileChanDedispB = fullfile(chanDir, "test_dedispB_chan");                % dedisperseChannels, blanked
fileChanPowerB  = fullfile(chanDir, "test_dedispB_chan_power.dat");      % detectChannels, blanked
fileChanWeight  = fullfile(chanDir, "test_dedispB_chan_weights.dat");    % blankingWeights
fileChanDedisp = fullfile(chanDir, "test_dedisp_chan");            % dedisperseChannels (base name)
fileChanPower  = fullfile(chanDir, "test_dedisp_chan_power.dat");  % detectChannels
fileFoldChan   = fullfile(dataDir, "test_fold_chan.mat");
