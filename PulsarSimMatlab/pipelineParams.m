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
rfiList   = 'illustrative';  % which list below: 'illustrative' (5 sources) or 'realistic'
                             % (B5e, 9 Oct 2026: 16 sources); both stay defined for runRFITest
rfiSelect = [];         % sources of the chosen list to switch on, e.g. 4 = carrier only
                        % (illustrative: 1 L1, 2 L2, 3 radar, 4 carrier, 5 impulses); [] = all

% Illustrative L-band RFI scenario; INR = power relative to all receiver
% noise in the 400 MHz band (while the source is on).
rfi = [ ...
    rfiSource('bpsk',    'Freq', 1575.42e6, 'ChipRate', 1.023e6, 'INRdB', -10, 'Label', 'GNSS L1 C/A-like'), ...
    rfiSource('bpsk',    'Freq', 1227.60e6, 'ChipRate', 10.23e6, 'INRdB', -15, 'Label', 'GNSS L2 P(Y)-like'), ...
    rfiSource('pulsed',  'Freq', 1300e6, 'PulseWidth', 2e-6, 'PRF', 373, 'ChirpBW', 1e6, ...
              'INRdB', 20, 'Label', 'L-band radar'), ...
    rfiSource('cw',      'Freq', 1351.3e6, 'INRdB', -5, 'Label', 'spurious carrier'), ... % not 1350 (ch 48/49 edge)
    rfiSource('impulse', 'Rate', 50, 'Duration', 200e-9, 'INRdB', 20, 'Label', 'broadband impulses')];
rfiScenario = rfi;      % the full illustrative list (runRFITest picks from it)

% Realistic L-band scenario (B5e, 9 Oct 2026) for the reference hardware: 3x3 array, main
% beam ~22 dBi (~15 deg), ~0-3 dBi to the rest of the sky, assumed -10 dBi to the horizon;
% T_sys 70 K -> receiver noise -124 dBW in 400 MHz. GNSS: ICD minimum received power
% + ~3 dB, summed over the visible satellites (L1 ~30: many satellites add up to Gaussian
% noise -> 'noise'; the weaker bands as one bpsk each). Other levels are link-budget
% estimates (+-10 dB). Radar 1 at +40 dB = some tens of km away with terrain in between;
% a nearby radar is much stronger and would clip a real ADC (not modelled). Frequencies
% verified 8 Oct (notes §10); radars and carrier off the channel grid. LTE: 90 % of each
% channel bandwidth occupied. Expected partly covered channels: 81 86 89 94 (LTE),
% 105 106 109 111 112 113 (Inmarsat), 120 121 (L1).
rfiRealistic = [ ...
    rfiSource('noise',   'Freq', 1575.42e6, 'Bandwidth', 2.046e6, 'INRdB', -14, 'Label', 'GNSS L1 band (~30 satellites)'), ...
    rfiSource('bpsk',    'Freq', 1561.098e6, 'ChipRate', 2.046e6, 'INRdB', -23, 'Label', 'BeiDou B1I'), ...
    rfiSource('bpsk',    'Freq', 1278.75e6, 'ChipRate', 5.115e6, 'INRdB', -17, 'Label', 'Galileo E6'), ...
    rfiSource('bpsk',    'Freq', 1268.52e6, 'ChipRate', 10.23e6, 'INRdB', -23, 'Label', 'BeiDou B3I'), ...
    rfiSource('bpsk',    'Freq', 1246e6,    'ChipRate', 5.11e6,  'INRdB', -29, 'Label', 'GLONASS G2'), ...
    rfiSource('bpsk',    'Freq', 1227.60e6, 'ChipRate', 10.23e6, 'INRdB', -22, 'Label', 'GPS L2'), ...
    rfiSource('pulsed',  'Freq', 1332.9e6, 'PulseWidth', 2e-6, 'PRF', 973, 'ChirpBW', 1e6, 'INRdB', 40, ...
              'ScanPeriod', 10, 'BeamTime', L/2, 'BeamWidth', 39e-3, 'SidelobeDB', -30, ...
              'Label', 'en-route radar, beam passage'), ...       % 1.4 deg beam, 6 rpm
    rfiSource('pulsed',  'Freq', 1257.3e6, 'PulseWidth', 2e-6, 'PRF', 351, 'ChirpBW', 1e6, 'INRdB', 20, ...
              'StartTime', 0.7e-3, 'ScanPeriod', 12, 'BeamTime', 6, 'BeamWidth', 47e-3, ...
              'SidelobeDB', -35, 'Label', 'radar 2, sidelobes only'), ...  % beam half a scan away
    rfiSource('noise',   'Freq', 1459.5e6, 'Bandwidth', 13.5e6, 'INRdB', 0,   'Label', 'LTE band 32, operator 1'), ...
    rfiSource('noise',   'Freq', 1472e6,   'Bandwidth', 9e6,    'INRdB', -6,  'Label', 'LTE band 32, operator 2'), ...
    rfiSource('noise',   'Freq', 1484.5e6, 'Bandwidth', 13.5e6, 'INRdB', -3,  'Label', 'LTE band 32, operator 3'), ...
    rfiSource('noise',   'Freq', 1528.0e6, 'Bandwidth', 4e6,    'INRdB', -17, 'Label', 'Inmarsat block 1'), ...
    rfiSource('noise',   'Freq', 1541.4e6, 'Bandwidth', 6.4e6,  'INRdB', -20, 'Label', 'Inmarsat block 2'), ...
    rfiSource('noise',   'Freq', 1550.5e6, 'Bandwidth', 2.8e6,  'INRdB', -23, 'Label', 'Inmarsat block 3'), ...
    rfiSource('impulse', 'Rate', 50, 'Duration', 200e-9, 'INRdB', 20, 'Label', 'broadband impulses'), ...
    rfiSource('cw',      'Freq', 1388.9e6, 'INRdB', -15, 'Label', 'local spurious carrier')];
switch rfiList
    case 'illustrative'
    case 'realistic', rfi = rfiRealistic;
    otherwise, error('pipelineParams:rfiList', 'rfiList must be ''illustrative'' or ''realistic''.');
end
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
periodicMask  = false;       % with excision: find periodic RFI (radars) in the detections and
                             % blank all its predicted pulses, also those too weak to be seen
                             % (periodicRFI, B6; 8 Oct); false = detectRFI's mask only.
                             % Default off since 9 Oct (Jasper): with dense false flags it
                             % accepts false emitters and blanked ~30 % of the data (B5e);
                             % only needed for radars locked to the pulsar. State and fix
                             % design: logs/2026-10-09_handover-B-closed.md
periodicArgs  = {};          % name-value options for periodicRFI

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
