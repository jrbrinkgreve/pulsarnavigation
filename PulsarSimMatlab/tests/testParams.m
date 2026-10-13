%TESTPARAMS  The frozen configuration of the unit tests (script; the tests' pipelineParams).
%{
The unit tests call testParams instead of pipelineParams, so editing pipelineParams (or
running main.m with other settings) does not change what the tests check. These are the
values of the seed-43 reference data the tests were validated on (13 Oct 2026, the same as
pipelineParams then); all files are under data/test/, made once by tests/makeTestData.m
(a copy of main.m's seed-43 files and the frozen references of data/mc/). Do not edit these
values: the bit-for-bit regressions (data/test/mc/*Ref*.mat) belong to this data.

Needs functions/ on the path and MATLAB's current folder = PulsarSimMatlab, like
pipelineParams. Defines the same variable names as pipelineParams (no RFI lists: the tests
define their own sources).
%}

% --- Pulsar (ground truth for generation)
T          = 10e-3;     % [s]  pulsar period
A          = 1;         % [V]  pulsar noise std at pulse peak
dutycycle  = 5;         % [%]  pulse FWHM as percent of T
genEnvMode = 'power';
DM         = 5;         % [pc cm^-3]

% --- Simulation
f_in = 4e9;             % [Hz] RF sampling rate
L    = 100e-3;          % [s]  total signal length
seed = 42;              % pulsar signal

% --- Receiver
fLow        = 1.2e9;
fHigh       = 1.6e9;
fLO         = 1.3e9;
fs          = 800e6;    % [Hz] IQ sample rate
filterOrder = 2048;

% --- Receiver noise and interference
snrDB     = -5;
noiseSeed = seed + 1;   % 43
rfiOn     = false;
rfi       = struct([]);

% --- Processing
refFreq       = fHigh;
f_out         = 1e6;
nBin          = 2048;
subintPeriods = 1;
frontEnd      = 'channels';
weighting     = 'optimal';
chanWidth     = 3.125e6;
excision      = true;
excisionArgs  = {};
periodicMask  = false;
periodicArgs  = {};

% --- Observer knowledge (ephemeris)
ephem.f0   = 1/T;
ephem.F1   = 0;
ephem.TRef = 0.5 * T;
ephem.DM   = DM;
ephem.profileFWHM = dutycycle/100;

% --- Files: the frozen copy (tests/makeTestData.m)
dataDir       = "data/test";
fileRaw       = fullfile(dataDir, "test.dat");             % info only (no data copied)
fileDispersed = fullfile(dataDir, "test_dispersed.dat");   % info only
fileRx        = fullfile(dataDir, "test_rx.dat");          % info only
fileIQ        = fullfile(dataDir, "test_rx_IQ.dat");
fileDedisp    = fullfile(dataDir, "test_IQ_dedispersed.dat");
fileEnvelope  = fullfile(dataDir, "test_envelope.dat");
fileFold      = fullfile(dataDir, "test_fold.mat");
chanDir        = fullfile(dataDir, "chan");                % the tests' channel files
fileChanIQ     = fullfile(chanDir, "test_rx_IQ_chan");
fileChanDedisp = fullfile(chanDir, "test_dedisp_chan");
fileChanPower  = fullfile(chanDir, "test_dedisp_chan_power.dat");
fileFoldChan   = fullfile(dataDir, "test_fold_chan.mat");
% (run() executes this script from tests/, so check the folder from here, not relative)
if ~exist('makingTestData', 'var') && ~isfolder(fullfile(fileparts(mfilename('fullpath')), '..', dataDir))
    error('testParams:noData', 'No frozen test data in %s: run run(''tests/makeTestData.m'') once.', dataDir);
end
