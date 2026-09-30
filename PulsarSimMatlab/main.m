%{
PULSAR NAVIGATION - synthetic data pipeline (main script)

Main thought: pulsar signals arrive at a known location at very accurate
times. A pulse arriving earlier means we are closer to the pulsar in that
direction; combining several pulsars gives our position.

Pipeline (status):
  Synthetic data generation
    1. generatePulsarSignal    noise-like pulses, ground truth  [done]
    2. applyDispersionStream   interstellar dispersion          [done]
    3. applyIQmodulation       downconversion to baseband       [done]
    -  receiver noise, Earth RFI, clock jitter, pulsar-Earth motion,
       polarisation?, 3x3 array element signals, multiple pulsars  [todo]
  Processing
    (  array beamforming, RFI excision - linear, before squaring  [todo] )
    4. applyInverseDispersion  coherent dedispersion            [done]
    5. detectPower             square-law detection             [done]
    6. foldProfile             folding with a phase model       [done]
    7. estimateTOA             FFT template matching (FFTFIT)   [done, to verify]
    -  noise normalization, NP detector                          [todo]
    -  barycentric / timing corrections, residuals               [todo]
    -  navigation solution (multi-pulsar)                        [todo]
  Validation
    - plot/check functions per stage against ground truth        [done]
    8. validateTOA             TOA vs ground truth              [done, to verify]
    - Monte Carlo over seeds / noise levels                      [todo]

Rule: the PROCESSING stages only use what an observer would know
(ephemeris + receiver settings). Ground truth (info_gen, info_disp) is
only used by the generation stages and by the check/plot functions.
%}

% Run control
runStage.generate = false;    % false: reuse files on disk (info loaded from *_info.mat)
runStage.process  = true;    % false: reuse dedispersed + detected files
runStage.fold     = true;
runStage.toa      = true;     % TOA estimation + validation (needs the fold)

plots.dispersion = true;   % raw vs dispersed (reads the big RF files)
plots.iq         = true;   % dispersed IQ vs dedispersed IQ
plots.detected   = true;
plots.fold       = true;
plots.toa        = true;
closeFigures     = true;

if closeFigures, close all; end

%add paths
scriptDir = fileparts(mfilename('fullpath')); 
addpath(fullfile(scriptDir, 'functions'));


% Parameters
% --- Pulsar (ground truth for generation)
T          = 10e-3;     % [s]  pulsar period
A          = 1;         % [V]  noise std dev at pulse peak
dutycycle  = 5;         % [%]  pulse FWHM as percent of T
genEnvMode = 'power';   % 'power': power profile has FWHM = dutycycle% of T
DM         = 5;         % [pc cm^-3] dispersion measure

% --- Simulation
f_in = 4e9;             % [Hz] RF generation sampling rate
L    = 100e-3;         % [s]  total signal length
seed = 42;              % RNG seed for reproducible runs

% --- Receiver (known to the observer)
fLow        = 1.2e9;    % [Hz] lower edge of observed RF band
fHigh       = 1.6e9;    % [Hz] upper edge of observed RF band
fLO         = 1.4e9;    % [Hz] local oscillator
fs          = 500e6;    % [Hz] IQ sample rate
filterOrder = 2048;     % IQ low-pass order

% --- Processing
refFreq       = fHigh;  % [Hz] dedispersion reference frequency (same as generation)
f_out         = 1e6;   % [Hz] detected-power bin rate
nBin          = 1024;   % phase bins in the fold
subintPeriods = 1;      % turns per sub-integration

% --- Files
dataDir       = "data";
fileRaw       = fullfile(dataDir, "test.dat");
fileDispersed = fullfile(dataDir, "test_dispersed.dat");
fileIQ        = fullfile(dataDir, "test_dispersed_IQ_modulated.dat");
fileDedisp    = fullfile(dataDir, "test_IQ_dedispersed.dat");
fileEnvelope  = fullfile(dataDir, "test_envelope.dat");
fileFold      = fullfile(dataDir, "test_fold.mat");

% Disk estimate for the streamed files
% raw + dispersed (real float32 at f_in), IQ + dedispersed (complex at fs), power
diskGB = L * (2*f_in*4 + 2*fs*8 + f_out*4) / 1e9;
fprintf('main: L = %.3g s -> about %.1f GB of data files in "%s"\n', L, diskGB, dataDir);

% Synthetic data generation
if runStage.generate
    % 1. Pulsar signal
    info_gen = generatePulsarSignal(fileRaw, T, f_in, A, L, dutycycle, ...
        'Seed', seed, 'EnvelopeMode', genEnvMode, 'Verbose', true);

    % 2. Interstellar dispersion
    info_disp = applyDispersionStream(info_gen.file, fileDispersed, DM, ...
        info_gen.actualFsOut, fLow, fHigh, 'RefFreq', refFreq);

    % 3. IQ downconversion
    info_IQ = applyIQmodulation(info_disp.file, fileIQ, info_disp.actualFsOut, fs, fLO, ...
        'FilterOrder', filterOrder, 'Band', [fLow fHigh]);
else
    info_gen  = loadInfo(fileRaw);
    info_disp = loadInfo(fileDispersed);
    info_IQ   = loadInfo(fileIQ);
    checkConsistency(info_gen, info_disp, info_IQ, T, f_in, L, DM, fLow, fHigh, fLO, fs);
end

if plots.dispersion
    plotDispersionCheck(info_gen, info_disp, 0, 20e-3);
end

% Observer knowledge (ephemeris) - in real use from a pulsar catalogue
% Synthetic case: spin frequency from T; phase 0 at the first pulse centre
% (T/2), referenced to refFreq.
ephem.f0   = 1/T;
ephem.F1   = 0;
ephem.TRef = 0.5 * T;
ephem.DM   = DM;
ephem.profileFWHM = dutycycle/100;   % [turns] power-profile FWHM, for the template
                                     % (genEnvMode 'amplitude' would need /sqrt(2))






% Processing
if runStage.process
    % 4. Coherent dedispersion
    info_dedisp = applyInverseDispersion(info_IQ.file, fileDedisp, info_IQ.actualFsOut, ...
        info_IQ.fLO, ephem.DM, fLow, fHigh, 'RefFreq', refFreq);

    % 5. Square-law detection
    info_det = detectPower(info_dedisp.file, fileEnvelope, info_dedisp.fs, info_dedisp.fLO, ...
        f_out, 'FullySupported', info_dedisp.fullySupported);
else
    info_dedisp = loadInfo(fileDedisp);
    info_det    = loadInfo(fileEnvelope);
end

if plots.iq
    [~, check_iq] = plotIQCheck(info_gen, info_disp, info_IQ, info_dedisp, 0, 20e-3);
end
if plots.detected
    [~, check_det] = plotDetectedPower(info_det, info_gen, 0, 20e-3, ...
        'InfoDisp', info_disp, 'InfoIQ', info_IQ, 'InfoDedisp', info_dedisp);
end

% 6. Folding
if runStage.fold
    [info_fold, fold] = foldProfile(info_det, fileFold, ephem.f0, ...
        'F1', ephem.F1, 'TRef', ephem.TRef, 'NBin', nBin, 'SubintPeriods', subintPeriods);
    if plots.fold
        [~, check_fold] = plotFoldCheck(fold, info_fold, info_gen, ...
            'InfoDisp', info_disp, 'InfoIQ', info_IQ, 'InfoDedisp', info_dedisp);
    end
elseif runStage.toa
    S = load(fileFold, 'fold', 'info');                 % reuse the saved fold
    fold = S.fold; info_fold = S.info;
end

% 7. TOA estimation
if runStage.toa
    template = gaussianTemplate(nBin, ephem.profileFWHM);
    % Noise-equivalent bandwidth of the detected band. The observer knows the
    % (calibrated) receiver bandpass and its own dedispersion taper; here both
    % tapers are taken from the stage info.
    Bnoise = noiseBandwidth(fLow, fHigh, [info_disp.edgeWidth, info_dedisp.edgeWidth]);
    [toa, info_toa] = estimateTOA(fold, info_fold, template, 'Bnoise', Bnoise);

    % 8. Validation against ground truth
    val = validateTOA(toa, info_gen, 'Plot', plots.toa);
end


% ======================================================================
%  Local functions
%  ======================================================================
function info = loadInfo(dataFile)
%LOADINFO  Load the <name>_info.mat saved next to a stage's output file.
[d, name] = fileparts(dataFile);
f = fullfile(d, name + "_info.mat");
if ~isfile(f)
    error('main:noInfo', 'No info file "%s"; run the stage that creates it first.', f);
end
s = load(f, 'info');
info = s.info;
end


function tmpl = gaussianTemplate(nBin, fwhmTurns)
%GAUSSIANTEMPLATE  Periodic Gaussian profile, peak 1 at phase 0 (bin 1).
ph  = (0:nBin-1).' / nBin;
sig = fwhmTurns / (2*sqrt(2*log(2)));
tmpl = zeros(nBin, 1);
for j = -2:2                                   % periodic images
    tmpl = tmpl + exp(-0.5*((ph - j)/sig).^2);
end
tmpl = tmpl / max(tmpl);
end


function Bn = noiseBandwidth(fLow, fHigh, edgeWidths)
%NOISEBANDWIDTH  (int W^2)^2 / int W^4 for the product of raised-cosine
% band-edge tapers (each of the given width, inside [fLow, fHigh]).
f = linspace(fLow, fHigh, 200001);
W = ones(size(f));
for e = edgeWidths
    w = ones(size(f));
    lo = f < fLow + e;   w(lo) = sin(pi/2 * (f(lo) - fLow) / e).^2;
    hi = f > fHigh - e;  w(hi) = sin(pi/2 * (fHigh - f(hi)) / e).^2;
    W = W .* w;
end
Bn = trapz(f, W.^2)^2 / trapz(f, W.^4);
end


function checkConsistency(info_gen, info_disp, info_IQ, T, f_in, L, DM, fLow, fHigh, fLO, fs)
%CHECKCONSISTENCY  Warn when reused files were made with other parameters.
d = {};
if abs(info_gen.T - T) > 1e-12,       d{end+1} = 'T';    end
if abs(info_gen.fs - f_in) > 1e-3,    d{end+1} = 'f_in'; end
if abs(info_gen.L - L) > 1e-12,       d{end+1} = 'L';    end
if abs(info_disp.DM - DM) > 1e-12,    d{end+1} = 'DM';   end
if abs(info_disp.fLow - fLow) > 1e-3 || abs(info_disp.fHigh - fHigh) > 1e-3
    d{end+1} = 'band';
end
if abs(info_IQ.fLO - fLO) > 1e-3,     d{end+1} = 'fLO';  end
if abs(info_IQ.actualFsOut - f_in/round(f_in/fs)) > 1e-3
    d{end+1} = 'fs';
end
if ~isempty(d)
    warning('main:stale', ['Files on disk were generated with different ' ...
        'parameters (%s). Set runStage.generate = true to regenerate.'], strjoin(d, ', '));
end
end