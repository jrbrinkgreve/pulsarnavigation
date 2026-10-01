%{
PULSAR NAVIGATION - synthetic data pipeline (main script)

Main thought: pulsar signals arrive at a known location at very accurate
times. A pulse arriving earlier means we are closer to the pulsar in that
direction; combining several pulsars gives our position.

Pipeline (status):
  Sky (synthetic)
    1. generatePulsarSignal    noise-like pulses, ground truth  [done]
    2. applyDispersionStream   interstellar dispersion          [done]
  Receiver (synthetic)
    3. addNoiseAndRFI          receiver noise + RFI at RF       [new, to verify]
    4. applyIQmodulation       downconversion to baseband       [done]
    -  clock jitter, pulsar-Earth motion, polarisation?,
       3x3 array element signals, multiple pulsars               [todo]
  Processing
    (  array beamforming, RFI excision - linear, before squaring  [todo] )
    5. applyInverseDispersion  coherent dedispersion            [done]
    6. detectPower             square-law detection             [done]
    7. foldProfile             folding with a phase model       [done]
    8. estimateTOA             FFT template matching (FFTFIT)   [done]
    -  noise normalization, NP detector                          [todo]
    -  barycentric / timing corrections, residuals               [todo]
    -  navigation solution (multi-pulsar)                        [todo]
  Validation
    - plot/check functions per stage against ground truth        [done]
    9. validateTOA             TOA vs ground truth              [done]
    - Monte Carlo over seeds / SNR                               [todo]

Rule: the PROCESSING stages only use what an observer would know
(ephemeris + receiver settings). Ground truth (info_gen, info_disp,
info_rx) is only used by the synthetic stages and the check/plot functions.
%}

% Run control
runStage.sky      = false;   % pulsar signal + dispersion (the slow part; independent of noise)
runStage.receiver = true;   % receiver noise/RFI + IQ (rerun this alone to change SNR or RFI)
runStage.process  = true;   % dedispersion + detection
runStage.fold     = true;
runStage.toa      = true;   % TOA estimation + validation (needs the fold)

plots.dispersion = false;   % raw vs dispersed (reads the big RF files)
plots.iq         = false;   % dispersed IQ vs dedispersed IQ
plots.detected   = true;
plots.fold       = true;
plots.toa        = true;
closeFigures     = true;

if closeFigures, close all; end

% Paths
scriptDir = fileparts(mfilename('fullpath'));
addpath(fullfile(scriptDir, 'functions'));

% Parameters
% --- Pulsar (ground truth for generation)
T          = 10e-3;     % [s]  pulsar period
A          = 1;         % [V]  pulsar noise std at pulse peak (sets the signal scale)
dutycycle  = 5;         % [%]  pulse FWHM as percent of T
genEnvMode = 'power';   % 'power': power profile has FWHM = dutycycle% of T
DM         = 5;         % [pc cm^-3] dispersion measure

% --- Simulation
f_in = 4e9;             % [Hz] RF sampling rate
L    = 100e-3;         % [s]  total signal length
seed = 42;              % RNG seed of the pulsar signal

% --- Receiver (known to the observer)
fLow        = 1.2e9;    % [Hz] lower edge of observed RF band
fHigh       = 1.6e9;    % [Hz] upper edge of observed RF band
fLO         = 1.4e9;    % [Hz] local oscillator
fs          = 500e6;    % [Hz] IQ sample rate (D = 8)
filterOrder = 2048;     % IQ low-pass order

% --- Receiver noise and interference (added at RF, after dispersion)
% snrDB: pulse-peak signal PSD / receiver-noise PSD in the band (S_peak/SEFD).
% -20 dB -> ~3.8 SNR per pulse (0.5 ms pulses, 400 MHz); Inf -> no noise.
% addNoiseAndRFI prints what a value means (SNR per pulse / folded, TOA error).
snrDB     = -20;
noiseSeed = seed + 1;   % different seed -> new noise, same pulsar realization
rfiOn     = false;      % validate noise-only first, then switch RFI on

% Illustrative L-band RFI scenario; INR = power relative to all receiver
% noise in the 400 MHz band (while the source is on).
rfi = [ ...
    rfiSource('bpsk',    'Freq', 1575.42e6, 'ChipRate', 1.023e6, 'INRdB', -10, 'Label', 'GNSS L1 C/A-like'), ...
    rfiSource('bpsk',    'Freq', 1227.60e6, 'ChipRate', 10.23e6, 'INRdB', -15, 'Label', 'GNSS L2 P(Y)-like'), ...
    rfiSource('pulsed',  'Freq', 1300e6, 'PulseWidth', 2e-6, 'PRF', 373, 'ChirpBW', 1e6, ...
              'INRdB', 20, 'Label', 'L-band radar'), ...
    rfiSource('cw',      'Freq', 1350e6, 'INRdB', -5, 'Label', 'spurious carrier'), ...
    rfiSource('impulse', 'Rate', 50, 'Duration', 200e-9, 'INRdB', 20, 'Label', 'broadband impulses')];
if ~rfiOn || isinf(snrDB), rfi = struct([]); end

% --- Processing
refFreq       = fHigh;  % [Hz] dedispersion reference frequency (same as generation)
f_out         = 1e6;    % [Hz] detected-power bin rate
nBin          = 1024;   % phase bins in the fold
subintPeriods = 10;     % turns per sub-integration (aim for >~10 SNR per sub-int)

% --- Files
dataDir       = "data";
fileRaw       = fullfile(dataDir, "test.dat");
fileDispersed = fullfile(dataDir, "test_dispersed.dat");
fileRx        = fullfile(dataDir, "test_rx.dat");
fileIQ        = fullfile(dataDir, "test_rx_IQ.dat");
fileDedisp    = fullfile(dataDir, "test_IQ_dedispersed.dat");
fileEnvelope  = fullfile(dataDir, "test_envelope.dat");
fileFold      = fullfile(dataDir, "test_fold.mat");

% Disk estimate for the streamed files
% raw + dispersed + receiver (real float32 at f_in), IQ + dedispersed (complex at fs), power
diskGB = L * (3*f_in*4 + 2*fs*8 + f_out*4) / 1e9;
fprintf('main: L = %.3g s -> about %.1f GB of data files in "%s"\n', L, diskGB, dataDir);

% Sky: pulsar signal and interstellar dispersion
if runStage.sky
    info_gen = generatePulsarSignal(fileRaw, T, f_in, A, L, dutycycle, ...
        'Seed', seed, 'EnvelopeMode', genEnvMode, 'Verbose', true);
    info_disp = applyDispersionStream(info_gen.file, fileDispersed, DM, ...
        info_gen.actualFsOut, fLow, fHigh, 'RefFreq', refFreq);
else
    info_gen  = loadInfo(fileRaw);
    info_disp = loadInfo(fileDispersed);
end

if plots.dispersion
    plotDispersionCheck(info_gen, info_disp, 0, 20e-3);
end

% Receiver: noise + interference at RF (not dispersed), then IQ conversion
if runStage.receiver
    info_rx = addNoiseAndRFI(info_disp.file, fileRx, info_disp.actualFsOut, ...
        'Band', [fLow fHigh], 'SNRdB', snrDB, 'SignalInfo', info_gen, ...
        'RFI', rfi, 'Seed', noiseSeed);
    info_IQ = applyIQmodulation(info_rx.file, fileIQ, info_rx.actualFsOut, fs, fLO, ...
        'FilterOrder', filterOrder, 'Band', [fLow fHigh]);
else
    info_rx = loadInfo(fileRx);
    info_IQ = loadInfo(fileIQ);
end
checkConsistency(info_gen, info_disp, info_rx, info_IQ, T, f_in, L, DM, fLow, fHigh, fLO, fs, snrDB);

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
    % Coherent dedispersion
    info_dedisp = applyInverseDispersion(info_IQ.file, fileDedisp, info_IQ.actualFsOut, ...
        info_IQ.fLO, ephem.DM, fLow, fHigh, 'RefFreq', refFreq);
    % Square-law detection
    info_det = detectPower(info_dedisp.file, fileEnvelope, info_dedisp.fs, info_dedisp.fLO, ...
        f_out, 'FullySupported', info_dedisp.fullySupported);
else
    info_dedisp = loadInfo(fileDedisp);
    info_det    = loadInfo(fileEnvelope);
end

truthArgs = {'InfoDisp', info_disp, 'InfoIQ', info_IQ, 'InfoDedisp', info_dedisp, 'InfoRx', info_rx};
if plots.iq
    [~, check_iq] = plotIQCheck(info_gen, info_disp, info_IQ, info_dedisp, 0, 20e-3, 'InfoRx', info_rx);
end
if plots.detected
    [~, check_det] = plotDetectedPower(info_det, info_gen, 0, 20e-3, truthArgs{:});
end

% Folding
if runStage.fold
    [info_fold, fold] = foldProfile(info_det, fileFold, ephem.f0, ...
        'F1', ephem.F1, 'TRef', ephem.TRef, 'NBin', nBin, 'SubintPeriods', subintPeriods);
    if plots.fold
        [~, check_fold] = plotFoldCheck(fold, info_fold, info_gen, truthArgs{:});
    end
elseif runStage.toa
    S = load(fileFold, 'fold', 'info');                 % reuse the saved fold
    fold = S.fold; info_fold = S.info;
end

% TOA estimation and validation
if runStage.toa
    template = gaussianTemplate(nBin, ephem.profileFWHM);
    % Noise-equivalent bandwidth for the radiometer noise model, from the
    % observer's own dedispersion taper (receiver noise only passes this one).
    % When the pulsar dominates (high SNR, in-pulse) the true value is ~0.5 %
    % lower, since the simulated signal also passed the dispersion-stage taper.
    Bnoise = noiseBandwidth(fLow, fHigh, info_dedisp.edgeWidth);
    [toa, info_toa] = estimateTOA(fold, info_fold, template, 'Bnoise', Bnoise);

    val = validateTOA(toa, info_gen, 'Plot', plots.toa);
    if isfield(info_rx.prediction, 'toaErrPulse')
        fprintf(['main: predicted best TOA error per sub-int (%d turns) %.3g us, ' ...
                 'SNR per sub-int %.1f; whole file %.3g us\n'], subintPeriods, ...
            info_rx.prediction.toaErrPulse/sqrt(subintPeriods)*1e6, ...
            info_rx.prediction.snrPulse*sqrt(subintPeriods), ...
            info_rx.prediction.toaErrFolded*1e6);
    end
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


function checkConsistency(info_gen, info_disp, info_rx, info_IQ, T, f_in, L, DM, fLow, fHigh, fLO, fs, snrDB)
%CHECKCONSISTENCY  Warn when reused files were made with other parameters.
d = {};
if abs(info_gen.T - T) > 1e-12,       d{end+1} = 'T';    end
if abs(info_gen.fs - f_in) > 1e-3,    d{end+1} = 'f_in'; end
if abs(info_gen.L - L) > 1e-12,       d{end+1} = 'L';    end
if abs(info_disp.DM - DM) > 1e-12,    d{end+1} = 'DM';   end
if abs(info_disp.fLow - fLow) > 1e-3 || abs(info_disp.fHigh - fHigh) > 1e-3
    d{end+1} = 'band';
end
if ~(isequal(info_rx.snrDB, snrDB) || (isinf(snrDB) && info_rx.noiseStd == 0))
    d{end+1} = 'snrDB';
end
if ~strcmp(info_rx.inFile, info_disp.file), d{end+1} = 'receiver input'; end
if abs(info_IQ.fLO - fLO) > 1e-3,     d{end+1} = 'fLO';  end
if abs(info_IQ.actualFsOut - f_in/round(f_in/fs)) > 1e-3
    d{end+1} = 'fs';
end
if ~isempty(d)
    warning('main:stale', ['Files on disk were generated with different ' ...
        'parameters (%s). Rerun the corresponding stage.'], strjoin(d, ', '));
end
end