function plotEnvelope(inFile, fs, tStart, tDur, maxPoints)
%PLOTENVELOPE Plot the real-valued envelope/power trace written by
% envelopeReconstruction.m.
%
%   plotEnvelope(inFile, fs, tStart, tDur, maxPoints)
%
% Inputs:
%   inFile    path to the real float32 .dat file from
%             envelopeReconstruction (single channel, NOT interleaved --
%             this is real-valued power/magnitude, not complex I/Q).
%   fs        [Hz] sample rate of inFile (= info.actualFout from
%             envelopeReconstruction, NOT the original fs/f_in of the
%             earlier RF/baseband files -- this file is at the reduced
%             envelope rate).
%   tStart    (opt) [s] start of the window to read. Default 0.
%   tDur      (opt) [s] duration of the window to read. Default: the
%             rest of the file from tStart.
%   maxPoints (opt) cap on points actually drawn. If the requested
%             window has more samples than this, the plotted trace is
%             strided down to maxPoints points -- for DISPLAY ONLY; the
%             read itself is not truncated for analysis, only for the
%             plot. Default 2e6 (MATLAB plots get sluggish well before
%             that many points anyway).
%
% This reads only the requested window (fseek + fread), never the whole
% file, following the same windowed-read pattern used throughout this
% pipeline -- appropriate since even the reduced-rate envelope file can
% still be large for a long observation.

if nargin < 3 || isempty(tStart), tStart = 0; end
if nargin < 5 || isempty(maxPoints), maxPoints = 2e6; end

fInfo = dir(inFile);
if isempty(fInfo)
    error('plotEnvelope: inFile "%s" not found.', inFile);
end
Nx = floor(fInfo.bytes / 4); % float32, real, single channel

offsetSamples = round(tStart * fs);
if offsetSamples < 0 || offsetSamples >= Nx
    error('plotEnvelope: tStart (%.6g s) is outside the file (%.6g s total).', ...
        tStart, Nx/fs);
end

if nargin < 4 || isempty(tDur)
    nSamples = Nx - offsetSamples; % rest of the file
else
    nSamples = min(round(tDur * fs), Nx - offsetSamples);
end

fid = fopen(inFile, 'r');
if fid == -1
    error('plotEnvelope: could not open inFile "%s".', inFile);
end
fseek(fid, offsetSamples * 4, 'bof');
y = fread(fid, nSamples, 'single=>single');
fclose(fid);

t = (offsetSamples : offsetSamples + numel(y) - 1) / fs;

% Downsample for display only, if the window is large. Simple stride
% (not a proper anti-alias decimation) -- fine for a quick look at
% envelope shape, not intended as an analysis step.
if numel(y) > maxPoints
    strideFactor = ceil(numel(y) / maxPoints);
    tPlot = t(1:strideFactor:end);
    yPlot = y(1:strideFactor:end);
    strideNote = sprintf(' (displayed at 1/%d stride for %d of %d samples)', ...
        strideFactor, numel(yPlot), numel(y));
else
    tPlot = t;
    yPlot = y;
    strideNote = '';
end

figure;
plot(tPlot, yPlot);
xlabel('t [s]');
ylabel('amplitude');
title(sprintf('Envelope: %s%s', strrep(inFile, '_', '\_'), strideNote));
grid on;

end