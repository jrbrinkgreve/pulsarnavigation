function fig = plotDispersionCheck(info_gen, info_disp, tStart, tSpan, opts)
%PLOTDISPERSIONCHECK  Side-by-side view of the raw and dispersed signals.
%{
Reads the same time window from the generator output and the dispersed
file and shows, per column (left: raw, right: dispersed):

  top    - dynamic spectrum (power vs time and frequency). A raw pulse is
           a vertical stripe; after dispersion it becomes a curve with LOW
           frequencies arriving LATER. The expected curve
           t = tc + tau(f) - tau(fRef) is overlaid (dashed), so the sign
           and size of the delay can be checked by eye.
  bottom - band-integrated power (sum over fLow..fHigh channels). The
           dispersed pulse is spread over the smear time, with lower peak.

  plotDispersionCheck(info_gen, info_disp)                  % first 20 ms
  plotDispersionCheck(info_gen, info_disp, tStart, tSpan)
  plotDispersionCheck(..., 'ChanBW', 2e6, 'TimeRes', 20e-6)

Inputs:
  info_gen   info struct from generatePulsarSignal
  info_disp  info struct from applyDispersionStream
  tStart     [s] window start (default 0)
  tSpan      [s] window length (default 20e-3)

Name-value options:
  'ChanBW'      [Hz] frequency resolution of the dynamic spectrum (default 2e6)
  'TimeRes'     [s]  time resolution after averaging (default 20e-6)
  'FreqMargin'  fraction of the band shown outside [fLow, fHigh] on each
                side, to show the band-limiting (default 0.25)
  'MaxSamples'  refuse windows longer than this many samples per file
                (default 2e8, ~0.8 GB single each)
%}

arguments
    info_gen  struct
    info_disp struct
    tStart    (1,1) double {mustBeNonnegative} = 0
    tSpan     (1,1) double {mustBePositive}    = 20e-3
    opts.ChanBW     (1,1) double {mustBePositive}    = 2e6
    opts.TimeRes    (1,1) double {mustBePositive}    = 20e-6
    opts.FreqMargin (1,1) double {mustBeNonnegative} = 0.25
    opts.MaxSamples (1,1) double {mustBePositive}    = 2e8
end

fs = info_gen.fs;
if abs(info_disp.fs - fs) > 1e-6*fs
    error('plotDispersionCheck:fs', 'Sample rates of the two files differ.');
end

% ---- Window in samples ------------------------------------------------------
i0 = round(tStart * fs);                         % 0-based first sample
n  = min(round(tSpan * fs), info_gen.N - i0);
if n <= 0
    error('plotDispersionCheck:window', 'Window starts beyond the end of the file.');
end
if n > opts.MaxSamples
    error('plotDispersionCheck:tooLong', ...
        'Window is %d samples per file; shorten tSpan or raise MaxSamples.', n);
end
t0 = i0 / fs;

xRaw  = readSegment(info_gen.file,  info_gen.byteOrder,  i0, n);
xDisp = readSegment(info_disp.file, info_disp.byteOrder, i0, n);

% ---- Dynamic spectra ----------------------------------------------------------
fLow  = info_disp.fLow;
fHigh = info_disp.fHigh;
B     = fHigh - fLow;
fMin  = max(0,      fLow  - opts.FreqMargin*B);
fMax  = min(fs/2,   fHigh + opts.FreqMargin*B);

frameLen = max(8, round(fs / opts.ChanBW));
nAvg     = max(1, round(opts.TimeRes * fs / frameLen));

[Sraw,  tAx, fAx] = dynSpec(xRaw,  fs, t0, frameLen, nAvg, fMin, fMax);
clear xRaw
[Sdisp, ~,   ~  ] = dynSpec(xDisp, fs, t0, frameLen, nAvg, fMin, fMax);
clear xDisp

inBand = fAx >= fLow & fAx <= fHigh;
Praw   = sum(Sraw(inBand, :), 1);
Pdisp  = sum(Sdisp(inBand, :), 1);

% ---- Expected dispersion sweep ------------------------------------------------
K    = info_disp.dispersionConst;                % s MHz^2 pc^-1 cm^3
DM   = info_disp.DM;
fRef = info_disp.refFreq;
tau  = @(f) K * DM ./ (f/1e6).^2;
dMin = tau(fHigh) - tau(fRef);                   % delay of top of band
dMax = tau(fLow)  - tau(fRef);                   % delay of bottom of band

tEnd = t0 + n/fs;
tc   = info_gen.pulseCenters;
tc   = tc(tc + dMax >= t0 & tc + dMin <= tEnd);  % pulses visible in window
fCurve = linspace(fLow, fHigh, 200);

% ---- Plot --------------------------------------------------------------------------
cMax = quantileSimple([Sraw(:); Sdisp(:)], 0.999);
pMax = 1.05 * max([Praw, Pdisp]);
tMs  = tAx * 1e3;

fig = figure('Name', 'Raw vs dispersed', 'Color', 'w');
tl  = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, sprintf(['DM = %.3g pc cm^{-3}, ref %.3f GHz, smear %.3f ms  ' ...
    '(%.2g MHz channels, %.3g \\mus bins)'], DM, fRef/1e9, ...
    (dMax - dMin)*1e3, fs/frameLen/1e6, nAvg*frameLen/fs*1e6));

% Top-left: raw dynamic spectrum
ax1 = nexttile(tl, 1);
imagesc(ax1, tMs, fAx/1e9, Sraw, [0 cMax]); axis(ax1, 'xy'); hold(ax1, 'on');
for k = 1:numel(tc)
    plot(ax1, [tc(k) tc(k)]*1e3, [fLow fHigh]/1e9, 'w--', 'LineWidth', 1);
end
yline(ax1, [fLow fHigh]/1e9, 'w:');
hold(ax1, 'off');
title(ax1, 'Raw (generator output)'); ylabel(ax1, 'Frequency [GHz]');

% Top-right: dispersed dynamic spectrum with expected sweep
ax2 = nexttile(tl, 2);
imagesc(ax2, tMs, fAx/1e9, Sdisp, [0 cMax]); axis(ax2, 'xy'); hold(ax2, 'on');
for k = 1:numel(tc)
    plot(ax2, (tc(k) + tau(fCurve) - tau(fRef))*1e3, fCurve/1e9, ...
        'w--', 'LineWidth', 1);
end
yline(ax2, [fLow fHigh]/1e9, 'w:');
hold(ax2, 'off');
title(ax2, 'Dispersed (dashed: expected t_c + \tau(f) - \tau(f_{ref}))');
cb = colorbar(ax2); cb.Label.String = 'Power [a.u.]';

% Bottom-left: raw band-integrated power
ax3 = nexttile(tl, 3);
plot(ax3, tMs, Praw, 'k'); hold(ax3, 'on');
for k = 1:numel(tc)
    xline(ax3, tc(k)*1e3, 'r--');
end
hold(ax3, 'off'); grid(ax3, 'on'); ylim(ax3, [0 pMax]);
xlabel(ax3, 'Time [ms]'); ylabel(ax3, 'In-band power [a.u.]');
title(ax3, 'Raw, summed over f_{low}..f_{high}');

% Bottom-right: dispersed band-integrated power with expected sweep span
ax4 = nexttile(tl, 4);
plot(ax4, tMs, Pdisp, 'k'); hold(ax4, 'on');
for k = 1:numel(tc)
    xline(ax4, (tc(k) + dMin)*1e3, 'r--');       % top of band arrives
    xline(ax4, (tc(k) + dMax)*1e3, 'b--');       % bottom of band arrives
end
hold(ax4, 'off'); grid(ax4, 'on'); ylim(ax4, [0 pMax]);
xlabel(ax4, 'Time [ms]');
title(ax4, 'Dispersed (red: f_{high} arrival, blue: f_{low} arrival)');

colormap(fig, 'parula');
linkaxes([ax1 ax2 ax3 ax4], 'x');
linkaxes([ax1 ax2], 'y');
xlim(ax1, [tMs(1) tMs(end)]);
end


% =================================================================================
function x = readSegment(file, byteOrder, i0, n)
[fid, msg] = fopen(file, 'r', byteOrder);
if fid == -1
    error('plotDispersionCheck:open', 'Could not open "%s": %s', file, msg);
end
c = onCleanup(@() fclose(fid)); %#ok<NASGU>
if fseek(fid, i0 * 4, 'bof') ~= 0
    error('plotDispersionCheck:seek', 'Could not seek in "%s".', file);
end
x = fread(fid, [1 n], 'single=>single');
if numel(x) < n
    error('plotDispersionCheck:short', 'Only %d of %d samples read from "%s".', ...
        numel(x), n, file);
end
end


% =================================================================================
function [S, tAx, fAx] = dynSpec(x, fs, t0, frameLen, nAvg, fMin, fMax)
% Hann-windowed, non-overlapping FFT frames; power averaged over nAvg
% frames per output pixel; only bins in [fMin, fMax] kept. Processed in
% chunks so the complex FFT buffer stays around 16M samples.
pix     = frameLen * nAvg;                       % samples per output column
nOut    = floor(numel(x) / pix);
if nOut < 1
    error('plotDispersionCheck:short', 'Window shorter than one time bin.');
end
kLo = ceil(fMin * frameLen / fs);
kHi = floor(fMax * frameLen / fs);
fAx = (kLo:kHi) * fs / frameLen;
nF  = numel(fAx);
w   = single(0.5 - 0.5*cos(2*pi*(0:frameLen-1).'/frameLen));

S      = zeros(nF, nOut, 'single');
perChk = max(1, floor(2^24 / pix));              % output columns per chunk
for c = 1:perChk:nOut
    cEnd = min(c + perChk - 1, nOut);
    seg  = x((c-1)*pix + 1 : cEnd*pix);
    X    = fft(reshape(seg, frameLen, []) .* w); % frames as columns
    P    = abs(X(kLo+1:kHi+1, :)).^2;
    P    = reshape(P, nF, nAvg, []);
    S(:, c:cEnd) = reshape(mean(P, 2), nF, []);
end
tAx = t0 + ((0:nOut-1) + 0.5) * pix / fs;
end


function q = quantileSimple(v, p)
% Quantile without the Statistics Toolbox.
v = sort(v(:));
q = v(max(1, round(p * numel(v))));
if q <= 0, q = max(v); end
if q <= 0, q = 1; end
end