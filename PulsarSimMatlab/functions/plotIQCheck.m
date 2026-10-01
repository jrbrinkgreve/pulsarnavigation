function [fig, check] = plotIQCheck(info_gen, info_disp, info_IQ, info_dedisp, tStart, tSpan, opts)
%PLOTIQCHECK  Visual and numerical check of the IQ and dedispersion stages.
%{
Reads the same time window from the dispersed IQ file and the dedispersed
IQ file and shows (left: dispersed IQ, right: dedispersed IQ):

  row 1  dynamic spectra on an RF frequency axis (f = fLO + f_baseband).
         Left: expected sweep tc + tau(f) - tau(fRef) dashed.
         Right: pulses should be vertical stripes on the true centres.
  row 2  power |z|^2 averaged per time bin. Right panel also shows the
         EXPECTED power profile, computed from the generator ground truth
         (A, envelope, pulse centres) and the band/gain conventions of the
         stages, plus the region that is not fully supported by the
         dedispersion (grey).
  row 3  left: mean spectrum of both files (band edges, tapers, LPF).
         right: I, Q and |z| of the dedispersed signal around one pulse
         peak (raw complex samples: should look like noise).

It also prints, per fully supported pulse, the timing offset of the
measured power centroid relative to the expected one, and the measured /
expected pulse energy. Expected: offsets of a few us or less (noise-
limited), energy ratios ~1 within the noise.

  plotIQCheck(info_gen, info_disp, info_IQ, info_dedisp)          % 0-20 ms
  [fig, check] = plotIQCheck(..., tStart, tSpan, 'TimeRes', 10e-6)

Name-value options:
  'ChanBW'      [Hz] dynamic-spectrum channel width (default 2e6)
  'TimeRes'     [s]  time bin (default 20e-6)
  'ZoomSpan'    [s]  length of the I/Q zoom (default 200e-9)
  'FreqMargin'  fraction of the band shown beyond each edge (default 0.1)
  'MaxSamples'  refuse windows longer than this per file (default 1e8)

Expected power (dedispersed, white generator noise, no receiver noise):
  E|z(t)|^2 = gain^2 * A^2 * p(t) * Beff / fsIn  (+ receiver-noise baseline
  when 'InfoRx' from addNoiseAndRFI is given; see expectedPowerModel)
  p(t) = G(t) ('power' envelope mode) or G(t)^2 ('amplitude' mode),
  Beff = integral over the band of (W_fwd(f) * W_inv(f))^2 df,
  gain = 2 for the 'envelope' IQ convention, fsIn = generator rate.
The IQ low-pass is assumed flat over the band (ripple ~1e-4).
%}

arguments
    info_gen    struct
    info_disp   struct
    info_IQ     struct
    info_dedisp struct
    tStart      (1,1) double {mustBeNonnegative} = 0
    tSpan       (1,1) double {mustBePositive}    = 20e-3
    opts.ChanBW     (1,1) double {mustBePositive}    = 2e6
    opts.TimeRes    (1,1) double {mustBePositive}    = 20e-6
    opts.ZoomSpan   (1,1) double {mustBePositive}    = 200e-9
    opts.FreqMargin (1,1) double {mustBeNonnegative} = 0.1
    opts.MaxSamples (1,1) double {mustBePositive}    = 1e8
    opts.InfoRx     struct = struct([])
end

fs  = info_dedisp.fs;
fLO = info_dedisp.fLO;
if abs(info_IQ.actualFsOut - fs) > 1e-6*fs || abs(info_IQ.fLO - fLO) > 1
    error('plotIQCheck:mismatch', 'IQ and dedispersed files use different fs or fLO.');
end

% ---- Window ------------------------------------------------------------------------
i0 = round(tStart * fs);
n  = min(round(tSpan * fs), info_dedisp.N - i0);
if n <= 0
    error('plotIQCheck:window', 'Window starts beyond the end of the file.');
end
if n > opts.MaxSamples
    error('plotIQCheck:tooLong', ...
        'Window is %d samples per file; shorten tSpan or raise MaxSamples.', n);
end
t0 = i0 / fs;                                     % IQ output sample k at t = k/fs

frameLen = max(8, round(fs / opts.ChanBW));
nAvg     = max(1, round(opts.TimeRes * fs / frameLen));
pix      = frameLen * nAvg;                       % samples per time bin
fLow  = info_dedisp.fLow;
fHigh = info_dedisp.fHigh;
B     = fHigh - fLow;
fMin  = max(fLO - fs/2, fLow  - opts.FreqMargin*B);
fMax  = min(fLO + fs/2, fHigh + opts.FreqMargin*B);

% ---- Read + reduce each file (one at a time to bound memory) ------------------------------
z = readComplex(info_IQ.file, info_IQ.byteOrder, i0, n);
[Sd, tAx, fAx] = dynSpecComplex(z, fs, fLO, t0, frameLen, nAvg, fMin, fMax);
Pd = binPower(z, pix);
clear z

z = readComplex(info_dedisp.file, info_dedisp.byteOrder, i0, n);
Sx = dynSpecComplex(z, fs, fLO, t0, frameLen, nAvg, fMin, fMax);
Px = binPower(z, pix);
zDedisp = z;                                      % kept for the zoom panel
clear z
nBins = min([numel(Pd), numel(Px), numel(tAx)]);
Pd = Pd(1:nBins); Px = Px(1:nBins); tAx = tAx(1:nBins);
Sd = Sd(:, 1:nBins); Sx = Sx(:, 1:nBins);

% ---- Ground truth ------------------------------------------------------------------------------
K    = info_disp.dispersionConst;
DM   = info_disp.DM;
fRef = info_dedisp.refFreq;
tau  = @(f) K * DM ./ (f/1e6).^2;
dMin = tau(fHigh) - tau(fRef);
dMax = tau(fLow)  - tau(fRef);
tEnd = t0 + n/fs;
T    = info_gen.T;
sigG = info_gen.sigma;

% Expected dedispersed power at the bin centres: pulsar + receiver noise
M    = expectedPowerModel(info_gen, info_disp, info_IQ, info_dedisp, opts.InfoRx);
PsX  = M.sigScale * M.envelope(tAx);
Pexp = PsX + M.Pn;

% Fully supported time range of the dedispersed file
fsup = (info_dedisp.fullySupported - 1) / fs;     % [s]

% Pulses to annotate
tcAll = info_gen.pulseCenters;
tcDisp = tcAll(tcAll + dMax >= t0 & tcAll + dMin <= tEnd);
tcIn   = tcAll(tcAll >= t0 & tcAll <= tEnd);

% ---- Numerical per-pulse check -------------------------------------------------------------------
if strcmpi(info_gen.envelopeMode, 'power')
    sigP = sigG;                                  % power profile ~ G
else
    sigP = sigG / sqrt(2);                        % power profile ~ G^2
end
halfW = 4 * sigP;
binDur = pix / fs;
check = struct('tc', {}, 'offset', {}, 'energyRatio', {}, 'offsetNoise', {});
offB = PsX < 1e-6 * M.sigScale;
if any(offB), baseX = median(Px(offB)); else, baseX = 0; end   % measured baseline
fprintf('plotIQCheck: dedispersed pulses vs ground truth (bins %.3g us)\n', binDur*1e6);
for k = 1:numel(tcIn)
    tc = tcIn(k);
    if tc - halfW < max(t0, fsup(1)) || tc + halfW > min(tEnd, fsup(2))
        continue                                  % window not fully usable
    end
    sel = tAx >= tc - halfW & tAx <= tc + halfW;
    Q   = Px(sel) - baseX;
    cm  = sum(tAx(sel) .* Q) / sum(Q);
    ce  = sum(tAx(sel) .* PsX(sel)) / sum(PsX(sel));
    er  = sum(Q) / sum(PsX(sel));
    % Noise-limited centroid scatter from the signal+noise variance per bin
    sdX  = sqrt(M.var(PsX(sel), binDur));
    sOff = sqrt(sum(((tAx(sel) - ce) .* sdX).^2)) / sum(PsX(sel));
    check(end+1) = struct('tc', tc, 'offset', cm - ce, ...
        'energyRatio', er, 'offsetNoise', sOff); %#ok<AGROW>
    fprintf('  pulse at %8.4f ms: centroid offset %+8.3f us (noise ~%.3f us), energy ratio %.4f\n', ...
        tc*1e3, (cm - ce)*1e6, sOff*1e6, er);
end
if isempty(check)
    fprintf('  no fully supported pulse in this window.\n');
end

% ---- Plot -----------------------------------------------------------------------------------------
tMs  = tAx * 1e3;
cMax = quantileSimple([Sd(:); Sx(:)], 0.999);
pMax = 1.05 * max([Pd, Px, Pexp]);

fig = figure('Name', 'IQ / dedispersion check', 'Color', 'w');
tl  = tiledlayout(fig, 3, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, sprintf(['fs = %.4g MHz, LO %.4g GHz, DM %.3g, ref %.3f GHz  ' ...
    '(%.2g MHz channels, %.3g \\mus bins)'], fs/1e6, fLO/1e9, DM, fRef/1e9, ...
    fs/frameLen/1e6, binDur*1e6));

% Row 1: dynamic spectra
ax1 = nexttile(tl, 1);
imagesc(ax1, tMs, fAx/1e9, Sd, [0 cMax]); axis(ax1, 'xy'); hold(ax1, 'on');
fC = linspace(fLow, fHigh, 200);
for k = 1:numel(tcDisp)
    plot(ax1, (tcDisp(k) + tau(fC) - tau(fRef))*1e3, fC/1e9, 'w--');
end
yline(ax1, [fLow fHigh]/1e9, 'w:'); hold(ax1, 'off');
title(ax1, 'Dispersed IQ (dashed: expected sweep)'); ylabel(ax1, 'RF frequency [GHz]');

ax2 = nexttile(tl, 2);
imagesc(ax2, tMs, fAx/1e9, Sx, [0 cMax]); axis(ax2, 'xy'); hold(ax2, 'on');
for k = 1:numel(tcIn)
    plot(ax2, [tcIn(k) tcIn(k)]*1e3, [fLow fHigh]/1e9, 'w--');
end
yline(ax2, [fLow fHigh]/1e9, 'w:'); hold(ax2, 'off');
title(ax2, 'Dedispersed IQ (dashed: true pulse centres)');
cb = colorbar(ax2); cb.Label.String = 'Power [a.u.]';

% Row 2: power profiles
ax3 = nexttile(tl, 3);
plot(ax3, tMs, Pd, 'k'); hold(ax3, 'on');
for k = 1:numel(tcDisp)
    xline(ax3, (tcDisp(k) + dMin)*1e3, 'r--');
    xline(ax3, (tcDisp(k) + dMax)*1e3, 'b--');
end
hold(ax3, 'off'); grid(ax3, 'on'); ylim(ax3, [0 pMax]);
ylabel(ax3, 'mean |z|^2'); title(ax3, 'Dispersed IQ power (red: f_{high}, blue: f_{low} arrival)');

ax4 = nexttile(tl, 4);
hold(ax4, 'on');
shadeUnsupported(ax4, tMs, fsup*1e3, pMax);
plot(ax4, tMs, Px, 'k');
plot(ax4, tMs, Pexp, 'r', 'LineWidth', 1.2);
for k = 1:numel(tcIn)
    xline(ax4, tcIn(k)*1e3, 'r:');
end
hold(ax4, 'off'); grid(ax4, 'on'); box(ax4, 'on'); ylim(ax4, [0 pMax]);
xlabel(ax4, 'Time [ms]');
legend(ax4, {'not fully supported', 'measured', 'expected (ground truth)'}, ...
    'Location', 'northeast');
title(ax4, 'Dedispersed IQ power vs expectation');
xlabel(ax3, 'Time [ms]');

% Row 3 left: mean spectra
ax5 = nexttile(tl, 5);
ref = max(mean(Sx, 2));
plot(ax5, fAx/1e9, 10*log10(mean(Sd, 2)/ref + eps), 'Color', [0.5 0.5 0.5]); hold(ax5, 'on');
plot(ax5, fAx/1e9, 10*log10(mean(Sx, 2)/ref + eps), 'k');
xline(ax5, [fLow fHigh]/1e9, 'r:'); hold(ax5, 'off'); grid(ax5, 'on');
ylim(ax5, [-60 5]); xlabel(ax5, 'RF frequency [GHz]'); ylabel(ax5, 'dB');
legend(ax5, {'dispersed IQ', 'dedispersed'}, 'Location', 'south');
title(ax5, 'Mean spectrum over the window');

% Row 3 right: raw I/Q zoom around a pulse peak
ax6 = nexttile(tl, 6);
tz = pickZoomCentre(tcIn, t0, tEnd, fsup);
nz = max(8, round(opts.ZoomSpan * fs));
iz = round((tz - t0) * fs) + (-floor(nz/2):ceil(nz/2)-1) + 1;
iz = iz(iz >= 1 & iz <= numel(zDedisp));
tz_ns = ((iz - 1)/fs + t0 - tz) * 1e9;
v = zDedisp(iz);
plot(ax6, tz_ns, real(v), 'b', tz_ns, imag(v), 'r', tz_ns, abs(v), 'k', 'LineWidth', 1);
grid(ax6, 'on'); xlabel(ax6, sprintf('Time - %.4f ms [ns]', tz*1e3));
legend(ax6, {'I', 'Q', '|z|'}, 'Location', 'northeast');
title(ax6, 'Dedispersed samples at a pulse peak');
clear zDedisp

colormap(fig, 'parula');
linkaxes([ax1 ax2 ax3 ax4], 'x');
linkaxes([ax1 ax2], 'y');
xlim(ax1, [tMs(1) tMs(end)]);
end


function z = readComplex(file, byteOrder, i0, n)
[fid, msg] = fopen(file, 'r', byteOrder);
if fid == -1
    error('plotIQCheck:open', 'Could not open "%s": %s', file, msg);
end
c = onCleanup(@() fclose(fid)); %#ok<NASGU>
if fseek(fid, i0 * 8, 'bof') ~= 0
    error('plotIQCheck:seek', 'Could not seek in "%s".', file);
end
raw = fread(fid, [2 n], 'single=>single');
if size(raw, 2) < n
    error('plotIQCheck:short', 'Only %d of %d samples read from "%s".', ...
        size(raw, 2), n, file);
end
z = complex(raw(1, :), raw(2, :));
end


function P = binPower(z, pix)
nOut = floor(numel(z) / pix);
P = mean(reshape(abs(z(1:nOut*pix)).^2, pix, nOut), 1);
P = double(P);
end


function [S, tAx, fAx] = dynSpecComplex(z, fs, fLO, t0, frameLen, nAvg, fMin, fMax)
% Hann-windowed complex FFT frames, fftshifted, power averaged over nAvg
% frames per column; RF axis fLO + f_baseband; chunked to bound memory.
pix  = frameLen * nAvg;
nOut = floor(numel(z) / pix);
if nOut < 1
    error('plotIQCheck:short', 'Window shorter than one time bin.');
end
fb   = ((0:frameLen-1) - floor(frameLen/2)) * fs / frameLen;
fRF  = fLO + fb;
keep = fRF >= fMin & fRF <= fMax;
fAx  = fRF(keep);
nF   = nnz(keep);
w    = single(0.5 - 0.5*cos(2*pi*(0:frameLen-1).'/frameLen));
S    = zeros(nF, nOut, 'single');
perChk = max(1, floor(2^23 / pix));
for c = 1:perChk:nOut
    cEnd = min(c + perChk - 1, nOut);
    seg  = z((c-1)*pix + 1 : cEnd*pix);
    X    = fftshift(fft(reshape(seg, frameLen, []) .* w), 1);
    P    = abs(X(keep, :)).^2;
    S(:, c:cEnd) = reshape(mean(reshape(P, nF, nAvg, []), 2), nF, []);
end
tAx = t0 + ((0:nOut-1) + 0.5) * pix / fs;
end


function shadeUnsupported(ax, tMs, supMs, yTop)
x0 = tMs(1); x1 = tMs(end);
c = [0.88 0.88 0.88];
if supMs(1) > x0
    patch(ax, [x0 supMs(1) supMs(1) x0], [0 0 yTop yTop], c, 'EdgeColor', 'none');
else
    patch(ax, nan(1,4), nan(1,4), c, 'EdgeColor', 'none');   % legend entry
end
if supMs(2) < x1
    patch(ax, [supMs(2) x1 x1 supMs(2)], [0 0 yTop yTop], c, ...
        'EdgeColor', 'none', 'HandleVisibility', 'off');
end
end


function tz = pickZoomCentre(tcIn, t0, tEnd, fsup)
ok = tcIn(tcIn >= max(t0, fsup(1)) & tcIn <= min(tEnd, fsup(2)));
if ~isempty(ok)
    tz = ok(1);
elseif ~isempty(tcIn)
    tz = tcIn(1);
else
    tz = (t0 + tEnd) / 2;
end
end


function q = quantileSimple(v, p)
v = sort(v(:));
q = v(max(1, round(p * numel(v))));
if q <= 0, q = max(v); end
if q <= 0, q = 1; end
end