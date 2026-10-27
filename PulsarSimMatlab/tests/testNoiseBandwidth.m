function testNoiseBandwidth()
%TESTNOISEBANDWIDTH  Unit tests for noiseBandwidth (the noise-equivalent bandwidth of the band tapers).
%{
Run from the PulsarSimMatlab folder: run('tests/testNoiseBandwidth.m'). < 1 s. Writes no files.
Ground truth: B_n = (int W^2 df)^2 / int W^4 df in closed form. W is a product of sin^2
band-edge tapers of width e inside a band B; with int_0^1 sin^(2n)(pi u/2) du = C(2n, n)/4^n,
m equal tapers give int W^2 = B - 2e + 2e C(4m, 2m)/4^(2m) and int W^4 = B - 2e +
2e C(8m, 4m)/4^(4m). One 8 MHz taper in 400 MHz: (B - 5e/4)^2 / (B - 93e/64) = 391.63180 MHz.
processIQ, foldToTOA and dedisperseChannels take their radiometer noise from B_n; its measured
meaning (bin variance 1/(B_n dt)) is testDetectPower test 4.

  1. Closed forms: one and two equal tapers, for the full band (1.2-1.6 GHz, 8 MHz edges) and
     for one channel (-dF/2..dF/2, dF 3.125 MHz, edges 2 % = 62.5 kHz), each within 1e-7.
  2. Scale invariance: B_n/B depends only on e/B, so 128 x the channel's B_n = the full band's
     (dedisperseChannels' BnoiseTotal relies on it), within 1e-7.
  3. Unequal widths (8 and 2 MHz; no simple closed form): against adaptive quadrature of the
     product taper (waypoints at the edge kinks) within 1e-7; the order of the widths does not
     matter.
  4. Edge cases: no taper gives B; the widest taper, e = B/2 (one hump), gives 0.375^2/0.2734375 B
     (closed form); the same channel at baseband (-dF/2..dF/2) and at RF (around 1.4 GHz) agree.
Not checked: edge widths above B/2 (the two edge regions then overlap and the second overwrites
the first; no caller does that) or <= 0. Errors at the end if any check fails.
%}

root = fullfile(fileparts(mfilename('fullpath')), '..');
addpath(fullfile(root, 'functions'));
fails = {};

fLow = 1.2e9; fHigh = 1.6e9; B = fHigh - fLow; e = 0.02 * B;    % 8 MHz (main.m's tapers)
dF = 3.125e6; eC = 0.02 * dF;                                    % one channel, 62.5 kHz edges
I  = @(n) nchoosek(2*n, n) / 4^n;                                % int_0^1 sin^(2n)(pi u/2) du
Bn = @(B, e, m) (B - 2*e + 2*e*I(2*m))^2 / (B - 2*e + 2*e*I(4*m));   % m equal tapers, closed form

% ---------------------------------------------------------------------------------
% 1. Closed forms
% ---------------------------------------------------------------------------------
cases1 = {'full band, one taper',   noiseBandwidth(fLow, fHigh, e),        Bn(B, e, 1)
          'full band, two tapers',  noiseBandwidth(fLow, fHigh, [e e]),    Bn(B, e, 2)
          'channel, one taper',     noiseBandwidth(-dF/2, dF/2, eC),       Bn(dF, eC, 1)
          'channel, two tapers',    noiseBandwidth(-dF/2, dF/2, [eC eC]),  Bn(dF, eC, 2)};
pass = true;
for c = 1:size(cases1, 1)
    [label, got, want] = cases1{c, :};
    d = got / want - 1;
    pass = pass && abs(d) < 1e-7;
    fprintf('1. closed form, %s: %.6f MHz vs %.6f MHz (%+.1e, < 1e-7)\n', label, got / 1e6, want / 1e6, d);
end
fprintf('1. closed forms: %s\n', passStr(pass));
if ~pass, fails{end+1} = 'closed forms'; end

% ---------------------------------------------------------------------------------
% 2. Scale invariance
% ---------------------------------------------------------------------------------
d2 = 128 * cases1{3, 2} / cases1{1, 2} - 1;
pass = abs(d2) < 1e-7;
fprintf('2. scale invariance: 128 x channel %.6f MHz = %.6f MHz vs full band %.6f MHz (%+.1e, < 1e-7): %s\n', ...
    cases1{3, 2} / 1e6, 128 * cases1{3, 2} / 1e6, cases1{1, 2} / 1e6, d2, passStr(pass));
if ~pass, fails{end+1} = 'scale invariance'; end

% ---------------------------------------------------------------------------------
% 3. Unequal widths
% ---------------------------------------------------------------------------------
e2 = [8e6 2e6];
W  = @(f) taper(f, fLow, fHigh, e2(1)) .* taper(f, fLow, fHigh, e2(2));
q  = @(g) integral(g, fLow, fHigh, 'Waypoints', sort([fLow + e2, fHigh - e2]), 'RelTol', 1e-12, 'AbsTol', 0);
ref3 = q(@(f) W(f).^2)^2 / q(@(f) W(f).^4);
got3 = noiseBandwidth(fLow, fHigh, e2);
rev3 = noiseBandwidth(fLow, fHigh, fliplr(e2));
d3 = got3 / ref3 - 1;
pass = abs(d3) < 1e-7 && abs(rev3 / got3 - 1) < 1e-12;
fprintf(['3. unequal widths (8 and 2 MHz): %.6f MHz vs quadrature %.6f MHz (%+.1e, < 1e-7); widths reversed ' ...
         '%+.1e: %s\n'], got3 / 1e6, ref3 / 1e6, d3, rev3 / got3 - 1, passStr(pass));
if ~pass, fails{end+1} = 'unequal widths'; end

% ---------------------------------------------------------------------------------
% 4. Edge cases
% ---------------------------------------------------------------------------------
dFlat = noiseBandwidth(fLow, fHigh, []) / B - 1;
dHump = noiseBandwidth(fLow, fHigh, B/2) / Bn(B, B/2, 1) - 1;
fc = 1.4e9 + dF/2;                               % a channel of the 128-channel grid, at RF
dRF = noiseBandwidth(fc - dF/2, fc + dF/2, eC) / cases1{3, 2} - 1;
pass = abs(dFlat) < 1e-12 && abs(dHump) < 1e-7 && abs(dRF) < 1e-9;
fprintf(['4. edge cases: no taper / B %+.1e (< 1e-12); e = B/2: %.6f B vs %.6f B (%+.1e, < 1e-7); the channel at ' ...
         '%.4f GHz vs at baseband %+.1e (< 1e-9): %s\n'], dFlat, noiseBandwidth(fLow, fHigh, B/2) / B, ...
    Bn(B, B/2, 1) / B, dHump, fc / 1e9, dRF, passStr(pass));
if ~pass, fails{end+1} = 'edge cases'; end

% ---------------------------------------------------------------------------------
if isempty(fails)
    fprintf('testNoiseBandwidth: ALL PASSED\n');
else
    error('testNoiseBandwidth:failed', 'FAILED: %s', strjoin(fails, ', '));
end
end


% =================================================================================
function w = taper(f, fLow, fHigh, e)
% One sin^2 band-edge taper of width e inside [fLow, fHigh] (the help text's W), any array f.
w = ones(size(f));
lo = f < fLow + e;   w(lo) = sin(pi/2 * (f(lo) - fLow) / e).^2;
hi = f > fHigh - e;  w(hi) = sin(pi/2 * (fHigh - f(hi)) / e).^2;
end

function s = passStr(p)
if p, s = 'PASS'; else, s = 'FAIL'; end
end
