function testApplyIQmodulation()
%TESTAPPLYIQMODULATION  Unit tests for applyIQmodulation (stage 4: RF to complex baseband).
%{
Run from the PulsarSimMatlab folder: run('tests/testApplyIQmodulation.m'). ~2 s.
Writes its own files to tempdir (~0.1 GB, deleted at the end); needs no data files.
Ground truth is the help text's complex envelope, computed here: x(t) = Re{z(t) e^(i 2 pi fLO t)},
so a tone A cos(2 pi f t + phi) becomes z = A e^(i phi) e^(i 2 pi (f - fLO) t) (no inversion,
factor 2: |z| = A, mean |z|^2 = 2 mean x^2), output sample k at t = k/fs (the group delay
G = 1024 input samples removed, LO phase 0 at input sample 0). The low-pass enters only
through its design limit delta = 10^(-80/20) = 1e-4 (Kaiser, 80 dB): passband ripple and
stopband both <= delta; the test never uses its coefficients. Settings as main.m and
testParams: f_in 4 GHz, fLO 1.3 GHz, fs 800 MHz (D = 5), FilterOrder 2048; 1 ms inputs.

  1. Format and info (4,000,003 input samples, not a multiple of D): N = ceil(N_in/D) =
     800,001 complex samples, 8*N bytes, interleaved I,Q; fs, D, 2049 taps, G = 1024, t0 = 0,
     cutoff 360 MHz, transition width (80 - 8)/(2.285 (M - 1)) f_in/(2 pi), gain factor 2;
     the saved _info.mat = the returned info; main's settings (with 'Band') give no warning.
  2. Tones, each alone (1.2, 1.25, 1.3 = fLO -> DC, 1.4, 1.55, 1.6 GHz; own amplitude and
     phase): every output sample with full filter support equals A e^(i phi) e^(i 2 pi
     (f - fLO) k/fs) within delta*A. Checks the mapping and its sign, the factor 2, the LO
     phase reference and sample k at k/fs at once (one input sample off: 0.16 rad at 100 MHz).
  3. Timing: a Gaussian burst (sigma 0.5 us) on a 1.4 GHz carrier centred at input sample
     5*k0: the energy centroid of |z|^2 is output sample k0 within 1e-3 sample (1.25 ps); the
     filter is linear-phase, so this is exact once G is removed (without: 204.8 samples late).
  4. Power ratio: white noise limited by FFT to 1.205-1.595 GHz and tapered to zero near the
     file ends: D * sum|z|^2 / sum x^2 = 2 within 2*delta (Parseval; no statistics).
  5. Suppression, each <= delta (-80 dB), measured by a Hann-weighted projection onto the
     frequency where it lands: the 1.4 GHz tone's real-signal mirror (at -300 MHz, about
     -134 dB), a tone at 1.68 GHz (+380 MHz, in the stopband, about -92 dB), a tone at
     1.75 GHz (+450 MHz, beyond the output Nyquist; would alias to -350 MHz, about -104 dB).
  6. Block length (the noise of check 4): BlockLen 2048 (= M - 1, the minimum, 1954 blocks),
     99,999 (odd, not a multiple of D), 1,000,003, against the default (one block here):
     the same length, relative rms difference <= 1e-6 (not bit for bit: other FFT sizes).
  7. Warnings and input checks: fs 700 MHz warns :rate (D = 6), Band 1.2-1.7 GHz warns :band,
     Cutoff 398 MHz warns :alias (shown silenced; lastwarn still records them); inFile =
     outFile, fLO >= f_in/2, an unknown Gain, a cutoff outside (0, f_in/2) (both sides), a
     missing file, an empty file and BlockLen < M - 1 are refused, each with its own id.
Not tested: 'Gain' 'unity' (unused; review R-I4 proposes removing it).
Errors at the end if any check fails.
%}

root = fullfile(fileparts(mfilename('fullpath')), '..');
addpath(fullfile(root, 'functions'));
oldDir = cd(root); restoreDir = onCleanup(@() cd(oldDir));
tmp = tempname(tempdir);                         % own folder per run
mkdir(tmp);
cleanTmp = onCleanup(@() rmdir(tmp, 's'));
fails = {};

f_in = 4e9; fs = 800e6; fLO = 1.3e9; band = [1.2e9 1.6e9];   % main.m / testParams
D = 5; M = 2049; G = 1024;                       % decimation, filter length, group delay
delta = 10^(-80/20);                             % Kaiser design: ripple and stopband <= 1e-4
fIn = fullfile(tmp, 'in.dat'); fOut = fullfile(tmp, 'out.dat');
iq = @(varargin) applyIQmodulation(fIn, fOut, f_in, fs, fLO, 'Verbose', false, 'SaveInfo', false, varargin{:});
rs = RandStream('mt19937ar', 'Seed', 9);
Nx = 4e6;                                        % 1 ms
n  = (0:Nx-1).';
N  = Nx / D;
k  = (0:N-1).';
kIn = (ceil(G / D) : floor((Nx - 1 - G) / D)).'; % output samples with full filter support (0-based)

% ---------------------------------------------------------------------------------
% 1. Format and info
% ---------------------------------------------------------------------------------
Nx1 = Nx + 3;
writeF32(fIn, single(randn(rs, Nx1, 1)));
lastwarn('');
i1 = applyIQmodulation(fIn, fOut, f_in, fs, fLO, 'Band', band, 'Verbose', false);
[~, wId1] = lastwarn;
f = dir(fOut);
S = load(fullfile(tmp, 'out_info.mat'));
N1 = ceil(Nx1 / D);
tw = (80 - 8) / (2.285 * (M - 1)) * f_in / (2*pi);
ok = [i1.N == N1, f.bytes == 8*N1, i1.samplesWritten == N1, i1.samplesRead == Nx1, ...
      strcmp(i1.precision, 'single'), strcmp(i1.byteOrder, 'ieee-le'), i1.isComplex, ...
      i1.fs == fs, i1.actualFsOut == fs, i1.fsIn == f_in, i1.decimationFactor == D, i1.fLO == fLO, ...
      i1.filterLen == M, i1.groupDelayRemoved == G, i1.t0 == 0, i1.gainFactor == 2, ...
      strcmp(i1.gainConvention, 'envelope'), abs(i1.cutoff / 360e6 - 1) < 1e-12, ...
      abs(i1.transitionWidth / tw - 1) < 1e-12, i1.stopbandDB == 80, isempty(wId1), ...
      isequal(rmfield(i1, 'infoFile'), S.info)];
pass = all(ok);
fprintf(['1. format and info: %d input -> %d complex samples (%d bytes), D %d, %d taps, G %d, cutoff %.4g MHz, ' ...
         'transition %.4f MHz, no warning at main''s settings %d: %d of %d checks: %s\n'], Nx1, i1.N, f.bytes, ...
    i1.decimationFactor, i1.filterLen, i1.groupDelayRemoved, i1.cutoff / 1e6, i1.transitionWidth / 1e6, ...
    isempty(wId1), nnz(ok), numel(ok), passStr(pass));
if ~pass, fails{end+1} = 'format / info'; end

% ---------------------------------------------------------------------------------
% 2. Tones
% ---------------------------------------------------------------------------------
fT = [1.2 1.25 1.3 1.4 1.55 1.6] * 1e9;
aT = [1.5 0.7 1 2 1.2 0.9];
pT = [0.3 -1.1 2.0 0.7 -2.5 1.4];
eT = zeros(size(fT));
for j = 1:numel(fT)
    writeF32(fIn, single(aT(j) * cos(2*pi * fT(j) * n / f_in + pT(j))));
    iq();
    z = readCF32(fOut);
    ze = aT(j) * exp(1i * (pT(j) + 2*pi * (fT(j) - fLO) * k / fs));
    eT(j) = max(abs(z(kIn + 1) - ze(kIn + 1))) / aT(j);
    if fT(j) == 1.4e9, z14 = z; end              % for the mirror in check 5
end
pass = all(eT < delta);
fprintf(['2. tones at 1.2, 1.25, 1.3 (DC), 1.4, 1.55, 1.6 GHz: largest |z - A e^(i phi) e^(i 2 pi (f - fLO) k/fs)| / A ' ...
         '%s (< %.0e): %s\n'], sprintf('%.1e ', eT), delta, passStr(pass));
if ~pass, fails{end+1} = 'tones'; end

% ---------------------------------------------------------------------------------
% 3. Timing: a Gaussian burst
% ---------------------------------------------------------------------------------
k0 = 400000; c = D * k0; sg = 0.5e-6 * f_in;     % centre (input sample 2,000,000), sigma 2000 samples
writeF32(fIn, single(exp(-0.5 * ((n - c) / sg).^2) .* cos(2*pi * 1.4e9 * (n - c) / f_in)));
iq();
e = abs(readCF32(fOut)).^2;
dT = sum(k .* e) / sum(e) - k0;
pass = abs(dT) < 1e-3;
fprintf(['3. timing: burst centred at input sample %d -> energy centroid at output sample %d %+.1e ' ...
         '(< 1e-3; without removing G it would be %+.1f): %s\n'], c, k0, dT, G / D, passStr(pass));
if ~pass, fails{end+1} = 'timing'; end

% ---------------------------------------------------------------------------------
% 4. Power ratio
% ---------------------------------------------------------------------------------
X  = fft(randn(rs, Nx, 1));
fa = abs(mod((0:Nx-1).' * (f_in / Nx) + f_in/2, f_in) - f_in/2);   % |signed frequency| of each bin
X(fa < 1.205e9 | fa > 1.595e9) = 0;
xb = real(ifft(X));
clear X fa
r  = 1e5;                                        % zero for r samples at each end, then r-sample ramps
w  = zeros(Nx, 1);
u  = (0:r-1).' / r;
w(r+1 : 2*r) = 0.5 * (1 - cos(pi * u));
w(2*r+1 : Nx-2*r) = 1;
w(Nx-2*r+1 : Nx-r) = flipud(w(r+1 : 2*r));
x4 = single(w .* xb);
clear xb w
writeF32(fIn, x4);
iq();
z4 = readCF32(fOut);
ratio = D * sum(abs(z4).^2) / sum(double(x4).^2);
pass = abs(ratio / 2 - 1) < 2*delta;
fprintf('4. power ratio: D sum|z|^2 / sum x^2 = %.6f (2 within %.0e): %s\n', ratio, 2*delta, passStr(pass));
if ~pass, fails{end+1} = 'power ratio'; end

% ---------------------------------------------------------------------------------
% 5. Suppression
% ---------------------------------------------------------------------------------
proj = @(z, fb) abs(sum(hann(numel(kIn)) .* z(kIn + 1) .* exp(-2i*pi * fb * kIn / fs)) / sum(hann(numel(kIn))));
sM = proj(z14, -300e6) / 2;                      % the 1.4 GHz tone (A = 2): mirror at -(1.4 + 1.3) GHz -> -300 MHz
writeF32(fIn, single(cos(2*pi * 1.68e9 * n / f_in))); iq();
sS = proj(readCF32(fOut), 380e6);                % in the stopband, below the output Nyquist
writeF32(fIn, single(cos(2*pi * 1.75e9 * n / f_in))); iq();
sA = proj(readCF32(fOut), 450e6 - fs);           % beyond the Nyquist: aliases to -350 MHz
s5 = [sM sS sA];
pass = all(s5 < delta);
fprintf(['5. suppression (relative to the tone): mirror of the 1.4 GHz tone at -300 MHz %.1f dB, 1.68 GHz ' ...
         '(+380 MHz) %.1f dB, 1.75 GHz aliased to -350 MHz %.1f dB (each < -80 dB): %s\n'], 20*log10(s5), passStr(pass));
if ~pass, fails{end+1} = 'suppression'; end

% ---------------------------------------------------------------------------------
% 6. Block length
% ---------------------------------------------------------------------------------
writeF32(fIn, x4);
iDef = iq();
zDef = readCF32(fOut);
bls = [M - 1, 99999, 1000003];
pass = true;
for b = 1:numel(bls)
    ib = iq('BlockLen', bls(b));
    zb = readCF32(fOut);
    okB = numel(zb) == numel(zDef);
    dR = NaN;
    if okB
        dR = sqrt(mean(abs(zb - zDef).^2) / mean(abs(zDef).^2));
        okB = dR <= 1e-6;
    end
    pass = pass && okB;
    fprintf('6. block length %d: %d blocks, relative rms difference from the default (%d block) %.1e (<= 1e-6): %s\n', ...
        bls(b), ib.nBlocks, iDef.nBlocks, dR, passStr(okB));
end
if ~pass, fails{end+1} = 'block length'; end
clear x4 z4 zDef zb z14 z e

% ---------------------------------------------------------------------------------
% 7. Warnings and input checks
% ---------------------------------------------------------------------------------
writeF32(fIn, single(randn(rs, 10000, 1)));
wIds = {warnId(@() applyIQmodulation(fIn, fOut, f_in, 700e6, fLO, 'Verbose', false, 'SaveInfo', false), ...
               'applyIQmodulation:rate'), ...
        warnId(@() iq('Band', [1.2e9 1.7e9]), 'applyIQmodulation:band'), ...
        warnId(@() iq('Cutoff', 398e6), 'applyIQmodulation:alias')};
wWant = strcat('applyIQmodulation:', {'rate', 'band', 'alias'});
fEmpty = fullfile(tmp, 'empty.dat');
fclose(fopen(fEmpty, 'w'));
bad = @(fi, fo, fl, varargin) errId(@() applyIQmodulation(fi, fo, f_in, fs, fl, 'Verbose', false, ...
    'SaveInfo', false, varargin{:}));
ids = {bad(fIn, fIn, fLO), bad(fIn, fOut, 2.1e9), bad(fIn, fOut, fLO, 'Gain', 'bogus'), ...
       bad(fIn, fOut, fLO, 'Cutoff', 2.5e9), bad(fIn, fOut, fLO, 'Cutoff', -1), ...
       bad(fullfile(tmp, 'none.dat'), fOut, fLO), bad(fEmpty, fOut, fLO), ...
       bad(fIn, fOut, fLO, 'BlockLen', M - 2)};
want = strcat('applyIQmodulation:', {'sameFile', 'lo', 'gain', 'cutoff', 'cutoff', 'notFound', 'empty', 'blockLen'});
pass = isequal(wIds, wWant) && isequal(ids, want);
fprintf(['7. warnings and input checks: %d of %d warnings with the right id (fs 700 MHz, Band to 1.7 GHz, ' ...
         'Cutoff 398 MHz); %d of %d bad inputs refused with the right id: %s\n'], nnz(strcmp(wIds, wWant)), ...
    numel(wWant), nnz(strcmp(ids, want)), numel(want), passStr(pass));
if ~pass, fails{end+1} = 'warnings / input checks'; end

% ---------------------------------------------------------------------------------
if isempty(fails)
    fprintf('testApplyIQmodulation: ALL PASSED\n');
else
    error('testApplyIQmodulation:failed', 'FAILED: %s', strjoin(fails, ', '));
end
end


% =================================================================================
function writeF32(f, x)
fid = fopen(f, 'w', 'ieee-le'); fwrite(fid, x, 'single'); fclose(fid);
end

function z = readCF32(f)
% Interleaved I,Q float32 -> complex double column.
fid = fopen(f, 'r', 'ieee-le'); v = fread(fid, [2 Inf], 'single=>double'); fclose(fid);
z = complex(v(1, :), v(2, :)).';
end

function w = hann(L)
% Hann window, column (base MATLAB; the Signal Processing Toolbox is not needed).
w = 0.5 * (1 - cos(2*pi * (0:L-1).' / (L - 1)));
end

function id = warnId(f, expected)
% The identifier of the last warning f() raises ('' if none), with the expected one silenced
% (lastwarn still records a disabled warning).
s = warning('off', expected);
restore = onCleanup(@() warning(s));
lastwarn('');
f();
[~, id] = lastwarn;
end

function id = errId(f)
% The identifier of the error f() raises ('' if none).
id = '';
try
    f();
catch e
    id = e.identifier;
end
end

function s = passStr(p)
if p, s = 'PASS'; else, s = 'FAIL'; end
end
