function testGeneratePulsarSignal()
%TESTGENERATEPULSARSIGNAL  Unit tests for generatePulsarSignal (stage 1: the synthetic pulsar).
%{
Run from the PulsarSimMatlab folder: run('tests/testGeneratePulsarSignal.m'). ~15 s.
Writes its own small files to tempdir (deleted at the end); needs no data files.
Ground truth here is independent of the generator's code: the envelope as an exact
periodic sum of 21 Gaussian images (not the generator's nearest pulse + neighbours),
normalized by its exact peak, on the phase grid m/P (f_in*T = P is an integer, so the
sample phases repeat every period). 100 MHz sampling instead of 4 GHz: the generator does
not depend on the rate, and the statistics need many pulses.

  1. Format and ground truth (T = 1/1024 s, f_in = 2^27 Hz, L = 12.5 T: binary-exact, so
     the last sample is exactly on a pulse centre): N = round(L*f_in) + 1, file of 4*N
     bytes, duration, pulseCenters = (k + 1/2)*T for k = 0..12 (13, the last one on the
     last sample), sigma, FWHM_power / FWHM_amplitude, format fields, and the saved
     _info.mat = the returned info.
  2. Power profile (5 % duty, 'power', A = 2, 200 pulses): the folded mean of x^2 in 1 us
     phase bins = A^2*G averaged over the same samples. Exact per-bin sigma (Gaussian x:
     var(x^2) = 2*(A^2*G)^2 per sample) -> normalized residuals on the pulse (G > 1e-3)
     mean 0 and std 1 within 4 sigma; pulse centre (power centroid) and width (second
     moment about phase 1/2) as the model within 4 sigma (linear error propagation);
     measured FWHM printed.
  3. 'amplitude' mode: power profile A^2*G^2, power FWHM = FWHM/sqrt(2); as 2.
  4. The widest real pulses (20 % duty; pulsars stay below ~20 %): as 2, over ALL bins.
     Halfway between pulses the two neighbouring Gaussians overlap (each ~3e-8 of the
     peak), so these bins check the generator's sum over neighbours (without it they would
     be off by a factor 2 there).
  5. Reproducibility: same seed twice identical; BlockSize 333,333 and 7,919 = one block,
     bit for bit (randn drawn in pieces gives the same numbers as in one piece); A = 4 is
     exactly 4 x A = 1 wherever the samples are normal float32 numbers (far off-pulse,
     |x| < 1.2e-38, float32 has fewer bits: differences of a few 1e-45 there); another seed
     differs; the global RNG is not touched; without 'Seed' the reported clock seed
     reproduces the file.
  6. Input checks: a bad EnvelopeMode and bad Seed values (-1, 2^32, 1.5) are refused;
     dutycycle 25 % warns (wider than real pulsars), 20 % does not.
  7. Components (21 Oct 2026): three Gaussian components (amplitudes 0.6 / 0.4 / 0.3, centres
     0 / -0.03 / +0.06 turns, FWHM 0.02 / 0.05 / 0.015 turns; 'power', A = 2, 200 pulses):
     the folded power against the exact periodic sum of the components (21 images each,
     normalized by its maximum on a 2^20-point phase grid - not the generator's peak search)
     as in 2 (residuals, centroid, width); one centred component [1 0 0.05] gives the
     single 5 % Gaussian bit for bit; Components with 'amplitude' mode are refused.
Uses its own temporary folder per run (21 Oct 2026).
Errors at the end if any check fails.
%}

root = fullfile(fileparts(mfilename('fullpath')), '..');
addpath(fullfile(root, 'functions'));
oldDir = cd(root); restoreDir = onCleanup(@() cd(oldDir));
tmp = tempname(tempdir);                         % own folder per run (21 Oct)
mkdir(tmp);
cleanTmp = onCleanup(@() rmdir(tmp, 's'));
fails = {};

% ---------------------------------------------------------------------------------
% 1. Format and ground truth
% ---------------------------------------------------------------------------------
T1 = 1/1024; f1 = 2^27; L1 = 12.5 * T1; A1 = 2; duty = 5;
[x1, i1] = gen(tmp, 'fmt.dat', T1, f1, A1, L1, duty, 'Seed', 1);
N1 = round(L1 * f1) + 1;
fwhm = duty/100 * T1;
pcTrue = ((0:12) + 0.5) * T1;
f = dir(fullfile(tmp, 'fmt.dat'));
S = load(fullfile(tmp, 'fmt_info.mat'));
ok = [i1.N == N1, numel(x1) == N1, f.bytes == 4*N1, i1.samplesWritten == N1, ...
      i1.bytesWritten == 4*N1, i1.duration == (N1 - 1) / f1, ...
      i1.nPulses == 13, isequal(i1.pulseCenters, pcTrue), i1.pulseCenters(end) == i1.duration, ...
      abs(i1.sigma / (fwhm / (2*sqrt(2*log(2)))) - 1) < 1e-12, ...
      abs(i1.FWHM_power / fwhm - 1) < 1e-12, abs(i1.FWHM_amplitude / (fwhm*sqrt(2)) - 1) < 1e-12, ...
      strcmp(i1.precision, 'single'), strcmp(i1.byteOrder, 'ieee-le'), ~i1.isComplex, ...
      i1.fs == f1, i1.actualFsOut == f1, i1.seed == 1, strcmp(i1.envelopeMode, 'power'), ...
      isequal(rmfield(i1, 'infoFile'), S.info)];
pass = all(ok);
fprintf(['1. format and ground truth: N %d (= round(L f_in) + 1), %d bytes, duration %.10g s, ' ...
         '%d pulse centres, the last at %.10g s = the last sample; sigma, FWHMs, info file: ' ...
         '%d of %d checks: %s\n'], i1.N, f.bytes, i1.duration, i1.nPulses, i1.pulseCenters(end), ...
    nnz(ok), numel(ok), passStr(pass));
if ~pass, fails{end+1} = 'format / ground truth'; end

% ---------------------------------------------------------------------------------
% 2.-4. Power profiles against the exact periodic envelope
% ---------------------------------------------------------------------------------
T = 1e-3; fIn = 100e6; P = round(T * fIn); A = 2; L = 0.2; nB = 1000;
cases = {'2. power mode, 5 %',      5, 'power',     1e-3
         '3. amplitude mode, 5 %',  5, 'amplitude', 1e-3
         '4. widest real pulses, 20 %', 20, 'power', 0};
for c = 1:size(cases, 1)
    [label, dc, mode, cut] = cases{c, :};
    [x, info] = gen(tmp, 'prof.dat', T, fIn, A, L, dc, 'Seed', 10 + c, 'EnvelopeMode', mode);
    r = profileCheck(x, envGrid(P, dc, mode), A, P, nB, cut);
    fwhmNom = info.FWHM_power;                       % what the power profile should have
    k = 2*sqrt(2*log(2)) * T;
    pass = abs(r.zMean) < 4 / sqrt(r.nOn) && abs(r.zStd - 1) < 4 / sqrt(2 * r.nOn) && ...
           abs(r.C - r.Cm) < 4 * r.sC && abs(r.V - r.Vm) < 4 * r.sV;
    if c < 3, pass = pass && abs(r.Cm - 0.5) < 1e-7; end   % model centre at 1/2 (0.1 ns)
    if strcmp(mode, 'amplitude'), pass = pass && abs(fwhmNom / (dc/100 * T / sqrt(2)) - 1) < 1e-12; end
    fprintf(['%s: %d bins, residual / sigma mean %+.3f (+- %.3f), std %.3f (1 +- %.3f); ' ...
             'centre %.6f turns (model %.6f +- %.1e); width %.5e vs model %.5e (+- %.1e), ' ...
             'FWHM %.4f us (nominal %.4f us): %s\n'], label, r.nOn, r.zMean, 1 / sqrt(r.nOn), ...
        r.zStd, 1 / sqrt(2 * r.nOn), r.C, r.Cm, r.sC, r.V, r.Vm, r.sV, k * sqrt(r.V) * 1e6, ...
        fwhmNom * 1e6, passStr(pass));
    if ~pass, fails{end+1} = label; end %#ok<AGROW>
end

% ---------------------------------------------------------------------------------
% 5. Reproducibility
% ---------------------------------------------------------------------------------
L5 = 0.01;
xa = gen(tmp, 'a.dat', T, fIn, 1, L5, 5, 'Seed', 7);
xb = gen(tmp, 'b.dat', T, fIn, 1, L5, 5, 'Seed', 7);
xc = gen(tmp, 'c.dat', T, fIn, 1, L5, 5, 'Seed', 7, 'BlockSize', 333333);
xd = gen(tmp, 'd.dat', T, fIn, 1, L5, 5, 'Seed', 7, 'BlockSize', 7919);
x4 = gen(tmp, 'e.dat', T, fIn, 4, L5, 5, 'Seed', 7);
xo = gen(tmp, 'f.dat', T, fIn, 1, L5, 5, 'Seed', 8);
gState = rng;                                        % the caller's global RNG, restored below
rng(123, 'twister'); r1 = rand(1, 5);
rng(123, 'twister'); gen(tmp, 'g.dat', T, fIn, 1, L5, 5, 'Seed', 7); r2 = rand(1, 5);
rng(gState);
[xn, iN] = gen(tmp, 'h.dat', T, fIn, 1, L5, 5);
xr = gen(tmp, 'i.dat', T, fIn, 1, L5, 5, 'Seed', iN.seed);
nrm = abs(xa) >= realmin('single');                 % normal float32 (not subnormal)
scaleOk = isequal(x4(nrm), 4 * xa(nrm)) && max(abs(double(x4) - 4 * double(xa))) < 1e-43;
ok = [isequal(xa, xb), isequal(xa, xc), isequal(xa, xd), scaleOk, ~isequal(xa, xo), ...
      isequal(r1, r2), iN.seed == floor(iN.seed) && iN.seed >= 0 && iN.seed < 2^32, isequal(xn, xr)];
pass = all(ok);
fprintf(['5. reproducibility: same seed %d; BlockSize 333333 / 7919 = one block %d / %d; ' ...
         'A = 4 exactly 4 x A = 1 (normal floats) %d; other seed differs %d; global RNG untouched %d; clock ' ...
         'seed %d reproduces the file %d: %s\n'], ok(1), ok(2), ok(3), ok(4), ok(5), ok(6), ...
    iN.seed, ok(8), passStr(pass));
if ~pass, fails{end+1} = 'reproducibility'; end

% ---------------------------------------------------------------------------------
% 6. Input checks
% ---------------------------------------------------------------------------------
fb = fullfile(tmp, 'bad.dat');
ids = {errId(@() generatePulsarSignal(fb, T, fIn, 1, L5, 5, 'EnvelopeMode', 'bogus')), ...
       errId(@() generatePulsarSignal(fb, T, fIn, 1, L5, 5, 'Seed', -1)), ...
       errId(@() generatePulsarSignal(fb, T, fIn, 1, L5, 5, 'Seed', 2^32)), ...
       errId(@() generatePulsarSignal(fb, T, fIn, 1, L5, 5, 'Seed', 1.5))};
want = {'generatePulsarSignal:badMode', 'generatePulsarSignal:badSeed', ...
        'generatePulsarSignal:badSeed', 'generatePulsarSignal:badSeed'};
duties = [25 20];                                  % wider than real pulsars; the widest real
warns = false(size(duties));
for i = 1:numel(duties)
    lastwarn('');
    generatePulsarSignal(fb, T, fIn, 1, 2e-3, duties(i), 'Seed', 1, 'SaveInfo', false);
    [~, wId] = lastwarn;
    warns(i) = strcmp(wId, 'generatePulsarSignal:overlap');
end
pass = isequal(ids, want) && isequal(warns, [true false]);
fprintf(['6. input checks: %d of 4 bad inputs refused with the right id; dutycycle 25 %% warns %d, ' ...
         '20 %% warns %d: %s\n'], nnz(strcmp(ids, want)), warns(1), warns(2), passStr(pass));
if ~pass, fails{end+1} = 'input checks'; end

% ---------------------------------------------------------------------------------
% 7. Components (multi-component power profile)
% ---------------------------------------------------------------------------------
cmp7 = [0.6 0 0.02; 0.4 -0.03 0.05; 0.3 0.06 0.015];
[x7, info7] = gen(tmp, 'comp.dat', T, fIn, A, L, 5, 'Seed', 21, 'Components', cmp7);
r7 = profileCheck(x7, envGridComp(P, cmp7), A, P, nB, 1e-3);
okP7 = abs(r7.zMean) < 4 / sqrt(r7.nOn) && abs(r7.zStd - 1) < 4 / sqrt(2 * r7.nOn) && ...
       abs(r7.C - r7.Cm) < 4 * r7.sC && abs(r7.V - r7.Vm) < 4 * r7.sV;
xa = gen(tmp, 'g1.dat', T, fIn, A, 0.02, 5, 'Seed', 22);
xb = gen(tmp, 'g2.dat', T, fIn, A, 0.02, 5, 'Seed', 22, 'Components', [1 0 0.05]);
okId7 = isequal(xa, xb);
okEr7 = strcmp(errId(@() generatePulsarSignal(fullfile(tmp, 'e.dat'), T, fIn, A, 0.01, 5, ...
    'Components', cmp7, 'EnvelopeMode', 'amplitude')), 'generatePulsarSignal:components');
pass = okP7 && okId7 && okEr7 && isequal(info7.components, cmp7);
fprintf(['7. components (3 Gaussians): %d bins, residual / sigma mean %+.3f (+- %.3f), std %.3f ' ...
         '(1 +- %.3f); centroid %.6f vs model %.6f (+- %.1e); width %.5e vs model %.5e (+- %.1e); ' ...
         '[1 0 0.05] = single Gaussian bit for bit %d; amplitude mode refused %d: %s\n'], r7.nOn, ...
    r7.zMean, 1 / sqrt(r7.nOn), r7.zStd, 1 / sqrt(2 * r7.nOn), r7.C, r7.Cm, r7.sC, r7.V, r7.Vm, ...
    r7.sV, okId7, okEr7, passStr(pass));
if ~pass, fails{end+1} = 'components'; end

% ---------------------------------------------------------------------------------
if isempty(fails)
    fprintf('testGeneratePulsarSignal: ALL PASSED\n');
else
    error('testGeneratePulsarSignal:failed', 'FAILED: %s', strjoin(fails, ', '));
end
end


% =================================================================================
function [x, info] = gen(tmp, name, T, fIn, A, L, duty, varargin)
% generatePulsarSignal into tmp/name; the samples as a single column.
f = fullfile(tmp, name);
info = generatePulsarSignal(f, T, fIn, A, L, duty, varargin{:});
fid = fopen(f, 'r', 'ieee-le'); x = fread(fid, Inf, 'single=>single'); fclose(fid);
end

function G = envGrid(P, duty, mode)
% Exact periodic POWER profile on the phase grid m/P (m = 0..P-1), peak 1 at phase 1/2:
% the sum of 21 Gaussian images of the generator's G (std from the FWHM in turns), divided
% by its exact peak; 'amplitude' mode: G is the amplitude envelope, power = G^2.
sig = duty/100 / (2*sqrt(2*log(2)));
ph  = (0:P-1).' / P;
j   = -10:10;
G   = sum(exp(-0.5 * ((ph - 0.5 - j) / sig).^2), 2) / sum(exp(-0.5 * (j / sig).^2));
if strcmp(mode, 'amplitude'), G = G.^2; end
end

function G = envGridComp(P, comps)
% Exact periodic power profile of the components on the phase grid m/P, their reference point
% at phase 1/2: 21 images of each Gaussian, divided by the maximum on a 2^20-point grid
% (independent of the generator's own peak search).
sg = comps(:, 3) / (2*sqrt(2*log(2)));
j  = -10:10;
f  = @(ph) sum(cell2mat(arrayfun(@(i) comps(i, 1) * sum(exp(-0.5 * ((ph - 0.5 - comps(i, 2) - j) ...
    / sg(i)).^2), 2), 1:size(comps, 1), 'UniformOutput', false)), 2);
G  = f((0:P-1).' / P) / max(f((0:2^20-1).' / 2^20));
end

function r = profileCheck(x, Gpow, A, P, nB, cut)
% Fold x^2 into nB phase bins and compare with A^2*Gpow over the same samples. Bins with
% model > cut*A^2 are used: normalized residuals, and the centroid and second moment about
% phase 1/2 of data and model with their sigma (linear propagation of the per-bin sigma).
n    = numel(x);
m    = mod((0:n-1).', P);                           % phase index of every sample
b    = floor(m * nB / P) + 1;                       % phase bin
mod1 = A^2 * Gpow(m + 1);                           % expected power of every sample
cnt  = accumarray(b, 1, [nB 1]);
meas = accumarray(b, double(x).^2, [nB 1]) ./ cnt;
modl = accumarray(b, mod1, [nB 1]) ./ cnt;
sig  = sqrt(accumarray(b, 2 * mod1.^2, [nB 1])) ./ cnt;
phB  = accumarray(b, m / P, [nB 1]) ./ cnt;          % mean phase of each bin
on   = modl > cut * A^2;
z    = (meas(on) - modl(on)) ./ sig(on);
w = phB(on); s = sig(on); M = meas(on); Mm = modl(on);
C  = sum(w .* M) / sum(M);   Cm = sum(w .* Mm) / sum(Mm);
d2 = (w - 0.5).^2;
V  = sum(d2 .* M) / sum(M);  Vm = sum(d2 .* Mm) / sum(Mm);
r = struct('nOn', nnz(on), 'zMean', mean(z), 'zStd', std(z), 'C', C, 'Cm', Cm, ...
    'sC', sqrt(sum(((w - C) / sum(M)).^2 .* s.^2)), 'V', V, 'Vm', Vm, ...
    'sV', sqrt(sum(((d2 - V) / sum(M)).^2 .* s.^2)));
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
