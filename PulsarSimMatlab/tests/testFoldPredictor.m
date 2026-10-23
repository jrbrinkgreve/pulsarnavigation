function testFoldPredictor()
%TESTFOLDPREDICTOR  Unit tests for foldProfile's 'Predictor' option (block D2b: folding at the
%predicted phase).
%{
Run from the PulsarSimMatlab folder: run('tests/testFoldPredictor.m'), or
runAllTests("testFoldPredictor"). ~10 s; small files in its own temporary folder.

Synthetic, noise-free detected power of one channel: 1 + 5 G(phi), G a periodic Gaussian (FWHM 0.02
turns) of the pulse phase phi. A 50 Hz pulsar seen from a probe flying along the line of sight at
300 km/s ('SSB' trajectory): Doppler factor 1e-3, so over the 20 s of data (20 us time bins) the
pulse drifts by one full turn against a fixed f0. The data's phase comes from delayModel directly
(phi = f0 (t + delta(t)) - N0, delta on a 10 ms grid, spline-interpolated), not from the predictor.

  1. Folded with the predictor (NBin 256, 100 turns per sub-int): the whole-data profile equals
     the fold made with the data's TRUE phases (the same linear assignment, computed here;
     1e-5 of the peak): the Doppler is removed. Shown: the difference from Doppler-free data
     folded with f0 (~2e-3: there every turn samples the same phases, here they drift through
     each bin - a sampling pattern, not an error).
  2. The same Doppler data folded with f0 alone: smeared over the turn (its peak < 30 % of 1).
  3. Sub-int references: pred.phase(tRef) = turnRef (1e-9 turns; the partial last sub-int's
     tRef lies just past the span, where the phase continues linearly with the edge
     frequency), fRef = pred.freq(tRef) (exactly), fRef / f0 - 1 = the Doppler factor (1e-6).
  4. A predictor without delays ('DelayArgs' {'Terms', []}) and TRef = t1 + TRef_file gives the
     fold of the old f0 / TRef model: profiles 1e-9, tRef 1e-9 s, fRef 1e-12 relative (TRef_file =
     12/1024 s, so that t1 + TRef_file is exact in double: an inexact sum would shift the phase by
     ~1e-7 s, far above 1e-9).
  5. A predictor that does not cover the data and a non-predictor struct are refused
     ('foldProfile:predictor').
  6. mergeSubints on the predictor fold (23 Oct 2026): the 100-turn sub-ints merged in pairs =
     the predictor fold with 200 turns per sub-int (sums 1e-12, turnRef equal, tRef 1e-9 s, fRef
     1e-9 relative). Shown: how far the f0 model's tRef would be off (~1e-3 x t, ~10 ms). A
     group whose reference turn falls in a gap (its middle sub-int left out) is refused
     ('mergeSubints:predictor').
Errors at the end if any check fails. (21 Oct 2026, block D2b.)
%}

root = fullfile(fileparts(mfilename('fullpath')), '..');
addpath(fullfile(root, 'functions'));
fails = {};
tmp = tempname(tempdir); mkdir(tmp);
cleanTmp = onCleanup(@() rmdir(tmp, 's'));
c = 299792458;
f0 = 50; dt = 20e-6; Ld = 20; N = round(Ld / dt); NB = 256; nPer = 100;
t1 = 9555.5 * 86400;                                  % absolute time of the file's t = 0
psr = struct('f0', f0, 'F1', 0, 'TRef', t1, 'Phi0', 0, 'raDeg', 30, 'decDeg', 10);
n = [cosd(10) * cosd(30); cosd(10) * sind(30); sind(10)];
r0 = 1.3 * 149597870700 * [0.6; 0.8; 0]; v = 3e5 * n;  % Doppler 1e-3
probe = makeTrajectory('custom', 'Frame', 'SSB', 'PosFcn', @(t) r0 + v * (t - t1), ...
    'VelFcn', @(t) repmat(v, 1, numel(t)));
pred = makePredictor(psr, probe, [t1, t1 + Ld]);
t = (0:N-1) * dt;                                    % file time = predictor local time
tg = 0:0.01:Ld;                                       % delays on a 10 ms grid, then spline
dg = delayModel(t1 + tg, probe, psr).total;
del = spline(tg, dg, t);
phiD = f0 * (t + del) - pred.N0;                      % the data's phase (Phi0 = 0, TRef = t1)
sg = 0.02 / (2 * sqrt(2 * log(2)));
pulse = @(ph) 1 + 5 * exp(-0.5 * ((ph - round(ph)) / sg).^2);
detD = writePower(fullfile(tmp, 'dop.dat'), pulse(phiD), dt);
detR = writePower(fullfile(tmp, 'ref.dat'), pulse(f0 * t), dt);
fa = {'NBin', NB, 'SubintPeriods', nPer, 'SaveFile', false, 'Verbose', false};

% 1. Doppler data folded with the predictor = Doppler-free data folded with f0
[iP, fP] = foldProfile(detD, 'x', f0, 'Predictor', pred, fa{:});
[~, fR]  = foldProfile(detR, 'x', f0, fa{:});
x  = (phiD - round(phiD)) * NB;                       % the fold with the true phases
j0 = floor(x); a = x - j0; pv = pulse(phiD);
Lb = [mod(j0, NB) + 1, mod(j0 + 1, NB) + 1].'; wv = [1 - a, a].'; pv2 = [pv, pv].';
pT = accumarray(Lb(:), wv(:) .* pv2(:), [NB 1]) ./ accumarray(Lb(:), wv(:), [NB 1]);
e1 = max(abs(fP.profTotal - pT)) / 5;
full = fP.subint.nTimeBins >= 0.99 * nPer / (f0 * dt);
pP = mean(fP.prof(:, full), 2); pR = mean(fR.prof(:, 2:end-1), 2);
e1r = max(abs(pP - pR)) / 5;
ok1 = e1 < 1e-5;
fprintf(['1. Doppler 1e-3, folded with the predictor vs the fold with the true phases: %.1e of the peak ' ...
         '(vs Doppler-free data folded with f0: %.1e, a sampling pattern): %s\n'], e1, e1r, passStr(ok1));
if ~ok1, fails{end+1} = 'predictor fold'; end

% 2. The Doppler data folded with f0 alone: smeared
[~, fF] = foldProfile(detD, 'x', f0, 'TRef', 0, fa{:});
pk2 = (max(fF.profTotal) - 1) / 5;
ok2 = pk2 < 0.3;
fprintf('2. the same data folded with f0 alone: peak %.3f of the true one (smeared over %.2f turns): %s\n', ...
    pk2, f0 * Ld * 1e-3, passStr(ok2));
if ~ok2, fails{end+1} = 'smeared without predictor'; end

% 3. Sub-int references
sr = fP.subint;
tc = min(max(sr.tRef.', pred.seg(1, 1)), pred.seg(end, 2));       % linear beyond the span
e3a = max(abs(pred.phase(tc) + pred.freq(tc) .* (sr.tRef.' - tc) - sr.turnRef.'));
e3b = max(abs(sr.fRef.' - pred.freq(tc)));
dopT = (n.' * v) / c;
e3c = max(abs(sr.fRef / f0 - 1 - dopT));
ok3 = e3a < 1e-9 && e3b == 0 && e3c < 1e-6 && ~isempty(iP.predictor);
fprintf(['3. sub-int references: phase(tRef) - turnRef %.1e turns, fRef - freq(tRef) %.1e Hz, ' ...
         'fRef/f0 - 1 vs v.n/c %.1e: %s\n'], e3a, e3b, e3c, passStr(ok3));
if ~ok3, fails{end+1} = 'sub-int references'; end

% 4. A predictor without delays = the old f0 / TRef model
TRefF = 12 / 1024;                                    % exact in binary: t1 + TRefF is exact (8e8 s)
p0 = makePredictor(struct('f0', f0, 'F1', 0, 'TRef', t1 + TRefF, 'Phi0', 0.2, 'raDeg', 30, 'decDeg', 10), ...
    probe, [t1, t1 + Ld], 'DelayArgs', {'Terms', string.empty(1, 0)});
[~, f4p] = foldProfile(detD, 'x', f0, 'Predictor', p0, fa{:});
[~, f4o] = foldProfile(detD, 'x', f0, 'TRef', TRefF, 'Phi0', 0.2, fa{:});
e4p = max(abs(f4p.prof(:) - f4o.prof(:)), [], 'omitnan');
e4t = max(abs(f4p.subint.tRef - f4o.subint.tRef));
e4f = max(abs(f4p.subint.fRef ./ f4o.subint.fRef - 1));
ok4 = isequal(size(f4p.prof), size(f4o.prof)) && e4p < 1e-9 && e4t < 1e-9 && e4f < 1e-12;
fprintf('4. predictor without delays vs the f0 / TRef model: profiles %.1e, tRef %.1e s, fRef %.1e: %s\n', ...
    e4p, e4t, e4f, passStr(ok4));
if ~ok4, fails{end+1} = 'equivalence'; end

% 5. Input checks
pShort = makePredictor(psr, probe, [t1, t1 + Ld / 2]);
ok5 = strcmp(errId(@() foldProfile(detD, 'x', f0, 'Predictor', pShort, fa{:})), 'foldProfile:predictor') && ...
      strcmp(errId(@() foldProfile(detD, 'x', f0, 'Predictor', struct('a', 1), fa{:})), 'foldProfile:predictor');
fprintf('5. predictor too short / not a predictor refused %d: %s\n', ok5, passStr(ok5));
if ~ok5, fails{end+1} = 'input checks'; end

% 6. mergeSubints on a predictor fold: pairs of sub-ints = the fold with twice the sub-int length
[~, fP2] = foldProfile(detD, 'x', f0, 'Predictor', pred, 'NBin', NB, 'SubintPeriods', 2 * nPer, ...
    'SaveFile', false, 'Verbose', false);
grp = floor((sr.turnFirst - sr.turnFirst(1)) / (2 * nPer)) + 1;
fM = mergeSubints(fP, iP, grp);
rel = @(a, b) max(abs(a(:) - b(:))) / max(abs(b(:)));
okS = isequal(size(fM.sum), size(fP2.sum)) && rel(fM.sum, fP2.sum) < 1e-12 && ...
    rel(fM.weight, fP2.weight) < 1e-12 && rel(fM.weight2, fP2.weight2) < 1e-12 && ...
    rel(fM.weightX, fP2.weightX) < 1e-12 && isequal(fM.subint.turnRef, fP2.subint.turnRef);
e6t = max(abs(fM.subint.tRef - fP2.subint.tRef));
e6f = max(abs(fM.subint.fRef ./ fP2.subint.fRef - 1));
e6o = max(abs(fM.subint.turnRef / f0 - fP2.subint.tRef));    % the f0 model (TRef 0, F1 0) instead
trAll = round((sr.turnFirst(1) + sr.turnLast(end)) / 2);     % one group of all, its middle sub-int left out
keep = find(~(sr.turnFirst <= trAll & sr.turnLast >= trAll));
okG = strcmp(errId(@() mergeSubints(subFold(fP, keep), iP, ones(numel(keep), 1))), 'mergeSubints:predictor');
ok6 = okS && e6t < 1e-9 && e6f < 1e-9 && okG;
fprintf(['6. mergeSubints on the predictor fold: pairs = the 200-turn fold (sums, weights, turnRef) %d, ' ...
         'tRef %.1e s, fRef %.1e relative (the f0 model would be off by %.1e s); turnRef in a gap ' ...
         'refused %d: %s\n'], okS, e6t, e6f, e6o, okG, passStr(ok6));
if ~ok6, fails{end+1} = 'mergeSubints predictor'; end

if isempty(fails)
    fprintf('testFoldPredictor: ALL PASSED\n');
else
    error('testFoldPredictor:failed', 'FAILED: %s', strjoin(fails, ', '));
end
end


% =================================================================================
function det = writePower(file, p, dt)
% One-channel float32 power file and the detectChannels-like info foldProfile needs.
fid = fopen(file, 'w', 'ieee-le'); fwrite(fid, single(p), 'float32'); fclose(fid);
det = struct('file', file, 'nChan', 1, 'N', numel(p), 'binTime0', 0, 'binDt', dt, ...
    'fullySupportedBins', [1 numel(p)], 'byteOrder', 'ieee-le', 'chanFreqs', 0);
end

function f = subFold(f, idx)
% Sub-ints idx of a fold (the arrays along dimension 2, the subint fields alike).
for fn = ["sum", "prof", "weight", "weight2", "weightX"], f.(fn) = f.(fn)(:, idx, :, :); end
for fn = string(fieldnames(f.subint)).', f.subint.(fn) = f.subint.(fn)(idx); end
end

function id = errId(f)
id = '';
try
    f();
catch err
    id = err.identifier;
end
end

function s = passStr(p)
if p, s = 'PASS'; else, s = 'FAIL'; end
end
