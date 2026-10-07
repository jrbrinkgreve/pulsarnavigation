function makeToaReference()
%MAKETOAREFERENCE  Save reference outputs of estimateTOA and detectPulsar BEFORE A2.
%{
Run once, from the PulsarSimMatlab folder, with estimateTOA / detectPulsar of
commit 918cae3 (before the A2 changes): run('tests/makeToaReference.m').
Writes data/mc/toaRef_preA2.mat. tests/testChannelTOA.m then checks that the
new functions reproduce these outputs bit for bit on single-channel folds
(the full-band path of main.m must not change). To regenerate: check out
functions/estimateTOA.m and functions/detectPulsar.m of 918cae3 first.

Input: the frozen full-band folds of data/mc/foldRef_pre3c.mat
(tests/makeFoldReference.m; seed-43 data, main.m's parameters), so the
reference does not depend on data files that main.m may rewrite.
Cases:
  A    linear assignment, radiometer noise (as main.m)
  Aoff linear, 'offpulse' noise model
  Acov linear, MinCoverage 0.3 (the partial last sub-int, 37 %, gap filling)
  B    nearest assignment
  C    3 turns per sub-integration
detectPulsar for A, Acov, B, C (and A with 'Phase' 0.01).
%}

root = fullfile(fileparts(mfilename('fullpath')), '..');
addpath(fullfile(root, 'functions'));
oldDir = cd(root); restoreDir = onCleanup(@() cd(oldDir));
pipelineParams;

R = load(fullfile(dataDir, 'mc', 'foldRef_pre3c.mat')); fr = R.ref;
template = gaussianTemplate(nBin, ephem.profileFWHM);
info_dedisp = loadInfo(fileDedisp);
Bnoise = noiseBandwidth(fLow, fHigh, info_dedisp.edgeWidth);
q = {'Bnoise', Bnoise, 'Verbose', false};

ref = struct();
[ref.toa.A,    ref.toaInfo.A]    = estimateTOA(fr.A.fold, fr.A.info, template, q{:});
[ref.toa.Aoff, ref.toaInfo.Aoff] = estimateTOA(fr.A.fold, fr.A.info, template, q{:}, 'NoiseModel', 'offpulse');
[ref.toa.Acov, ref.toaInfo.Acov] = estimateTOA(fr.A.fold, fr.A.info, template, q{:}, 'MinCoverage', 0.3);
[ref.toa.B,    ref.toaInfo.B]    = estimateTOA(fr.B.fold, fr.B.info, template, q{:});
[ref.toa.C,    ref.toaInfo.C]    = estimateTOA(fr.C.fold, fr.C.info, template, q{:});
[ref.det.A,    ref.detInfo.A]    = detectPulsar(fr.A.fold, fr.A.info, template, q{:});
[ref.det.Aph,  ref.detInfo.Aph]  = detectPulsar(fr.A.fold, fr.A.info, template, q{:}, 'Phase', 0.01);
[ref.det.Acov, ref.detInfo.Acov] = detectPulsar(fr.A.fold, fr.A.info, template, q{:}, 'MinCoverage', 0.3);
[ref.det.B,    ref.detInfo.B]    = detectPulsar(fr.B.fold, fr.B.info, template, q{:});
[ref.det.C,    ref.detInfo.C]    = detectPulsar(fr.C.fold, fr.C.info, template, q{:});
ref.Bnoise = Bnoise;
ref.made = sprintf('%s, estimateTOA / detectPulsar before A2 (commit 918cae3)', ...
    char(datetime('now')));
save(fullfile(dataDir, 'mc', 'toaRef_preA2.mat'), 'ref');
fprintf('makeToaReference: saved 5 estimateTOA and 5 detectPulsar reference outputs to %s\n', ...
    fullfile(dataDir, 'mc', 'toaRef_preA2.mat'));
fprintf('  A: %d valid TOAs, median error %.3f us; Acov: %d valid\n', nnz(ref.toa.A.valid), ...
    median(ref.toa.A.toaErr(ref.toa.A.valid))*1e6, nnz(ref.toa.Acov.valid));
end
