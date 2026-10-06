function makeFoldReference()
%MAKEFOLDREFERENCE  Save reference folds made with foldProfile BEFORE unit 3c.
%{
Run once, from the PulsarSimMatlab folder, with the foldProfile of commit
e852da6 (before the unit-3c changes): run('tests/makeFoldReference.m').
Writes data/mc/foldRef_pre3c.mat. tests/testFoldWeights.m then checks that
the new foldProfile reproduces these folds bit for bit when no noise
statistics or data weights are given (the full-band path of main.m must not
change). To regenerate: check out functions/foldProfile.m of e852da6 first.

Cases (main.m's parameters, seed-43 data):
  A  full-band power (data/test_envelope.dat), linear assignment
  B  the same, nearest assignment
  C  the same, 3 turns per sub-integration (phase wraps inside a sub-int)
  D  128-channel power (data/chan/test_dedisp_chan_power.dat), linear
%}

root = fullfile(fileparts(mfilename('fullpath')), '..');
addpath(fullfile(root, 'functions'));
oldDir = cd(root); restoreDir = onCleanup(@() cd(oldDir));
pipelineParams;

detA = loadInfo(fileEnvelope);
detD = loadInfo(fullfile(dataDir, 'chan', 'test_dedisp_chan_power.dat'));
base = {'F1', ephem.F1, 'TRef', ephem.TRef, 'NBin', nBin, 'SaveFile', false, 'Verbose', false};
ref = struct();
[ref.A.info, ref.A.fold] = foldProfile(detA, 'x', ephem.f0, base{:}, 'SubintPeriods', subintPeriods);
[ref.B.info, ref.B.fold] = foldProfile(detA, 'x', ephem.f0, base{:}, 'SubintPeriods', subintPeriods, ...
    'Assign', 'nearest');
[ref.C.info, ref.C.fold] = foldProfile(detA, 'x', ephem.f0, base{:}, 'SubintPeriods', 3);
[ref.D.info, ref.D.fold] = foldProfile(detD, 'x', ephem.f0, base{:}, 'SubintPeriods', subintPeriods);
ref.made = sprintf('%s, foldProfile before unit 3c (commit e852da6)', datestr(now));
ref.inputs = {detA.file, detD.file};
save(fullfile(dataDir, 'mc', 'foldRef_pre3c.mat'), 'ref');
fprintf('makeFoldReference: saved 4 reference folds to %s\n', fullfile(dataDir, 'mc', 'foldRef_pre3c.mat'));
end
