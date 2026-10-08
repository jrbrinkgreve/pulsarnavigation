function makeRfiReference()
%MAKERFIREFERENCE  Save a reference output of addNoiseAndRFI BEFORE B5a (radar rise time).
%{
Run once, from the PulsarSimMatlab folder, with addNoiseAndRFI / rfiSource of
commit f4181f5 (before B5a): run('tests/makeRfiReference.m').
Writes data/mc/rfiRef_preB5a.mat: 20 us of receiver noise + all five RFI
types at RF (80,000 float32 samples at 4 GHz), so tests/testRadarEdges.m can
check that the new code with 'RiseTime' 0 reproduces it bit for bit. To
regenerate: check out functions/addNoiseAndRFI.m and functions/rfiSource.m
of f4181f5 first.

Sources (strong, so every one is visible in 20 us): GNSS L1- and L2-like
bpsk, a radar pulse 5..7 us (1300 MHz, 1 MHz chirp), a drifting carrier,
impulses at 2e5/s (about 4 bursts).
%}

root = fullfile(fileparts(mfilename('fullpath')), '..');
addpath(fullfile(root, 'functions'));
oldDir = cd(root); restoreDir = onCleanup(@() cd(oldDir));
pipelineParams;

[src, N, seed] = rfiReferenceCase();
tmp = fullfile(tempdir, 'makeRfiReference');
if ~isfolder(tmp), mkdir(tmp); end
zeroFile = fullfile(tmp, 'zero.dat'); outFile = fullfile(tmp, 'rx.dat');
fid = fopen(zeroFile, 'w', 'ieee-le'); fwrite(fid, zeros(1, N, 'single'), 'single'); fclose(fid);
info = addNoiseAndRFI(zeroFile, outFile, f_in, 'Band', [fLow fHigh], 'NoiseStd', 1, ...
    'RFI', src, 'Seed', seed, 'SaveInfo', false, 'Verbose', false);
fid = fopen(outFile, 'r', 'ieee-le'); y = fread(fid, Inf, 'single=>single'); fclose(fid);
delete(zeroFile); delete(outFile);
ref = struct('y', y, 'N', N, 'seed', seed, 'fs', f_in, 'band', [fLow fHigh], ...
    'nBursts', numel(info.rfi(5).evStart), 'commit', 'f4181f5');
save(fullfile(dataDir, 'mc', 'rfiRef_preB5a.mat'), 'ref');
fprintf('makeRfiReference: %d samples, %d impulse bursts -> data/mc/rfiRef_preB5a.mat\n', N, ref.nBursts);
end


function [src, N, seed] = rfiReferenceCase()
% The reference case (the same list is rebuilt in tests/testRadarEdges.m).
N = 80000; seed = 3;
src = [ ...
    rfiSource('bpsk',    'Freq', 1575.42e6, 'ChipRate', 1.023e6, 'INRdB', 0), ...
    rfiSource('bpsk',    'Freq', 1227.60e6, 'ChipRate', 10.23e6, 'INRdB', 0), ...
    rfiSource('pulsed',  'Freq', 1300e6, 'PulseWidth', 2e-6, 'PRF', 373, 'ChirpBW', 1e6, ...
              'StartTime', 5e-6, 'INRdB', 20), ...
    rfiSource('cw',      'Freq', 1351.3e6, 'Drift', 1e9, 'INRdB', 0), ...
    rfiSource('impulse', 'Rate', 2e5, 'Duration', 200e-9, 'INRdB', 20)];
end
