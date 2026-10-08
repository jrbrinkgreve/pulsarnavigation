function testBlankChannels()
%TESTBLANKCHANNELS  Unit tests for blankChannels (B1): apply a blanking mask to the channel files.
%{
Run from the PulsarSimMatlab folder: run('tests/testBlankChannels.m').
Uses the seed-43 channel files of tests/testDetectChannels.m (data/chan/
test_rx_IQ_chan_*, test_dedisp_chan_*); makes them first if they are
missing. Writes ~1 GB of temporary files to tempdir (deleted on the way).

  1. Exactness, all 128 channels: masked samples exactly 0 (I and Q), every
     other sample bit-identical to the raw file; blanked counts per channel.
     Mask: random overlapping blanks of 1-3000 samples (every 3rd channel,
     some reaching outside the file), short overlapping blanks across the
     block edges of check 3 (channel 2), a whole channel plus a duplicate
     inside it (channel 5), only an empty and an outside row (channel 7:
     nothing blanked), touching blanks at both file ends (channel 128).
     The expected samples come from a plain loop over the mask rows.
  2. Normalized mask: the explicit cases of channels 2, 5, 7, 128 merge and
     clip as expected; the normalized mask fed back in gives identical
     files and the same normalized mask.
  3. Block size: BlockSize 4099 (blanks across block edges) gives files
     identical to the default (one block).
  4. Empty mask: byte-identical copies of the raw files, nothing blanked,
     and info = info_chan apart from the file names and the added fields.
  5. Interface: loadInfo(outBase) returns the info; dedisperseChannels runs
     on it with the same grid as the unblanked run; every channel without
     blanks gives a dedispersed file bit-identical to data/chan/
     test_dedisp_chan_*; in channel 128 the blanks at the file ends change
     the dedispersed output only within the filter's reach (nPast,
     nFuture) and leave it unchanged (rounding) everywhere else.
Errors at the end if any check fails.
%}

root = fullfile(fileparts(mfilename('fullpath')), '..');
addpath(fullfile(root, 'functions'));
oldDir = cd(root); restoreDir = onCleanup(@() cd(oldDir));
pipelineParams;
chanDir = fullfile(dataDir, 'chan');
tmp = fullfile(tempdir, 'testBlankChannels');
if ~isfolder(tmp), mkdir(tmp); end
cleanTmp = onCleanup(@() cleanFolder(tmp));
fails = {};

baseCh = char(fullfile(chanDir, 'test_rx_IQ_chan'));
baseDC = char(fullfile(chanDir, 'test_dedisp_chan'));
if isfile([baseCh '_info.mat']) && isfile([baseDC '_info.mat'])
    info_chan = loadInfo(baseCh);
    info_dc   = loadInfo(baseDC);
else
    info_IQ = loadInfo(fileIQ);
    info_chan = channelizeIQ(info_IQ.file, baseCh, info_IQ.actualFsOut, info_IQ.fLO, ...
        fLow, fHigh, 'T0', info_IQ.t0);
    info_dc = dedisperseChannels(info_chan, baseDC, ephem.DM, 'RefFreq', refFreq);
end
nCh = info_chan.nChan; Nc = info_chan.N;

% ---------------------------------------------------------------------------------
% 1. Exactness
% ---------------------------------------------------------------------------------
t1 = tic;
rng(5);
randCh = setdiff(1:3:nCh, 7);
rows = {};
for c = randCh                                       % random overlapping blanks
    s = randi([-100, Nc + 100], 40, 1);
    rows{end+1} = [c*ones(40, 1), s, s + randi([0, 2999], 40, 1)]; %#ok<AGROW>
end
k = (1:20).' * 4099;                                 % channel 2: across the BlockSize-4099 edges
rows{end+1} = [2*ones(20, 1), k, k + 1; 2*ones(20, 1), k + 1, k + 5];
rows{end+1} = [5, 1, Nc; 5, 100, 200];               % channel 5: whole, plus a duplicate inside
rows{end+1} = [7, 50, 40; 7, Nc + 5, Nc + 10];       % channel 7: empty row, outside the file
rows{end+1} = [128, -20, 10; 128, 11, 11; 128, Nc - 9, Nc + 50];   % touching, both ends
mask1 = vertcat(rows{:});
mask1 = mask1(randperm(size(mask1, 1)), :);          % unsorted
base1 = fullfile(tmp, 'b1');
info1 = blankChannels(info_chan, mask1, base1, 'Verbose', false);

maxBlank = 0; same = true; countOk = true; nBl = 0;
for c = 1:nCh
    keep = keepOf(mask1, c, Nc);
    x = readRaw(info_chan.chanFiles(c));
    y = readRaw(info1.chanFiles(c));
    maxBlank = max([maxBlank, max(abs(y(:, ~keep)), [], 'all')]);
    same = same && isequal(y(:, keep), x(:, keep));
    countOk = countOk && info1.blankedSamples(c) == nnz(~keep);
    nBl = nBl + nnz(~keep);
end
pass = maxBlank == 0 && same && countOk;
fprintf(['1. %d mask rows, %d blanked samples (%.3f %%): masked samples max |value| %g, other ' ...
         'samples bit-identical %d, counts per channel %d (%.1f s): %s\n'], size(mask1, 1), nBl, ...
    100*nBl/(nCh*Nc), maxBlank, same, countOk, toc(t1), passStr(pass));
if ~pass, fails{end+1} = 'exactness'; end

% ---------------------------------------------------------------------------------
% 2. Normalized mask
% ---------------------------------------------------------------------------------
ivOf = @(info, c) info.mask(info.mask(:, 1) == c, 2:3);
exp2 = [k, k + 5];
okCases = isequal(ivOf(info1, 2), exp2) && isequal(ivOf(info1, 5), [1, Nc]) && ...
    isempty(ivOf(info1, 7)) && isequal(ivOf(info1, 128), [1, 11; Nc - 9, Nc]);
sortedOk = issortedrows(info1.mask) && all(info1.mask(:, 2) <= info1.mask(:, 3));
base2 = fullfile(tmp, 'b2');
info2 = blankChannels(info_chan, info1.mask, base2, 'Verbose', false);
sameFiles = filesEqual(info1.chanFiles, info2.chanFiles);
idem = isequal(info2.mask, info1.mask);
deleteFiles(cellstr(info2.chanFiles));
pass = okCases && sortedOk && sameFiles && idem;
fprintf(['2. explicit cases (channels 2, 5, 7, 128) %d, sorted %d; normalized mask fed back: ' ...
         'identical files %d, same mask %d: %s\n'], okCases, sortedOk, sameFiles, idem, passStr(pass));
if ~pass, fails{end+1} = 'normalized mask'; end

% ---------------------------------------------------------------------------------
% 3. Block size
% ---------------------------------------------------------------------------------
base3 = fullfile(tmp, 'b3');
info3 = blankChannels(info_chan, mask1, base3, 'BlockSize', 4099, 'Verbose', false);
pass = filesEqual(info1.chanFiles, info3.chanFiles) && isequal(info3.mask, info1.mask);
deleteFiles(cellstr(info3.chanFiles));
fprintf('3. BlockSize 4099 vs one block: identical files and mask: %s\n', passStr(pass));
if ~pass, fails{end+1} = 'block size'; end

% ---------------------------------------------------------------------------------
% 4. Empty mask
% ---------------------------------------------------------------------------------
base4 = fullfile(tmp, 'b4');
info4 = blankChannels(info_chan, [], base4, 'Verbose', false);
sameFiles = filesEqual(info_chan.chanFiles, info4.chanFiles);
deleteFiles(cellstr(info4.chanFiles));
added = {'file', 'chanFiles', 'rawFile', 'mask', 'nMaskRows', 'blankedSamples', ...
         'blankedFraction', 'blockSize', 'elapsed', 'infoFile'};
a = rmfield(info4, intersect(fieldnames(info4), added));
b = rmfield(info_chan, intersect(fieldnames(info_chan), added));
infoOk = isequal(a, b) && isempty(info4.mask) && all(info4.blankedSamples == 0) && ...
    strcmp(info4.rawFile, info_chan.file);
pass = sameFiles && infoOk;
fprintf('4. empty mask: byte-identical copies %d, info = info_chan apart from file names %d: %s\n', ...
    sameFiles, infoOk, passStr(pass));
if ~pass, fails{end+1} = 'empty mask'; end

% ---------------------------------------------------------------------------------
% 5. Interface: dedisperseChannels on the blanked files
% ---------------------------------------------------------------------------------
t5 = tic;
L1 = loadInfo(base1);                                % saved before infoFile was added
loadOk = isequal(L1, rmfield(info1, 'infoFile'));
info_dcB = dedisperseChannels(info1, fullfile(tmp, 'dcB'), ephem.DM, 'RefFreq', refFreq, ...
    'Verbose', false);
gridOk = info_dcB.N == info_dc.N && isequal(info_dcB.Nfft, info_dc.Nfft) && ...
    isequal(info_dcB.nPast, info_dc.nPast) && isequal(info_dcB.nFuture, info_dc.nFuture) && ...
    isequal(info_dcB.fullySupported, info_dc.fullySupported) && strcmp(info_dcB.inFile, info1.file);
clean = find(info1.blankedSamples == 0).';
cleanOk = filesEqual(info_dc.chanFiles(clean), info_dcB.chanFiles(clean));
% channel 128: blanks [1, 11] and [Nc-9, Nc]; output a uses inputs a-nP .. a+nF
c = 128; nP = info_dc.nPast(c); nF = info_dc.nFuture(c);
yU = readRaw(info_dc.chanFiles(c)); yB = readRaw(info_dcB.chanFiles(c));
dY = sqrt(sum((double(yB) - double(yU)).^2, 1));
rmsU = sqrt(mean(sum(double(yU).^2, 1)));
Nd   = size(yU, 2);                                  % same grid: Nd = Nc
far  = 12 + nP : Nd - 10 - nF;                       % beyond the reach of both blanks
near = [1 : 11 + nP, Nd - 9 - nF : Nd];
farMax = max(dY(far)) / rmsU; nearMax = max(dY(near)) / rmsU;
reachOk = farMax < 1e-5 && nearMax > 1e-2;
deleteFiles(cellstr(info_dcB.chanFiles));
deleteFiles(cellstr(info1.chanFiles));
pass = loadOk && gridOk && cleanOk && reachOk;
fprintf(['5. loadInfo %d; dedisperseChannels: same grid %d, %d unblanked channels bit-identical ' ...
         '%d; channel 128 (reach -%d..+%d samples): max change beyond the reach %.1e, within ' ...
         '%.2f of the rms (%.0f s): %s\n'], loadOk, gridOk, numel(clean), cleanOk, nF, nP, ...
    farMax, nearMax, toc(t5), passStr(pass));
if ~pass, fails{end+1} = 'interface'; end

% ---------------------------------------------------------------------------------
if isempty(fails)
    fprintf('testBlankChannels: ALL PASSED\n');
else
    error('testBlankChannels:failed', 'FAILED: %s', strjoin(fails, ', '));
end
end


% =================================================================================
function keep = keepOf(mask, c, Nc)
% Samples of channel c not covered by any mask row (plain loop, the reference).
keep = true(1, Nc);
r = mask(mask(:, 1) == c, 2:3);
for i = 1:size(r, 1)
    a = max(1, r(i, 1)); b = min(Nc, r(i, 2));
    if a <= b, keep(a:b) = false; end
end
end

function x = readRaw(file)
% cf32 file as [2 x N] single (rows I, Q), for bit-exact comparisons.
fid = fopen(file, 'r', 'ieee-le');
c = onCleanup(@() fclose(fid));
x = fread(fid, [2 Inf], 'single=>single');
end

function eq = filesEqual(filesA, filesB)
% Byte-identical, file by file.
eq = numel(filesA) == numel(filesB);
for i = 1:numel(filesA)
    if ~eq, return; end
    eq = isequal(readBytes(filesA(i)), readBytes(filesB(i)));
end
end

function b = readBytes(file)
fid = fopen(file, 'r');
c = onCleanup(@() fclose(fid));
b = fread(fid, Inf, 'uint8=>uint8');
end

function s = passStr(p)
if p, s = 'PASS'; else, s = 'FAIL'; end
end

function deleteFiles(files)
for i = 1:numel(files)
    if isfile(files{i}), delete(files{i}); end
end
end

function cleanFolder(d)
f = dir(fullfile(d, '*'));
for i = 1:numel(f)
    if ~f(i).isdir, delete(fullfile(d, f(i).name)); end
end
end
