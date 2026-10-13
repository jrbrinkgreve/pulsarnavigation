function makeTestData(opts)
%MAKETESTDATA  Freeze main.m's seed-43 files as the unit tests' own data (data/test/), once.
%{
Run once, from the PulsarSimMatlab folder: run('tests/makeTestData.m'). ~10 s, ~2.2 GB.
Not a test (like the make*Reference scripts): it copies, it does not compute.

Why: the unit tests used main.m's files in data/ (data/test_rx_IQ.dat, ...) and the channel
files in data/chan/, which main.m rewrites; a main run with other settings changed what the
tests checked, or broke their bit-for-bit regressions. Now the tests read only data/test/,
with their own frozen settings (tests/testParams.m).

What it copies (paths in the info structs are rewritten to data/test/):
  - data files + info: test_rx_IQ.dat, test_IQ_dedispersed.dat, test_envelope.dat;
  - info only (the tests read the ground truth, not the samples): test.dat, test_dispersed.dat,
    test_rx.dat;
  - the channel files of the tests: chan/test_rx_IQ_chan, chan/test_dedisp_chan,
    chan/test_dedisp_chan_power.dat (the tests remake them, but single tests need them);
  - the frozen references: mc/foldRef_pre3c.mat, mc/toaRef_preA2.mat, mc/rfiRef_preB5a.mat.
First it checks that the files in data/ are the reference configuration of testParams
(pulsar, dispersion, -5 dB noise with seed 43 and no RFI, receiver, f_out, channels, and
that each file was made from the one before it); otherwise it stops and lists what
differs. It does not check that the contents belong together beyond that: the tests do
(e.g. testDedisperseChannels' bit-for-bit regression), so run run('tests/runAllTests.m')
afterwards.

  makeTestData                      % stops if data/test exists
  makeTestData('Overwrite', true)   % replace it (only to re-freeze on purpose)
%}

arguments
    opts.Overwrite (1,1) logical = false
end

root = fullfile(fileparts(mfilename('fullpath')), '..');
addpath(fullfile(root, 'functions'));
oldDir = cd(root); restoreDir = onCleanup(@() cd(oldDir));
P = readTestParams(root);                           % the reference values and dataDir = data/test
dst = char(P.dataDir);
src = 'data';
if isfolder(dst) && ~opts.Overwrite
    error('makeTestData:exists', ['%s exists already (frozen test data). Use ' ...
        'makeTestData(''Overwrite'', true) only to re-freeze on purpose.'], dst);
end

% ---- 1. the source must be the reference configuration ---------------------------------
in = @(f) loadInfo(fullfile(src, f));
g  = in("test.dat");             d  = in("test_dispersed.dat");  r  = in("test_rx.dat");
q  = in("test_rx_IQ.dat");       dd = in("test_IQ_dedispersed.dat");
e  = in("test_envelope.dat");
c  = in("chan/test_rx_IQ_chan"); dc = in("chan/test_dedisp_chan");
p  = in("chan/test_dedisp_chan_power.dat");
same = @(a, b) strcmp(string(a), string(b));
checks = {
    'pulsar (T, f_in, L, A, duty, mode, seed)', g.T == P.T && g.fs == P.f_in && g.L == P.L && ...
        g.A == P.A && g.dutycycle == P.dutycycle && strcmp(g.envelopeMode, P.genEnvMode) && g.seed == P.seed
    'dispersion (DM, band, refFreq)', d.DM == P.DM && d.fLow == P.fLow && d.fHigh == P.fHigh && ...
        d.refFreq == P.refFreq && same(d.inFile, g.file)
    'receiver noise (SNR, seed, no RFI)', isequal(r.snrDB, P.snrDB) && r.seed == P.noiseSeed && ...
        isempty(r.rfi) && same(r.inFile, d.file)
    'IQ (fLO, fs, filter)', q.fLO == P.fLO && q.actualFsOut == P.f_in / round(P.f_in / P.fs) && ...
        q.filterLen == P.filterOrder + 1 && same(q.inFile, r.file)
    'full-band dedispersion (DM, refFreq, band)', dd.DM == P.ephem.DM && dd.refFreq == P.refFreq && ...
        dd.fLow == P.fLow && dd.fHigh == P.fHigh && same(dd.inFile, q.file)
    'full-band power (f_out)', e.binLen == round(dd.fs / P.f_out) && same(e.inFile, dd.file)
    'channels (width, from the IQ file)', c.chanWidth == P.chanWidth && same(c.inFile, q.file)
    'channel dedispersion (DM, refFreq)', dc.DM == P.ephem.DM && dc.refFreq == P.refFreq && ...
        same(dc.inFile, c.file)
    'channel power (f_out)', p.binLen == round(dc.fs / P.f_out) && same(p.inFile, dc.file)};
bad = checks(~[checks{:, 2}], 1);
if ~isempty(bad)
    error('makeTestData:source', ['The files in data/ are not the seed-43 reference of ' ...
        'testParams: %s. Run main.m with these settings first (and the tests that make the ' ...
        'channel files).'], strjoin(bad, '; '));
end
fprintf('makeTestData: data/ is the seed-43 reference (%d checks)\n', size(checks, 1));

% ---- 2. copy -----------------------------------------------------------------------------
if isfolder(dst), rmdir(dst, 's'); end
mkdir(fullfile(dst, 'chan')); mkdir(fullfile(dst, 'mc'));
for f = ["test.dat" "test_dispersed.dat" "test_rx.dat"]                 % info only
    copyInfo(src, dst, f);
end
for f = ["test_rx_IQ.dat" "test_IQ_dedispersed.dat" "test_envelope.dat"] % data + info
    copyfile(fullfile(src, f), fullfile(dst, f));
    copyInfo(src, dst, f);
end
for base = ["chan/test_rx_IQ_chan" "chan/test_dedisp_chan"]              % channel sets
    copyfile(fullfile(src, base + "_ch*.dat"), fullfile(dst, 'chan'));
    copyInfo(src, dst, base);
end
copyfile(fullfile(src, 'chan', 'test_dedisp_chan_power.dat'), fullfile(dst, 'chan'));
copyInfo(src, dst, "chan/test_dedisp_chan_power.dat");
for f = ["foldRef_pre3c.mat" "toaRef_preA2.mat" "rfiRef_preB5a.mat"]   % frozen references
    copyfile(fullfile(src, 'mc', f), fullfile(dst, 'mc', f));
end
s = dir(fullfile(dst, '**', '*.*'));
fprintf('makeTestData: %d files, %.2f GB in %s\n', nnz(~[s.isdir]), sum([s.bytes]) / 1e9, dst);
fprintf('makeTestData: now run run(''tests/runAllTests.m'')\n');
end


% =========================================================================================
function P = readTestParams(root)
% testParams' variables as the fields of P. Read through a struct because a name such as
% dutycycle is also a function on the path (Signal Processing Toolbox): used directly in
% makeTestData it would be parsed as that function, not as the script's variable.
makingTestData = true; %#ok<NASGU>                  % testParams: data/test may not exist yet
run(fullfile(root, 'tests', 'testParams.m'));
P = struct();
for v = setdiff(who, {'P', 'root', 'makingTestData'}).'
    P.(v{1}) = eval(v{1});
end
end

function copyInfo(src, dst, f)
% <f>_info.mat from src to dst, with every top-level path field starting with 'data/'
% (file, inFile, chanFiles, ...) pointing into dst instead.
[d, name] = fileparts(f);
S = load(fullfile(src, d, name + "_info.mat"), 'info');
info = S.info;
for fn = fieldnames(info).'
    v = info.(fn{1});
    if (ischar(v) || isstring(v)) && any(startsWith(string(v), "data/"))
        nv = replace(string(v), regexpPattern("^data/"), string(dst) + "/");
        if ischar(v), nv = char(nv); end
        info.(fn{1}) = nv;
    end
end
save(fullfile(dst, d, name + "_info.mat"), 'info');
end
