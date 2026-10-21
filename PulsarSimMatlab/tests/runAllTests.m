function results = runAllTests(names)
%RUNALLTESTS  Run the unit tests in tests/ one after another in this session; summary at the end.
%{
Run from the PulsarSimMatlab folder: run('tests/runAllTests.m') (all), or
runAllTests(["testRfiNoise" "testRadarEdges"]) with tests/ on the path (a subset). 6-9 min
for all 20 (testOptimalTOA and testBlankingWeights 2-4 min each; testSimulateFold ~1 min,
testRFIModel ~20 s).

Each test runs in its own try/catch, so one failure does not stop the others; its own
output is printed as usual. At the end a table: PASS / FAIL per test, run time, and the
first line of the error. Errors at the end if any test failed (like the tests themselves),
so a batch run (matlab -batch) exits with an error.

Data: the tests read only their frozen copy in data/test/ with their own settings
(tests/testParams.m, made once by tests/makeTestData.m), not main.m's files, so main.m
runs and pipelineParams edits do not change them.

Order: pipeline order, and the tests that make shared files before the ones that read them:
testGeneratePulsarSignal (own files only) first; testChannelizeIQ and testDedisperseChannels
rewrite the seed-43 channel files in data/test/chan/, testDetectChannels the 128-channel power
file; the rest read them (or make their own data in tempdir). Not run: demo*, exp*
(demonstrations, experiments), make* (they freeze data and references) and testParams (the
settings). A test*.m file that is not in the list below runs at the end, with a note.
Slow tests (testRFIChain: a 1 s voltage-chain run per case, ~75 min per case, ~2.5 h for its two
default cases) are skipped in the default run and run only when named: runAllTests("testRFIChain").

Never run two test sessions at once: they share data/test/chan/ and their tempdir folders
(notes §8). main.m at the same time is fine now (it does not touch data/test/).
%}

arguments
    names (1,:) string = strings(1, 0)
end

order = ["testGeneratePulsarSignal", "testChannelizeIQ", "testDedisperseChannels", "testDetectChannels", ...
         "testFoldWeights", "testChannelTOA", "testOptimalTOA", "testBlankingWeights", ...
         "testBlankChannels", "testDetectRFI", "testPeriodicRFI", "testRadarEdges", ...
         "testRfiGating", "testRfiNoise", "testSimulateFold", "testRFIModel", "testMakeTrajectory", "testDelayModel", "testMakePredictor", "testFoldPredictor"];
testDir = fileparts(mfilename('fullpath'));
found = string({dir(fullfile(testDir, 'test*.m')).name});
found = setdiff(erase(found, ".m"), "testParams");  % the settings script, not a test
slow  = "testRFIChain";                             % only when named (~2.5 h)
extra = setdiff(found, [order, slow]);              % new tests not yet in the list
order = [order(ismember(order, found)), extra, slow(ismember(slow, found))];
if ~isempty(names)
    bad = setdiff(names, order);
    if ~isempty(bad)
        error('runAllTests:name', 'Unknown test(s): %s.', strjoin(bad, ', '));
    end
    order = order(ismember(order, names));          % the listed ones, in dependency order
else
    skipped = order(ismember(order, slow));
    order = order(~ismember(order, slow));
    for k = 1:numel(skipped)
        fprintf('runAllTests: skipped (slow): %s - run runAllTests("%s")\n', skipped(k), skipped(k));
    end
end

oldDir = cd(fullfile(testDir, '..'));               % the tests expect the PulsarSimMatlab folder
restoreDir = onCleanup(@() cd(oldDir));

n = numel(order);
results = struct('name', cellstr(order), 'passed', false, 'seconds', NaN, 'message', '');
tAll = tic;
for i = 1:n
    fprintf('\n===== %s (%d of %d) =====\n', order(i), i, n);
    if ismember(order(i), extra)
        fprintf('(not in the list of runAllTests: run last)\n');
    end
    t = tic;
    try
        run(char(fullfile(testDir, order(i) + ".m")));
        results(i).passed = true;
    catch e
        results(i).message = strtok(e.message, newline);
        fprintf(2, '%s: %s\n', order(i), e.message);
    end
    results(i).seconds = toc(t);
end

fprintf('\n===== runAllTests: %d of %d passed (%.1f min) =====\n', nnz([results.passed]), n, toc(tAll) / 60);
for i = 1:n
    status = 'PASS';
    if ~results(i).passed, status = 'FAIL'; end
    fprintf('  %-24s %s  %6.1f s  %s\n', results(i).name, status, results(i).seconds, results(i).message);
end
if ~all([results.passed])
    error('runAllTests:failed', '%d test(s) failed: %s.', nnz(~[results.passed]), ...
        strjoin({results(~[results.passed]).name}, ', '));
end
end
