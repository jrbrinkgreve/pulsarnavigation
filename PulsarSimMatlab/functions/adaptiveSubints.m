function [foldM, infoM, sel] = adaptiveSubints(fold, info_fold, template, Bnoise, opts)
%ADAPTIVESUBINTS  Choose the sub-int length from the data: merge short sub-ints until SNR ~7.
%{
Observer side (block C3c, notes §10 item 6). The SNR of a sub-int follows the pulsar's flux,
which scintillation changes by large factors from hour to hour (physics.html#scint). A fixed
sub-int length is then too long when the pulsar is bright (TOAs wasted) and too short when it
is faint (no TOA). So the data are folded into short sub-ints, and consecutive ones are merged
until the merged profile reaches the target SNR - the SNR per TOA of the design rule (~7,
detection.html), where every TOA is good.

  [foldM, infoM, sel] = adaptiveSubints(fold, info_fold, template, Bnoise)
  ... = adaptiveSubints(..., 'TargetSNR', 7, 'Weighting', 'optimal')
  ... = adaptiveSubints(..., 'AcrossGaps', true)     % groups may span gaps between passes

Steps:
  1. Runs: consecutive sub-ints whose turns follow without a gap (turnFirst(s+1) =
     turnLast(s) + 1); a gap between passes ends a run - nothing is merged across it.
     With 'AcrossGaps' all sub-ints form one run: a group continues in the next pass.
  2. Within a run, from its first sub-int on: add sub-ints one at a time to a running sum
     and compute the matched-filter statistic T0 of the merged profile at the predicted
     phase (detectPulsar, known phase: the SNR as the noise sees it). As soon as T0 reaches
     TargetSNR (or MaxSubints sub-ints are merged) the group is closed and the next starts.
  3. The sub-ints left at the end of a run form a last group below the target (kept;
     detectPulsar and the good-TOA rule decide about it; sel.complete is false for it).
  4. mergeSubints makes the merged fold.

Inputs:
  fold, info_fold  short sub-ints (foldProfile or simulateFold)
  template         [NBin x 1] (as estimateTOA)
  Bnoise           [Hz] as estimateTOA / detectPulsar (scalar or per channel)

Name-value options:
  'TargetSNR'  T0 at which a group is closed (default 7)
  'Weighting'  'equal' (default) or 'optimal', passed to detectPulsar
  'MaxSubints' at most this many sub-ints per group (default Inf); with 'TargetSNR' Inf:
               fixed groups of MaxSubints sub-ints (for comparison)
  'AcrossGaps' false (default): groups end at a gap between passes; true: a group may span
               gaps (20 Oct 2026) - needed when one TOA takes longer than a pass (J0437 in
               the worst case: ~15 h against 8-15 h passes). Physically this needs a phase
               model that stays coherent over the gap (motion, barycentric correction:
               notes §10 item 4); the fast simulator's ephemeris equals the truth, so there
               it is exact.
  'Verbose'    print a summary (default true)

Outputs:
  foldM, infoM  the merged fold (mergeSubints)
  sel           group [nSub x 1], nMerged [nG x 1], T0 [nG x 1] (at closing), complete
                [nG x 1] (reached the target), run [nG x 1] (all 1 with 'AcrossGaps')

Note on the stopping rule: a group closes when its noise happens to push T0 up, so the
SNR of the closed groups is biased up a little; the arrival time is not (with a symmetric
template the phase estimate is nearly uncorrelated with the amplitude) - checked by the TOA
pulls in testSimulateFold 12.
%}

arguments
    fold      struct
    info_fold struct
    template  double
    Bnoise    double {mustBePositive}
    opts.TargetSNR  (1,1) double {mustBePositive} = 7
    opts.Weighting  {mustBeTextScalar} = 'equal'
    opts.MaxSubints (1,1) double {mustBePositive} = Inf
    opts.AcrossGaps (1,1) logical = false
    opts.Verbose    (1,1) logical = true
end

nSub = size(fold.weight, 2);
S = fold.subint;
newRun = [true; S.turnFirst(2:end) ~= S.turnLast(1:end-1) + 1];
if opts.AcrossGaps, newRun(2:end) = false; end      % one run: groups continue over gaps
runId = cumsum(newRun);
nW = size(fold.weight, 3); D = size(fold.weightX, 4); nChan = size(fold.sum, 3);
N = size(fold.weight, 1);
group = zeros(nSub, 1);
G = 0; T0g = []; nMg = []; compl = []; runG = [];
info1 = info_fold; info1.nSub = 1;
s = 1;
while s <= nSub
    G = G + 1;
    acc = struct('sum', zeros(N, 1, nChan), 'weight', zeros(N, 1, nW), ...
        'weight2', zeros(N, 1, nW), 'weightX', zeros(N, 1, nW, D));
    n = 0; T0 = -Inf; done = false;
    while s <= nSub && ~done
        acc.sum     = acc.sum     + fold.sum(:, s, :);
        acc.weight  = acc.weight  + fold.weight(:, s, :);
        acc.weight2 = acc.weight2 + fold.weight2(:, s, :);
        acc.weightX = acc.weightX + fold.weightX(:, s, :, :);
        group(s) = G; n = n + 1;
        f1 = oneSubint(acc, S, s, nChan, nW);
        d1 = detectPulsar(f1, info1, template, 'Bnoise', Bnoise, 'Weighting', opts.Weighting, ...
            'Verbose', false);
        T0 = d1.T0(1);
        lastOfRun = s == nSub || runId(s + 1) ~= runId(s);
        done = T0 >= opts.TargetSNR || n >= opts.MaxSubints || lastOfRun;
        s = s + 1;
    end
    T0g(G, 1) = T0; nMg(G, 1) = n; compl(G, 1) = T0 >= opts.TargetSNR; runG(G, 1) = runId(s - 1); %#ok<AGROW>
end
[foldM, infoM] = mergeSubints(fold, info_fold, group);
sel = struct('group', group, 'nMerged', nMg, 'T0', T0g, 'complete', logical(compl), 'run', runG);
if opts.Verbose
    fprintf(['adaptiveSubints: %d sub-ints in %d run(s) -> %d groups (%d reached T0 >= %g); ' ...
             'sub-ints per group median %g, range %d-%d\n'], nSub, max(runId), G, nnz(compl), ...
        opts.TargetSNR, median(nMg), min(nMg), max(nMg));
end
end


function f1 = oneSubint(acc, S, s, nChan, nW)
%ONESUBINT  A one-sub-int fold struct from the running sums (for detectPulsar).
f1 = acc;
wFull = repmat(acc.weight, 1, 1, nChan / nW);
f1.prof = acc.sum ./ wFull;
f1.prof(wFull == 0) = NaN;
f1.profTotal = reshape(f1.prof, size(f1.prof, 1), nChan);
f1.subint = struct('turnFirst', S.turnFirst(s), 'turnLast', S.turnLast(s), 'turnRef', S.turnRef(s), ...
    'tRef', S.tRef(s), 'fRef', S.fRef(s), 'tMean', S.tMean(s), 'nTimeBins', S.nTimeBins(s));
end
