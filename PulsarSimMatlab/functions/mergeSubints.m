function [foldM, infoM] = mergeSubints(fold, info_fold, group)
%MERGESUBINTS  Merge groups of consecutive sub-integrations of a fold exactly (observer side).
%{
Sums and weights of a fold add, so sub-ints can be combined after folding without loss:
the merged sub-int's sum, weight, weight2 and weightX are the sums over its members, and its
profile is sum ./ weight. Used to build longer sub-ints from short ones (adaptiveSubints,
block C3c: sub-int length chosen from the data). Works for foldProfile's and simulateFold's
folds alike (one or more channels, shared or per-channel weights).

  [foldM, infoM] = mergeSubints(fold, info_fold, group)

Inputs:
  fold, info_fold  from foldProfile or simulateFold
  group            [nSub x 1] group number of every sub-int: 1, 1, 1, 2, 2, 3, ... - each
                   group a run of consecutive sub-ints, numbered in order (a sub-int may be
                   left out with 0 or NaN)

Outputs:
  foldM  the fold struct with one sub-int per group: sum, weight, weight2, weightX summed;
         prof = sum ./ weight (NaN where weight = 0); profTotal and phase unchanged;
         subint: turnFirst = first member's, turnLast = last member's, turnRef = their
         middle (round((turnFirst + turnLast)/2), as foldProfile), tRef / fRef from the
         phase model at turnRef, tMean weighted by the time bins, nTimeBins summed.
         A predictor fold (info_fold.predictor, foldProfile 'Predictor'; 23 Oct 2026): its
         phase is not the f0 / F1 model (f0 only labels it), so tRef / fRef come from the
         member sub-int that holds turnRef, whose tRef / fRef foldProfile solved on the
         predictor: tRef = tRef_m + (turnRef - turnRef_m) / fRef_m, fRef = fRef_m (exact
         but for the frequency change within half a member, (fdot/f) dt^2 / 2: ~45 ns for
         60 s sub-ints on the Earth). The f0 / F1 model would be off by ~(v.n/c) t, 0.3 s
         per hour of recording for B0329+54.
  infoM  info_fold with nSub, subintPeriods = NaN (variable), subintGroup = group

Exact except for one term foldProfile does not store either: the noise covariance between
time bins on both sides of a sub-int boundary (with correlated time bins, NoiseCoeffs lags
> 0, a merged sub-int lacks those few pairs at each internal boundary; independent time bins
are exact). Merging across a gap in the data (between passes) is allowed here; whether it
makes sense (phase coherence over the gap) is the caller's choice. For a predictor fold a
group whose turnRef falls in such a gap is refused ('mergeSubints:predictor'): its tRef
would need the predictor there (one recording's sub-ints are contiguous, so this does not
arise within one foldProfile fold).
%}

arguments
    fold      struct
    info_fold struct
    group     double {mustBeVector}
end

nSub = size(fold.weight, 2);
group = group(:);
if numel(group) ~= nSub
    error('mergeSubints:group', 'group must have one value per sub-int (%d).', nSub);
end
use = isfinite(group) & group > 0;
g = group(use);
idx = find(use);
if any(diff(g) < 0) || any(diff(g) > 1) || g(1) ~= 1
    error('mergeSubints:group', 'Groups must be numbered 1, 2, ... in sub-int order.');
end
for k = 1:max(g)
    m = idx(g == k);
    if any(diff(m) ~= 1)
        error('mergeSubints:group', 'Group %d is not a run of consecutive sub-ints.', k);
    end
end
nG = max(g);

foldM = fold;
foldM.sum     = groupSum(fold.sum, idx, g, nG);
foldM.weight  = groupSum(fold.weight, idx, g, nG);
foldM.weight2 = groupSum(fold.weight2, idx, g, nG);
foldM.weightX = groupSum(fold.weightX, idx, g, nG);
nChan = size(fold.sum, 3);
wFull = repmat(foldM.weight, 1, 1, nChan / size(foldM.weight, 3));
foldM.prof = foldM.sum ./ wFull;
foldM.prof(wFull == 0) = NaN;

% sub-int bookkeeping
S  = fold.subint;
tf = zeros(nG, 1); tl = zeros(nG, 1); nTB = zeros(nG, 1); tM = zeros(nG, 1);
for k = 1:nG
    m = idx(g == k);
    tf(k)  = S.turnFirst(m(1));
    tl(k)  = S.turnLast(m(end));
    nTB(k) = sum(S.nTimeBins(m));
    tM(k)  = sum(S.tMean(m) .* S.nTimeBins(m)) / max(nTB(k), 1);
end
tr = round((tf + tl) / 2);
if isfield(info_fold, 'predictor') && ~isempty(info_fold.predictor)
    tRef = zeros(nG, 1); fRef = zeros(nG, 1);       % predictor fold: from the member holding turnRef
    for k = 1:nG
        m = idx(g == k);
        j = m(S.turnFirst(m) <= tr(k) & S.turnLast(m) >= tr(k));
        if isempty(j)
            error('mergeSubints:predictor', ['Group %d''s reference turn %d lies in a gap between ' ...
                'its sub-ints: a predictor fold''s tRef needs the predictor there.'], k, tr(k));
        end
        tRef(k) = S.tRef(j) + (tr(k) - S.turnRef(j)) / S.fRef(j);
        fRef(k) = S.fRef(j);
    end
else
    f0 = info_fold.f0; F1 = info_fold.F1; Phi0 = 0;
    if isfield(info_fold, 'Phi0'), Phi0 = info_fold.Phi0; end
    d = tr - Phi0;
    dtRef = 2*d ./ (f0 + sqrt(f0^2 + 2*F1*d));      % stable root of phi(t) = turnRef
    tRef = info_fold.TRef + dtRef; fRef = f0 + F1*dtRef;
end
foldM.subint = struct('turnFirst', tf, 'turnLast', tl, 'turnRef', tr, ...
    'tRef', tRef, 'fRef', fRef, 'tMean', tM, 'nTimeBins', nTB);

infoM = info_fold;
infoM.nSub = nG;
infoM.subintPeriods = NaN;
infoM.subintGroup = group;
end


function Y = groupSum(X, idx, g, nG)
%GROUPSUM  Sum X(:, idx, ...) over the sub-ints of each group (dimension 2).
sz = size(X); sz(end+1:4) = 1;
Y = zeros([sz(1), nG, sz(3:end)]);
for k = 1:nG
    Y(:, k, :, :) = sum(X(:, idx(g == k), :, :), 2);
end
end
