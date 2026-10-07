function ch = combineChannels(prof, W, W2, WX, B)
%COMBINECHANNELS  Channels used for one folded profile, their sum and their weights.
%{
Shared by estimateTOA and detectPulsar. For one profile (a sub-integration
or the total fold) of a fold with several channels:

Exclusion rule: a channel is used only if it has data (W > 0) in every phase
bin where any channel has data, so every bin of the sum holds the same
channels (a channel missing in some bins would leave a dip in the sum and
pull the TOA). Channels with partial data (0 < W) stay in. When all channels
share one set of weights (nW = 1), all channels are used.

  ch = combineChannels(prof, W, W2, WX, B)

Inputs:
  prof  [N x nChan]        channel profiles (fold.prof of one sub-int, or
                           fold.profTotal)
  W, W2 [N x nW]           fold.weight, fold.weight2 of that profile
  WX    [N x nW x nLag]    fold.weightX of that profile
  B     Bnoise: scalar, [] or one value per channel

Output ch: p [N x 1] sum of the used channel profiles; Pc [N x nUsed] the
  used profiles; W, W2, WX, B of the used channels (unchanged when nW = 1);
  have [N x 1] phase bins with data; nUsed; coverage (fraction of bins with
  data; 0 when no channel is used).
%}

nChan = size(prof, 2);
have  = any(W > 0, 2);
if size(W, 2) == 1
    use = true(1, nChan);
else
    use = all(W(have, :) > 0, 1);
    W = W(:, use); W2 = W2(:, use); WX = WX(:, use, :);
end
if numel(B) > 1, B = B(use); end
Pc = prof(:, use);
ok = any(use) && any(have);
ch = struct('p', sum(Pc, 2), 'Pc', Pc, 'W', W, 'W2', W2, 'WX', WX, 'B', B, ...
    'have', have & ok, 'nUsed', nnz(use) * ok, 'coverage', mean(have) * ok);
end
