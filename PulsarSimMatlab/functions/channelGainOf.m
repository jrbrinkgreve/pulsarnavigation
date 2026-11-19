function g = channelGainOf(opt, fold)
%CHANNELGAINOF  A 'ChannelGain' option as the row vector that estimateTOA and detectPulsar take.
%{
The 'optimal' fit and detection model the pulse in channel c as b * s_c * T, with s_c the
relative pulse amplitude of the channel ('ChannelGain', default 1 in every channel). foldToTOA
and adaptiveSubints take the option in three forms; this turns each into what they pass on.

  g = channelGainOf(opt, fold)

Inputs:
  opt   [] : not given. g = [], and the caller passes nothing on (s_c = 1, as before
             17 Nov 2026).
        [1 x nChan] values: gains the observer knows, e.g. the pulsar's spectrum, or a
             band-edge response that shapes the pulse but not the noise. Returned as a
             row; estimateTOA / detectPulsar check the values (>= 0, finite, one > 0).
        "measured": each channel's mean power, the mean of fold.profTotal over the phase
             bins, normalized to mean 1 over the channels with data. Right for a receiver
             bandpass, which scales the pulse and the noise of a channel alike (risk
             register FM-41). A channel without data gets 0 (left out of the fit). It
             cannot see a gain on the pulse alone, and leftover RFI raises a channel's
             mean power.
  fold  the fold struct (foldProfile, simulateFold): prof [NBin x nSub x nChan] for the
        channel count, profTotal [NBin x nChan] for "measured".

Output:
  g     [] or [1 x nChan].

Observer information only: the fold itself and what the caller knows.
(17 Nov 2026 as a local function of foldToTOA; its own file 19 Nov 2026, for adaptiveSubints.)
%}

g = [];
if isempty(opt), return, end
nChan = size(fold.prof, 3);
if ischar(opt) || isstring(opt)
    if ~(isscalar(string(opt)) && strcmpi(opt, 'measured'))
        error('channelGainOf:gain', 'ChannelGain must be one value >= 0 per channel, or "measured".');
    end
    g = mean(fold.profTotal, 1, 'omitnan');         % [1 x nChan] mean power over the phase bins
    g(~isfinite(g) | g < 0) = 0;                    % a channel without data: left out of the fit
    if ~any(g > 0)
        error('channelGainOf:gain', 'ChannelGain "measured": no channel has data.');
    end
    g = g / mean(g(g > 0));
elseif isnumeric(opt)
    g = reshape(double(opt), 1, []);
    if numel(g) ~= nChan
        error('channelGainOf:gain', 'ChannelGain must have one value per channel (%d).', nChan);
    end
else
    error('channelGainOf:gain', 'ChannelGain must be one value >= 0 per channel, or "measured".');
end
end
