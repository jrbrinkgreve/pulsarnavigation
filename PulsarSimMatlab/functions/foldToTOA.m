function [toa, detection, info] = foldToTOA(fold, info_fold, profileFWHM, Bnoise, opts)
%FOLDTOTOA  TOAs, detection and the good-TOA rule from a fold (the back end of the processing).
%{
The last processing step: template matching and detection on a fold, and the
quality rule that decides which TOAs are good. main.m calls it after
processIQ (IQ -> fold; out.Bnoise, out.weighting are the inputs here); kept
as its own function so that a fold from elsewhere (the
fast simulator of block C, simulateFold, which draws the fold struct
directly) goes through exactly the same code. Uses only observer
information: the fold, the template width from the ephemeris and the noise
bandwidth of the observer's own filters.

  [toa, detection, info] = foldToTOA(fold, info_fold, profileFWHM, Bnoise, ...
                                     'Weighting', 'optimal')
  ... = foldToTOA(..., 'Weighting', 'optimal', 'ChannelGain', "measured")

Inputs:
  fold, info_fold  from foldProfile (or the same fields from a simulator)
  profileFWHM      [turns] power-profile FWHM (ephem.profileFWHM), or [nComp x 3]
                   components [amplitude, centre, FWHM] (20 Oct 2026): the
                   template is gaussianTemplate(info_fold.NBin, profileFWHM)
  Bnoise           [Hz] noise-equivalent bandwidth for the radiometer noise
                   model (estimateTOA, detectPulsar): full band
                   noiseBandwidth(fLow, fHigh, edgeWidth) of the
                   dedispersion; channels BnoiseChan of one channel
                   (detectChannels' info.noise.Bnoise)

Name-value options:
  'Weighting'  how channels are combined, 'equal' (default) or 'optimal';
               passed to estimateTOA and detectPulsar (one channel: no effect)
  'ChannelGain' relative pulse amplitude per channel, s_c in the 'optimal' fit and
               detection (prof_cj = a_c + b s_c T + noise); 'equal' does not use it.
               [] (default): not passed on, so s_c = 1 in every channel, as before
               17 Nov 2026 (results unchanged bit for bit).
               [1 x nChan] values >= 0: gains the observer knows, e.g. the pulsar's
               spectrum, or a band-edge response that shapes the pulse but not the
               noise (a channel with gain 0 is left out).
               "measured": each channel's mean power (mean of fold.profTotal over the
               phase bins), normalized to mean 1. Right for a receiver bandpass, which
               scales pulse and noise of a channel alike (risk register FM-41; a channel
               without data gets 0). It cannot see a gain on the pulse alone, and
               leftover RFI raises a channel's mean power.
  'Verbose'    print the stages' summaries and the good-TOA line (default true)

Steps:
  1. template = gaussianTemplate(NBin, profileFWHM)
  2. estimateTOA: TOA, error bar, red. chi^2 (flagChi2), pChi2 per sub-int
  3. detectPulsar: Neyman-Pearson detection, known and unknown phase
  4. toa.good = valid & detected (unknown phase) & not chi^2-flagged: the
     pulsar is found without the ephemeris phase, and the profile is
     template + noise in the gross sense (MaxRedChi2 2; pChi2 is reported,
     it does not decide; notes §5.8)

Outputs:
  toa        estimateTOA's struct plus toa.good [nSub x 1 logical]
  detection  detectPulsar's struct
  info       template, Bnoise, weighting, channelGain ([] or the [1 x nChan] values
             passed on), toa (estimateTOA's info), detect (detectPulsar's info)

(19 Oct 2026: moved out of main.m, next to processIQ; same calls. 17 Nov 2026:
'ChannelGain'.)
%}

arguments
    fold        struct
    info_fold   struct
    profileFWHM double {mustBeNonempty, mustBeFinite}
    Bnoise      double {mustBePositive}
    opts.Weighting {mustBeTextScalar} = 'equal'
    opts.ChannelGain = []
    opts.Verbose   (1,1) logical = true
end

template = gaussianTemplate(info_fold.NBin, profileFWHM);
toaArgs  = {'Bnoise', Bnoise, 'Weighting', opts.Weighting, 'Verbose', opts.Verbose};
gain = channelGainOf(opts.ChannelGain, fold);       % [] = not given: the fit's own default (1)
if ~isempty(gain)
    toaArgs = [toaArgs, {'ChannelGain', gain}];
    if opts.Verbose
        fprintf('foldToTOA: channel gains %.3g ... %.3g passed to the fit and the detection\n', ...
            min(gain), max(gain));
    end
end
[toa, info_toa] = estimateTOA(fold, info_fold, template, toaArgs{:});

% Detection (Neyman-Pearson, P_FA 1e-3 per sub-int) and TOA quality.
% Good TOA: fitted, pulsar detected without using the ephemeris phase
% (unknown phase), and the fit residuals consistent with the noise (chi^2).
[detection, info_detect] = detectPulsar(fold, info_fold, template, toaArgs{:});
toa.good = toa.valid & detection.detectedUnknown & ~toa.flagChi2;
if opts.Verbose
    fprintf('foldToTOA: %d of %d fitted TOAs good (%d not detected, %d chi2-flagged)\n', ...
        nnz(toa.good), nnz(toa.valid), nnz(toa.valid & ~detection.detectedUnknown), ...
        nnz(toa.valid & toa.flagChi2));
end

info = struct('template', template, 'Bnoise', Bnoise, 'weighting', char(opts.Weighting), ...
    'channelGain', gain, 'toa', info_toa, 'detect', info_detect);
end

% =====================================================================================
function g = channelGainOf(opt, fold)
%CHANNELGAINOF  The 'ChannelGain' option as a row vector for estimateTOA and detectPulsar;
% [] when it was not given. "measured": each channel's mean power, normalized to mean 1.
g = [];
if isempty(opt), return, end
nChan = size(fold.prof, 3);
if ischar(opt) || isstring(opt)
    if ~(isscalar(string(opt)) && strcmpi(opt, 'measured'))
        error('foldToTOA:gain', 'ChannelGain must be one value >= 0 per channel, or "measured".');
    end
    g = mean(fold.profTotal, 1, 'omitnan');         % [1 x nChan] mean power over the phase bins
    g(~isfinite(g) | g < 0) = 0;                    % a channel without data: left out of the fit
    if ~any(g > 0)
        error('foldToTOA:gain', 'ChannelGain "measured": no channel has data.');
    end
    g = g / mean(g(g > 0));
elseif isnumeric(opt)
    g = reshape(double(opt), 1, []);
    if numel(g) ~= nChan
        error('foldToTOA:gain', 'ChannelGain must have one value per channel (%d).', nChan);
    end
else
    error('foldToTOA:gain', 'ChannelGain must be one value >= 0 per channel, or "measured".');
end
end
