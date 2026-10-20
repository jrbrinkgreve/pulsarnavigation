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
  info       template, Bnoise, weighting, toa (estimateTOA's info),
             detect (detectPulsar's info)

(19 Oct 2026: moved out of main.m, next to processIQ; same calls.)
%}

arguments
    fold        struct
    info_fold   struct
    profileFWHM double {mustBeNonempty, mustBeFinite}
    Bnoise      double {mustBePositive}
    opts.Weighting {mustBeTextScalar} = 'equal'
    opts.Verbose   (1,1) logical = true
end

template = gaussianTemplate(info_fold.NBin, profileFWHM);
toaArgs  = {'Bnoise', Bnoise, 'Weighting', opts.Weighting, 'Verbose', opts.Verbose};
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
    'toa', info_toa, 'detect', info_detect);
end
