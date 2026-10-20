function tmpl = gaussianTemplate(nBin, profile)
%GAUSSIANTEMPLATE  Periodic Gaussian profile (one or several components), peak 1, phase 0 at bin 1.
%{
Template for estimateTOA, on the phase grid of foldProfile.

  tmpl = gaussianTemplate(nBin, fwhmTurns)
  tmpl = gaussianTemplate(nBin, components)        % [amplitude, centre, FWHM] per row

Inputs:
  nBin       number of phase bins (as in foldProfile 'NBin')
  profile    [turns] FWHM of a single Gaussian power profile (e.g. ephem.profileFWHM),
             or [nComp x 3] components [amplitude, centre (turns), FWHM (turns)] (20 Oct
             2026; e.g. B0329+54's five components): the template is their sum. Phase 0 is
             the components' reference point (centre 0), which is not necessarily the peak:
             a TOA is the arrival of that reference point (the same convention as
             simulateFold 'Components').

Output:
  tmpl  [nBin x 1], bin j (1-based) at phase (j-1)/nBin, maximum 1 (at bin 1 for a single
        Gaussian). Periodic: the images at -2..+2 turns are summed, so the half at
        negative phase wraps to the end of the array, as in the fold.

estimateTOA accepts any [nBin x 1] profile with phase 0 at bin 1, e.g. a high-SNR observed
profile for real data.
%}
ph  = (0:nBin-1).' / nBin;
if isscalar(profile)                           % one Gaussian at phase 0
    sig = profile / (2*sqrt(2*log(2)));
    tmpl = zeros(nBin, 1);
    for j = -2:2                               % periodic images
        tmpl = tmpl + exp(-0.5*((ph - j)/sig).^2);
    end
else
    if size(profile, 2) ~= 3 || any(profile(:, 3) <= 0) || ~all(isfinite(profile), 'all')
        error('gaussianTemplate:components', ...
            'Components must be [nComp x 3] rows [amplitude, centre, FWHM] (turns, FWHM > 0).');
    end
    tmpl = zeros(nBin, 1);
    for i = 1:size(profile, 1)
        sig = profile(i, 3) / (2*sqrt(2*log(2)));
        for j = -2:2
            tmpl = tmpl + profile(i, 1) * exp(-0.5*((ph - profile(i, 2) - j)/sig).^2);
        end
    end
end
tmpl = tmpl / max(tmpl);
end
