function tmpl = gaussianTemplate(nBin, fwhmTurns)
%GAUSSIANTEMPLATE  Periodic Gaussian profile, peak 1 at phase 0 (bin 1).
%{
Template for estimateTOA, on the phase grid of foldProfile.

  tmpl = gaussianTemplate(nBin, fwhmTurns)

Inputs:
  nBin       number of phase bins (as in foldProfile 'NBin')
  fwhmTurns  [turns] FWHM of the power profile (e.g. ephem.profileFWHM)

Output:
  tmpl  [nBin x 1], bin j (1-based) at phase (j-1)/nBin, peak 1 at bin 1.
        Periodic: the images at -2..+2 turns are summed, so the half at
        negative phase wraps to the end of the array, as in the fold.

This is one choice of template (the synthetic pulsar is Gaussian).
estimateTOA accepts any [nBin x 1] profile with phase 0 at bin 1, e.g. a
multi-component model or a high-SNR observed profile for real data.
%}
ph  = (0:nBin-1).' / nBin;
sig = fwhmTurns / (2*sqrt(2*log(2)));
tmpl = zeros(nBin, 1);
for j = -2:2                                   % periodic images
    tmpl = tmpl + exp(-0.5*((ph - j)/sig).^2);
end
tmpl = tmpl / max(tmpl);
end
