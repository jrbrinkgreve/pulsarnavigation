function Bn = noiseBandwidth(fLow, fHigh, edgeWidths)
%NOISEBANDWIDTH  Noise-equivalent bandwidth of raised-cosine band tapers.
%{
  Bn = noiseBandwidth(fLow, fHigh, edgeWidths)

Bn = (int W^2 df)^2 / int W^4 df, for W the product of raised-cosine
(sin^2) band-edge tapers, each of the given width, inside [fLow, fHigh]
(the taper used by applyDispersionStream / applyInverseDispersion).

Why this quantity: after square-law detection, the mean of |z|^2 over a
time bin dt has relative variance 1/(Bn*dt) for Gaussian noise with power
spectrum ~ W^2. For a flat band Bn = fHigh - fLow; tapered edges lower it
slightly (8 MHz edges in 400 MHz: 391.6 MHz for one taper, 389.6 MHz for
two). estimateTOA uses it in the radiometer noise model.

Inputs:
  fLow, fHigh  [Hz] band edges
  edgeWidths   [Hz] taper width(s); several values = product of tapers
               (e.g. info_dedisp.edgeWidth for the observer's own taper)

For real data, a measured receiver bandpass can be multiplied into W in
the same way.
%}
f = linspace(fLow, fHigh, 200001);
W = ones(size(f));
for e = edgeWidths
    w = ones(size(f));
    lo = f < fLow + e;   w(lo) = sin(pi/2 * (f(lo) - fLow) / e).^2;
    hi = f > fHigh - e;  w(hi) = sin(pi/2 * (fHigh - f(hi)) / e).^2;
    W = W .* w;
end
Bn = trapz(f, W.^2)^2 / trapz(f, W.^4);
end
