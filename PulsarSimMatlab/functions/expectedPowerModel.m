function M = expectedPowerModel(info_gen, info_disp, info_IQ, info_dedisp, info_rx)
%EXPECTEDPOWERMODEL  Ground-truth mean and variance of the detected power.
%{
Used by the check/plot functions only (it uses ground truth). After IQ
conversion and dedispersion, the complex baseband power spectral density is

  S(f,t) = s_s * p(t) * Wf(f)^2 * Wi(f)^2  +  s_n * Wi(f)^2

  signal: s_s = g^2 A^2 / fsIn, shaped by the forward (dispersion) and the
          inverse (dedispersion) band tapers Wf, Wi; p(t) = unit-peak power
          profile of the generator
  noise:  s_n = g^2 sigma_n^2 / fsIn; receiver noise was added AFTER
          dispersion, so it only sees the dedispersion taper Wi
  (g = IQ gain factor, fsIn = RF sample rate; the IQ low-pass is flat
   over the band)

Mean power per sample:   Ps(t) = M.sigScale * p(t),  Pn = M.Pn,
                          E[P] = Ps + Pn
Variance of a bin average over dt (Gaussian fields, dt >> 1/B):
  var = int S^2 df / dt
      = (css*Ps^2 + 2*csn*Ps*Pn + cnn*Pn^2) / dt       -> M.var(Ps, dt)
with css = int (Wf Wi)^4 / Js^2, csn = int Wf^2 Wi^4 / (Js Jn),
     cnn = int Wi^4 / Jn^2, Js = int (Wf Wi)^2, Jn = int Wi^2.
Without noise this reduces to the earlier self-noise model
(relative std 1/sqrt(dt * Bnoise)). RFI is NOT included: in the checks it
shows up as a deviation from this model.

  M = expectedPowerModel(info_gen, info_disp, info_IQ, info_dedisp)            % no noise
  M = expectedPowerModel(info_gen, info_disp, info_IQ, info_dedisp, info_rx)   % with noise

Best achievable per-pulse SNR and TOA error (Fisher information of the
detected power under this model; total-power detection, optimal time
weighting), with v = css*Ps^2 + 2*csn*Ps*Pn + cnn*Pn^2 (dt cancels):
  SNR^2 = int Ps^2 / v dt,     1/sigma_TOA^2 = int Ps'^2 / v dt
For a flat band (css = csn = cnn = 1/B) this is the flat-band prediction of
addNoiseAndRFI; with the real tapers it is ~1.5 % less optimistic.
Only with receiver noise (noise-free, the integrals diverge in the tails).

Fields: sigScale, Pn, css, csn, cnn, Bsignal (=1/css), Bnoise (=1/cnn),
        var(Ps, dt), envelope(t) (unit-peak power profile p(t)),
        snrPulse, toaErrPulse [s] (NaN without noise), hasNoise, hasRFI.
%}

arguments
    info_gen    struct
    info_disp   struct
    info_IQ     struct
    info_dedisp struct
    info_rx     struct = struct([])
end

fLow = info_dedisp.fLow; fHigh = info_dedisp.fHigh;
f  = linspace(fLow, fHigh, 200001);
Wf = taperW(f, info_disp.fLow, info_disp.fHigh, info_disp.edgeWidth);
Wi = taperW(f, fLow, fHigh, info_dedisp.edgeWidth);
Js  = trapz(f, (Wf.*Wi).^2);
Jn  = trapz(f, Wi.^2);
Iss = trapz(f, (Wf.*Wi).^4);
Isn = trapz(f, Wf.^2 .* Wi.^4);
Inn = trapz(f, Wi.^4);

g2   = info_IQ.gainFactor^2;
fsIn = info_IQ.fsIn;
sigN = 0;
if ~isempty(info_rx) && isfield(info_rx, 'noiseStd'), sigN = info_rx.noiseStd; end

M = struct();
M.sigScale = g2 * info_gen.A^2 * Js / fsIn;
M.Pn       = g2 * sigN^2 * Jn / fsIn;
M.css = Iss / Js^2;
M.csn = Isn / (Js * Jn);
M.cnn = Inn / Jn^2;
M.Bsignal = 1 / M.css;
M.Bnoise  = 1 / M.cnn;
M.hasNoise = M.Pn > 0;
M.hasRFI   = ~isempty(info_rx) && isfield(info_rx, 'rfi') && ~isempty(info_rx.rfi);
css = M.css; csn = M.csn; cnn = M.cnn; Pn = M.Pn;
M.var      = @(Ps, dt) (css*Ps.^2 + 2*csn*Ps.*Pn + cnn*Pn.^2) ./ dt;
M.envelope = @(t) envelopePower(t, info_gen);

% Best achievable per-pulse SNR and TOA error (see header)
M.snrPulse = NaN; M.toaErrPulse = NaN;
if M.hasNoise
    sp = info_gen.sigma;                            % std of the power profile
    if ~strcmpi(info_gen.envelopeMode, 'power'), sp = sp / sqrt(2); end
    t   = linspace(-8*sp, 8*sp, 40001);
    Ps  = M.sigScale * exp(-0.5*(t/sp).^2);
    dPs = -t/sp^2 .* Ps;
    v   = css*Ps.^2 + 2*csn*Ps*Pn + cnn*Pn^2;
    M.snrPulse    = sqrt(trapz(t, Ps.^2 ./ v));
    M.toaErrPulse = 1 / sqrt(trapz(t, dPs.^2 ./ v));
end
end


function p = envelopePower(t, info_gen)
%ENVELOPEPOWER  Generator ground truth: unit-peak periodic power profile.
T = info_gen.T; sig = info_gen.sigma;
nN = max(1, ceil(6*sig/T + 0.5));
peakNorm = sum(exp(-0.5*((-nN:nN)*T/sig).^2));
k0 = round(t/T - 0.5);
G = zeros(size(t));
for j = -nN:nN
    G = G + exp(-0.5*((t - (k0 + j + 0.5)*T)/sig).^2);
end
G = G / peakNorm;
if strcmpi(info_gen.envelopeMode, 'power')
    p = G;
else
    p = G.^2;
end
end


function W = taperW(f, fLow, fHigh, e)
W = zeros(size(f));
ib = f >= fLow & f <= fHigh;
W(ib) = 1;
lo = ib & f < fLow + e;   W(lo) = sin(pi/2 * (f(lo) - fLow) / e).^2;
hi = ib & f > fHigh - e;  W(hi) = sin(pi/2 * (fHigh - f(hi)) / e).^2;
end
