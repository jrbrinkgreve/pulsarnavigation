function src = rfiSource(type, opts)
%RFISOURCE  Define one radio-frequency interference source for addNoiseAndRFI.
%{
All sources share the same fields, so several can be concatenated:
  rfi = [rfiSource('cw', 'Freq', 1.35e9, 'INRdB', -5), ...
         rfiSource('impulse', 'Rate', 50, 'Duration', 200e-9, 'INRdB', 20)];

Power: 'INRdB' is the interference-to-noise ratio: the RFI power (while it
is on) relative to the receiver-noise power in the analysis band
[fLow, fHigh] given to addNoiseAndRFI. 0 dB = as much power as all the
receiver noise in the band.

Types and their parameters:
  'cw'       continuous carrier            Freq [Hz], Drift [Hz/s], Phase [rad]
  'bpsk'     spread-spectrum (GNSS-like)   Freq [Hz], ChipRate [chip/s], Phase [rad]
             random +-1 chips, rectangular; main lobe +-ChipRate around Freq
  'pulsed'   radar pulses                  Freq [Hz], PulseWidth [s], PRF [Hz],
                                           ChirpBW [Hz] (linear chirp within the
                                           pulse, 0 = plain carrier), StartTime [s]
  'impulse'  broadband bursts              Rate [1/s] (Poisson), Duration [s];
             white noise bursts over the whole sampled band
  'Label'    free text for logs/plots
%}

arguments
    type {mustBeTextScalar}
    opts.Freq       (1,1) double = NaN
    opts.INRdB      (1,1) double = 0
    opts.Label             = ''
    opts.Phase      (1,1) double = 0
    opts.Drift      (1,1) double = 0
    opts.ChipRate   (1,1) double = NaN
    opts.PulseWidth (1,1) double = NaN
    opts.PRF        (1,1) double = NaN
    opts.ChirpBW    (1,1) double = 0
    opts.StartTime  (1,1) double = 0
    opts.Rate       (1,1) double = NaN
    opts.Duration   (1,1) double = NaN
end

type = lower(char(type));
need = struct('cw', {{'Freq'}}, 'bpsk', {{'Freq', 'ChipRate'}}, ...
    'pulsed', {{'Freq', 'PulseWidth', 'PRF'}}, ...
    'impulse', {{'Rate', 'Duration'}});
if ~isfield(need, type)
    error('rfiSource:type', 'Unknown RFI type "%s" (cw, bpsk, pulsed, impulse).', type);
end
for f = need.(type)
    v = opts.(f{1});
    if ~isfinite(v) || v <= 0
        error('rfiSource:param', 'RFI type "%s" needs a positive ''%s''.', type, f{1});
    end
end
if strcmp(type, 'pulsed') && opts.PulseWidth * opts.PRF >= 1
    error('rfiSource:duty', 'PulseWidth * PRF must be < 1.');
end

label = char(opts.Label);
if isempty(label), label = type; end
src = struct('type', type, 'label', label, 'freq', opts.Freq, 'INRdB', opts.INRdB, ...
    'phase', opts.Phase, 'drift', opts.Drift, 'chipRate', opts.ChipRate, ...
    'pulseWidth', opts.PulseWidth, 'prf', opts.PRF, 'chirpBW', opts.ChirpBW, ...
    'startTime', opts.StartTime, 'rate', opts.Rate, 'duration', opts.Duration);
end