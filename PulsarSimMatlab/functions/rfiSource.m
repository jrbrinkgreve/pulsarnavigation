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
                                           pulse, 0 = plain carrier), StartTime [s],
                                           RiseTime [s] (default 0.1e-6, see below)
  'impulse'  broadband bursts              Rate [1/s] (Poisson), Duration [s];
             white noise bursts over the whole sampled band
  'noise'    band-limited Gaussian noise   Freq [Hz] (centre), Bandwidth [Hz]:
             e.g. an LTE / OFDM downlink (many subcarriers with random data
             -> Gaussian by the central limit theorem); flat over the band,
             -6 dB at Freq +- Bandwidth/2, transition 5 % of Bandwidth,
             60 dB stopband (B5d, 8 Oct 2026)
  'Label'    free text for logs/plots
  any type:  rotating antenna              ScanPeriod [s], BeamTime [s],
                                           BeamWidth [s], SidelobeDB (see below)

Radar edges ('pulsed', 'RiseTime'): a real transmitter ramps its pulse up
and down. The envelope rises and falls as a raised cosine of duration
RiseTime (0 -> 100 %; the 10-90 % rise time is 0.59 x RiseTime), centred on
the nominal edges, so PulseWidth stays the width between the 50 % points
(the usual radar definition) and the energy is (PulseWidth - RiseTime/4)
x the power while on. Spectrum: with instant edges (RiseTime 0, the
behaviour before 8 Oct 2026) the power falls only as 1/f^2 away from the
carrier, sidelobes over tens of MHz; with raised-cosine edges it falls as
1/f^6 beyond ~1/(2*RiseTime) (5 MHz for 0.1 us). Default 0.1e-6: an
assumed realistic value for a solid-state L-band surveillance radar (to be
checked against ITU emission masks); RiseTime <= PulseWidth.

Rotating antenna (any type; 'ScanPeriod' finite): a surveillance radar's
antenna turns (5-6 rpm: ScanPeriod 10-12 s), so its main beam points at us
once per scan and only its sidelobes the rest of the time. The source's
power is multiplied by the antenna gain seen at time t,
  g(t) = max( exp(-4 ln2 (dt/BeamWidth)^2), 10^(SidelobeDB/10) ),
dt = t - BeamTime wrapped to (-ScanPeriod/2, ScanPeriod/2]: a Gaussian main
beam with -3 dB width BeamWidth in time (= azimuth beamwidth / 360 deg x
ScanPeriod; 1.4 deg at 10 s -> 39 ms), centred at BeamTime + k ScanPeriod,
on a flat sidelobe floor (default -30 dB). INRdB is the main-beam peak.
Default ScanPeriod Inf: no rotation (g = 1), the behaviour before 8 Oct 2026.
Beamwidth and sidelobe levels of real L-band radars are assumptions here.
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
    opts.RiseTime   (1,1) double = NaN
    opts.Rate       (1,1) double = NaN
    opts.Duration   (1,1) double = NaN
    opts.Bandwidth  (1,1) double = NaN
    opts.ScanPeriod (1,1) double = Inf
    opts.BeamTime   (1,1) double = 0
    opts.BeamWidth  (1,1) double = NaN
    opts.SidelobeDB (1,1) double = -30
end

type = lower(char(type));
need = struct('cw', {{'Freq'}}, 'bpsk', {{'Freq', 'ChipRate'}}, ...
    'pulsed', {{'Freq', 'PulseWidth', 'PRF'}}, ...
    'impulse', {{'Rate', 'Duration'}}, 'noise', {{'Freq', 'Bandwidth'}});
if ~isfield(need, type)
    error('rfiSource:type', 'Unknown RFI type "%s" (cw, bpsk, pulsed, impulse, noise).', type);
end
for f = need.(type)
    v = opts.(f{1});
    if ~isfinite(v) || v <= 0
        error('rfiSource:param', 'RFI type "%s" needs a positive ''%s''.', type, f{1});
    end
end
riseTime = 0;                                       % only 'pulsed' has edges
if strcmp(type, 'pulsed')
    riseTime = opts.RiseTime;
    if isnan(riseTime), riseTime = 0.1e-6; end      % realistic default (8 Oct 2026)
    if riseTime < 0 || riseTime > opts.PulseWidth
        error('rfiSource:rise', 'RiseTime must lie in [0, PulseWidth].');
    end
    if (opts.PulseWidth + riseTime) * opts.PRF >= 1
        error('rfiSource:duty', '(PulseWidth + RiseTime) * PRF must be < 1.');
    end
end

if isfinite(opts.ScanPeriod)
    if opts.ScanPeriod <= 0 || ~(opts.BeamWidth > 0 && opts.BeamWidth < opts.ScanPeriod)
        error('rfiSource:scan', ['A rotating antenna needs ScanPeriod > 0 and ' ...
            '0 < BeamWidth < ScanPeriod.']);
    end
    if opts.SidelobeDB > 0
        error('rfiSource:scan', 'SidelobeDB must be <= 0 (relative to the main beam).');
    end
elseif opts.ScanPeriod ~= Inf
    error('rfiSource:scan', 'ScanPeriod must be positive (Inf = no rotation).');
end

label = char(opts.Label);
if isempty(label), label = type; end
src = struct('type', type, 'label', label, 'freq', opts.Freq, 'INRdB', opts.INRdB, ...
    'phase', opts.Phase, 'drift', opts.Drift, 'chipRate', opts.ChipRate, ...
    'pulseWidth', opts.PulseWidth, 'prf', opts.PRF, 'chirpBW', opts.ChirpBW, ...
    'startTime', opts.StartTime, 'riseTime', riseTime, 'rate', opts.Rate, ...
    'duration', opts.Duration, 'scanPeriod', opts.ScanPeriod, 'beamTime', opts.BeamTime, ...
    'beamWidth', opts.BeamWidth, 'sidelobeDB', opts.SidelobeDB, 'bandwidth', opts.Bandwidth);
end