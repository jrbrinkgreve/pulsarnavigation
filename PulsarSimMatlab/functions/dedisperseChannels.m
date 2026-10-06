function info = dedisperseChannels(info_chan, outBase, DM, opts)
%DEDISPERSECHANNELS  Coherent dedispersion of every channel of channelizeIQ.
%{
Second stage of the channelized front end. Each channel file is an
ordinary IQ file (fLO = its centre, fs = channel rate), so the validated
applyInverseDispersion is applied per channel with:
  - band = the channel's useful width, centre +- ChanWidth/2: the chirp
    filter is zero outside it, which removes the channelizer's transition
    region (oversampling) and makes the channels tile the band exactly;
  - RefFreq = one common frequency for all channels (default the top of
    the whole band, as in the full-band path), with AllowRefOutsideBand:
    the chirp of channel j then also contains the delay
    tau(f) - tau(RefFreq) between this channel and the reference, so all
    channels come out aligned at RefFreq. No separate shift step.

Why per channel: the chirp kernel only has to span the sweep from the
channel to RefFreq at the channel rate (DM 100, bottom channel: 126 ms x
4.17 MHz ~ 5e5 samples) instead of the full-band sweep at 800 MHz (~1e8
samples): this solves the processing-side memory problem at high DM.
Later, RFI excision acts on the channel files before this stage.

  info = dedisperseChannels(info_chan, outBase, DM)
  info = dedisperseChannels(info_chan, outBase, DM, 'RefFreq', 1.6e9, 'EdgeFrac', 0.02)

Inputs:
  info_chan  info struct of channelizeIQ (chanFiles, chanFreqs, chanWidth,
             fs, N, fullySupported, t0, fHigh).
  outBase    path without extension; channel j -> <outBase>_ch001.dat ...
             (cf32), info -> <outBase>_info.mat (loadInfo(outBase) works).
  DM         [pc cm^-3] from the ephemeris.

Name-value options:
  'RefFreq'      [Hz] common reference (default info_chan.fHigh).
  'EdgeFrac'     sin^2 band-edge taper per channel, fraction of ChanWidth
                 (default 0.02, inside the channel). Neighbouring channels
                 do not overlap, so their noise is independent; the taper
                 costs the same fraction of the band as the full-band path
                 with EdgeFrac 0.02 of 400 MHz (int W^2 390.0 MHz in both).
  'GuardTime'    [s] passed on (default [] = 20/edge width).
  'MaxMemoryGB'  per channel, passed on (default 4).
  'SaveInfo'     save info to <outBase>_info.mat (default true).
  'Verbose'      print a summary (default true).

Output info: chanFiles, chanFreqs, chanWidth, fs, N, t0, DM, refFreq,
  edgeWidth, per-channel nPast / nFuture / Nfft / leakage,
  fullySupported (1-based sample range valid in EVERY channel: the
  channelizer's supported range shrunk by each channel's nPast/nFuture),
  BnoiseChan (noise bandwidth of one channel's taper, noiseBandwidth) and
  BnoiseTotal = nChan*BnoiseChan (all channels summed with equal weight;
  independent channels), bulkDelayRef.
%}

arguments
    info_chan         struct
    outBase                 {mustBeTextScalar}
    DM                (1,1) double {mustBeNonnegative, mustBeFinite}
    opts.RefFreq            double = []
    opts.EdgeFrac     (1,1) double {mustBePositive} = 0.02
    opts.GuardTime          double = []
    opts.MaxMemoryGB  (1,1) double {mustBePositive} = 4
    opts.SaveInfo     (1,1) logical = true
    opts.Verbose      (1,1) logical = true
end

tStart  = tic;
outBase = char(outBase);
fRef    = opts.RefFreq;
if isempty(fRef), fRef = info_chan.fHigh; end
nChan = info_chan.nChan;
dF    = info_chan.chanWidth;

outDir = fileparts(outBase);
if ~isempty(outDir) && ~isfolder(outDir)
    mkdir(outDir);
end

chanFiles = strings(nChan, 1);
nPast = zeros(nChan, 1); nFuture = zeros(nChan, 1);
Nfft  = zeros(nChan, 1); leakage = zeros(nChan, 1);
for j = 1:nChan
    chanFiles(j) = sprintf('%s_ch%03d.dat', outBase, j);
    fc = info_chan.chanFreqs(j);
    d = applyInverseDispersion(info_chan.chanFiles(j), chanFiles(j), info_chan.fs, ...
        fc, DM, fc - dF/2, fc + dF/2, 'RefFreq', fRef, 'AllowRefOutsideBand', true, ...
        'EdgeFrac', opts.EdgeFrac, 'GuardTime', opts.GuardTime, ...
        'MaxMemoryGB', opts.MaxMemoryGB, 'SaveInfo', false, 'Verbose', false);
    nPast(j) = d.nPast; nFuture(j) = d.nFuture;
    Nfft(j)  = d.Nfft;  leakage(j) = d.leakage;
    if j == 1
        edgeW = d.edgeWidth;
        bulk  = d.bulkDelayRef;
        Nout  = d.N;
    end
    if opts.Verbose && (j == nChan || mod(j, 32) == 0)
        fprintf('  dedisperseChannels: %d/%d channels\n', j, nChan);
    end
end

% Valid in every channel: channelizer support shrunk by each filter's reach
s = info_chan.fullySupported;
fullySupported = [max(s(1) + nPast), min(s(2) - nFuture)];
BnoiseChan = noiseBandwidth(-dF/2, dF/2, edgeW);

info = struct();
info.file            = string(outBase);
info.chanFiles       = chanFiles;
info.inFile          = info_chan.file;
info.format          = 'cf32 interleaved (I0,Q0,I1,Q1,...), one file per channel';
info.precision       = 'single';
info.byteOrder       = 'ieee-le';
info.isComplex       = true;
info.nChan           = nChan;
info.N               = Nout;
info.fs              = info_chan.fs;
info.actualFsOut     = info_chan.fs;
info.chanFreqs       = info_chan.chanFreqs;
info.chanWidth       = dF;
info.fLow            = info_chan.fLow;
info.fHigh           = info_chan.fHigh;
info.t0              = info_chan.t0;
info.DM              = DM;
info.refFreq         = fRef;
info.timingConvention = 'same grid as the channel files; every channel aligned at refFreq';
info.bulkDelayRef    = bulk;                        % tau(refFreq), not removed
info.edgeWidth       = edgeW;
info.nPast           = nPast;
info.nFuture         = nFuture;
info.Nfft            = Nfft;
info.leakage         = leakage;
info.fullySupported  = fullySupported;              % 1-based, all channels
info.BnoiseChan      = BnoiseChan;
info.BnoiseTotal     = nChan * BnoiseChan;
info.elapsed         = toc(tStart);

if opts.SaveInfo
    infoFile = [outBase '_info.mat'];
    save(infoFile, 'info');
    info.infoFile = infoFile;
end
if opts.Verbose
    fprintf(['dedisperseChannels: %d channels, DM %.4g, ref %.4f GHz, edge %.4g kHz, ' ...
             'Nfft 2^%d..2^%d, max leakage %.1e, supported %d..%d of %d, %.1f s\n'], ...
        nChan, DM, fRef/1e9, edgeW/1e3, log2(min(Nfft)), log2(max(Nfft)), max(leakage), ...
        fullySupported(1), fullySupported(2), Nout, info.elapsed);
end
end
