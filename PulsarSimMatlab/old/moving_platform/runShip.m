%{
RUNSHIP - A ship navigating by pulsars: the closed loop of block D4d (fast simulator + navFilter)

A ship sails a route that it knows only by dead reckoning (DR). Its array follows B0329+54 and
B0950+08 in simultaneous digital beams (Jasper, 22 Oct 2026), its SDR clock drifts, and every day
the observer runs navFilter on all its TOAs so far and folds the next day with the corrected track
and clock. The script reports how far the observer's track is from the truth, day by day, against
pure dead reckoning.

The truth (simulator side only):
  - the DR track: a rhumb line from 'route.start' with 'route.heading' and 'route.speed';
  - the true track: the DR track plus a drift in the ship's own east/north frame - an unknown mean
    current (route.current m/s per axis, drawn) plus a wander around it that reverts to it
    (route.currentWander m/s per axis, correlation time route.currentTime: a first-order
    Gauss-Markov process, so it stays bounded, unlike a random walk), integrated on an hourly grid
    (at sea level to first order: the drop d^2/2R below the tangent plane is left out, and the loop
    keeps d small). The filter models the wander as a random walk with the same change over a day:
    2 currentWander^2 / currentTime = 1.16e-7 m^2/s^3, its default;
  - the clock: the SDR's crystal noise (makeClock "crystal") and a rate drawn within +-2 ppm.
The observer:
  - starts from port: position known to 1 km, clock calibrated to 1 ms and 1e-9 (1 sigma, as
    runScenario; the errors drawn, portClock);
  - each pulsar is folded whenever it is above 20 deg elevation at the ship's DR position: passes
    centred on its transits there, recomputed every day as the ship moves;
  - folds with its current track (DR + the filter's predicted correction) and clock estimate:
    simulateFold's phase offset per sub-int = trueResiduals(true track and clock, its track and
    clock estimate); TOAs from adaptiveSubints (SNR ~7) and foldToTOA, as runScenario;
  - each sub-int keeps what the observer used for it: delta_DR - delta_used (delayModel) and its clock
    estimate - the port one, at the sub-int's own time and with the predictor it was folded with; a
    TOA merges sub-ints, possibly from two passes a day apart, so it is a time-average: its
    residual is referred back to the DR track and the port clock with the average of its sub-ints'
    values, and it is tagged at the time nearest its mean time whose sidereal phase is the circular
    mean of its sub-ints' (for two passes at the same hour angle: that hour angle; the midpoint
    between them, which mergeSubints uses, lies near anti-transit, where the daily position signal
    has the opposite sign - found 22 Oct: B0950+08's TOAs then carried none of it);
  - after each day: refers all TOAs so far back to the DR track and the port clock (as above),
    runs navFilter ('Clock' = its oscillator's h; 'Start' "self", or "calibrated" with the port's
    clock errors as 'ClockPrior'), and predicts the next day: position and velocity correction from
    the last filtered state (constant velocity), the clock from its offset and rate, plus B0329+54's
    constant with "self" (the clock offset at the start sits in the constants there; B0950+08's
    differs from it only by noise). "calibrated" knows the constants (0), so both pulsars share one
    clock, and their different declinations separate the polar error from the clock's wander.
Pulsars: B0329+54 as runScenario (five components, NBin 2048, 120 s sub-ints); B0950+08 a single
Gaussian with the catalogue values of scenario.html's table (P 253.1 ms, W50 8.6 ms, S1400 100 mJy;
NBin 512, 300 s sub-ints), no scintillation; both clean (receiver noise only), 'NPol' 2.

Output data/scenario/ship_<nDays>d_<flux>.mat (R: settings, TOAs with their truth, the daily
observer errors, the final smoothed fit) and a printed daily table: TOAs per pulsar, the DR error
|true - DR| and the observer's error |true - its track| at the end of each day [km], and the
clock it had not estimated [us]. Then the final smoother over all TOAs against the truth, and its
twin: the same smoother on ideal TOAs (the truth averaged over each TOA's sub-ints plus Gaussian
noise of their stated sigma; evaluation only), so that the TOA chain's errors and the filter's own
show apart. The filter takes each TOA's sub-int times (navFilter toa.tSub): its rows are averaged
over them, as the TOA is.
Other runs: set shipOpts before running the script, with any of the fields nDays, flux, seed and
start ("self" | "calibrated"), e.g. shipOpts = struct('seed', 3, 'start', "calibrated"); runShip.
The seed draws the current, the clock rate, the port's clock errors, the clock's noise and every
fold's noise; "calibrated" adds _cal to the file name, seeds other than 1 add _seed<n>.
(22 Oct 2026, block D4d-2; design agreed with Jasper the same day.)
%}

nDays = 14; flux = "nominal"; startMode = "self";
portClock = [1e-3 1e-9];                             % the port calibration's clock errors, 1 sigma: offset [s], rate
route = struct('start', [50 -5], 'heading', 250, 'speed', 5 * 1852 / 3600, ...   % 5 kn
    'current', 0.3, 'currentWander', 0.1, 'currentTime', 2 * 86400, 'seed', 1);
if exist('shipOpts', 'var')                          % a caller's settings (a dedicated name: runScenario also has nDays)
    if isfield(shipOpts, 'nDays'), nDays = shipOpts.nDays; end
    if isfield(shipOpts, 'flux'), flux = string(shipOpts.flux); end
    if isfield(shipOpts, 'seed'), route.seed = shipOpts.seed; end
    if isfield(shipOpts, 'start'), startMode = string(shipOpts.start); end
end
epoch = 9555.5 * 86400;                              % 2026-03-01 00:00 TT

scriptDir = fileparts(mfilename('fullpath'));
addpath(fullfile(scriptDir, 'functions'));
pipelineParams;
outDir = fullfile(dataDir, "scenario");
if ~isfolder(outDir), mkdir(outDir); end
P0 = load(fullfile(dataDir, "rfi", "profile", "rfiProfile_noiseOnly_1s.mat")).prof;   % channel settings
nCh = P0.nChan; BnC = P0.BnoiseChan; dtC = P0.binDt; ncC = P0.modelCoeffs;
gc = chanGain(fLow, fHigh, chanWidth, nCh);
nPol = 2; elMin = 20; tSid = 86164.0905; c = 299792458;
sefd = struct('nominal', 3.2e5, 'worst', 7.4e5);
cB = [0.159 -11.427 1.479; 0.439 -1.123 3.281; 0.585 0 1.055; 0.152 7.719 2.751; 0.170 9.482 0.847];
cB = [cB(:, 1), cB(:, 2) / 360, cB(:, 3) * 2*sqrt(2*log(2)) / 360];
psr = struct('name', {'B0329+54', 'B0950+08'}, 'f0', {1/0.714519699726, 1/0.2530651649482}, ...
    'raDeg', {53.2474, 15 * (9 + 53/60 + 9.3/3600)}, 'decDeg', {54.5787, 7 + 55/60 + 36/3600}, ...
    'comps', {cB, []}, 'fwhm', {0.0092, 8.6 / 253.065}, 'S1400', {0.203, 0.100}, 'nBin', {2048, 512}, ...
    'subintTime', {120, 300}, 'chunkHours', {2, 4});
nP = numel(psr);
for k = 1:nP
    prof = gaussianTemplate(2^16, profileOf(psr(k)));
    psr(k).rho = 10*log10(psr(k).S1400 / (sefd.(flux) * sum(prof) / 2^16));
    psr(k).nPer = max(1, round(psr(k).subintTime * psr(k).f0));
    psr(k).tmpl = gaussianTemplate(psr(k).nBin, profileOf(psr(k)));
end
pTime = @(k) struct('f0', psr(k).f0, 'F1', 0, 'TRef', epoch, 'raDeg', psr(k).raDeg, 'decDeg', psr(k).decDeg);
pDir = struct('raDeg', {psr.raDeg}, 'decDeg', {psr.decDeg});

% ---- the truth: DR track, true track, clock ---------------------------------------------------------
rs = RandStream('mt19937ar', 'Seed', route.seed);
tg = epoch + (0:3600:(nDays + 1) * 86400);
[latDR, lonDR] = rhumb(route.start, route.heading, route.speed * (tg - epoch));
trajDR = makeTrajectory('track', 'Time', tg, 'Lat', latDR, 'Lon', lonDR, 'Interp', 'pchip');
vMean = route.current * randn(rs, 2, 1);
u = route.currentWander * randn(rs, 2, numel(tg)); aGM = exp(-3600 / route.currentTime);
for i = 2:numel(tg), u(:, i) = aGM * u(:, i-1) + sqrt(1 - aGM^2) * u(:, i); end   % Gauss-Markov, stationary from the start
v = vMean + u;
qV = 2 * route.currentWander^2 / route.currentTime;  % the filter's random walk with the same daily change [m^2/s^3]
d = [zeros(2, 1), cumsum(3600 * (v(:, 1:end-1) + v(:, 2:end)) / 2, 2)];   % east, north [m]
rDR = trajDR.pos(tg); [Eg, Ng] = enFrame(rDR);
rT = rDR + Eg .* d(1, :) + Ng .* d(2, :);
trajT = makeTrajectory('custom', 'Frame', 'ECEF', 'PosFcn', @(t) interp1(tg, rT.', t, 'pchip').');
offT = @(t) [interp1(tg, d(1, :), t, 'pchip'); interp1(tg, d(2, :), t, 'pchip')];   % true - DR, local E N
clkRate = 2e-6 * (2 * rand(rs) - 1);
clk = makeClock([epoch, epoch + (nDays + 1) * 86400], 'Preset', "crystal", 'Rate', clkRate, 'Seed', 10 + route.seed);
portErr = portClock .* randn(rs, 1, 2);
ehat0 = @(cc) portErr(1) + (clkRate + portErr(2)) * (cc - epoch);   % the port calibration
fArgs = {'Clock', clk.h, 'Start', startMode, 'ClockPrior', portClock, 'PositionPrior', 1e3, ...   % the observer's filter
    'VelocityPrior', hypot(route.current, route.currentWander), 'VelocityNoise', qV};
fprintf('runShip: %d days (%s flux, seed %d, start %s) from %s at %.1f kn, heading %d deg; mean current %s m/s; clock rate %.2e\n', nDays, ...
    flux, route.seed, startMode, mat2str(route.start), route.speed * 3600 / 1852, route.heading, mat2str(round(vMean.', 3)), clkRate);

% ---- the observer's day-by-day loop ------------------------------------------------------------------
corr = struct('tL', epoch, 'd', [0; 0], 'v', [0; 0]);  % its track correction (none on day 1)
ehat = ehat0; trajA = trajDR;
G = struct('k', {}, 'day', {}, 'good', {}, 'tAbs', {}, 'resid', {}, 'sigma', {}, 'snr', {}, 'rRef', {}, ...
    'trueEN', {}, 'trueClock', {}, 'tSub', {});
carry = cell(1, nP); daily = zeros(nDays, 6);
for day = 1:nDays
    t0d = epoch + (day - 1) * 86400; t1d = t0d + 86400;
    for k = 1:nP
        ps = psr(k);
        win = passesOfDay(t0d, t1d, ps, trajDR, elMin, tSid);       % absolute [start end] rows
        if isempty(win), continue, end
        ca = {'Clock', clk, 'ClockEstimate', ehat};
        truthOf = @(s) trueResiduals(epoch + s, pTime(k), trajT, trajA, ca{:});
        simArgs = {'NChan', nCh, 'ChannelGain', gc, 'NPol', nPol, 'NBin', ps.nBin, 'SubintPeriods', ps.nPer, ...
            'Bnoise', BnC, 'BinDt', dtC, 'NoiseCoeffs', ncC, 'Verbose', false, 'PhaseOffset', @(s) offsetAt(truthOf, s)};
        if isempty(ps.comps), simArgs = [simArgs, {'FWHM', ps.fwhm}]; else, simArgs = [simArgs, {'Components', ps.comps}]; end %#ok<AGROW>
        TRef = 0.5 / ps.f0; tSub = ps.nPer / ps.f0;
        for iw = 1:size(win, 1)
            s1 = win(iw, 1) - epoch; s2 = win(iw, 2) - epoch;
            k1 = ceil(ps.f0 * (s1 - TRef) + 0.5); k2 = floor(ps.f0 * (s2 - TRef) - 0.5);
            nSub = floor((k2 - k1 + 1) / ps.nPer);
            mChunk = max(1, floor(ps.chunkHours * 3600 / tSub));
            for c0 = 0:mChunk:nSub - 1
                kA = k1 + c0 * ps.nPer; kB = k1 + min(c0 + mChunk, nSub) * ps.nPer - 1;
                w = TRef + ([kA - 0.5 - 1e-6, kB + 0.5 + 1e-6]) / ps.f0;
                [iC, fC] = simulateFold(ps.f0, ps.rho, [], 'Passes', w, simArgs{:}, 'Seed', 1e7 * (route.seed - 1) + 1e6 * k + 1000 * day + 10 * iw + c0);
                tN = epoch + fC.subint.tRef.';                % this chunk's sub-ints: the observer's reference values
                qN = delayModel(tN, trajDR, pTime(k)).total - delayModel(tN, trajA, pTime(k)).total + ehat(tN) - ehat0(tN);
                rN = truthOf(tN - epoch); cN = rN.clock - rN.clockEstimate;
                if ~isempty(carry{k})
                    [fC, iC] = catFolds(carry{k}.fold, carry{k}.info, fC, iC);
                    tN = [carry{k}.t, tN]; qN = [carry{k}.q, qN]; cN = [carry{k}.c, cN];
                end
                [fM, iM, sel] = adaptiveSubints(fC, iC, ps.tmpl, BnC * nPol, 'TargetSNR', 7, 'Weighting', 'optimal', ...
                    'AcrossGaps', true, 'Verbose', false);
                nG = numel(sel.nMerged); nClosed = nG - ~sel.complete(end);
                if ~sel.complete(end)
                    carry{k} = struct(); [carry{k}.fold, carry{k}.info] = subFold(fC, iC, find(sel.group == nG));
                    isOpen = sel.group(:).' == nG;
                    carry{k}.t = tN(isOpen); carry{k}.q = qN(isOpen); carry{k}.c = cN(isOpen);
                else
                    carry{k} = [];
                end
                if nClosed == 0, continue, end
                [fT, iT] = subFold(fM, iM, 1:nClosed);
                toa = foldToTOA(fT, iT, profileOf(ps), BnC * nPol, 'Weighting', 'optimal', 'Verbose', false);
                [tA, qT, cT] = deal(zeros(1, nClosed));
                for g = 1:nClosed                              % each TOA: its sub-ints' average, and its tag
                    m = sel.group(:).' == g;
                    tA(g) = sideTag(tN(m), tSid); qT(g) = mean(qN(m)); cT(g) = mean(cN(m));
                end
                % the observer refers its residuals back to the DR track and the port clock
                resid = toa.toa.' - fT.subint.tRef.';
                rRef = resid + qT;
                en = offT(tA);
                for g = 1:nClosed
                    G(end+1) = struct('k', k, 'day', day, 'good', toa.good(g), 'tAbs', tA(g), 'resid', resid(g), ...
                        'sigma', toa.toaErr(g), 'snr', toa.snr(g), 'rRef', rRef(g), 'trueEN', en(:, g), ...
                        'trueClock', cT(g), 'tSub', tN(sel.group(:).' == g)); %#ok<SAGROW>
                end
            end
        end
    end
    % end of the day: the filter on all good TOAs so far, and tomorrow's track and clock
    g = [G.good];
    tEnd = t1d; enTrue = offT(tEnd);
    obsErr = norm(enTrue - (corr.d + (tEnd - corr.tL) * corr.v));
    if nnz(g) >= 10
        tq = struct('t', [G(g).tAbs], 'r', [G(g).rRef], 'sigma', [G(g).sigma], 'k', [G(g).k]); tq.tSub = {G(g).tSub};
        ws = warning('off', 'navFilter:converge');
        f = navFilter(tq, pDir, trajDR, fArgs{:});
        warning(ws);
        xL = f.xF(:, end); ix = f.index;
        corr = struct('tL', f.t(end), 'd', xL(ix.P), 'v', xL(ix.V));
        b0 = xL(ix.B) / c; r0 = xL(ix.D) / c; tL = f.t(end);
        J1 = 0; if startMode == "self", J1 = xL(ix.J(1)) / c; end   % "calibrated": no constants
        ehat = @(cc) ehat0(cc) + J1 + b0 + r0 * (cc - tL);
        trajA = correctedTrack(trajDR, corr);
    end
    gd = [G.day] == day & [G.good];
    daily(day, :) = [day, nnz(gd & [G.k] == 1), nnz(gd & [G.k] == 2), norm(enTrue) / 1e3, obsErr / 1e3, ...
        1e6 * rms([G(gd).trueClock])];
    fprintf('runShip: day %2d: TOAs %2d + %d; DR error %6.1f km, the observer''s %6.1f km; clock not estimated %6.0f us\n', ...
        daily(day, :));
end

% ---- the final smoother over all TOAs against the truth ---------------------------------------------
g = [G.good];
tq = struct('t', [G(g).tAbs], 'r', [G(g).rRef], 'sigma', [G(g).sigma], 'k', [G(g).k]); tq.tSub = {G(g).tSub};
fS = navFilter(tq, pDir, trajDR, fArgs{:});
enT = offT(fS.t);
errS = fS.enu(1:2, :) - enT;
pull = errS ./ fS.enuErr(1:2, :);
fprintf(['runShip: final smoother, %d TOAs: track error rms east %.1f km, north %.1f km (DR rms %.1f km); stated error median ' ...
         '%.1f / %.1f km; pulls rms %.2f / %.2f (one realization: its errors are correlated along the track, so this is a ' ...
         'size, not a distribution); at the end %.1f km (DR %.1f km)\n'], numel(fS.t), 1e-3 * rms(errS(1, :)), ...
    1e-3 * rms(errS(2, :)), 1e-3 * rms(vecnorm(enT)), 1e-3 * median(fS.enuErr(1, :)), 1e-3 * median(fS.enuErr(2, :)), ...
    rms(pull(1, :)), rms(pull(2, :)), 1e-3 * norm(errS(:, end)), 1e-3 * norm(enT(:, end)));
% the twin: the same smoother on ideal TOAs, averaged over the same sub-ints (evaluation only)
kq = [G(g).k]; tI = {G(g).tSub}; rI = zeros(size(kq));
for i = 1:numel(rI)                                  % the truth averaged over each TOA's sub-ints, as the data
    rI(i) = mean(trueResiduals(tI{i}, pTime(kq(i)), trajT, trajDR, 'Clock', clk, 'ClockEstimate', ehat0).residual);
end
rsTwin = RandStream('mt19937ar', 'Seed', 1000 + route.seed);
fI = navFilter(setfield(tq, 'r', rI + tq.sigma .* randn(rsTwin, size(rI))), pDir, trajDR, fArgs{:}); %#ok<SFLD>
errI = fI.enu(1:2, :) - enT;
fprintf('runShip: the twin (ideal TOAs): track error rms east %.1f km, north %.1f km; pulls rms %.2f / %.2f\n', ...
    1e-3 * rms(errI(1, :)), 1e-3 * rms(errI(2, :)), rms(errI(1, :) ./ fI.enuErr(1, :)), rms(errI(2, :) ./ fI.enuErr(2, :)));
R = struct('route', route, 'nDays', nDays, 'flux', flux, 'start', startMode, 'portErr', portErr, ...
    'psr', rmfield(psr, 'tmpl'), 'clockRate', clkRate, ...
    'G', G, 'daily', daily, 'final', struct('t', fS.t, 'enu', fS.enu, 'enuErr', fS.enuErr, 'truth', enT), ...
    'twin', struct('enu', fI.enu, 'enuErr', fI.enuErr), ...
    'made', char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm')));
tag = ""; if startMode == "calibrated", tag = "_cal"; end
if route.seed ~= 1, tag = tag + sprintf("_seed%d", route.seed); end
save(fullfile(outDir, sprintf('ship_%dd_%s%s.mat', nDays, flux, tag)), 'R');


% =================================================================================
% Helpers. catFolds, subFold, profileOf and chanGain are copies of runScenario's (a script cannot share
% local functions); the rest is the ship's.
function [lat, lon] = rhumb(start, heading, dist)
%RHUMB  Positions along a rhumb line from start [lat lon] deg with heading [deg] after distances dist [m].
R = 6371e3; lat0 = deg2rad(start(1)); hd = deg2rad(heading);
lat = lat0 + dist * cos(hd) / R;
dq = log(tan(pi/4 + lat / 2) / tan(pi/4 + lat0 / 2));
qf = cos(lat0) * ones(size(dist)); m = abs(dq) > 1e-12; qf(m) = (lat(m) - lat0) ./ dq(m);
lon = deg2rad(start(2)) + dist * sin(hd) ./ (R * qf);
lat = rad2deg(lat); lon = rad2deg(lon);
end

function [E, N] = enFrame(r)
%ENFRAME  Local east and north unit vectors at Earth-fixed r [3 x M] (WGS84 normal, as navFilter).
a = 6378137; f = 1 / 298.257223563; e2 = f * (2 - f); b = a * (1 - f); ep2 = (a^2 - b^2) / b^2;
p = hypot(r(1, :), r(2, :)); lon = atan2(r(2, :), r(1, :));
th = atan2(r(3, :) * a, p * b);
lat = atan2(r(3, :) + ep2 * b * sin(th).^3, p - e2 * a * cos(th).^3);
E = [-sin(lon); cos(lon); zeros(size(lon))];
N = [-sin(lat) .* cos(lon); -sin(lat) .* sin(lon); cos(lat)];
end

function traj = correctedTrack(trajDR, corr)
%CORRECTEDTRACK  The observer's track: DR plus the filter's predicted correction (constant velocity from
% its last filtered state), in the local east/north frame of the DR track.
pos = @(t) shifted(trajDR.pos(t), t, corr);
traj = makeTrajectory('custom', 'Frame', 'ECEF', 'PosFcn', pos);
end

function r = shifted(r0, t, corr)
[E, N] = enFrame(r0);
dd = corr.d + (t - corr.tL) .* corr.v;
r = r0 + E .* dd(1, :) + N .* dd(2, :);
end

function W = passesOfDay(t0, t1, ps, trajDR, elMin, tSid)
%PASSESOFDAY  [start end] of the pulsar's time above elMin at the ship's DR position within [t0 t1]:
% passes centred on its transits (hour angle 0 from the Earth rotation angle, as delayModel).
r = trajDR.pos((t0 + t1) / 2);
lat = atan2d(r(3), hypot(r(1), r(2))); lon = atan2d(r(2), r(1));
cH = (sind(elMin) - sind(lat) * sind(ps.decDeg)) / (cosd(lat) * cosd(ps.decDeg));
W = zeros(0, 2);
if cH >= 1, return, end
half = acos(max(cH, -1)) / (2*pi) * tSid;
if cH <= -1, W = [t0 t1]; return, end
era = @(t) 2*pi * (0.7790572732640 + 1.00273781191135448 * (t - 69.184) / 86400);
H0 = mod(era(t0) + deg2rad(lon - ps.raDeg), 2*pi);
tTr = t0 + mod(-H0, 2*pi) / (2*pi) * tSid + (-1:1) * tSid;
W = [max(tTr.' - half, t0), min(tTr.' + half, t1)];
W = W(W(:, 2) > W(:, 1), :);
end

function v = offsetAt(truthOf, s)
r = truthOf(s); v = r.phaseOffset;
end

function tg = sideTag(t, tSid)
%SIDETAG  A TOA's time tag: the time nearest the mean of its sub-ints' times t whose sidereal phase is the
% circular mean of theirs (the hour angle its data represent; equal weights, as the sub-ints are equal).
w = 2*pi / tSid; tm = mean(t);
tg = tm + angle(mean(exp(1i * w * (t - tm)))) / w;
tg = min(max(tg, min(t)), max(t));                   % within its data (a partial pass can put that phase just outside)
end

function [f, info] = catFolds(fA, iA, fB, iB)
f = fB;
for fn = ["sum", "prof", "weight", "weight2", "weightX"], f.(fn) = cat(2, fA.(fn), fB.(fn)); end
for fn = string(fieldnames(fB.subint)).', f.subint.(fn) = [fA.subint.(fn); fB.subint.(fn)]; end
nChan = size(f.sum, 3);
f.profTotal = reshape(sum(f.sum, 2), [], nChan) ./ reshape(sum(f.weight, 2), size(f.weight, 1), []);
info = iB; info.nSub = iA.nSub + iB.nSub;
end

function [f, info] = subFold(f0, i0, idx)
f = f0;
for fn = ["sum", "prof", "weight", "weight2", "weightX"], f.(fn) = f0.(fn)(:, idx, :, :); end
for fn = string(fieldnames(f0.subint)).', f.subint.(fn) = f0.subint.(fn)(idx); end
nChan = size(f.sum, 3);
f.profTotal = reshape(sum(f.sum, 2), [], nChan) ./ reshape(sum(f.weight, 2), size(f.weight, 1), []);
info = i0; info.nSub = numel(idx);
end

function prof = profileOf(ps)
prof = ps.comps; if isempty(prof), prof = ps.fwhm; end
end

function g = chanGain(fLow, fHigh, chanWidth, nCh)
eF = 0.02 * (fHigh - fLow); eC = 0.02 * chanWidth;
g = zeros(1, nCh);
for c = 1:nCh
    lo = fLow + (c - 1) * chanWidth; hi = lo + chanWidth;
    f = linspace(lo, hi, 20001);
    Wi = ones(size(f)); a = f < lo + eC; Wi(a) = sin(pi/2 * (f(a) - lo) / eC).^2;
    b = f > hi - eC; Wi(b) = sin(pi/2 * (hi - f(b)) / eC).^2;
    Wf = ones(size(f)); a = f < fLow + eF; Wf(a) = sin(pi/2 * (f(a) - fLow) / eF).^2;
    b = f > fHigh - eF; Wf(b) = sin(pi/2 * (fHigh - f(b)) / eF).^2;
    g(c) = trapz(f, Wf.^2 .* Wi.^2) / trapz(f, Wi.^2);
end
end
