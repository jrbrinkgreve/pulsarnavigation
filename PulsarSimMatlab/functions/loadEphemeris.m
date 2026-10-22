function eph = loadEphemeris(file)
%LOADEPHEMERIS  JPL's planetary ephemeris from a binary SPK file (DE440): where the Earth, the Moon
%and the Sun are relative to the solar-system barycentre (block D5a; ephemeris).
%{
delayModel's analytic Earth orbit is good to about 1e4 km (~30 ms of Roemer delay). In a simulation
that error cancels, because the truth and the observer share the model; real TOAs need the Earth's
true position. JPL's DE440 gives it to well below a metre for the present era. This reads JPL's
binary file (de440s.bsp from NAIF: 1849-2150, 32.7 MB; data/ephem/) once and returns functions
that evaluate it:

  eph = loadEphemeris(fullfile(dataDir, "ephem", "de440s.bsp"));
  [r, v] = eph.earth(t);              % the geocentre relative to the SSB [3 x N, m and m/s]
  [r, v] = eph.sun(t);                % the Sun;  eph.moon(t), eph.emb(t) (Earth-Moon barycentre)
  [r, v] = eph.state(target, center, t)   % any NAIF ids in the file (0 SSB, 3 EMB, 10 Sun,
                                          % 301 Moon, 399 Earth, 5 Jupiter barycentre, ...)
t: seconds since J2000 in TDB (the ephemeris' own time scale), a row. Frame: ICRF, which is the
J2000 equator and equinox to 0.02" (delayModel's frame).

The file format (NAIF's "DAF" and "SPK Required Reading"): 1024-byte records of little-endian
doubles. The first record names the format and points to the summary records; each summary gives
a segment's time span, target, center, frame, data type and the addresses of its data. DE440 uses
type 2: for each body, consecutive time intervals of equal length, each with Chebyshev polynomial
coefficients of x, y and z (km) in the scaled time s = (t - mid) / radius in [-1, 1]:
  x(t) = sum_k c_k T_k(s),    dx/dt = sum_k c_k T_k'(s) / radius,
with T_0 = 1, T_1 = s, T_k = 2 s T_(k-1) - T_(k-2) (and its derivative by the same recurrence). A
segment's last four numbers are the first interval's start, the interval length, the record size
and the number of records. Bodies chain to the SSB through their centers: the Earth (399) is given
relative to the EMB (3), the EMB relative to the SSB (0).

TDB - TT at the geocentre (block D5c), from the same file: a clock on the Earth runs slow by its
speed and the Sun's and planets' potential, relative to the barycentre's time, so
  d(TDB - TT)/dt = (v_E^2 / 2 + U) / c^2 - L_C,   U = sum GM_i / |r_E - r_i| (Sun, Moon, planets),
with L_C = 1.48082686741e-8 the IAU's mean rate (Resolution B3, 2006), so that TDB keeps pace with
TT on average. Integrated with the trapezoid rule on a 6-hour grid over the whole file (1 s; 12 ns
from a 12-hour grid) and evaluated by a cubic spline. GM values: DE440's (Park et al. 2021). The
constant: the mean of the USNO series (Kaplan 2005, USNO Circular 179, eq. 2.6, good to ~10 us)
over 2000-2050, which the integral follows to 2.2 us rms there. (Fairhead & Bretagnon's series is
this integral expanded analytically; delayModel's analytic model keeps only its two largest terms,
up to 60 us off.)

Output eph: earth, sun, moon, emb (@(t) [r, v]), state (@(target, center, t) [r, v]), tdbMinusTT
(@(t) TDB - TT at the geocentre [s]), span [t1 t2] (s TDB since J2000, covered by all segments),
segments (struct array: target, center, frame, type, span), file.
Errors: 'loadEphemeris:file' (no file), 'loadEphemeris:format' (not a little-endian SPK file, or a
segment type other than 2), 'loadEphemeris:body' (an id not in the file), 'loadEphemeris:range'
(a time outside a segment).
%}

arguments
    file (1,1) string
end

if ~isfile(file)
    error('loadEphemeris:file', 'No ephemeris file %s (de440s.bsp from naif.jpl.nasa.gov, in data/ephem/).', file);
end
fid = fopen(file, 'r', 'ieee-le');
closer = onCleanup(@() fclose(fid));
rec = fread(fid, 1024, '*uint8').';
if numel(rec) < 1024 || ~strcmp(char(rec(1:7)), 'DAF/SPK') || ~strcmp(strtrim(char(rec(89:96))), 'LTL-IEEE')
    error('loadEphemeris:format', '%s is not a little-endian SPK file (DAF/SPK, LTL-IEEE).', file);
end
nd = double(typecast(rec(9:12), 'int32')); ni = double(typecast(rec(13:16), 'int32'));
ss = nd + ceil(ni / 2);                              % a summary's size in doubles

% the summaries: a chain of summary records
sm = zeros(0, ni + nd);                              % target center frame type begin end | t1 t2
r = double(typecast(rec(77:80), 'int32'));
while r > 0
    fseek(fid, (r - 1) * 1024, 'bof');
    d = fread(fid, 128, 'double');
    for j = 1:d(3)
        o = 3 + (j - 1) * ss;
        ic = typecast(d(o + nd + (1:ceil(ni / 2))), 'int32');   % the integers packed in the doubles
        sm(end+1, :) = [double(ic(1:ni)).', d(o + (1:nd)).']; %#ok<AGROW>
    end
    r = d(1);
end

% the type-2 segments' data
segs = struct('target', {}, 'center', {}, 'frame', {}, 'type', {}, 'span', {}, 'init', {}, ...
    'intlen', {}, 'n', {}, 'coef', {});
for i = 1:size(sm, 1)
    if sm(i, 4) ~= 2
        error('loadEphemeris:format', 'Segment %d (target %d) has SPK type %d; only type 2 is read.', i, sm(i, 1), sm(i, 4));
    end
    fseek(fid, (sm(i, 6) - 4) * 8, 'bof');          % the trailer: init, intlen, rsize, n
    tr = fread(fid, 4, 'double');
    fseek(fid, (sm(i, 5) - 1) * 8, 'bof');
    coef = reshape(fread(fid, tr(3) * tr(4), 'double'), tr(3), tr(4));
    segs(end+1) = struct('target', sm(i, 1), 'center', sm(i, 2), 'frame', sm(i, 3), 'type', sm(i, 4), ...
        'span', sm(i, ni + (1:2)), 'init', tr(1), 'intlen', tr(2), 'n', tr(4), 'coef', coef); %#ok<AGROW>
end

eph = struct();
eph.state = @(target, center, t) relState(segs, target, center, t);
eph.earth = @(t) relState(segs, 399, 0, t);
eph.sun = @(t) relState(segs, 10, 0, t);
eph.moon = @(t) relState(segs, 301, 0, t);
eph.emb = @(t) relState(segs, 3, 0, t);
eph.span = [max(arrayfun(@(s) s.span(1), segs)), min(arrayfun(@(s) s.span(2), segs))];
eph.tdbMinusTT = tdbTable(segs, eph.span);
eph.segments = rmfield(segs, {'init', 'intlen', 'n', 'coef'});
eph.file = file;
end


% =================================================================================
function f = tdbTable(segs, span)
% TDB - TT at the geocentre: the integral of (v_E^2/2 + U)/c^2 - L_C on a 6-hour grid, 0 at J2000 plus
% the USNO series' mean over 2000-2050, as a spline.
c = 299792458; LC = 1.48082686741e-8; h = 6 * 3600;
emrat = 81.3005682214972154;                         % DE440's Earth/Moon mass ratio
ids = [10 301 1 2 4 5 6 7 8 9];                      % Sun, Moon, Mercury ... Pluto (system barycentres)
gm = 1e9 * [132712440041.279419, 403503.235502 / (1 + emrat), 22031.868551, 324858.592, 42828.375816, ...
            126712764.1, 37940584.8418, 5794556.4, 6836527.10058, 975.5];   % [m^3/s^2]
tg = (ceil(span(1) / h):floor(span(2) / h)) * h;
[rE, vE] = toSSB(segs, 399, tg);
U = zeros(size(tg));
for j = 1:numel(ids)
    U = U + gm(j) ./ vecnorm(rE - toSSB(segs, ids(j), tg));
end
rate = (sum(vE.^2, 1) / 2 + U) / c^2 - LC;
F = cumsum([0, (rate(1:end-1) + rate(2:end)) / 2 * h]);
F = F - F(tg == 0);
w = tg >= 0 & tg <= 50 * 365.25 * 86400;
F = F + mean(usnoSeries(tg(w)) - F(w));
pp = spline(tg, F);
f = @(t) ppval(pp, t);
end

function dt = usnoSeries(t)
% TDB - TT, USNO Circular 179 (Kaplan 2005) eq. 2.6 [s]; T in Julian centuries since J2000.
T = t / (36525 * 86400);
dt = 0.001657 * sin(628.3076 * T + 6.2401) + 0.000022 * sin(575.3385 * T + 4.2970) + ...
     0.000014 * sin(1256.6152 * T + 6.1969) + 0.000005 * sin(606.9777 * T + 4.0212) + ...
     0.000005 * sin(52.9691 * T + 0.4444) + 0.000002 * sin(21.3299 * T + 5.5431) + ...
     0.000010 * T .* sin(628.3076 * T + 4.2490);
end

function [r, v] = relState(segs, target, center, t)
% The state of target relative to center: both chained to the SSB (id 0).
[rT, vT] = toSSB(segs, target, t);
[rC, vC] = toSSB(segs, center, t);
r = rT - rC; v = vT - vC;
end

function [r, v] = toSSB(segs, id, t)
% A body's state relative to the SSB, following its segments' centers.
t = t(:).';
r = zeros(3, numel(t)); v = r;
while id ~= 0
    j = find([segs.target] == id, 1, 'last');       % a later segment takes precedence (SPK rule)
    if isempty(j), error('loadEphemeris:body', 'Body %d is not in the ephemeris file.', id); end
    [ri, vi] = chebyshev(segs(j), t);
    r = r + ri; v = v + vi;
    id = segs(j).center;
end
end

function [r, v] = chebyshev(S, t)
% One type-2 segment at times t: position [m] and velocity [m/s].
if any(t < S.span(1) | t > S.span(2))
    error('loadEphemeris:range', 'Times outside the ephemeris segment of body %d (%.0f to %.0f s TDB since J2000).', ...
        S.target, S.span(1), S.span(2));
end
idx = min(floor((t - S.init) / S.intlen) + 1, S.n);  % the interval of each time (the last one includes its end)
nc = (size(S.coef, 1) - 2) / 3;
rad = S.coef(2, idx);
s = (t - S.coef(1, idx)) ./ rad;
T = zeros(nc, numel(t)); dT = T;
T(1, :) = 1;
if nc > 1, T(2, :) = s; dT(2, :) = 1; end
for k = 3:nc
    T(k, :) = 2 * s .* T(k-1, :) - T(k-2, :);
    dT(k, :) = 2 * T(k-1, :) + 2 * s .* dT(k-1, :) - dT(k-2, :);
end
r = zeros(3, numel(t)); v = r;
for a = 1:3
    C = S.coef(2 + (a - 1) * nc + (1:nc), idx);
    r(a, :) = sum(C .* T, 1);
    v(a, :) = sum(C .* dT, 1) ./ rad;
end
r = 1e3 * r; v = 1e3 * v;                            % km -> m
end
