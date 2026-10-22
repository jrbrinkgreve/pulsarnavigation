%{
EXPNAVPLAN - Which extra pulsars help the navigation fit? A planning study with navFit 'Predict' (block D4c-0)

Experiment, not a unit test; no simulation. For a fixed site at 52 N 5 E it builds TOA schedules of
B0329+54 and the weaker pulsars visible from there, and asks navFit for the errors those schedules
would give ('Predict': the covariance of the fit, Gaussian errors, no data needed):
  - the east and radial (equatorial) components, and the polar one with a calibrated start;
  - with the crystal clock as correlated noise: the same, and how well the clock is predicted one day
    after the last TOA (what limits a millisecond pulsar's fold).
The extra pulsars observe in their own simultaneous beams (Jasper, 22 Oct), so they cost B0329+54
nothing.

Schedules. Each pulsar is observed whenever it is above 20 deg elevation: passes centred on its
transits at the site (hour angle 0, from the Earth rotation angle as in delayModel), half-width
from the elevation mask. Every 'hours per TOA' of accumulated pass time gives one TOA, placed at the
middle of that stretch (adaptive groups also span the gap between passes). Numbers per pulsar:
  B0329+54   the simulated 30-day runs: 10.07 TOAs/day at 738 us (nominal), 1.97 at 744 us (worst),
             with 19.6 h visible a day: 1.95 h and 9.95 h per TOA;
  others     scenario.html's table (the radiometer estimate at SNR 7 with the catalogue W50 and S1400):
             B0950+08 12 / 67 h per TOA, 738 us; B1642-03 52 / 278 h, 326 us; B1929+10 106 / 571 h,
             480 us (nominal / worst). B1749-28 (dec -28) is never above 20 deg from 52 N.
Positions (J2000): B0329+54 53.2474 / +54.5787 deg; B0950+08 09h53m09.3s / +07d55m36s (ATNF);
B1929+10 and B1642-03 from their J2000 names (J1932+1059, J1645-0317), to an arcminute - more than
enough for a planning study.

Run from the PulsarSimMatlab folder: run('tests/expNavPlan.m'). ~10 s. (22 Oct 2026, block D4c-0.)
%}

scriptDir = fileparts(mfilename('fullpath'));
addpath(fullfile(scriptDir, '..', 'functions'));
t1 = 9555.5 * 86400;                                  % 2026-03-01 00:00 TT
lat = 52; lon = 5; elMin = 20;
site = makeTrajectory('fixed', 'Lat', lat, 'Lon', lon);
hCr = makeClock([0 1], 'Preset', "crystal").h;
P = struct('name', {'B0329+54', 'B0950+08', 'B1642-03', 'B1929+10'}, ...
    'raDeg', {53.2474, 15 * (9 + 53/60 + 9.3/3600), 15 * (16 + 45/60), 15 * (19 + 32/60)}, ...
    'decDeg', {54.5787, 7 + 55/60 + 36/3600, -(3 + 17/60), 10 + 59/60}, ...
    'hNom', {1.95, 12, 52, 106}, 'hWorst', {9.95, 67, 278, 571}, ...
    'sigNom', {738e-6, 738e-6, 326e-6, 480e-6}, 'sigWorst', {744e-6, 738e-6, 326e-6, 480e-6});
sets = {1, [1 2], [1 2 3 4], [2 3 4]};
setNames = {'B0329+54 alone', '+ B0950+08', '+ all three', 'the three without B0329+54'};

for nDays = [30 100]
    for fl = ["nominal", "worst"]
        fprintf('\nexpNavPlan: %d days, %s flux, 52 N; predicted 1-sigma errors [km], clock one day ahead [us]\n', nDays, fl);
        fprintf('  %-28s | TOAs per pulsar      | perfect clock: east radial polar(cal) | crystal clock: east radial polar(cal) clock(self, cal)\n', 'pulsars');
        for is = 1:numel(sets)
            idx = sets{is};
            [tt, ss, kk, nPer] = deal([], [], [], zeros(1, numel(idx)));
            for j = 1:numel(idx)
                p = P(idx(j));
                if fl == "nominal", h = p.hNom; s = p.sigNom; else, h = p.hWorst; s = p.sigWorst; end
                tj = schedule(t1, nDays, lat, lon, elMin, p.raDeg, p.decDeg, h * 3600);
                tt = [tt, tj]; ss = [ss, s * ones(size(tj))]; kk = [kk, j * ones(size(tj))]; %#ok<AGROW>
                nPer(j) = numel(tj);
            end
            if numel(tt) < 4
                fprintf('  %-28s | %-20s | too few TOAs\n', setNames{is}, mat2str(nPer)); continue
            end
            [tt, o] = sort(tt); ss = ss(o); kk = kk(o);
            toa = struct('t', tt, 'sigma', ss, 'k', kk);
            ps = P(idx);
            fS = navFit(toa, ps, site, 'Predict', true);
            fC = navFit(toa, ps, site, 'Predict', true, 'Start', "calibrated");
            gS = navFit(toa, ps, site, 'Predict', true, 'Clock', hCr, 'ClockStart', t1);
            gC = navFit(toa, ps, site, 'Predict', true, 'Start', "calibrated", 'Clock', hCr, 'ClockStart', t1);
            tq = tt(end) + 86400;
            fprintf('  %-28s | %-20s | %6.1f %6.1f %6.1f              | %6.1f %6.1f %6.1f   %6.0f %6.0f\n', setNames{is}, ...
                mat2str(nPer), fS.eastErr / 1e3, fS.radialErr / 1e3, fC.polarErr / 1e3, gS.eastErr / 1e3, ...
                gS.radialErr / 1e3, gC.polarErr / 1e3, 1e6 * gS.clockAtErr(tq), 1e6 * gC.clockAtErr(tq));
        end
    end
end


% =================================================================================
function tq = schedule(t1, nDays, lat, lon, elMin, ra, dec, dtTOA)
%SCHEDULE  TOA times of one pulsar: one per dtTOA seconds of accumulated time above elMin, at the middle
% of that stretch; passes centred on the transits at the site (hour angle 0).
c = (sind(elMin) - sind(lat) * sind(dec)) / (cosd(lat) * cosd(dec));
if c >= 1, tq = []; return, end
half = acos(max(c, -1)) / (2*pi) * 86164.0905;        % half a pass [s]
if c <= -1, half = 86164.0905 / 2 - 60; end
era = @(t) 2*pi * (0.7790572732640 + 1.00273781191135448 * (t - 69.184) / 86400);
H0 = mod(era(t1) + deg2rad(lon - ra), 2*pi);          % hour angle at t1
tTr = t1 + mod(-H0, 2*pi) / (2*pi) * 86164.0905 + (-1:nDays + 1) * 86164.0905;   % transits
P = [tTr - half; tTr + half].';
P(:, 1) = max(P(:, 1), t1); P(:, 2) = min(P(:, 2), t1 + nDays * 86400);
P = P(P(:, 2) > P(:, 1), :);
tq = []; acc = 0; tStart = P(1, 1);
for i = 1:size(P, 1)
    a = P(i, 1);
    if acc == 0, tStart = a; end                      % a new stretch starts with this pass
    while true
        need = dtTOA - acc;
        if a + need <= P(i, 2)                         % the stretch completes in this pass
            tEnd = a + need; tq(end+1) = (tStart + tEnd) / 2; %#ok<AGROW>
            acc = 0; a = tEnd; tStart = a;
        else                                           % carried into the next pass
            acc = acc + (P(i, 2) - a);
            break
        end
    end
end
end
