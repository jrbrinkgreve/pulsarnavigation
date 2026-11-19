# 19 Nov 2026 — idea to explore later: the Earth's motion as a moving receiver

Jasper: "can you note down the idea of using the movement of the earth + measurements for
simultaneous estimation of more than just a range to a pulsar, but also infer information about 3d
position? ... its like having a moving receiver. just note this idea, i want to explore it later".
A note only. Nothing was computed or run for it, and no code or docs changed.

## The idea

One TOA of one pulsar gives one number: the site's position along that pulsar's line of sight (a
range). A snapshot fix therefore needs several pulsars in different directions. But the receiver is
not standing still in space: the Earth carries it around. Measurements spread over time see the
same pulsar from a site vector that has turned, so they sample different components of the
position. Time takes the place of extra pulsars, as a moving receiver would: a synthetic baseline.

## The geometry in one line

Delay of the site against the Earth's centre, for a pulsar at right ascension alpha, declination
delta, a site at distance rho from the spin axis, height z above the equatorial plane, longitude
lambda, Earth rotation angle theta(t):

  dt(t) = [ rho cos(delta) cos(theta(t) + lambda - alpha) + z sin(delta) ] / c

- The daily sine's phase gives the longitude, its amplitude gives rho: two of the three position
  components from one pulsar.
- The constant z sin(delta) / c cannot be told from that pulsar's unknown phase offset and the
  clock offset. z needs something more: the site on the Earth's surface (then rho gives the
  latitude, north error = radial error / sin(lat), with a hemisphere choice and a weak band near the
  equator), or absolute pulse phases, or pulsars at different declinations.

## What is already in the project

- `logs/2026-11-04_static-redesign.md` builds on exactly this for a static site: item 3 (Z as a
  constant per pulsar, the two routes to it), the expected 30 / 47 km east / north after 30 days
  at 52 N, and the cold-start search.
- Block D showed it in simulation for a known rough site: B0329+54 alone finds a site 100 km off
  to about ±30 km after 30 days.
- The amplitude at 52 N: rho = 3930 km = 13 ms of light time; ×cos(delta) = 7.6 ms for B0329+54,
  against a TOA error of about 0.74 ms.

## What to explore (questions, not results)

1. **Observability, written out.** Which of X, Y, Z, clock offset, clock rate and the per-pulsar
   phase constants are separable, for 1, 2, 3 pulsars and for 1 day, several days, a year. A
   covariance study in the style of the 4 Nov Python toy.
2. **More than the daily sine.** The Earth's orbit moves the geocentre, which is known, so it does
   not locate the site by itself. Does anything of the annual motion help: the parallax-like
   curvature term, the Doppler-shifted period, the changing direction of the site's velocity?
3. **Rate as a second measurement.** Besides the arrival time, the pulse frequency (Doppler) and
   its change through the day carry the same geometry. Is anything gained by fitting them, or is it
   all in the TOAs already? (The satellite Doppler systems, Transit and Argos, fix a position from
   one moving source this way; here the receiver moves and the source is fixed.)
4. **The clock.** A daily clock term (temperature) has the same period as the signal (FM-11). How
   many pulsars, at which right ascensions, separate it? What does a better oscillator buy?
5. **A receiver that also moves by itself** (the ship, later). Known own motion adds baseline
   directions that the Earth's rotation does not give, such as along the spin axis. How much
   north-south travel fixes z without absolute phases?
6. **Time to a fix.** How does the 3-D precision grow with observing time, compared with the
   snapshot picture of four pulsars at once?
7. **What must be known for it to work:** the Earth's orientation (UT1: 1 s = 286 m east at 52 N),
   the pulsar's position and spin model, the height or a height model.

## Where it connects

The open decisions on the D6 plan (`logs/2026-11-13_handover-static-nav.md`) come first: the
redesign is this idea applied to the static site. This note is for the wider question behind it.
