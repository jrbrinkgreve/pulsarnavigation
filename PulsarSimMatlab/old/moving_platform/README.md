# Moving-platform code (archived 5 Nov 2026)

Jasper, 5 Nov 2026: "to prevent bloat, I want to work towards a static earth position, with code later developed if we want
to move to a moving platform", and "you can archive the old code". Nothing in this folder is on the MATLAB path.

Moved here, unchanged:
- `navFilter.m`, `testNavFilter.m`: the dead-reckoning filter and its test.
- `runShip.m`: the ship scenario.

Copies of the last versions with moving-platform support (commit 4433443). The live versions in `functions/` and `tests/` are simpler:
- `makeTrajectory.m`: also the `'track'` and `'barycentric'` types.
- `makePredictor.m`: also splits a segment that misses its check, around a ship's turn.
- `testMakeTrajectory.m`, `testMakePredictor.m`, `testDelayModel.m`: their tests of those features.
- `gen_predseg.py`: the generator of the old Figure 7a.3 (it needs `docs/tools/figures/`).
- `fold_before_static.html`: the stage 7 docs page with its ship material, for reference (it needs `docs/style.css` to look right).

To run the ship code again, put this folder on the path ahead of `functions/`, or restore the versions above.
