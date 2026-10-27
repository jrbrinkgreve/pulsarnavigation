# Which kinds of check each unit test contains (tests.html), read from the test descriptions on that page:
# R = regression (default output equal to a frozen reference, bit for bit), A = algebra (code equal to an independently computed
# formula or matrix), P = physics (a measured quantity, such as noise, false-flag rates or a found radar period, equal to its
# prediction within the predicted scatter); with the run times.
tests = [('testGeneratePulsarSignal', 'generatePulsarSignal', '', 'A', 'P', 5),   # run times: the table on tests.html
         ('testApplyDispersionStream', 'applyDispersionStream', '', 'A', 'P', 3),   # 27 Oct 2026 (review R-D1)
         ('testAddNoiseAndRFI', 'addNoiseAndRFI', '', 'A', 'P', 4),   # 27 Oct 2026 (review R-N1)
         ('testChannelizeIQ', 'channelizeIQ', '', 'A', 'P', 10),
         ('testDedisperseChannels', 'dedisperseChannels', 'R', 'A', 'P', 15),
         ('expBlankingVariance', 'experiment 3a', '', 'A', 'P', 15),
         ('testDetectChannels', 'detectChannels, powerCovariance', '', 'A', 'P', 1),
         ('testFoldWeights', 'foldProfile weights', 'R', 'A', 'P', 7),
         ('testChannelTOA', 'estimateTOA, detectPulsar', 'R', 'A', 'P', 15),
         ('testOptimalTOA', 'optimal weighting', 'R', 'A', 'P', 130),
         ('testBlankingWeights', 'blankingWeights', '', 'A', 'P', 110),
         ('testBlankChannels', 'blankChannels', '', 'A', '', 10),
         ('testDetectRFI', 'detectRFI', '', '', 'P', 150),
         ('testPeriodicRFI', 'periodicRFI', '', '', 'P', 35),
         ('testRadarEdges', 'rfiSource rise time', 'R', 'A', '', 5),
         ('testRfiGating', 'rfiSource rotation', 'R', 'A', '', 10),
         ('testRfiNoise', "rfiSource 'noise'", 'R', '', 'P', 5),
         # blocks C and D (20-22 Oct 2026), kinds of check from each test file's header
         ('testSimulateFold', 'simulateFold (fast simulator)', 'R', 'A', 'P', 60),
         ('testRFIModel', 'rfiProfile, simulateFold RFI', 'R', 'A', 'P', 30),
         ('testMakeTrajectory', 'makeTrajectory', '', 'A', '', 1),
         ('testDelayModel', 'delayModel (incl. DE440)', '', 'A', '', 2),
         ('testMakePredictor', 'makePredictor', 'R', 'A', '', 5),
         ('testFoldPredictor', "foldProfile 'Predictor'", '', 'A', '', 2),
         ('testTrueResiduals', 'trueResiduals', '', 'A', '', 2),
         ('testMakeClock', 'makeClock', '', 'A', 'P', 2),
         ('testNavFit', 'navFit', '', 'A', 'P', 120),
         ('testNavFilter', 'navFilter', '', 'A', 'P', 10),
         ('testLoadEphemeris', 'loadEphemeris (DE440)', '', 'A', '', 2)]
x0, y0, rh = 10, 52, 22
cx = [420, 510, 600]
bx0, bw = 670, 140
o = []
for x, lab in zip(cx, ['regression', 'algebra', 'physics']):
    o.append(f'<text class="ttl" x="{x}" y="{y0 - 22}" text-anchor="middle">{lab}</text>')
o.append(f'<text class="sub" x="{cx[0]}" y="{y0 - 8}" text-anchor="middle">bit for bit</text>')
o.append(f'<text class="sub" x="{cx[1]}" y="{y0 - 8}" text-anchor="middle">code = formula</text>')
o.append(f'<text class="sub" x="{cx[2]}" y="{y0 - 8}" text-anchor="middle">data = prediction</text>')
o.append(f'<text class="ttl" x="{bx0}" y="{y0 - 22}">run time</text>')
for i, (name, unit, r, a, p, sec) in enumerate(tests):
    y = y0 + i*rh
    if i % 2 == 0:
        o.append(f'<rect class="fillfg" x="{x0}" y="{y}" width="{bx0 + bw - x0}" height="{rh}" opacity="0.25"/>')
    o.append(f'<text class="mono" x="{x0 + 4}" y="{y + 15}">{name}</text>')
    o.append(f'<text class="sub" x="{x0 + 200}" y="{y + 15}">{unit}</text>')
    for x, flag, cls in zip(cx, (r, a, p), ('dot', 'dot', 'dot')):
        if flag:
            o.append(f'<circle cx="{x}" cy="{y + 11}" r="5.5" class="{cls}"/>')
        else:
            o.append(f'<circle cx="{x}" cy="{y + 11}" r="5" fill="none" stroke="var(--fg-muted)" stroke-width="1.2" opacity="0.6"/>')
    w = bw*sec/max(t[5] for t in tests)
    o.append(f'<rect class="box proc" x="{bx0}" y="{y + 5}" width="{max(w, 2):.1f}" height="{rh - 10}" rx="2"/>')
    o.append(f'<text class="sub" x="{bx0 + max(w, 2) + 5:.1f}" y="{y + 15}">{sec} s</text>')
H = y0 + len(tests)*rh + 14
svg = f'''<figure class="scroll">
<svg class="wide" viewBox="0 0 860 {H}" role="img" aria-label="Which kinds of check each unit test contains: regression, algebra and physics, with the run time of each test">
  ''' + '\n  '.join(o) + '''
</svg>
<figcaption>The 27 unit tests that <code>runAllTests</code> runs, in its order, and experiment 3a (sixth row; not run by <code>runAllTests</code>), with the kinds of check each one
contains (filled dot) and its approximate run time (the table below). The slow <code>testRFIChain</code>, ~2.5 h, is not shown. Regression checks exist where a function gained an option after it was validated; algebra checks
compare the code with an independent formula or matrix; physics checks compare a measured quantity (noise, false-flag rates, a
found radar period) with its prediction, within the predicted scatter.</figcaption>
</figure>'''
open('testmatrix.svg', 'w').write(svg)
