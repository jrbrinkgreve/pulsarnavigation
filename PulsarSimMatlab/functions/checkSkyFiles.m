function d = checkSkyFiles(info_gen, info_disp, T, f_in, L, A, dutycycle, genEnvMode, seed, DM, fLow, fHigh, refFreq)
%CHECKSKYFILES  Parameters in which reused sky files differ from the current settings.
%{
The sky stage (generatePulsarSignal + applyDispersionStream) is the slow part, so main.m and
the run scripts reuse its files. This compares their info structs with the current
pipelineParams values and returns the names of the parameters that differ ({} when all
match). The caller decides what to do: main.m warns (main:stale), the run scripts that read
main's sky files stop.

  d = checkSkyFiles(info_gen, info_disp, T, f_in, L, A, dutycycle, genEnvMode, seed, ...
                    DM, fLow, fHigh, refFreq)

Compared: the generator's T, f_in, L, A, dutycycle, envelope mode and seed; the dispersion's
DM, band and reference frequency, and that the dispersed file was made from this generator
file. Why each matters, e.g.: a different refFreq dedisperses to another reference than the
data were dispersed to, so every TOA shifts by tau(old) - tau(new) (1.6 -> 1.4 GHz at DM 5:
-2.48 ms, below validateTOA's T/4 warning); a different dutycycle leaves the template unlike
the pulses on disk; a different seed is another pulsar realization.
(13 Oct 2026; before, main.m checked only T, f_in, L, DM and band of the sky files.)
%}

d = {};
if abs(info_gen.T - T) > 1e-12,                     d{end+1} = 'T';          end
if abs(info_gen.fs - f_in) > 1e-3,                  d{end+1} = 'f_in';       end
if abs(info_gen.L - L) > 1e-12,                     d{end+1} = 'L';          end
if abs(info_gen.A - A) > 1e-12 * max(1, abs(A)),    d{end+1} = 'A';          end
if abs(info_gen.dutycycle - dutycycle) > 1e-12,     d{end+1} = 'dutycycle';  end
if ~strcmpi(info_gen.envelopeMode, genEnvMode),     d{end+1} = 'genEnvMode'; end
if info_gen.seed ~= seed,                           d{end+1} = 'seed';       end
if abs(info_disp.DM - DM) > 1e-12,                  d{end+1} = 'DM';         end
if abs(info_disp.fLow - fLow) > 1e-3 || abs(info_disp.fHigh - fHigh) > 1e-3
    d{end+1} = 'band';
end
if abs(info_disp.refFreq - refFreq) > 1e-3,         d{end+1} = 'refFreq';    end
if ~strcmp(string(info_disp.inFile), string(info_gen.file))
    d{end+1} = 'dispersion input';
end
end
