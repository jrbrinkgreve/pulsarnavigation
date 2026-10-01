function info = loadInfo(dataFile)
%LOADINFO  Load the <name>_info.mat saved next to a stage's output file.
%{
Every stage saves its info struct as <outputname>_info.mat next to its data
file. This reloads it, so a script can skip a stage and continue from the
files on disk.

  info = loadInfo(dataFile)     % e.g. loadInfo("data/test_dispersed.dat")
%}
[d, name] = fileparts(dataFile);
f = fullfile(d, name + "_info.mat");
if ~isfile(f)
    error('loadInfo:noInfo', 'No info file "%s"; run the stage that creates it first.', f);
end
s = load(f, 'info');
info = s.info;
end
