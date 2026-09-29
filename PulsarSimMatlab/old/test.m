f_in = 3.2e9;

info = dir('data/test.dat');
Nx = info.bytes / 4;
half = floor(Nx/2);

% original, first half
fid = fopen('data/test.dat', 'r');
xOrig = fread(fid, half, 'single=>single');
fclose(fid);

% dispersed, first half
fid = fopen('data/test_dispersed.dat', 'r');
yDisp = fread(fid, half, 'single=>single');
fclose(fid);

t = (0:half-1)/f_in;

figure;
subplot(2,1,1); plot(t, xOrig); title('Original, first half'); ylabel('amp');
subplot(2,1,2); plot(t, yDisp); title('Dispersed, first half'); ylabel('amp'); xlabel('t [s]');