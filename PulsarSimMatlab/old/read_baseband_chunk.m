function v_complex = read_baseband_chunk(filename, start_sample, n_samples)
% Reads a specific segment of complex Int8 voltage data from disk
fileID = fopen(filename, 'r');

% Seek to start position (2 bytes per complex sample)
fseek(fileID, (start_sample - 1) * 2, 'bof');

% Read interleaved Int8 bytes [Real1, Imag1, Real2, Imag2, ...]
raw_bytes = fread(fileID, 2 * n_samples, 'int8=>double');
fclose(fileID);

% De-interleave into complex array
v_real = raw_bytes(1:2:end);
v_imag = raw_bytes(2:2:end);
v_complex = v_real + 1i * v_imag;
end