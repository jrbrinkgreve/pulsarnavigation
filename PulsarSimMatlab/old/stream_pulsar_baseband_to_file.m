function stream_pulsar_baseband_to_file(filename, varargin)
% STREAM_PULSAR_BASEBAND_TO_FILE
% Streams 1.0 GS/s complex baseband voltage data directly to binary file
% without causing Memory/OOM errors.

% Default Parameters
p = inputParser;
addParameter(p, 'duration_sec', 20.0);    % Total simulation duration (s)
addParameter(p, 'chunk_ms', 10.0);        % Chunk size in milliseconds (1 ms = 1e6 samples)
addParameter(p, 'fs', 1.0e9);            % Sampling rate: 1.0 GS/s
addParameter(p, 'f_center', 1.4e9);      % 1.4 GHz center frequency
addParameter(p, 'bw', 400.0e6);          % 400 MHz bandwidth (1.2 - 1.6 GHz)
addParameter(p, 'dm', 50.0);             % Dispersion Measure (pc cm^-3)
addParameter(p, 'period', 0.033);        % Pulsar period (e.g., Crab = 33 ms)
addParameter(p, 'pulse_width', 0.001);   % Pulse width (1 ms)
addParameter(p, 'snr_v', 0.05);          % Single-sample voltage SNR (< 0 dB)
addParameter(p, 'scale_factor', 30.0);   % Quantization scaling to fill int8 dynamic range
parse(p, varargin{:});
opts = p.Results;

% Derived Parameters
samples_per_chunk = round((opts.chunk_ms / 1000.0) * opts.fs);
total_chunks = round(opts.duration_sec / (opts.chunk_ms / 1000.0));
dt = 1.0 / opts.fs;
D_const = 4.148808e3; % MHz^2 pc^-1 cm^3 s

% 1. Pre-compute Dispersion Transfer Function H(f) ONCE
% FFT frequency vector (-fs/2 to +fs/2)
freqs = (-samples_per_chunk/2 : (samples_per_chunk/2 - 1)) * (opts.fs / samples_per_chunk);
freqs = fftshift(freqs); % Match MATLAB fft indexing

H_f = ones(1, samples_per_chunk, 'like', complex(0));
for i = 1:samples_per_chunk
    f_offset = freqs(i);
    if abs(f_offset) <= (opts.bw / 2.0)
        % Quadratic phase delay for dispersion
        phase_shift = 2 * pi * D_const * opts.dm * (f_offset / 1e6)^2 / ((opts.f_center / 1e6)^3);
        H_f(i) = exp(1i * phase_shift);
    else
        H_f(i) = 0.0 + 0.0i; % Sharp bandpass filter (1.2 - 1.6 GHz)
    end
end

% 2. Pre-allocate fixed working buffers (RAM usage stays constant)
v_out_bytes = zeros(1, 2 * samples_per_chunk, 'int8');

% Open binary file for writing ('w' mode)
fileID = fopen(filename, 'w');
if fileID == -1
    error('Failed to open file %s for writing.', filename);
end
cleanup = onCleanup(@() fclose(fileID)); % Guarantees file closure on exit/error

fprintf('Streaming %.2f sec of baseband data to %s...\n', opts.duration_sec, filename);
fprintf('Chunk size: %d samples (%.2f MB RAM working memory)\n', ...
    samples_per_chunk, (samples_per_chunk * 16) / 1e6);

t_global = 0.0;
pulse_sigma = (opts.pulse_width / opts.period) / 2.355;

% 3. Main Streaming Loop
for chunk = 1:total_chunks
    % Time vector for the current chunk
    t_vec = t_global + (0 : samples_per_chunk - 1) * dt;

    % Calculate global pulsar phase phi(t)
    phase = mod(t_vec, opts.period) ./ opts.period;
    A_t = exp(-0.5 * ((phase - 0.5) ./ pulse_sigma).^2);

    % Synthesize complex Gaussian noise modulated by pulse envelope sqrt(A(t))
    v_intrinsic = sqrt(A_t) .* (randn(1, samples_per_chunk) + 1i * randn(1, samples_per_chunk));

    % Apply dispersion in Fourier Domain
    V_freq = fft(v_intrinsic);
    v_dispersed = ifft(V_freq .* H_f);

    % Add uncorrelated system/receiver noise
    v_noise = randn(1, samples_per_chunk) + 1i * randn(1, samples_per_chunk);
    v_total = (opts.snr_v * v_dispersed) + v_noise;

    % Quantize Float -> Int8 (Interleaved Real and Imaginary: [R1, I1, R2, I2, ...])
    r_int = int8(min(max(round(real(v_total) * opts.scale_factor), -128), 127));
    i_int = int8(min(max(round(imag(v_total) * opts.scale_factor), -128), 127));

    v_out_bytes(1:2:end) = r_int;
    v_out_bytes(2:2:end) = i_int;

    % Stream raw byte vector directly to disk
    fwrite(fileID, v_out_bytes, 'int8');

    t_global = t_global + samples_per_chunk * dt;
end

fprintf('Done! Total written file size: %.2f GB\n', ...
    (total_chunks * 2 * samples_per_chunk) / 1e9);
end