function process_large_baseband_file(filename, n_channels)
    fileID = fopen(filename, 'r');
    cleanup = onCleanup(@() fclose(fileID));

    % Block size parameters
    samples_per_block = 1e6; % 1 ms block at 1 GS/s
    bytes_per_sample = 2;    % Int8 Real + Int8 Imag
    block_bytes = samples_per_block * bytes_per_sample;

    % Dynamic Spectrum Accumulator (pre-allocate dynamic growth buffer)
    spectrogram_accumulated = [];

    % Setup Figure for Live Sanity Check Display
    fig = figure('Color', 'w', 'Name', 'Spectrogram Sanity Check', 'Position', [100, 100, 800, 450]);
    ax = axes('Parent', fig);

    block_count = 0;

    while ~feof(fileID)
        % 1. Read block from disk (Uses ~2 MB of RAM)
        raw_bytes = fread(fileID, block_bytes, 'int8=>double');
        if isempty(raw_bytes)
            break;
        end

        % 2. De-interleave into complex voltage
        v_complex = raw_bytes(1:2:end) + 1i * raw_bytes(2:2:end);

        % 3. Run DSP (Channelize via FFT)
        n_bins = floor(length(v_complex) / n_channels);
        v_matrix = reshape(v_complex(1 : n_channels * n_bins), n_channels, n_bins);

        % Compute power spectrum and shift 0 Hz center
        spectrogram_block = fftshift(abs(fft(v_matrix, n_channels, 1)).^2, 1); 

        % 4. Accumulate & Display Sanity Check
        % Append the processed block time-bins to the dynamic spectrum
        spectrogram_accumulated = [spectrogram_accumulated, spectrogram_block]; 
        
        block_count = block_count + 1;

        % --- LIVE SANITY CHECK PLOT ---
        imagesc(ax, spectrogram_accumulated);
        axis(ax, 'xy');
        colormap(ax, flipud(hot));
        colorbar(ax);
        title(ax, sprintf('Sanity Check: Streaming Block %d (Total Time Bins: %d)', ...
              block_count, size(spectrogram_accumulated, 2)));
        xlabel(ax, 'Time Bins (1000 samples / bin)');
        ylabel(ax, 'Frequency Channels');
        grid(ax, 'on');

        % Refresh figure window without interrupting execution
        drawnow limitrate;
    end

    fprintf('Successfully processed file in streaming mode.\n');
end