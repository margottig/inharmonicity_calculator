function plotSpectra(audio_files, note_names)
    % Initialize figure
    figure;
    hold on;

    % Loop through each note
    for i = 1:length(audio_files)
        [audio, fs] = audioread(audio_files{i});
        
        % Extract a segment (for example, the first 1 second)
        duration = 1; % duration in seconds
        segment_length = min(round(duration * fs), length(audio));
        segment = audio(1:segment_length);
        
        % Apply a window function
        window = hamming(length(segment));
        windowedSegment = segment .* window;
        
        % Compute the FFT
        N = length(windowedSegment);
        Y = fft(windowedSegment);
        
        % Compute the magnitude of the FFT and normalize
        Y_mag = abs(Y) / N;
        Y_mag = Y_mag(1:floor(N/2)); % Only need the positive frequencies
        f = (0:floor(N/2)-1) * (fs / N); % Frequency axis
        
        % Identify the fundamental frequency
        [~, idx] = max(Y_mag);
        fundamentalFreq = f(idx);
        disp(['Fundamental Frequency of ', note_names{i}, ': ', num2str(fundamentalFreq), ' Hz']);
        
        % Plot the spectrum
        plot(f, Y_mag, 'DisplayName', note_names{i},'LineWidth',2);
    end

    % Add plot details
    title('Magnitude Spectrum of Different Notes');
    xlabel('Frequency (Hz)');
    ylabel('Magnitude');
    legend show;
    grid on;
    hold off;
end
