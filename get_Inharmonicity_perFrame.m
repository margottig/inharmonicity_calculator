function inharmonicity_vector = get_Inharmonicity_perFrame(espectro, T, F)

    % Convert the spectrogram to decibels
    S_dB = 20*log10(abs(espectro));
    
    % Initialize an array to store the fundamental frequencies
    fundamental_freqs = zeros(1, length(T));
    
    % Loop through each time slice to find the peak frequency
    for ii = 1:length(T)
        % Get the magnitude spectrum for the current time slice
        spectrum_slice = S_dB(:, ii);
        % Find the frequency index with the maximum magnitude
        [~, maxIdx] = max(spectrum_slice);
        % Store the corresponding frequency
        fundamental_freqs(ii) = F(maxIdx);
    end
    
    % Calculate the most consistent fundamental frequency (mode)
    fundamental_freq = mode(fundamental_freqs);
    f1 = fundamental_freq;
    % Display the fundamental frequency
    disp(['Fundamental Frequency: ', num2str(fundamental_freq), ' Hz']);
    
    inharmonicity_vector = zeros(size(espectro));
    
    for ii = 1:width(espectro)           
           k = [1,2,3,4,5]; %constants
            n_num = 1:1:length(F);
            Ak = espectro(1:5,ii);  % amplitude of the kth harmonic
            DFk = F(1:5) - k'.*f1; %Dfn
            DFk = DFk';
            numerator = sum(Ak.*(DFk./(k'*f1)));
            numerator = sum(numerator);
            denominator = sum(Ak);
            fc = f1 * (1 + (numerator/denominator));
     
            % Calculate inharmonicity
            fn = F;
            nfc = n_num.*fc;
    
            % Calculate inharmonicity for each harmonic
            inh = ((fn' -nfc ) ./ (n_num.*f1))/10;
            inharmonicity_vector(:,ii) = inh';
        end   
    
    % Plot the spectrogram with the fundamental frequency overlaid
    % figure;
    % surf(T, F, S_dB, 'EdgeColor', 'none');
    % axis xy; axis tight; colormap(jet); view(0, 90);
    % xlabel('Time (s)');
    % ylabel('Frequency (Hz)');
    % title('Spectrogram with Fundamental Frequency');
    % hold on;
    % plot(T, fundamental_freqs, 'w', 'LineWidth', 2);
    % hold off;
    
    end
    
    
    
    
    
    
    