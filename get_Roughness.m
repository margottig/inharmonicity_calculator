function roughness_value = get_Roughness(signal,fs)


% Compute the FFT
n = length(signal);               % Length of the signal
f = (0:n-1)*(fs/n);               % Frequency vector
y = fft(signal);                  % Compute the FFT
magnitude = abs(y);               % Magnitude of the FFT


% Identify the peaks
threshold = 0.01 * max(magnitude); %  peaks  with a height greater than 10% of the highest peak will be considered.
distance = 50;                    % minimum distance between adjacent peaks (minimum number of samples between consecutive peaks.)
[pks, locs] = findpeaks(magnitude, 'MinPeakHeight', threshold, 'MinPeakDistance', distance);

% Get the corresponding frequencies of the peaks
peakFrequencies = f(locs);
%concatenate vector from findpeaks matlab function
combinedMatrix = [peakFrequencies', pks];

roughness_value = roughness(combinedMatrix);

%% PLOT Peaks
% figure;
% plot(f, magnitude);
% hold on;
% plot(peakFrequencies, pks, 'ro'); % Mark peaks with red circles
% title('Frequency Spectrum with Peaks');
% xlabel('Frequency (Hz)');
% ylabel('Magnitude');
% grid on;

end


%% -----------------  COMPLEMENT FUNCTIONS -----------

function roughnessScalar = roughness(peaks)
    %% calculates the roughness of a single spectrogram frame
    %% assumes first bin is at 0 Hz, last bin at fs/2 Hz
    %% based on formula from MacCallum and Einbond
    
    %peaks = cell2mat(peaks);
    
    % first, check for empty case
    if (isempty(peaks))
        roughnessScalar = 0;
    else
        % first, extract peak frequencies and amplitudes from the cell
        freqs = peaks(:,1);
        mags = peaks(:,2);
    
        % now, start adding things together
        numerator = 0;
        denominator = 1;
    
        for j=1:length(freqs)
            for k=1:(j-1)
                a_1 = mags(j);
                a_2 = mags(k);
                f_1 = freqs(j);
                f_2 = freqs(k);
                numerator = numerator + (a_1 * a_2 * e_parncuttG(f_1, f_2));
                denominator = denominator + (a_1^2);
            end
        end
    
        roughnessScalar = numerator / denominator;
    end

end

function g = e_parncuttG(f_1, f_2)
    %% calculate "g(f_cb)" as cited from Parncutt in MacCallum
    f_m = mean([f_1 f_2]);
    f_cb = 1.72 * (f_m^0.65);
    
    % "distance between them in critical bands"
    f_cb = abs(f_1 - f_2) / f_cb;
    
    % f_cb's relation to the space between the pitches is unclear...
    % I should also consult Parncutt directly
    
    g = exp(1) * f_cb * 4;
    g = g * exp(-f_cb * 4);
    g = g^2;
    
end        