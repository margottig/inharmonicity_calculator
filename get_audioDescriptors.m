clear all;
close all;  

%% Automatized running of all Audiofiles
inputDir = 'audiofiles';
outputDir = 'results';
audioFiles = dir(fullfile(inputDir, '*.wav'));

for k = 1:length(audioFiles)
    
    filepath = fullfile(audioFiles(k).folder,audioFiles(k).name);    
    [~, baseName, ~] = fileparts(filepath);
    filename_plot = fullfile(outputDir,[baseName '_plot.png']);
    filename_mat = fullfile(outputDir,[baseName '_results.mat']);
    sample_results = struct(); % initialized empty struct
    
    %% Start Code
    window = hamming(1024); %fftSize = 1024; hamming?
    noverlap = 512; % hopSize = 512;
    nfft = 1024; % number of fft points try 256
    
    %% AUDIOREAD
    [fullSignal, fs] = audioread(filepath);
    fullSignal = fullSignal(:,1); % mono signal
    auInfo = audioinfo(filepath);
    [s,f,t] = spectrogram(fullSignal, window, noverlap, nfft, fs, 'yaxis');
    %[s,f,t] = spectrogram(fullSignal, 256, 250, 256, fs, 'yaxis');
    
    % make magnitude spectrogram
    % espectrograma = abs(s).^2; %when raised to the square, envelope is obtained?
    espectrograma = abs(s);
    
    %% Calculate Spectral Descriptors
    fullCentroid = b_spectralCentroid(espectrograma, f); % Spectral Centroid
    spectralFlatness = get_spectralFlatness(espectrograma);
    spectralSpread = get_spectralSpread(espectrograma, fullCentroid,f);
    
    %% Calculate Harmonic Descriptors
    % inharmonicity (given a spectrogram, calculates inharmonicity for each frame)
    inharmonicity_matrix  = get_Inharmonicity_perFrame(espectrograma, t,f);
    inharmonicity_vector = mean(inharmonicity_matrix, 2);
    % Tristimulus calculation
    tristimulus_vector = harmonicEnergy(espectrograma);
    % harmonic odd to even ratio / even harmonic ratios (debug size vector)
    % [harmonicOddToEvenRatio_vector, even_harmonic_results] = harmonicOddToEvenRatio(espectrograma); 

    %% Calculate Roughness
    roughness_value = get_Roughness(fullSignal,fs); 

        
    %% now, plot!
    figure;
    
    % timepoints for waveform (in minutes!)
    waveformX = (1 / fs / 60) .* (1:length(fullSignal));
    envolvente = envelope(fullSignal,300,'peak');
    subplot(3,1,1);
    plot(waveformX, fullSignal);
    hold on;
    plot(waveformX,envolvente, 'r', 'LineWidth', 1);
    title("waveform");
    xlabel("time (minutes)");
    xlim([0 waveformX(end)]);
    ylabel("sample value");
    legend('Original Signal', 'Envelope', 'Location', 'southeast');
    
    % timepoints for spectral measures (in minutes!)
    value_of_interest = mean(fullCentroid); % Example value in seconds
    range = 2000; % Range around the value of interest in seconds
    x_max = min(max(f), value_of_interest + range);
    subplot(3,1,2);
    imagesc(t*0.0166667, f, 10*log10(abs(s))); % Plot in decibels
    axis xy; % Set axes to normal orientation
    %colorbar; % Add color bar to indicate power levels
    xlabel('Time (minutes)');
    ylabel('Frequency (Hz)');
    ylim([0 x_max]);
    title('Spectrogram');
    
    
    % timepoints for spectral measures (in minutes!)
    spectralX = (noverlap / fs / 60) .* (1:width(espectrograma));
    subplot(3,1,3);
    plot(spectralX, fullCentroid);
    title("spectral centroid");
    xlabel("time (minutes)");
    xlim([0 spectralX(end)]);
    ylabel("spectral centroid (Hz)");
    
    sgtitle(baseName);

    % Save the plot as a PNG file
    saveas(gcf, filename_plot,'png');

    %% SAVE IN STRUCT
    sample_results.espectrograma = espectrograma;
    sample_results.fullCentroid = fullCentroid;
    sample_results.spectralFlatness = spectralFlatness;
    sample_results.spectralSpread = spectralSpread;
    sample_results.roughness_value = roughness_value;
    sample_results.envolvente = envolvente;
    sample_results.f = f;
    sample_results.t = t;
    sample_results.inharmonicity_vector = inharmonicity_vector;
    sample_results.tristimulus_vector = tristimulus_vector;

    save(filename_mat, 'sample_results');
end




%% ----------------- CALL FUNCTIONS ---------------------------------------
% given a spectrogram, calculates spectral centroid, and so.. for each frame
% assumes columns are frames, rows cover exactly the positive frequency range

function centroidVector = b_spectralCentroid(spectro, f)

centroidVector = zeros(1, width(spectro));
for ii = 1:length(centroidVector)
    centroidVector(ii) = spectral_centroid(f,spectro(:,ii));
end
end

function [harmonicOddToEvenRatio_vector, even_harmonic_results]  = harmonicOddToEvenRatio(spectro)
    %[even_harmonics, odd_harmonics]
    harmonicOddToEvenRatio_vector = zeros(2, width(spectro));
    even_harmonic_results = zeros(1,width(spectro));
    
    for ii = 1:width(harmonicOddToEvenRatio_vector)
        total_power = sum(spectro(:,ii).^2);
        % Calculate the even harmonic ratio
        even_indices = 2:2:length(spectro); % Indices of even harmonics
        even_harmonics = sum(spectro(even_indices,ii).^2) / total_power;
        % Calculate the odd harmonic ratio
        odd_indices = 1:2:length(spectro); % Indices of odd harmonics
        odd_harmonics = sum(spectro(odd_indices,ii).^2) / total_power;
        
        harmonicOddToEvenRatio_vector(:,ii) = [even_harmonics odd_harmonics];
        even_harmonic_results(:,ii) = even_harmonics;
    end
end

function tristimulus_vector = harmonicEnergy(spectro)

    tristimulus_vector = zeros(3, width(spectro));
    for ii = 1:length(tristimulus_vector)
        [T1, T2, T3] = tristimulus(spectro(:,ii));
        tristimulus_vector(:,ii) = [T1 T2 T3];
    end  
end