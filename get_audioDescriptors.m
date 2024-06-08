clear all;
close all;  

filepath = 'bassoon_A3.mp3';
window = hamming(1024); %fftSize = 1024; hamming?
noverlap = 512; % hopSize = 512;
nfft = 1024; % number of fft points try 256


%% AUDIOREAD
[fullSignal, fs] = audioread(filepath);
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
% inharmonicity

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


%% ----------------- FUNCTIONS ------------------

function centroidVector = b_spectralCentroid(spectro, f)
%% given a spectrogram, calculates spectral centroid for each frame
%% assumes columns are frames, rows cover exactly the positive frequency range
centroidVector = zeros(1, width(spectro));
for i = 1:length(centroidVector)
    centroidVector(i) = spectral_centroid(f,spectro(:,i));
end
end


