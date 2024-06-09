close all;
clear all;
clc;
% cd('inharmonicity_calculator'); % change into working directory

%% MEASUREMENTS OF SIGNALS
% Data Flute A5
nota2_Flute = [673.5802	673.5802	1	1
1.3472e+03	1.3480e+03	0.1445	1.0006
2.0207e+03	2.0176e+03	0.4193	0.9985
2.6943e+03	2.6943e+03	0.0641	1
3.3679e+03	3.3608e+03	0.0327	0.9979
4.0415e+03	4.0350e+03	0.0169	0.9984
4.7151e+03	4.7055e+03	0.0295	0.9980
5.3886e+03	5.3780e+03	0.0218	0.9980
6.0622e+03	6.0505e+03	0.0301	0.9981
6.7358e+03	6.7231e+03	0.0094	0.9981] ;  

% Data Flute C5
nota1_Flute = [528.9587	528.9587	1	1
1.0579e+03	1.0587e+03	0.0379	1.0007
1.5869e+03	1.5883e+03	0.0638	1.0009
2.1158e+03	2.1180e+03	0.0108	1.0010
2.6448e+03	2.6529e+03	0.0212	1.0031
3.1738e+03	3.1774e+03	0.0068	1.0012
3.7027e+03	3.7074e+03	0.0096	1.0013
4.2317e+03	4.2371e+03	0.0048	1.0013
4.7606e+03	4.7660e+03	0.0155	1.0011
5.2896e+03	5.3063e+03	0.0039	1.0032 ];

% Data Guitar A5
nota1_Git = [876.3475	876.3475	1	1
1.7527e+03	1.7507e+03	0.9266	0.9988
2.6290e+03	2.6270e+03	0.1293	0.9992
3.5054e+03	3.5095e+03	0.2761	1.0012
4.3817e+03	4.3919e+03	0.1512	1.0023
5.2581e+03	5.2784e+03	0.0614	1.0039
6.1344e+03	6.1710e+03	0.0782	1.0060
7.0108e+03	7.0169e+03	0.0055	1.0009
7.8871e+03	7.9725e+03	0.0106	1.0108
8.7635e+03	8.7858e+03	0.0050	1.0026 ];

% Data Guitar C5
nota2_Git = [522.9747	522.9747	1	1
1.0459e+03	1.0459e+03	0.2030	1
1.5689e+03	1.5658e+03	0.3156	0.9980
2.0919e+03	2.0919e+03	0.1253	1
2.6149e+03	2.6149e+03	0.2176	1
3.1378e+03	3.1410e+03	0.3375	1.0010
3.6608e+03	3.6671e+03	0.0519	1.0017
4.1838e+03	4.1932e+03	0.0398	1.0023
4.7068e+03	4.7068e+03	0.0061	1
5.2297e+03	5.2518e+03	0.0271	1.0042 ];


%% TEST VALIDATION
% load nota_A3.mat
% Amp_Flute_A5 =out.ModArm(1:8,3); 
% F_real_Flute_A5 =out.ModArm(1:8,2);
% instruments(1).frequency = freq_cent(F_real_Flute_A5);
% instruments(1).amplitude = Amp_Flute_A5;

%% TASK 1:  SPECTRUM PLOT
% Define the audio files and note names
subfolder = 'bassoon_inharmonicity';  % Replace with your subfolder name

audio_files = {fullfile(subfolder,'A5_Guitar.wav'), fullfile(subfolder,'C5_Guitar.wav'), ...
               fullfile(subfolder,'A5_Flute.wav'), fullfile(subfolder,'C5_Flute.wav')};
note_names = {'A5 Guitar', 'C5 Guitar', 'A5 Flute', 'C5 Flute'};

% Call the function to plot the spectra
plotSpectra(audio_files, note_names);

%%  FLUTE eigenfrequencies fn, amplitudes An 
Amp_Flute_A5 =(nota1_Flute(1:8,3)); 
F_real_Flute_A5 =(nota1_Flute(1:8,2));
Amp_Flute_C5 = (nota2_Flute(1:8,3));
F_real_Flute_C5 =(nota2_Flute(1:8,2));

% GUITAR eigenfrequencies fn, amplitudes An
Amp_Git_A5 =(nota1_Git(1:8,3));
F_real_Git_A5 =(nota1_Git(1:8,2));
Amp_Git_C5 =(nota2_Git(1:8,3));
F_real_Git_C5 =(nota2_Git(1:8,2));

%% ----------------------------- Calculation of the spectral parameters----------------------
% Pre-allocate output array for calculations
instruments(1).frequency = F_real_Flute_A5;
instruments(1).amplitude = Amp_Flute_A5;
instruments(2).frequency = F_real_Flute_C5;
instruments(2).amplitude = Amp_Flute_C5;
instruments(3).frequency = F_real_Git_A5;
instruments(3).amplitude = Amp_Git_A5;
instruments(4).frequency = F_real_Git_C5;
instruments(4).amplitude = Amp_Git_C5;


%% Spectral centroid brightness
notes_centroids = zeros(1, length(instruments));
spectral_content_oddh_evenh = zeros(2, length(instruments));
tristimulus_results = zeros(3, length(instruments));
inh_results = zeros(8,length(instruments));
even_h_results = zeros(1,length(instruments));

%Loop through each instrument
for ii = 1:length(instruments)
  % Extract frequency and amplitude vectors
  f = instruments(ii).frequency;
  A = instruments(ii).amplitude;
  
  % Call function to calculate centroid
  notes_centroids(ii) = spectral_centroid(f, A);
  
  % Call function to calculate spectral content (odd and even harmonics)
  % for each note
  [even_h,odd_h] = harmonic_ratios(instruments(ii).amplitude)
  even_h_results(:,ii)=even_h;
  spectral_content_oddh_evenh(:,ii) = [even_h odd_h];
  
  % Tristimulus calculation
  [T1, T2, T3] = tristimulus(instruments(ii).amplitude);
  tristimulus_results(:,ii) = [T1 T2 T3];

  % Inharmonicity calculation 
  inh_results(:,ii) = inharmonicity(instruments(ii).frequency, instruments(ii).amplitude);
  
end


%% TEST
% figure();
% plot(1:8, inh_results(:,1));
% legend('Test note bassoon');
% ylabel('Inharmonicity');
% xlabel('Harmonic');
% axis tight  % Fit the plot snugly to the data
% ylim([-0.4 0.8])
% axis([-0.0002 0.0008 ymin-(ymax-ymin)*0.1 ymax+(ymax-ymin)*0.1])  % Add 10% padding to top and bottom
% grid on;


figure();
plot(1:8, inh_results(:,1), 1:8, inh_results(:,2), 1:8, inh_results(:,3),1:8, inh_results(:,4));
legend('Flute A5','Flute C5','Guitar A5','Guitar C5');
ylabel('Inharmonicity');
xlabel('Harmonic');
axis tight  % Fit the plot snugly to the data
ylim([-0.004 0.008])
grid on;


%% --------------  3.1 COMPARISON AND ANALYSIS OF THE SPECTRAL PARAMETERS ----------------
% a) difference in B greater betweem different notes played on the same
% instrument

% Define the spectral centroid values for each condition
BC_flute_A5 = notes_centroids(1); 
BC_flute_C5 = notes_centroids(2);
BC_guitar_A5 = notes_centroids(3); 
BC_guitar_C5 = notes_centroids(4);

% Calculate the differences for the same instrument with different notes
Delta_BC_flute = abs(BC_flute_A5 - BC_flute_C5);
Delta_BC_guitar = abs(BC_guitar_A5 - BC_guitar_C5);

% Calculate the differences for different instruments playing the same note
Delta_BC_FluteGuitar_A5 = abs(BC_flute_A5 - BC_guitar_A5);
Delta_BC_FluteGuitar_C5 = abs(BC_flute_C5 - BC_guitar_C5);

% Display the results
fprintf('Difference in B for Flute (different notes): %f\n', Delta_BC_flute);
fprintf('Difference in B for Guitar (different notes): %f\n', Delta_BC_guitar);
fprintf('Difference in B for A5 (different instruments): %f\n', Delta_BC_FluteGuitar_A5);
fprintf('Difference in B for C5 (different instruments): %f\n', Delta_BC_FluteGuitar_C5);

% Analyze which difference is greater
if Delta_BC_flute > Delta_BC_FluteGuitar_A5 && Delta_BC_flute > Delta_BC_FluteGuitar_C5
    fprintf('The difference in B is greater between different notes played on the same Flute.\n');
elseif Delta_BC_guitar > Delta_BC_FluteGuitar_A5 && Delta_BC_guitar > Delta_BC_FluteGuitar_C5
    fprintf('The difference in B is greater between different notes played on the same Guitar.\n');
else
    fprintf('The difference in B is greater between different instruments playing the same note.\n');
end

% Create a bar plot to visualize the differences
figure;
categories = {'Flute (A5-C5)', 'Guitar (A5-C5)', 'A5 (Diff Instruments)','C5 (Diff Instruments)'};
values = [Delta_BC_flute, Delta_BC_guitar, Delta_BC_FluteGuitar_A5, Delta_BC_FluteGuitar_C5];

bar(values);
set(gca, 'XTickLabel', categories, 'XTickLabelRotation', 45);
ylabel('Spectral Centroid Difference');
title('Comparison of Spectral Centroid Differences');
legend('Flute A5 vs C5', 'Guitar A5 vs C5', 'Flute vs Guitar A5', 'Flute vs Guitar C5', 'Location', 'northwest');
grid on;



% B)
% Calculate the average spectral centroid for the chordophone
B_flute_avg = (BC_flute_A5 + BC_flute_C5) / 2;

% Calculate the average spectral centroid for the aerophone
B_guitar_avg = (BC_guitar_A5 + BC_guitar_C5) / 2;

% Display the average spectral centroid values
fprintf('Average spectral centroid for the flute: %f\n', B_flute_avg);
fprintf('Average spectral centroid for the guitar: %f\n', B_guitar_avg);

% Determine which instrument has a higher spectral centroid value
if B_flute_avg > B_guitar_avg
    fprintf('The flute has a higher average spectral centroid value.\n');
else
    fprintf('The guitar has a higher average spectral centroid value.\n');
end

% Create a bar plot to visualize the average spectral centroid values
figure;
categories = {'Flute', 'Guitar'};
values = [B_flute_avg, B_guitar_avg];

bar(values);
set(gca, 'XTickLabel', categories);
ylabel('Average Spectral Centroid');
title('Comparison of Average Spectral Centroid Values');
grid on;

%% C Which of the two instruments has a greater weight in the spectrum of even  harmonics?
heven_C1=even_h_results(1,3); heven_A1=even_h_results(1,1);
heven_C2=even_h_results(1,4); heven_A2=even_h_results(1,2);

spectrum_Weight(heven_C1,heven_A1,heven_C2,heven_A2);

%% D Define the Tristimulus values for each condition
% tristimulus = each row one amplitude, each column one note
% columns = flute_A5, Flute_C5, Guitar_A5, Guitar_C5
T1_C1=tristimulus_results(1,1); T1_A1=tristimulus_results(1,2); T1_C2=tristimulus_results(1,3); T1_A2=tristimulus_results(1,4); 
T2_C1=tristimulus_results(2,1); T2_A1=tristimulus_results(2,2); T2_C2=tristimulus_results(2,3); T2_A2=tristimulus_results(2,4); 
T3_C1=tristimulus_results(3,1); T3_A1=tristimulus_results(3,2); T3_C2=tristimulus_results(3,3); T3_A2=tristimulus_results(3,4);

Ts = [T1_C1,T1_A1,T1_C2,T1_A2, T2_C1,T2_A1,T2_C2,T2_A2, T3_C1,T3_A1,T3_C2,T3_A2];
plot_Tristimulus(Ts);

