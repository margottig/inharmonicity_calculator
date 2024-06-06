close all;
clear all;
clc;


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


%  FLUTE eigenfrequencies fn, amplitudes An 
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

%Loop through each instrument
for ii = 1:length(instruments)
  % Extract frequency and amplitude vectors
  f = instruments(ii).frequency;
  A = instruments(ii).amplitude;
  
  % Call function to calculate centroid
  notes_centroids(ii) = spectral_centroid(f, A);
  
  % Call function to calculate spectral content (odd and even harmonics)
  % for each note
  [even_h,odd_h] = harmonic_ratios(instruments(ii).amplitude);
  spectral_content_oddh_evenh(:,ii) = [even_h odd_h];
  
  % Tristimulus calculation
  [T1, T2, T3] = tristimulus(instruments(ii).amplitude);
  tristimulus_results(:,ii) = [T1 T2 T3];

  % Inharmonicity calculation 
  inh_results(:,ii) = inharmonicity(instruments(ii).frequency, instruments(ii).amplitude);
  
end


%% TEST
% figure(2);
% plot(1:8, inh_results(:,1));
% legend('Test note bassoon');
% ylabel('Inharmonicity');
% xlabel('Harmonic');
% axis tight  % Fit the plot snugly to the data
% ylim([-0.4 0.8])
% axis([-0.0002 0.0008 ymin-(ymax-ymin)*0.1 ymax+(ymax-ymin)*0.1])  % Add 10% padding to top and bottom
% grid on;


figure(1);
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
BA_guitar_C5 = notes_centroids(4);

% Calculate the differences for the same instrument with different notes
Delta_BC_flute = abs(BC_flute_A5 - BC_flute_C5);
Delta_BC_guitar = abs(BC_guitar_A5 - BA_guitar_C5);

% Calculate the differences for different instruments playing the same note
Delta_BC_FluteGuitar_A5 = abs(BC_flute_A5 - BC_guitar_A5);
Delta_BC_FluteGuitar_C5 = abs(BC_flute_C5 - BA_guitar_C5);

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

% Define the spectral centroid values for each condition
BC1 = % Your calculated value for chordophone-note 1
BA1 = % Your calculated value for aerophone-note 1
BC2 = % Your calculated value for chordophone-note 2
BA2 = % Your calculated value for aerophone-note 2

% Calculate the average spectral centroid for the chordophone
B_C_avg = (BC1 + BC2) / 2;

% Calculate the average spectral centroid for the aerophone
B_A_avg = (BA1 + BA2) / 2;

% Display the average spectral centroid values
fprintf('Average spectral centroid for chordophone: %f\n', B_C_avg);
fprintf('Average spectral centroid for aerophone: %f\n', B_A_avg);

% Determine which instrument has a higher spectral centroid value
if B_C_avg > B_A_avg
    fprintf('The chordophone has a higher average spectral centroid value.\n');
else
    fprintf('The aerophone has a higher average spectral centroid value.\n');
end

% Create a bar plot to visualize the average spectral centroid values
figure;
categories = {'Chordophone', 'Aerophone'};
values = [B_C_avg, B_A_avg];

bar(values);
set(gca, 'XTickLabel', categories);
ylabel('Average Spectral Centroid');
title('Comparison of Average Spectral Centroid Values');
grid on;













%% -------------------------------FUNCTIONS--------------------
%SPECTRAL CENTROID
function B = spectral_centroid(frequencies, amplitudes)
    % This function calculates the spectral centroid B
    % Inputs:
    %   frequencies - array of eigenfrequencies (f1, f2, ..., fN)
    %   amplitudes - array of amplitudes (A1, A2, ..., AN)
    % Output:
    %   B - spectral centroid

    % Check if inputs are vectors and have the same length
    if ~isvector(frequencies) || ~isvector(amplitudes) || length(frequencies) ~= length(amplitudes)
        error('Inputs must be vectors of the same length');
    end

    % Calculate the numerator and denominator of the spectral centroid formula
    numerator = sum(frequencies .* amplitudes);
    denominator = sum(amplitudes);

    % Calculate the spectral centroid
    B = numerator / denominator;
end

% Spectral content of odd and even harmonics
function [even_h, odd_h] = harmonic_ratios(amplitudes)
    % This function calculates the even and odd harmonic ratios
    % Input:
    %   amplitudes - array of amplitudes (A1, A2, ..., AN)
    % Outputs:
    %   even_h - even harmonic ratio
    %   odd_h - odd harmonic ratio

    % Number of harmonics
    N = 8;

    % Sum of squares of all amplitudes
    total_power = sum(amplitudes.^2);

    % Calculate the even harmonic ratio
    even_indices = 2:2:N; % Indices of even harmonics
    even_h = sum(amplitudes(even_indices).^2) / total_power;

    % Calculate the odd harmonic ratio
    odd_indices = 1:2:N; % Indices of odd harmonics
    odd_h = sum(amplitudes(odd_indices).^2) / total_power;

end

function [T1, T2, T3] = tristimulus(amplitudes)
    % This function calculates the tristimulus values T1, T2, T3
    % Input:
    %   amplitudes - array of amplitudes (A1, A2, ..., AN)
    % Outputs:
    %   T1 - Tristimulus value T1
    %   T2 - Tristimulus value T2
    %   T3 - Tristimulus value T3

    % Check if input is a vector
    if ~isvector(amplitudes)
        error('Input must be a vector');
    end

    % Number of harmonics
    N = length(amplitudes);

    % Sum of squares of all amplitudes
    total_power = sum(amplitudes.^2);

    % Calculate T1
    T1 = amplitudes(1)^2 / total_power;

    % Calculate T2
    if N > 1
        T2 = sum(amplitudes(2:3).^2) / total_power;
    else
        T2 = 0;
    end

    % Calculate T3
    if N > 3
        T3 = sum(amplitudes(4:end).^2) / total_power;
    else
        T3 = 0;
    end

    % Ensure that T1 + T2 + T3 = 1
    if abs(T1 + T2 + T3 - 1) > 1e-6
        error('Tristimulus values do not sum to 1');
    end
end

function [inh] = inharmonicity(harmonics, amplitudes)
    % This function calculates the inharmonicity values
    % Inputs:
    %   frequencies - array of eigenfrequencies (f1, f2, ..., fN)
    %   fundamental_frequency - the fundamental frequency f1
    % Output:
    %   inharmonicity - array of inharmonicity values for each harmonic

    % Check if inputs are valid
    if ~isvector(harmonics) || ~isvector(amplitudes)
        error('Inputs must be a vector of frequencies and amplitudes');
    end

    k = [1,2,3,4,5]; %constants
    n_num = 1:1:8;

    %Calculate fc
    f1 = harmonics(1); %fundamental frequency for each note
    %k_th = k(jj); %kth harmonic
    Ak = amplitudes(1:5);  % amplitude of the kth harmonic
    DFk = harmonics(1:5) - k'.*f1; %Dfn
    numerator = sum(Ak.*(DFk./(k'*f1)));
    denominator = sum(Ak);
    fc = f1 * (1 + (numerator/denominator));

    % Calculate inharmonicity
    fn = harmonics;
    nfc = n_num*fc;

    % Calculate inharmonicity for each harmonic
    inh = (fn' -nfc ) ./ (n_num.*f1);
end

% CENTS Calculation
function cent = freq_cent(f)
    cent = 3986.*log10(f);
end







