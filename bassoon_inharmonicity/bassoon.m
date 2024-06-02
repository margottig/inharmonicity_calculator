clear all;
close all;
nota_G2=load('nota_G2.mat')
nota_A3=load('nota_A3.mat')
nota_C2=load('nota_C2.mat')
nota_D3=load('nota_D3.mat')

A_ii_G2 =(nota_G2.out.ModArm(1:end,3));
F_real_G2=(nota_G2.out.ModArm(1:end,2));
A_ii_A3 =(nota_A3.out.ModArm(1:end,3));
F_real_A3=(nota_A3.out.ModArm(1:end,2));

A_ii_C2 =(nota_C2.out.ModArm(1:end,3));
F_real_C2=(nota_C2.out.ModArm(1:end,2));
A_ii_D3 =(nota_D3.out.ModArm(1:end,3));
F_real_D3=(nota_D3.out.ModArm(1:end,2));

A_matrix = [A_ii_A3 A_ii_D3 A_ii_G2 A_ii_C2]; %amplitude matrix

%% SPECTRUM PLOT
% figure;
% plot(F_real_A3,A_ii_A3)
% hold on
% plot(F_real_G2,A_ii_G2)
% plot(F_real_C2,A_ii_C2)
% plot(F_real_D3,A_ii_D3)
% grid on
% title('Spectrum of the notes')
% xlabel('Frequency (Hz)')
% ylabel('Amplitude')
% legend('A3','G2','C2','D3')

%% 2. CALCULATIONS
% D3 is a descending 5th in relation to A3,  - - 
% G2 is a 2 descending 5th in relation to A3, reduced to the original 8ve and  
% C2 is a 3 descending 5th in relation to A3, reduced to the original 8ve.  
f1 = 220;
fifth_Pt = 3/2;
fifth_ET = 2^(7/12);
%PYTHAGOREAN
D3_Pt = f1/fifth_Pt;
G2_Pt = f1/(fifth_Pt^2);
C2_Pt = f1/(fifth_Pt^3);
% TEMPERAMENTAL
C2_ET = f1/(fifth_ET^3);
G2_ET = f1/(fifth_ET^2);
D3_ET = f1/fifth_ET;
%% 3. DEVIATION
% Calculate deviations in Hz
D3_dev_Hz = abs(D3_Pt - D3_ET);
G2_dev_Hz = abs(G2_Pt - G2_ET);
C2_dev_Hz = abs(C2_Pt - C2_ET);

% Calculate deviations in cents
D3_dev_cents = abs(1200 * log2(D3_Pt / D3_ET));
G2_dev_cents = abs(1200 * log2(G2_Pt / G2_ET));
C2_dev_cents = abs(1200 * log2(C2_Pt / C2_ET));

% 3b Mean deviation
Fund_Freq = [F_real_A3(1) F_real_D3(1) F_real_G2(1) F_real_C2(1)];
Pt_Vec = [f1 D3_Pt G2_Pt C2_Pt];
ET_Vec = [f1 D3_ET G2_ET C2_ET];

% Pt - real deviation
real_dev_PT_Hz = abs(Pt_Vec - Fund_Freq);
real_dev_PT_ct = abs(3986*log10(Pt_Vec) - (3986*log10(Fund_Freq)));

% ET - real deviation
real_dev_ET_Hz = abs(ET_Vec - Fund_Freq);
real_dev_ET_ct = abs(3986*log10(ET_Vec) - (3986*log10(Fund_Freq)));

%% 3c Find Tuning System
if mean(real_dev_PT_ct) < mean(real_dev_ET_ct)
    tuningSystem = [1 2 (3/2)*2 2^2 ((3/2)^4) (3/2)*(2^2) (2^4)*((3/2)^-2) 2^3 (2^2)*((3/2)^2) (2)*((3/2)^4)];
else 
    tuningSystem = [1 2 2*2^(7/12) 2^2 2^(28/12) 2^(7/12)*2^2 2^(34/12) 2^3 2^(7/12)^2*2^2 2^(7/12)^4*2];
end

%% 4  inharmonicity of the first 10 harmonics
harm_matrix = zeros(length(tuningSystem),length(Fund_Freq));
inharm_matrix = harm_matrix;
delta_fn = zeros(length(tuningSystem),length(Fund_Freq));


% MATRIX FOR THE 4 NOTES
denominator = zeros(1,width(A_matrix));
numerator = zeros(1,width(A_matrix));

% Calculate theoretical harmonics from using the chosen tuningSystem
for n = 1:length(tuningSystem)
   harm_matrix(n,:) = (tuningSystem(n)).*Fund_Freq;
   harm_matrix(n,:) = 3986.*log10(harm_matrix(n,:)); % make to cents

    delta_fn(n,:) = harm_matrix(n,:)-n.*(3986.*log10(Fund_Freq));
    
end


for k = 1:5 % fixed in that formula
        denominator = denominator + A_matrix(k, :);
        numerator = numerator + ( A_matrix(k,:).*( delta_fn(k,:)./( k.*(3986.*log10(Fund_Freq)))));  
       
end

% Previous calculations
term = 1+(numerator./denominator);
fc = (3986.*log10(Fund_Freq)).*term;

% Calculation of inharmonicity
for ii=1:10 
        inharm_matrix(ii,:) = (harm_matrix(ii,:)-(ii.*fc))./(ii.*3986.*log10(Fund_Freq));
end


figure
plot(round(tuningSystem),inharm_matrix(:,1))
hold on
plot(round(tuningSystem),inharm_matrix(:,2))
plot(round(tuningSystem),inharm_matrix(:,3))
plot(round(tuningSystem),inharm_matrix(:,4))
legend('A3','D3','G2','C2')
title('Inharmonicity of each harmonic.')
xlabel('harmonic')
ylabel('inharmonicity (cents)')
grid on