% Which of the two instruments has a greater weight in the spectrum of even  harmonics?
% % Define the even harmonic ratio values for each condition
% heven_C1 = % Your calculated value for chordophone-note 1
% heven_A1 = % Your calculated value for aerophone-note 1
% heven_C2 = % Your calculated value for chordophone-note 2
% heven_A2 = % Your calculated value for aerophone-note 2

function spectrum_Weight(heven_C1,heven_A1,heven_C2,heven_A2)

    % Calculate the average even harmonic ratio for the chordophone
    heven_C_avg = (heven_C1 + heven_C2) / 2;
    
    % Calculate the average even harmonic ratio for the aerophone
    heven_A_avg = (heven_A1 + heven_A2) / 2;
    
    % Display the average even harmonic ratio values
    fprintf('Average even harmonic ratio for guitar: %f\n', heven_C_avg);
    fprintf('Average even harmonic ratio for flute: %f\n', heven_A_avg);
    
    % Determine which instrument has a greater weight in the spectrum of even harmonics
    if heven_C_avg > heven_A_avg
        fprintf('The guitar has a greater weight in the spectrum of even harmonics.\n');
    else
        fprintf('The flute has a greater weight in the spectrum of even harmonics.\n');
    end
    
    % Create a bar plot to visualize the average even harmonic ratios
    figure;
    categories = {'Guitar', 'Flute'};
    values = [heven_C_avg, heven_A_avg];
    
    bar(values);
    set(gca, 'XTickLabel', categories);
    ylabel('Average Even Harmonic Ratio');
    grid on;
    title('Comparison of Average Even Harmonic Ratios');
    
end