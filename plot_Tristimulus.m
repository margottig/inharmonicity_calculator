function plot_Tristimulus(Ts)
    % Define the Tristimulus values for each condition
    % T1_C1 = % Your calculated value for chordophone-note 1
    % T1_A1 = % Your calculated value for aerophone-note 1
    % T1_C2 = % Your calculated value for chordophone-note 2
    % T1_A2 = % Your calculated value for aerophone-note 2
    % 
    % T2_C1 = % Your calculated value for chordophone-note 1
    % T2_A1 = % Your calculated value for aerophone-note 1
    % T2_C2 = % Your calculated value for chordophone-note 2
    % T2_A2 = % Your calculated value for aerophone-note 2
    % 
    % T3_C1 = % Your calculated value for chordophone-note 1
    % T3_A1 = % Your calculated value for aerophone-note 1
    % T3_C2 = % Your calculated value for chordophone-note 2
    % T3_A2 = % Your calculated value for aerophone-note 2
    

    % Ts should be a 1x4 array where:
    % Ts(1) = T1_C1, Ts(2) = T1_A1, Ts(3) = T1_C2, Ts(4) = T1_A2, 
    % Ts(5) = T2_C1, Ts(6) = T2_A1, Ts(7) = T2_C2, Ts(8) = T2_A2, 
    % Ts(9) = T3_C1, Ts(10) = T3_A1, Ts(11) = T3_C2, Ts(12) = T3_A2, 
    
     % Extract Tristimulus values from the input array
    T1_C1 = Ts(1); T1_A1 = Ts(2); T1_C2 = Ts(3); T1_A2 = Ts(4);
    T2_C1 = Ts(5); T2_A1 = Ts(6); T2_C2 = Ts(7); T2_A2 = Ts(8);
    T3_C1 = Ts(9); T3_A1 = Ts(10); T3_C2 = Ts(11); T3_A2 = Ts(12);
    
    % Combine values into matrices for bar plot
    T1_values = [T1_C1, T1_A1, T1_C2, T1_A2];
    T2_values = [T2_C1, T2_A1, T2_C2, T2_A2];
    T3_values = [T3_C1, T3_A1, T3_C2, T3_A2];
    
    % Combine all T values for stacking
    all_values = [T1_values; T2_values; T3_values];
    
    % Create the bar plot
    figure;
    bar(all_values', 'stacked');
    
    % Set the x-axis labels
    set(gca, 'XTickLabel', {'Guitar A5', 'Flute A5', 'Guitar C5', 'Flute C5'});
    
    % Add title and labels
    title('Tristimulus Values for Different Notes');
    xlabel('Notes');
    ylabel('Tristimulus Values');
    
    % Add legend
    legend('T1', 'T2', 'T3');
    grid on;
end