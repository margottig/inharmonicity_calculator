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