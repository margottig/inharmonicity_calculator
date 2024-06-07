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