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