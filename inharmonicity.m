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
