function flatnessVector = get_spectralFlatness(espectrograma)
% calculates the noisiness of the spectrum.
% High values indicate a flat (and therefore possibly noisy) spectrum which often sounds "harsh”.
flatnessVector = zeros(1, width(espectrograma));

for i = 1:length(flatnessVector)
    this_window = espectrograma(:,i);
    %geomean requires Statistics and Machine Learning Toolbox 
    flatnessVector(i) = geomean(this_window) / mean(this_window);
end
