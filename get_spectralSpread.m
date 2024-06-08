function sdVector = get_spectralSpread(espectro, centroidVector,f)
%% given a spectrogram, pre-calculated spectral centroid, and sample rate,
%% calculate the spectral standard deviation at each frame

% initialize
sdVector = zeros(size(centroidVector));
for i = 1:length(sdVector)
    sdVector(i) = spectralSD(espectro(:,i), centroidVector(i), f);
end

end


function spectralSDScalar = spectralSD(frame, centroidFrame, f)
% calculates the standard deviation of the spectrum around the spectral centroid. 
% High values indicate a rich spectrum.
    numerator = 0;
    denominator = 0;
    
    for i = 1:length(frame)
        thisBin = frame(i);
        numerator = numerator + (((f(i) - centroidFrame)^2) * thisBin);
       denominator = denominator + thisBin;
    end
    spectralSDScalar = sqrt(numerator)/denominator;
end


