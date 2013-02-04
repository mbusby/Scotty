function [ samplesN ] = getNormalizedSamples( samples )
    %Normalizes samples to a common sequencing depth based on a
    %comparison of the median values in both samples
    
    x=median(samples,2);
    samplesN=samples;
    
    for j=1:size(samples,2)
        samplesN(:,j)=samples(:,j)./getSlope(x,samples(:,j));
    end

end