function [ npStd controlSampleIds testSampleIds  ] = getDispersionAllComparisons( samples )
    %get the nonpoisson variance for all pairs of samples within an experiment
    
       
    npStd=0;
    ctr=1;
    
    samplesN=getNormalizedSamples(samples);   
    
    
    for i=1:size(samples,2)
        for j=i+1:(size(samples,2))
            s=[samplesN(:,i) samplesN(:,j)];            
            m=mean(s,2);
            v=transpose(var(transpose(s)))-m;
            npSt(:,1)=sqrt(v)./m;
            npStd(ctr,1)=sqrt(median(npSt(m>100).^2));            
            controlSampleIds(ctr)=i;
            testSampleIds(ctr)=j;
            ctr=ctr+1;
        end    
    end

    
end

