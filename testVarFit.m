function [h p k ] = testVarFit( relNPStd, mSample, nSamples, m, v, iterations )
    
    %This tests the measured relative non poissson standard deviation
    %agains the model
    %relNPStd=the relative non poisson standard deviation of the data
    %being fitted
    %mSample= the mean of the sample values
    %nSamples=the number of samples in the origingal dataset (usually 2)
    %m=mean tested
    %v=variance tested
   
    
    [mu sigma]=getMuSigmaLognormal(m, v) ;
        
    varModel=lognrnd(mu,  sigma, size(mSample,1),1);
    varModel=(varModel.*mSample).^2; %Make var model the non normalized value    
        
    mu1 = log((mSample.^2)./sqrt(varModel+mSample.^2));
    sigma1 = sqrt(log(varModel./(mSample.^2)+1));
              
       
    for i=1:iterations
        %Make the samples by using the paramaters above
        %The number of samples should be equal to the
        %number in the original dataset        
        
        for j=1:nSamples
            s(:,j)=poissrnd(lognrnd(mu1, sigma1));
        end
   
        mSamplet=mean(s,2);
        
        relNPStdTest=sqrt((transpose(var(transpose(s)))-mSamplet))./mSamplet;  
        [h p k]=kstest2(relNPStd(mSample>100).^2, relNPStdTest(mSamplet>100).^2);
        
        hs(i)=h;
        ps(i)=p;
        ks(i)=k;
    end
    
    h=median(hs);
    p=median(ps);
    k=median(ks);
    
end
