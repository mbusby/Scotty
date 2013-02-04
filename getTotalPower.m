function [ totalPwr ] = getTotalPower( powerByReadDepth, readDepths, meanReadDepth, varianceReadDepth )
    %{
    Get the total power by read depth based on 
    the power at each read depth
    and the fraction of each datapoint that is expected to be found
    at that read depth
    Power varies both by read depth and gene variance.
    Variance by read depth is highest at values below 100.
    This function assumes all variances are the same when greater than
    the top read depth.  This is because the top read depth is usually 100
    and beyond this point Poisson variance is usually a trivial source of 
    variance at these high counts and power at levels greater than 100
    can be approximated by the power at 100.
    Assumes readDepths correspond to values in powerByReadDepth and 
    that readDepths is sorted 0->Inf.
    Assumes both samples sequenced to approximately the same read depth
    (that is, deep enough for a 2X fold change to appear)
    %}
    
    [mu sigma]=getMuSigmaLognormal(meanReadDepth, varianceReadDepth);   
    
    meanPowerByReadDepth=mean(powerByReadDepth,2);
    
    cumProbEach=logncdf(readDepths, mu, sigma);
    
    probEach=diff([0; cumProbEach;]);
    
    totalPwr0to100=sum(probEach.*meanPowerByReadDepth);
    
    finalPower=meanPowerByReadDepth(length(meanPowerByReadDepth));
    
    totalPwrGreaterThan100=(1-logncdf(max(readDepths), mu, sigma))*finalPower;   
    
    totalPwr=totalPwr0to100+totalPwrGreaterThan100;
    
    

end

