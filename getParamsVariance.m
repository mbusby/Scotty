function [ mVControl, vVControl, pVControl, corrControl,  mVTest, vVTest, pVTest, corrTest] = getParamsVariance(data, nControlSamples, nTestSamples )
    %{
    Gets the mean and variance for the underlying 
    distribution of variances for the 

    Uses the fitVariance function which assumes
    the true underlying variance is itself distributed in a 
    lognormal distribution having mean m and variance v.

    Assumes that one (either the test or the control) will be entered
    with more than one sample (measurable variance)
    This should be checked on the PHP form.
    If the Control or Test is not entered then the variance reported
    is from the one that is entered.

    Also tests for a correlation between the non-Poisson variance and the
    read depth.  This is used to warn users if there is a correlation,
    which is not expected in the statistical model and may cause errors 
    when variance is extrapolated.
    %}

    out='Running Sequence Depth Parameters'
    
    tic 
    if nControlSamples>1
        Ci=data(:,1:nControlSamples);
        [ mVControl, vVControl, pVControl , corrControl] = fitVariance( Ci );
    end
    toc
    
    if nTestSamples>1       
        Ti=data(:,nControlSamples+1:nControlSamples+nTestSamples);
        [ mVTest, vVTest, pVTest, corrTest] = fitVariance( Ti );
    else
        mVTest=mVControl;
        vVTest=vVControl;
        pVTest=pVControl;  
        corrTest=0;
    end

    
    if nControlSamples<2
        mVControl=mVTest;
        vVControl=vVTest;
        pVControl=pVTest;    
        corrControl=0;
    end
end

