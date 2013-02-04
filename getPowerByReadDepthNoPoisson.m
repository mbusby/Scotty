function [ powerByReadDepth readDepths] = getPowerByReadDepthNoPoisson( mVControl, vVControl, mVTest, vVTest, mReadDepth, vReadDepth, nRepsControl, nRepsTest, fc, pCut )
  %{
    Gets the power for a dataset given a cert

    Creates a table with 99 columns and 99 rows.
    Power is based on sequencing depth (poisson noise)
    And variance is drwawn 
    
    Each column represents 1/100th of the distribution of variances 
    (i.e. the variances at the 0.05, 0.15, 0.25...0.95 point)
    Each row represents 1/100th of the sequencing depth of the control 
    sample (1:100).  Each matrix is decided on by using the inverse 
    function so that each block represents 1/9801 of the likely values.
    Total power can thus be calculated based on a straight mean of the
    resulting powerByReadDepth matrix.    

    mVControl, vVControl, mVTest, vVTest=the mean and variance of the
    lognormal distributions that describe the variance
    meanSeqDepth, varSeqDepth =the mean and variance of the reads per gene    
    
   %}

    [muVarianceC sigmaVarianceC]=getMuSigmaLognormal(mVControl, vVControl);    
    [muVarianceT sigmaVarianceT]=getMuSigmaLognormal(mVTest, vVTest);
    [muReadDepth sigmaReadDepth]=getMuSigmaLognormal(mReadDepth, vReadDepth);
    
    varSteps=transpose(0.02:0.02:.98);
    
    ods1=logninv(varSteps, muVarianceC, sigmaVarianceC);
    ods2=logninv(varSteps, muVarianceT, sigmaVarianceT);
    
    ods1(isnan(ods1)==1)=0;
    ods2(isnan(ods2)==1)=0; 
    
    readDepths=transpose(logninv(0.02:0.02:.98, muReadDepth, sigmaReadDepth));
   
    for i=1:length(readDepths)
                   
        mx1=readDepths(i);
        mx2=mx1*fc;
            
        vx1=((mx1).*ods1).^2;
        vx2=((mx2).*ods2).^2;
            
        [ pwr ] = getPowerTTestLogged( mx1, vx1, mx2, vx2, nRepsControl, nRepsTest, pCut );
        powerByReadDepth(i,:)=transpose(pwr);
    end
    


end

