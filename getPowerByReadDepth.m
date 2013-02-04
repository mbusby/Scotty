function [ powerByReadDepth ] = getPowerByReadDepth( mVControl, vVControl, mVTest, vVTest, seqDepths, nRepsControl, nRepsTest, fc, pCut, includePoisson )
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
    
    varSteps=transpose(0.02:0.02:.98);
    
    ods1All=logninv(varSteps, muVarianceC, sigmaVarianceC);
    ods2All=logninv(varSteps, muVarianceT, sigmaVarianceT);
    
    ods1All(isnan(ods1All)==1)=0;
    ods2All(isnan(ods2All)==1)=0; 
    
    ctr=1;      
    
    for i=1:length(ods1All)
                 
        mx1=seqDepths;
        mx2=mx1*fc;
        
        ods1=ods1All(i);
        ods2=ods2All(i);
            
        if includePoisson==1            
            vx1=mx1+((mx1).*ods1).^2;
            vx2=mx2+((mx2).*ods2).^2;
        else
            vx1=((mx1).*ods1).^2;
            vx2=((mx2).*ods2).^2;
        end 
        
        bad=(vx1+vx2==0);        
        vx1(bad)=1; %put in dummy variables so it won't error
        
        if mx1==0 & mx2==0
            pwr=zeros(size(vx1));        
        else
            [ pwr ] = getPowerTTest( mx1, vx1, mx2, vx2, nRepsControl, nRepsTest, pCut );
        end
                
        pwr(bad==1)=0; %Correct the dummy variables
        
        powerByReadDepth(i,:)=pwr;
    end
    


end

