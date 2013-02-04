function [nReps optReadsControl optReadsTest insuffFunds ]=getOptimizationByBudget( nControlSamples, nTestSamples,  ...
   readsRequired, costPerRepControl, costPerRepTest, costPerMillionReads, totalBudget)
    
    %{
    This function finds the maximum number of replicates for the test
    and control conditions that can be obtained given the budget
    constraints.  

    %}

    mReadsRequiredC=mean(readsRequired(1:nControlSamples));
    mReadsRequiredT=mean(readsRequired(nControlSamples+1:nControlSamples+nTestSamples));
    costPerRead=costPerMillionReads/(10^6);
    insuffFunds=0;

    %Is it possible to have at least three replicates of each condition?
    nReps=0;
    totalCost=0;
    while totalCost<totalBudget
        nReps=nReps+1;       
        costControl=(nReps*costPerRepControl)+(nReps*mReadsRequiredC*costPerRead);
        costTest=(nReps*costPerRepTest)+(nReps*mReadsRequiredT*costPerRead);    
        totalCost=costControl+costTest;
    end

    if nReps<3
        insuffFunds=1;
        nReps=3;        
    end
            
    budgetReps=(nReps*costPerRepControl)+(nReps*costPerRepTest);
    budgetRemaining=totalBudget-budgetReps ;
        
    allocationRatio=mReadsRequiredC/(mReadsRequiredC+mReadsRequiredT);     
    budgetReadsControl=allocationRatio*budgetRemaining;
    optReadsControl=(budgetReadsControl/nReps)/costPerRead ;  
        
    budgetReadsTest=(1/allocationRatio)*budgetRemaining;
    optReadsTest=(budgetReadsTest/nReps)/costPerRead;       
      

end

