function [ nReps ] = getOptimizationByGoals( mVControl, vVControl, mVTest, vVTest, fc, pCut, percFound, meanReadDepthControl, varianceReadDepthControl )
    %{
    This optimizes the experiment by the experimental goals set by the user
    
    %}

    nReps=2;
    
    pwr=0;
    
    mVControl
    mVTest
    
    while pwr<percFound
        [ powerByReadDepth readDepths]=getPowerByReadDepth( mVControl, vVControl, mVTest, vVTest, nReps, nReps, fc, pCut );
        pwr=getTotalPower( powerByReadDepth, readDepths, meanReadDepthControl, varianceReadDepthControl );
        nReps=nReps+1;        
    end
    

end

