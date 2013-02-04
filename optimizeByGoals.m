function [ nReps ] = optimizeByGoals( mVControl, vVControl, mVTest, vVTest, fc, fp, percFound )
    %{
    This optimizes the experiment by the experimental goals
    %}

    nReps=2;
    
    while pwr<percFound
        [ powerByReadDepth readDepths]=getPowerByReadDepth( mVControl, vVControl, mVTest, vVTest, nReps, nReps, fc, fp );
        pwr=getTotalPower( powerByReadDepth, readDepths, meanReadDepth, varianceReadDepth );
        nReps=nReps+1;
        
    end
    

end

