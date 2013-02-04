function [ readDepthsSampled ] = getReadDepthsRun( probSequenced, nReads, nSteps , nRepsOriginal)
%Gets evenly spaced readDepths
    
    readDepths=probSequenced.*nReads;
    readDepths=readDepths(readDepths>(1/(nRepsOriginal)) );
    
    nPoints=size(readDepths,1);
    
    increment=nPoints/nSteps;
    
    minInc=1;    
    maxInc=nPoints;    
    
    vals=minInc:increment:maxInc;
    
    %{
    Wiggle it, just a little bit
    
    For datasets with many low count reads the probabitiles might be the same
    at various points.  This wiggles the points a little so the 
    plot points don't exactly overlap.  The plot is in log space so 
    it won't matter much at higher read depths.
    
    %}
    
    for i=1:length(vals)
       
       f=floor(vals(i));
       c=ceil(vals(i));
       r=vals(i)-f;
       
       c(c>nPoints)=nPoints;
       f(f==0)=1;       
       thisVal=readDepths(f)+r*(readDepths(c)-readDepths(f));     
       thisVal(thisVal<0)=0;
       readDepthsSampled(i)=thisVal;
    end
    
    %makes sure plot goes from 0 to maximum read depth
    readDepthsSampled(i+1)=max(readDepths);
    
    max(readDepths)
    
    readDepthsSampled=sort(readDepthsSampled);
    
end

