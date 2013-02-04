function [i] = getCumProbLocation(cumProbSequenced, target)
    %cumProbSequenced is a vector sorted in descending order
    %cumulative so all values are unique
    %finds the smallest cumProbSequenced that is greater than x
    
    %Starts searching at top
    i=1;
    
    %while cumProbSequenced(i)>target       
    %    i=i+1;
    %end    
    %
    %binary search for closest
    
    outVal=floor(length(cumProbSequenced)/2);
    
    l=0;
    r=length(cumProbSequenced);
    
    while l<r-1        
        t=cumProbSequenced(outVal);
        if t>target           
            r=outVal;
            outVal=r-floor((r-l)/2);            
        elseif t<target        
            l=outVal;
            outVal=l+floor((r-l)/2);
        elseif t==target
            r=l;              
        end        
        
    end 
     
    if target<cumProbSequenced(outVal) 
       outVal=outVal-1; 
    end
    
    i=outVal;
    
    
end