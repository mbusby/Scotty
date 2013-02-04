function [ slope ] = getSlope( Ci, Ti )
    %This gets the slope from two vectors of reads per gene
    %based on normalizing them to a common median value
    
    total=Ci+Ti;   
    CiNonZero=Ci(total>0);      
    TiNonZero=Ti(total>0); 
    
    slope=median(TiNonZero./CiNonZero);
        
end

