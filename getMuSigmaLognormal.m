function [ mu sigma ] = getMuSigmaLognormal( m, v )
        
    mu = log((m.^2)./sqrt(v+m.^2));
    sigma = sqrt(log(v./(m.^2)+1));    
        
end

