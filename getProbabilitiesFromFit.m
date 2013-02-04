function [  probSequenced observedGenes totalGenesExpressed m v] = getProbabilitiesFromFit( data, totalReads )
    
    %{
    This function returns the probability of genes being captured per read
`   sequenced as well as the total number expressed.

    The probablility that each gene will be expressed is estimated as
    the number of reads found for the gene divided by the total number
    of genes.

    However, not all expressed genes will be observed if the pilot data
    is not sequenced deeply enough.  The probability for these sequences
    are filled in by fitting the data to a zero truncated Poisson-
    Lognormal distribution using maximum likelihood estimation.

    data = read counts per gene for the pilot data data
    totalReads = the number of reads in the pilot data
    %}
  
    %seed the random number generators so we aren't all over the place,
    %people    
    s = RandStream('mt19937ar','seed',4503);
    RandStream.setDefaultStream(s);       
    
    poissLogNormal_cdf = @(x,mu, sigma) zeroTruncatedPoissonLognormalCDF(x,mu, sigma);
    
    dataX=data(data>0);
    [params,lambdaCI] = mle(dataX, 'cdf',poissLogNormal_cdf, 'start',[0 var(data)], 'lower',[0 0]);

    [mu sigma]=getMuSigmaLognormal(  params(1), params(2)^2);
   
        
    %how many zeros are there? estimate fraction by bootstrapping - this 
    %workds better than the cdf
    observedGenes=sum(data>0);
    x=poissrnd(random('lognormal', mu, sigma, 100000,1));
    unobservedGenePortion=sum(x==0)/100000;
    totalGenesExpressed=round(observedGenes/(1-unobservedGenePortion));  
            
    modelDist=lognrnd(mu, sigma, totalGenesExpressed,1);

    if length(data)>length(modelDist)
        modelDist=[zeros(length(data)-length(modelDist), 1); modelDist];
    end
    
    modelDist=lognrnd(mu, sigma, size(data,1),1);

    modelProb=modelDist/totalReads;
    modelProb=sort(modelProb, 'descend');
   
    modelProb=modelProb(1:size(data,1));
    
    modelProb=sort(modelProb, 'ascend');
    probSequenced=sort(data/totalReads, 'ascend');
    
    probSequenced(1:sum(data<=1))=modelProb(1:sum(data<=1));         
  
    
end