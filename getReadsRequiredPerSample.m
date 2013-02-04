	 
function [  probSequenced totalGenesExpressed] = getProbabilitiesFromFit( dataModelled, readsModelled )
    
    %{
    This function returns the probability of genes being captured per read
`   sequenced as well as the total number expressed.

    The probablility that each gene will

    However, not all expressed genes will be observed if the pilot data
    is not sequenced deeply enough.  The probability for these sequences
    are filled in by fitting the sample to a truncated lognormal poisson
    distribution.  
    %}
  
    %seed the random number generators so we aren't all over the place,
    %people    
    s = RandStream('mt19937ar','seed',4503);
    RandStream.setDefaultStream(s);
       
    
    poissLogNormal_cdf = @(x,mu, sigma) zeroTruncatedPoissonLognormalCDF(x,mu, sigma);

    [params,lambdaCI] = mle(dataModelled, 'cdf',poissLogNormal_cdf, 'start',[0 var(dataModelled)], 'lower',[0 0]);

    m=params(1);
    s=params(2);
    v=s^2;
    s=sqrt(v);
    
    [mu sigma]=getMuSigmaLognormal(m, v);
    
    %how many zeros are there?    
    %estimate by bootstrapping, as the CDF doesn't work so good for 0
    observedGenes=sum(dataModelled>0);
    x=poissrnd(random('lognormal', mu, sigma, 10000,1));
    unobservedGenePortion=sum(x==0)/10000;
    totalGenesExpressed=round(observedGenes/(1-unobservedGenePortion));  
    
   %{
    It is possible that you may not be able to measure nTarget genes at
    the cutOff depth becasue that number is higher than the number of
    expressed genes.
    
    In these cases, read depths are calculated for the maximum number of
    genes that can be detected and an informational  method (error)
    is displayed for the user.    
    
    In cases where the number of genes that the user would like to detect
    is very close to the total number of available genes, a great high
    number of reads will be required to calculate all genes.
    
    To inform user that this may not be practical, there is a second message
    when the number requested is >90% of the total predicted to be
    expressed.  The user can then examine the rarefaction plots to
    determine if this is a practical goal.    
    %}
   
    if totalGenesExpressed<nTarget
        nTarget=totalGenesExpressed;
    end    
            
    modelDist=lognrnd(mu, sigma, totalGenesExpressed,1);

    if length(dataModelled)>length(modelDist)
        modelDist=[zeros(length(dataModelled)-length(modelDist), 1); modelDist];
    end

    modelProb=modelDist/readsModelled;
    modelProb=sort(modelProb, 'descend');
    modelProb=modelProb(1:size(dataModelled,1));
    
    modelProb=sort(modelProb, 'ascend');
    probSequenced=sort(dataModelled/readsModelled, 'ascend');
    
    probSequenced(1:sum(dataModelled<=10))=modelProb(1:sum(dataModelled<=10));    
    
    %{
    Calculate how many reads do you need at the end
    So if you want to sequence nTarget genes, sort the
    genes by probaibilty of sequencing. The final
    gene you will need to sequence is the one
    that is nTarget from the most frequently sequenced gene.
    r.g.
    You have a teeny tiny organism with only four genes.
    They each have the following probabilities of being sequenced:
    0.50 0.40 0.09 0.01
    
    To identify 3 genes by 5 reads you need to sequence
    5/0.09 = 55.6 = 60 reads
    
    %}
    len=length(probSequenced);
   
    nthGenesProbSequenced=probSequenced(len-nTarget);    
    
    readsRequired=ceil(cutOff/nthGenesProbSequenced);
    
    

end