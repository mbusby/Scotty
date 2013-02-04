function [observedGenes totalGenesExpressed mFinal vFinal readsFinal lognFitPlot probSequenced] = getSequenceDepthParameters( data, totalReads, dataTitle)
    %Get the probability of being sequenced by brute force
    %same as  getProbSequencedMonteCarlo but uses the user defined number of reads
    %instead of the target genes
   
    %s = RandStream('mt19937ar','seed',4503);
    %RandStream.setDefaultStream(s);
    
    [  probSequenced observedGenes totalGenesExpressed ] = getProbabilitiesFromFit( data, totalReads );
    
    probSequenced=sort(probSequenced, 'ascend');
    
    cumProbSequenced=[1;  cumsum(probSequenced); 0]; %0 is not expressed, 1 is not aligned
    cumProbSequenced=sort(cumProbSequenced, 'ascend');
    
    geneCounts=zeros(length(cumProbSequenced),1);
    
    %Set up variables for loop
    readsRequired=0;
    readsSequenced=0;

    readsUsed=0;
    totalReadsUsed=0;
    topSpot=length(cumProbSequenced);        
    genesFoundNow=0;
    ctr=1;%for saving the genes found
    ctr2=1;%For saving how many reads were modelled
      
    top=max(cumsum(probSequenced));   
  
    breakPt=2000; %does a thousand reads at a time.  Improves performance
                  %and gives answer within 1000 reads which is 
                  %only introduces a trivial error because the total
                  %read scale is in the millions
                  %usually Reads required
    randVars=random('unif', 0,top,breakPt,1);
    
    %Don't have to keep checking genes that are already found
    %top can get smaller with each iteration
    %and reads can be allocated at the other
    cumProbSequencedUsed=cumProbSequenced;
    
    readsFinal=10^5;
    
    while readsUsed<readsFinal
            x=randVars(ctr2);                          
            place=getCumProbLocation(cumProbSequencedUsed,x);
            geneCounts(place)=geneCounts(place)+1;
            readsUsed=readsUsed+1;
            ctr2=ctr2+1;
            if mod(readsUsed,breakPt)==0  
               readsSequenced(ctr,1)=readsUsed;
               meanExp(ctr,1)=mean(geneCounts);
               stdExp(ctr, 1)=std(geneCounts);
               readsSequenced(ctr,1)=readsUsed;
               ctr=ctr+1;                
               ctr2=1; 
               randVars=random('unif', 0,top,breakPt,1);
            end
    end
    
    
    p = polyfit( readsSequenced,meanExp,1);
    
    mFinal=p(1)*totalReads;
    
    p = polyfit( readsSequenced,stdExp,1);
    vFinal=(p(1)*totalReads)^2;    
    
    [mu sigma]=getMuSigmaLognormal(mFinal, vFinal);%mu and sigma are the mean and var of the corresponding lognormal
    
    lb=logninv(0.05, mu, sigma);
    lb(lb<10)=1;%just start from one if it's less than 10
    ub=logninv(0.95, mu, sigma);
    
 	lognFitPlot=figure();
    [h x]=hist(data(data>0),lb:ub+1);
    h=h/sum(h);
    bar(x(x<ub),h(x<ub));
    hold on  
    plot(x(x<ub),lognpdf(x(x<ub), mu, sigma),'r');
    titleText=strcat({'Fit of the '}, dataTitle);
    legend('Data', 'Model')
    title(titleText);

    
    
end

