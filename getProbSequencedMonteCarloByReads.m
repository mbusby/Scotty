function [ readsRequired observedGenes totalGenesExpressed readsSequenced genesFound finalMeanSample finalVarSample] = getProbSequencedMonteCarloByReads( data, totalReads, cutOff , targetReads )
    %Get the probability of being sequenced by brute force
    %same as  getProbSequencedMonteCarlo but uses the user defined number of reads
    %instead of the target genes
   
    s = RandStream('mt19937ar','seed',4503);
    RandStream.setDefaultStream(s);
    
    [  probSequenced observedGenes totalGenesExpressed ] = getProbabilitiesFromFit( data, totalReads );
    
    probSequenced=sort(probSequenced, 'ascend');
    
    cumProbSequenced=[1;  cumsum(probSequenced); 0]; %0 is not expressed, 1 is not aligned
    cumProbSequenced=sort(cumProbSequenced, 'ascend');
    
    geneCounts=zeros(length(cumProbSequenced)-2,1);
    
    %Set up variables for loop
    readsRequired=0;
    readsSequenced=0;
    genesFound=0;

    readsUsed=0;
    totalReadsUsed=0;
    topSpot=length(geneCounts);       
    genesFoundNow=0;
    ctr=1;%for saving the genes found
    ctr2=1;%For saving how many reads were modelled
      
    top=max(cumsum(probSequenced));   
 
    breakPt=1000; %does a hundred thousand reads at a time.  Improves performance
                  %and gives answer within 1000000 reads which is 
                  %only introduces a trivial error because the total
                  %read scale is in the millions
                  %usually Reads required
    randVars=random('unif', 0,top,breakPt,1);
    randVars(randVars==1)=0.99999;
    %Don't have to keep checking genes that are already found
    %top can get smaller with each iteration
    %and reads can be allocated at the other
    cumProbSequencedUsed=cumProbSequenced;
    lg=length(geneCounts);
    while readsRequired<targetReads && top>0
            x=randVars(ctr2);            
            place=getCumProbLocation(cumProbSequencedUsed,x);
            place(place>lg)=lg;
            geneCounts(place)=geneCounts(place)+1;
            readsUsed=readsUsed+1;
            ctr2=ctr2+1;
            
            if readsUsed==breakPt                
                totalReadsUsed=totalReadsUsed+readsUsed*(1/top);
                genesFoundNow=sum(geneCounts>=cutOff);
                readsSequenced(ctr,1)=totalReadsUsed;
                genesFound(ctr,1)=genesFoundNow;
                ctr=ctr+1;                
                ctr2=1;
                readsUsed=0;                  
                %This finds the final value where geneCounts is bigger than
                %the cutOff
                while topSpot>1 && geneCounts(topSpot-1)>cutOff+10
                    topSpot=topSpot-1;
                end        
                top=cumProbSequenced(topSpot);    
                cumProbSequencedUsed=cumProbSequencedUsed(1:topSpot+1);
                cumProbSequencedUsed=[cumProbSequencedUsed; 1;];            
                randVars=random('unif', 0,top,breakPt,1);     
                randVars(randVars==1)=0.99999;
                readsRequired=totalReadsUsed;
                
            end
    end
    
    simFinalData=poissrnd(readsRequired.*probSequenced);
    finalMeanSample=mean(simFinalData);
    finalVarSample=var(simFinalData);
end

