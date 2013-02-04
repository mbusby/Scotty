function [powerPlot, excludedPlot, biasPlot, costPlot, allowedPlot, cheapestExperiment, mostPowerfulExperiment, experimentsEvaluated, pwrsCalc] = getOptimizationCharts(...
    maxReps, probSequencedC, mVControl, vVControl,  mVTest, vVTest, fc, pCut, costPerRepControl, costPerRepTest, ...
    costPerMillionReads, totalBudget, minReadsPerRep, maxReadsPerRep, minPercUnBiasedGenes, minPercDetected, pwrBiasCutoff, nRepsOriginal )

    %{
    Gets the optimization charts
    The program creates an optimization table that has the power
    based on the variance and number of sequencing depths.   
    The probability of being sequenced is based on
    taking a random sample of 250 probailities.  This gives a
    sample at a fairly evenly distributed interval without any
    assumptions about the underlying distribution.
    Because the table is representative of the total experiment the
    mean of that table gives the power for the total experiment.
    %}
    
    out='Running optimization charts.'
    
    costPerRead=costPerMillionReads/10^6;
    nReps=2:maxReps;
    nReadsAllowedPerRep=Inf;
   
    
    if minReadsPerRep>maxReadsPerRep
        minReadsPerRep=maxReadsPerRep/10;
    end
    
    inc=(maxReadsPerRep-minReadsPerRep)/9;

    readsS=round(minReadsPerRep:inc:maxReadsPerRep);
    
    %Set up dummy 
    pwrsCalcA=zeros(length(nReps), length(readsS));
    pwrBiasesA=zeros(length(nReps), length(readsS));  
    
    %Sample the probabilities of being sequenced 
    %These are used to calculate a probability randomly distributed
    %
    
    seqProbs=randsample(probSequencedC,250); %same prob for each rep, perpetuates random sampling effects but prevents inconsistent results.
        
    for j=1:length(nReps)
        nReplicates=nReps(j);        
        nRepsControl=nReplicates;
        nRepsTest=nReplicates;   
        
        parfor k=1:length(readsS)
            
            readsPerRep=readsS(k);
            
            seqDepths=seqProbs.*readsPerRep; %The random sequencing depths to test            
            
            seqDepths=seqDepths(seqDepths>(1/(nRepsOriginal)) ); %Only include seqDepths that would have been observed in the original dataset
                                
            [ powerByReadDepth] = getPowerByReadDepth( mVControl, vVControl,  mVTest, vVTest, seqDepths, nRepsControl, nRepsTest, fc, pCut, 1 );             
       
            [ powerByReadDepthNoPoisson ] = getPowerByReadDepth( mVControl, vVControl,  mVTest, vVTest,  seqDepths, nRepsControl, nRepsTest, fc, pCut, 0 );          
            
            [ totalPwr ] = mean(mean(powerByReadDepth));
            
            pwrBias=getPwrBias( powerByReadDepth,powerByReadDepthNoPoisson, pwrBiasCutoff  );  
            pwrsCalcA(j,k)=totalPwr;
            pwrBiasesA(j, k)=pwrBias;             
            
        end       
    end
    
    
    pwrsCalc=pwrsCalcA*100;    
    pwrBiases=pwrBiasesA*100;
       
    
    %Create cost matrix
   
    for i=1:length(nReps)
            for j=1:length(readsS)
                costMatrix(i,j)=(costPerRepControl*nReps(i))+(costPerRepTest*nReps(i)) + readsS(j)*nReps(i)*2*costPerRead;                
           end
    end
    
    readsSX=readsS/10^6;    
    
    %======================================================
    %Make exclusion plot
    %======================================================
  
    
   	excludedPlot=figure;    
    optimizationExclusions=ones(size(pwrsCalc));  
    imagesc(readsSX,nReps,pwrsCalc)
    colorbar
    xlabel('Millions of Reads Aligned to Genes Per Replicate', 'FontSize', 14);
    ylabel('Number Reps in Each Condition', 'FontSize', 14);
    zlabel('Power')
    set(gca,'XTick',readsSX)
    colormap bone
    hold on
    
    %Add Optimization by total budget
   % plot(nReadsAllowedPerRep/10^6, nReps, 'Color', [0 .5 0], 'LineWidth', 3, 'HandleVisibility', 'Off')  
   
    %Add PowerBias
    plot(  -1, -1, 's', 'MarkerFaceColor', [ 0.9137    0.5686    0.3647] , 'MarkerEdgeColor', [ 0.9137    0.5686    0.3647] , 'MarkerSize', 12);    
    if minPercUnBiasedGenes<Inf
       failedPwr=0;      
       for i=1:length(nReps)
           for j=1:length(readsS)
                if pwrBiases(i,j)<minPercUnBiasedGenes %if too few genes are measured to good enough power                
                   plot(  readsSX(j),nReps(i), 's', 'MarkerFaceColor', [ 0.9137    0.5686    0.3647],  'MarkerEdgeColor', [ 0.9137    0.5686    0.3647] ,'MarkerSize', 12, 'HandleVisibility', 'Off');     
                   optimizationExclusions(i,j)=0;
                end
           end
       end
    end
   
    
    %Add BudgetConstraint
     plot(  -1, -1, '^','MarkerFaceColor', [0 .5 0], 'MarkerEdgeColor', [0 .5 0], 'MarkerSize', 10);
    if totalBudget<Inf
       for i=1:length(nReps)
            for j=1:length(readsS)            
                if costMatrix(i,j)>totalBudget                   
                  plot(readsSX(j),nReps(i), '^', 'MarkerFaceColor', [0 .5 0],'MarkerEdgeColor', [0 .5 0], 'MarkerSize', 10,'HandleVisibility', 'Off');  
                  optimizationExclusions(i,j)=0;
                end
           end
       end
    end
    
    %Add Min Percent Detected
    plot(  -1, -1, 'o', 'MarkerFaceColor', 'r','MarkerEdgeColor', 'k'); 
    
    if minPercDetected>0    
       for i=1:length(nReps)
           for j=1:length(readsS)
               if pwrsCalc(i,j)<=minPercDetected %if power isn't high enough
                    plot(  readsSX(j),nReps(i), 'o', 'MarkerFaceColor', 'r','MarkerEdgeColor', 'k','HandleVisibility', 'Off');   
                    optimizationExclusions(i,j)=0;
               end
           end
       end
    end
    
    titleText=strcat({'% of Genes with a '}, num2str(fc), 'x Fold Change Detected (p<', num2str(pCut), ')');
    title(titleText, 'FontSize', 14)
    legend( 'Measurements Too Biased','Too Expensive', 'Insufficient Power', 'Location', 'SouthOutside')
    set(gca,'YDir', 'normal')        
    
    
  	powerPlot=figure;
    imagesc(readsSX,nReps,pwrsCalc)
    xlabel('Millions of Reads Aligned to Genes Per Replicate', 'FontSize', 14);
    ylabel('Number Reps in Each Condition', 'FontSize', 14);
    zlabel('Power')
    colormap bone;
    set(gca,'YDir', 'normal')    
    set(gca,'XTick',readsSX)
    titleText=strcat({'% Genes with a '}, num2str(fc), 'x Fold Change Detected (p<', num2str(pCut), ')');
    title(titleText, 'FontSize', 14)
    colorbar
    
    
    
    biasPlot=figure;
    imagesc(readsSX,nReps,pwrBiases)
    xlabel('Millions of Reads Aligned to Genes Per Replicate', 'FontSize', 14);
    ylabel('Number Reps in Each Condition', 'FontSize', 14);
    zlabel('Detection Bias')
    colormap copper
    set(gca,'YDir', 'normal')    
    set(gca,'XTick',readsSX)
    titleText=strcat({'% of Genes with at Least '}, num2str(pwrBiasCutoff*100),'% of Maximum Power');
    title(titleText, 'FontSize', 14)
    colorbar

     
    costPlot=figure;
    imagesc(readsSX,nReps,costMatrix)
    hold on
    xlabel('Millions of Reads Aligned to Genes Per Replicate', 'FontSize', 14);
    ylabel('Number Reps in Each Condition', 'FontSize', 14);
    zlabel('Cost')
    colormap(moneyscale)
    set(gca,'YDir', 'normal')    
    set(gca,'XTick',readsSX)
    title('Cost of Each Experimental Configuration', 'FontSize', 14)   
    cb=colorbar;
    
    %=======================================================
    %Get final optimization chart
    %=======================================================
    
    %Get best power
    pwrsCalcSize=size(pwrsCalc)
    
    bestPwr=max(pwrsCalc(optimizationExclusions==1));
   
    if isempty(bestPwr)==0
        [rp rd]=find(pwrsCalc==bestPwr);
        bestReadDepth=readsS(rd(1));
        bestNReps=nReps(rp(1));    
    else
        bestReadDepth=[];
        bestNReps=[];
    end
   
    
    mostPowerfulExperiment=[];
    if isempty(bestReadDepth)==0
        mostPowerfulExperiment(1)=bestReadDepth;
        mostPowerfulExperiment(2)=bestNReps;
    end
    
    %Get cheapestAllowed
    bestPwr=min(costMatrix(optimizationExclusions==1));
   
    if isempty(bestPwr)==0
        [rp rd]=find(costMatrix==bestPwr);
        bestReadDepth=readsS(rd(1));
        bestNReps=nReps(rp(1));    
    else
        bestReadDepth=[];
        bestNReps=[];
    end
        
    cheapestExperiment=[];
   	if isempty(bestReadDepth)==0
        cheapestExperiment(1)=bestReadDepth;
        cheapestExperiment(2)=bestNReps;
    end

    %Plot final allowed optimization
    
    allowedPlot=figure;
    
       
    imagesc(readsSX,nReps,optimizationExclusions)    
    hold on
    if length(cheapestExperiment)==2 && length(mostPowerfulExperiment)==2 
        plot(readsSX(readsS==cheapestExperiment(1)), cheapestExperiment(2), '^', 'MarkerFaceColor',[0 .5 0],'MarkerEdgeColor',[0 .5 0], 'MarkerSize', 12 )
        plot(readsSX(readsS==mostPowerfulExperiment(1)), mostPowerfulExperiment(2), 'o', 'MarkerFaceColor', [0 0 .9], 'MarkerEdgeColor',[0 0 .9], 'MarkerSize', 12)    
        legend('Cheapest Allowed Experiment', 'Most Powerful Allowed Experiment', 'Location', 'SouthOutside' );
    end
    xlabel('Millions of Reads Aligned to Genes Per Replicate', 'FontSize', 14);
    ylabel('Number Reps in Each Condition', 'FontSize', 14);
    zlabel('Cost')
    colormap([ 1 0 0; 1 1 1])  
    if sum(sum(optimizationExclusions))==0
         colormap([1 0 0]) ;%Paint it black if nothing is allowed
    end
    
    minX=min(readsSX)-.5*inc;
    maxX=max(readsSX)+.5*inc;
    minY=min(nReps)-.5;
    maxY=max(nReps)+.5;
    
    for i=1:length(readsSX)-1
        x=(readsSX(i)+readsSX(i+1))/2;        
        plot([x x], [minY maxY],'Color',[0.8 0.8 0.8]);
    end
    
    for i=1:length(nReps)-1
        y=nReps(i)+.5;        
        plot([minX maxX] ,[y y],'Color',[0.8 0.8 0.8]);
    end
    
    set(gca,'XTick',readsSX)
    set(gca,'YDir', 'normal')
    title('Allowed (white) and Excluded (red) Experimental Designs', 'FontSize', 14)
    
   experimentsEvaluated=size(pwrsCalc,1)*size(pwrsCalc,2);
   
   out='Experiments Evaluated'
   experimentsEvaluated
   

end

