function [ powerPlot ] = getPowerPlot( powerByReadDepth, allReadDepths, readDepths, fcsPlot, pCut , nReps, nReadsPerRep)
    %{

    This returns a power plot for the data showing the power that will
    be achieved for each fold change in the plot.
    powerByReadDepths is a matrix of power by row=seq depth, column=variance
    layer=fold change
    %}

    %getReadDepth coords
    
    increment=length(allReadDepths)/4;
    
    for i=1:3
        quartVals(i)=allReadDepths(round(increment*i));
    end

    
    %Make it a percentage
    powerByReadDepth1=powerByReadDepth*100;   
    
    colorSpecs={'b' 'c' 'g'};
    colorSpecsMarker={'bo' 'co' 'go'};
    powerPlot=figure();    
   
    %Add the quantiles
    semilogx([quartVals(1) quartVals(1)], [0 100], 'Color', [.8 .8 .8],'HandleVisibility', 'On', 'LineWidth', 2)
    hold on
    semilogx([quartVals(2) quartVals(2)], [0 100], 'Color', [.7 .7 .7],'HandleVisibility', 'On', 'LineWidth', 3)
    semilogx([quartVals(3) quartVals(3)], [0 100], 'Color', [.8 .8 .8],'HandleVisibility', 'Off', 'LineWidth', 2)   
    legendLabels{1}='Read Depth Quartiles';
    legendLabels{2}='Median Read Depth';
    
    for k=1:size(powerByReadDepth1,3)  
        powerByReadDepthThis=powerByReadDepth1(:,:,k);
        meanPower=mean(powerByReadDepthThis,1);
        semilogx(readDepths(:,k), meanPower, colorSpecs{k}, 'LineWidth', 3)
        
        legendLabels{k+2}=strcat(num2str(fcsPlot(k)), 'X Fold Change');
    end
    
    
    titleStr=strcat({'% of Genes Detected at p<='}, ...
           num2str(pCut), {', '},num2str(nReps), {' Reps, '}, num2str(nReadsPerRep/10^6, '%.1f'), {' Million Reads'}   );          
    title(titleStr, 'FontSize', 12);
    xlabel('Read Depth (Control Sample)', 'FontSize', 12);
    ylabel('Percentage of Genes Detected', 'FontSize', 12);
    
    legend(legendLabels, 'Location', 'Best');     
    axis([0 max(max(readDepths)) 0 100]);  

end

