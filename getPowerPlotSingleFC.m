function [ powerPlotSingleFC ] = getPowerPlotSingleFC( powerByReadDepth0, readDepths, fc, pCut  )
    %{
    This returns a power plot for the data showing the power that will
    be achieved for high and low variance genes.
    %}
    
    %Make it a percentage
    powerByReadDepth=powerByReadDepth0*100;   

    meanPower=mean(powerByReadDepth,2);
    
    powerPlotSingleFC=figure();
    hold on 
    plot(readDepths, powerByReadDepth, '--','Color', [0.8    0.8     0.8 ], 'HandleVisibility', 'Off')
    plot(readDepths, meanPower, 'b' );    
    plot(readDepths, powerByReadDepth(:,1), '--','Color', [0.8  0.8  0.8 ]) %Just for the legend's sake! 
    titleStr=strcat({'Power to Detect a '}, num2str(fc), 'X Fold Change Among Genes of Different Variances, p<', num2str(pCut));          
    title(titleStr);
    xlabel('Read Depth (Control Sample)');
    ylabel('Percentage of Genes Detected');
    legend('Overall Power', 'Power At Each Variance Decile', 'Location', 'SouthEast'); 
    axis([0 max(readDepths) 0 100]);    


end

