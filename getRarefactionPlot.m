function [ rarefactionPlot ] = getRarefactionPlot( readsSequenced, genesFound, sampleIds, sampleNames, nControlSamples, nTestSamples, cutOff )
    %{
    Generates a plot from the reads sequenced and the number of genes
    found at each sequencing depth.  (A rarefaction plot)
    Because readsSequenced and genesExpressed differ in length,
    they are saved as long vectors with identifiers in the
    sampleIds vector.
    The legend for each sample is in the sampleNamesFile
    %}

    [ colorOfMarker lineStyles ]=getLineMarkers(nControlSamples, nTestSamples);
    
    totalSamples=nControlSamples+nTestSamples;
    
    rarefactionPlot=figure('Visible','on');
        
    hold on;
    for i=1:totalSamples
        r=readsSequenced(sampleIds==i);
        g=genesFound(sampleIds==i);
        plot(r,g, 'Color', colorOfMarker(i,:), 'LineStyle', lineStyles{i} );
    end
    xlabel('Reads Aligned to Genes', 'FontSize', 14)
    ylabelText=strcat({'Projected # of Genes Detected by '}, num2str(cutOff) ,{' or more reads'});
    ylabel(ylabelText, 'FontSize', 14);
    legend(sampleNames, 'Location', 'SouthEast');
    
    
end

