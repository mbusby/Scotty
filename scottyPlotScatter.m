function [ thePlot ] = scottyPlotScatter( xData, yData, genesToHighlight, colorsToHighlight )
    
    
    thePlot=figure('Visible','on');
    loglog(xData, yData, '.', 'Color', [0.7 0.7 0.7]);
    hold on
    for i=1:size(genesToHighlight,2)
        %loglog(xData(genesToHighlight(:, 1)==1), yData(genesToHighlight(:, 1)==1), 'Color', colorsToHighlight(i));
    end
    slope=getSlope(xData, yData);
    loglog(xData, xData*slope, 'k');
    xlabel('Control Data (Reads per Gene)', 'FontSize', 12);
    ylabel('Test Data (Reads per Gene)', 'FontSize', 12);
    title('Test Versus Control Data', 'FontSize', 12);

end

