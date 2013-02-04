function [ comparativeScatter ] = scottyPlotScatterComparingSamples( data, sampleNames )
    %{
        This makes a multi panelled plot comparing each sample with the 
        first sample to detect outliers

        Assumes data has 2 or more columns
    %}

    totalSamples=size(data,2);
    nPanelCols=ceil(sqrt(totalSamples-1));
    
    if nPanelCols*(nPanelCols-1)>=totalSamples-1
        nPanelRows=nPanelCols-1;
    else
        nPanelRows=nPanelCols;
    end
    
   comparativeScatter=figure(); 
    
   if nPanelRows>1
       ms=3;
   else
       ms=6;       
   end
   
   out='About to start plot'  
   
   %If there is enough room to make a 3x3 cover do that, or else
   %just do the best and worse
   if size(data,2)<=5
   
       for i=2:size(data,2)

            x=data(:,1);
            y=data(:,i);
            
            [ax]=subplot(nPanelRows, nPanelCols, i-1);

            loglog(x,y, '.', 'MarkerSize', ms);

            xlabel(sampleNames{1});
            ylabel(sampleNames{i});

            if nPanelRows>2
                set(ax,'XTickLabel',[]); %remove y1 ticks from the right side
                set(ax,'YTickLabel',[]); %remove y1 ticks from the right side
            end
        end
        % Extract axes handles of all subplots from the figure

        axesHandles = findobj(get(comparativeScatter,'Children'), 'flat','Type','axes');
        % Set the axis property to square
        axis(axesHandles,'square')
       
        
   else
        [ npStd controlSampleIds testSampleIds  ] = getDispersionAllComparisons(data);
       
        subplot(1,2,1)
        cMin=controlSampleIds(npStd==min(npStd));
        tMin=testSampleIds(npStd==min(npStd));
             
        Ci=data(:,cMin);
        Ti=data(:,tMin);        
        loglog(Ci, Ti, '.');
        xlabel(sampleNames(cMin), 'FontSize',12);
        ylabel(sampleNames(tMin), 'FontSize',12);
        title('Most Similar Replicates', 'FontSize', 12);
        m=max(max(Ci), max(Ti))+1000;
        axis([1 m 1 m]);
        
        
        subplot(1,2,2)        
        cMax=controlSampleIds(npStd==max(npStd));
        tMax=testSampleIds(npStd==max(npStd));    
        
        Ci=data(:,cMax);
        Ti=data(:,tMax);        
        loglog(Ci, Ti, '.');
        xlabel(sampleNames(cMax), 'FontSize',12);
        ylabel(sampleNames(tMax), 'FontSize',12);
        title('Least Similar Replicates', 'FontSize', 12);
        m=max(max(Ci), max(Ti))+1000;
        axis([1 m 1 m]);
       
        
        axesHandles = findobj(get(comparativeScatter,'Children'), 'flat','Type','axes');
         
        % Set the axis property to square
        axis(axesHandles,'square')
       
       
       
   end
    
end

