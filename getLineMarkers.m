function [ colorOfMarker lineStyles ] = getLineMarkers(nControlSamples, nTestSamples)
        %{
        This gets a set of line markers for plotting data on line graphs.
        The control samples are set to cool colors.
        The test samples are set to hot colors.
        %}

        potentialControlColors=[  0 0 1;
            0 1 1;  
            0 1 0;  
            0.6902  0.7686    0.8706;
            0 .5 0; ]; 

        potentialTestColors=[ 1 0 0; 1 0 1; 1.0000    0.7529    0.7961; 1 .6 0;    1.0000    0.4980    0.3137; ] ;

        potentialLineStyles={'-', ':', '--', '-.'};
        
        colorCtr=1;
        lineCtr=1;
        colorOfMarker=zeros(nControlSamples+nTestSamples,3);
        lineStyles={};
        
        for i=1:nControlSamples
           colorOfMarker(i,:)=potentialControlColors(colorCtr,:);           
           lineStyles{i}=potentialLineStyles{lineCtr};     
           %Inefficient but quick and dirty counters  
           if lineCtr==4
               lineCtr=1;
               colorCtr=colorCtr+1;               
           else
               lineCtr=lineCtr+1;
           end
           
           if colorCtr==5
               colorCtr=1;
           end
           
        end
        
        colorCtr=1;
        lineCtr=1;
        for i=nControlSamples+1:nControlSamples+nTestSamples
           colorOfMarker(i,:)=potentialTestColors(colorCtr,:);           
           lineStyles{i}=potentialLineStyles{lineCtr};         
           
           if lineCtr==4
               lineCtr=1;
               colorCtr=colorCtr+1;               
           else
               lineCtr=lineCtr+1;
           end         
           
        end
   

end

