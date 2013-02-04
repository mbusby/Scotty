function [ status ] = printResultFile( resultString, dataDirectory, outputTag )
    %Prints the result summary file, whether at the end or if it crashes
    
    filename=strcat(dataDirectory,'resultsSummaryFile',outputTag,'.txt');
    delete(filename);   
    fileID = fopen(filename,'w');
    fprintf(fileID,'%s',resultString);
    fclose(fileID);   
    
    status=0;
    
    
end

