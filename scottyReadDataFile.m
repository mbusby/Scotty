function [ statusReport errorMessage geneNames data sampleNames] = scottyReadDataFile( fileName, controlCols, testCols  )
    %{
    This function reads a data file set up as a tab-delimitted text file 
    in the following format:

    Gene_Name	Control_Rep_1	Control_Rep_2	Test_Rep_1	Test_Rep_2
    Gene_A	123	154	223	102
    Gene_B	12	7	8	4
    
    Mixing dos and unix files can give problems with line characters.
    Because Scotty runs on unix it executes a dos2unix script on 
    every file that is uploaded to ensure that the characters conform 
    to expected unix new lines.
    %}
    
    %initialize all returned variables
   
    statusReport=0;
    errorMessage='0';
    
    try
        %Set up file format
        fileFormat='%s';

        for i=1:controlCols+testCols
            fileFormat=strcat(fileFormat, '\t%f');
        end

        %Open file to get sample names
        fid=fopen(fileName, 'r');   

        columnNames = textscan(fid, '%s', 1+controlCols+testCols, 'delimiter', '\t');        
      
        %Transpose the column names into a sampleNames array of strings
        for i=2:1+controlCols+testCols
            sampleNames(i-1,1)=columnNames{1}(i,1);
        end

        %Now get the rest of the file
        %The control and test data is concatenated together
        %into one big data array.  Samples are the columns, genes are the
        %rows

        importedData = textscan(fid, fileFormat); 
       
        geneNames=importedData{1}(:,1);

        for i=2:controlCols+testCols+1
            data(:,i-1)=importedData{i}(:,1);
        end

        fclose(fid);
        
    catch errorMessage
        statusReport=1;
        %intialize returned variables
        geneNames=0;       
        sampleNames=0;
        data=0;
    end
    
end

