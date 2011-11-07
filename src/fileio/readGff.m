function intervalData = readGff(iFile)
% Reads gff or gtf file
% chr,source,feature,start,stop,score,strand

tmpInputFile = tmpCopyFile(iFile); % Make tmp copy of file if required

fp = fopen(tmpInputFile,'r');
formatString = '%s%s%s%n%n%s%c%*[^\n]';
try
    iFileData = textscan(fp,formatString,'CommentStyle','#', 'ReturnOnError', false);
catch ME
    fclose(fp);
    if ~strcmp(iFile,tmpInputFile)
        delete(tmpInputFile);
    end    
    error( ME.identifier , ME.message , 'ERROR: Parse error in GFF file %s' , iFile );
end
fclose(fp);

if ~strcmp(iFile,tmpInputFile)
    delete(tmpInputFile);
end

intervalData = dataset();
intervalData.Properties.UserData.fileName = iFile;
intervalData.Properties.UserData.fileFormat = 'gff';

i = 1;
intervalData.chr = nominal(iFileData{i});
i = i+1;
intervalData.source = iFileData{i};
i = i+1;
intervalData.feature = iFileData{i};
i = i+1;
intervalData.start = iFileData{i}; % already 1-based
i = i+1;
intervalData.stop = iFileData{i};
i = i+1;
intervalData.score = iFileData{i};
i = i+1;
intervalData.strand = iFileData{i};
intervalData.strand( ~ismember(intervalData.strand,{'+','-'}) ) = '+';
intervalData.strand = nominal(intervalData.strand);

end