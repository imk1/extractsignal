function intervalData = readNarrowPeak(iFile)
% Reads narrowPeak file
% chr,start,stop,name,score,strand,signalValue,pValue,qValue,summit

tmpInputFile = tmpCopyFile(iFile); % Make tmp copy of file if required

fp = fopen(tmpInputFile,'r');
formatString = '%s%n%n%s%s%c%n%n%n%n%*[^\n]'; % 9 columns
try
    iFileData = textscan(fp,formatString,'CommentStyle','#', 'ReturnOnError', false);
catch ME
    fclose(fp);
    if ~strcmp(iFile,tmpInputFile)
        delete(tmpInputFile);
    end
    error( ME.identifier , ME.message , 'ERROR: Parse error in narrowPeak file %s' , iFile );
end
fclose(fp);

if ~strcmp(iFile,tmpInputFile)
    delete(tmpInputFile);
end

intervalData = dataset();
intervalData.Properties.UserData.fileName = iFile;
intervalData.Properties.UserData.fileFormat = 'narrowpeak';

i = 1;
intervalData.chr = nominal(iFileData{i});
i = i+1;
intervalData.start = iFileData{i}+1;
i = i+1;
intervalData.stop = iFileData{i};
i = i+1;
intervalData.name = iFileData{i};
i = i+1;
intervalData.score = iFileData{i};
i = i+1;
intervalData.strand = iFileData{i};
intervalData.strand( ~ismember(intervalData.strand,{'+','-'}) ) = '+';
intervalData.strand = nominal(intervalData.strand);
i = i+1;
intervalData.signalValue = iFileData{i};
i = i+1;
intervalData.pValue = iFileData{i};
i = i+1;
intervalData.qValue = iFileData{i};
i = i+1;
intervalData.summit = intervalData.start + iFileData{i};    

end
