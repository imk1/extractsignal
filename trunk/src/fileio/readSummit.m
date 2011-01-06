function intervalData = readSummit(iFile)
% Read summit format
% chr,summit,strand

tmpInputFile = tmpCopyFile(iFile); % Make tmp copy of file if required

fp = fopen(tmpInputFile,'r');
% Estimate number of columns in the data
firstLine = fgetl(fp);
splitLine = textscan(firstLine,'%s');
nCols = numel(splitLine{1});
fseek(fp,0,-1); % Rewind the file

fullFormatString = '%s%n%c'; % 3 columns
formatString = [ fullFormatString(1:2*min(3,nCols)) , '%*[^\n]' ];

try
    iFileData = textscan(fp,formatString,'CommentStyle','#', 'ReturnOnError', false);
catch ME
    fclose(fp);
    if ~strcmp(iFile,tmpInputFile)
        delete(tmpInputFile);
    end    
    error( ME.identifier , ME.message , 'ERROR: Parse error in summit file %s' , iFile );
end
fclose(fp);

if ~strcmp(iFile,tmpInputFile)
    delete(tmpInputFile);
end

intervalData = dataset();
intervalData.Properties.UserData.fileName = iFile;
intervalData.Properties.UserData.fileFormat = 'summit';

i = 1;
intervalData.chr = nominal(iFileData{i});
i = i+1;
intervalData.start = iFileData{i};
intervalData.strand = nominal( repmat( '+' , size(intervalData.chr) ) ); % initialize strand
if nCols > 2
    i = i+1;
    intervalData.strand = iFileData{i};
    intervalData.strand( ~ismember(intervalData.strand,{'+','-'}) ) = '+';
    intervalData.strand = nominal(intervalData.strand);
end

end