function intervalData = readBed(iFile)
% Reads BedFile
% chr,start,stop, name,score,strand, thickStart,thickEnd,itemRgb,blockCount,blockSizes,blockStarts

tmpInputFile = tmpCopyFile(iFile); % Make tmp copy of file if required

fp = fopen(tmpInputFile,'r');
% Estimate number of columns in the data
firstLine = fgetl(fp);
splitLine = textscan(firstLine,'%s');
nCols = numel(splitLine{1});
fseek(fp,0,-1); % Rewind the file

fullFormatString = '%s%n%n%s%s%c%n%n%s%n%s%s'; % 12 columns

if nCols > 12
    nCols = 12;
end

try % Try reading all 12 BED columns
    formatString = [ fullFormatString(1:2*nCols) , '%*[^\n]' ]; 
    iFileData = textscan(fp,formatString,'CommentStyle','#', 'ReturnOnError', false);
catch
    if nCols >=6 % Try reading first 6 columns
        try
            nCols = 6;
            fseek(fp,0,-1); % Rewind the file
            formatString = [ fullFormatString(1:2*nCols) , '%*[^\n]' ];
            iFileData = textscan(fp,formatString,'CommentStyle','#', 'ReturnOnError', false);
        catch
            try % Try reading first 3 columns
                nCols = 3;
                fseek(fp,0,-1); % Rewind the file
                formatString = [ fullFormatString(1:2*nCols) , '%*[^\n]' ];
                iFileData = textscan(fp,formatString,'CommentStyle','#', 'ReturnOnError', false);
            catch ME
                fclose(fp);
                if ~strcmp(iFile,tmpInputFile)
                    delete(tmpInputFile);
                end
                error( ME.identifier , ME.message , 'ERROR: Parse error in BED file %s' , iFile );
            end
        end
    else
        try % Try reading first 3 columns
            nCols = 3;
            fseek(fp,0,-1); % Rewind the file
            formatString = [ fullFormatString(1:2*nCols) , '%*[^\n]' ];
            iFileData = textscan(fp,formatString,'CommentStyle','#', 'ReturnOnError', false);
        catch ME
            fclose(fp);
            if ~strcmp(iFile,tmpInputFile)
                delete(tmpInputFile);
            end            
            error( ME.identifier , ME.message , 'ERROR: Parse error in BED file %s' , iFile );
        end                
    end
end
fclose(fp);

if ~strcmp(iFile,tmpInputFile)
    delete(tmpInputFile);
end


intervalData = dataset();
intervalData.Properties.UserData.fileName = iFile;
intervalData.Properties.UserData.fileFormat = 'bed';

i = 1;
intervalData.chr = nominal(iFileData{i});
i = i+1;
intervalData.start = iFileData{i}+1;
i = i+1;
intervalData.stop = iFileData{i};
intervalData.strand = nominal( repmat( '+' , size(intervalData.chr) ) ); % initialize strand
if nCols > 3
    i = i+1;
    intervalData.name = iFileData{i};
    if nCols > 4
        i = i+1;
        intervalData.score = iFileData{i};
        if nCols > 5
            i = i+1;
            intervalData.strand = iFileData{i};
            intervalData.strand( ~ismember(intervalData.strand,{'+','-'}) ) = '+';
            intervalData.strand = nominal(intervalData.strand);
            if nCols > 6
                i = i+1;
                intervalData.thickStart = iFileData{i};
                i = i+1;
                intervalData.thickEnd = iFileData{i};
                i = i+1;
                intervalData.itemRgb = iFileData{i};
                i = i+1;
                intervalData.blockCount = iFileData{i};
                i = i+1;
                intervalData.blockSizes = iFileData{i};
                i = i+1;
                intervalData.blockStarts = iFileData{i};                
            end
        end
    end
end

end
