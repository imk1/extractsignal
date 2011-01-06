function intervalData = readIntervalFileMat(iFile)
% Reads interval data mat file
var = who( '-file' , iFile , 'intervalData' );
assert( ~isempty(var) , 'ERROR: Interval File %s does not contain variable named intervalData' , iFile );
load( iFile , 'intervalData' );
assert(isa(intervalData,'dataset') , 'ERROR: Interval File %s contains variable intervalData but it is not a dataset array', iFile);
intervalData.Properties.UserData.fileName = iFile;
intervalData.Properties.UserData.fileFormat = 'mat';
% Check for mandatory fields
assert( all( ismember( {'chr','start','stop','strand'} , get(intervalData,'VarNames') ) ) , 'ERROR: interval dataset array in file %s does not contain all mandatory fields chr,start,stop,strand' , iFile );
end
