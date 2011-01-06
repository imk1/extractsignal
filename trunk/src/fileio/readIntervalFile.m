function intervalData = readIntervalFile( iFile , iFormat , readMode )
% Reads various types of genomic interval files
% function intervalData = readIntervalFile( iFile , iFormat , readMode )
% ===============
% INPUT ARGUMENTS
% ===============
% iFile : input Interval File (can be gzipped with a .gz extension)
% iFormat: input file format {narrowpeak','bed','mat','summit','gff','gtf'}
% readMode: (OPTIONAL) if set to 'minimal', returned fields are chr,start,stop,strand,summit. Else all fields are returned
% ===============
% OUTPUT ARGUMENTS
% ===============
% intervalData<dataset>
% .chr[nominal] : chromosome names
% .start[double]: start positions (1-based)
% .stop[double]: stop positions (1-based)
% .strand[nominal]: +/-
% .summit[double]: absolute position of summit (1-based) [OPTIONAL]

if ~exist('readMode','var')
    readMode = 'full';
end

%% Parse interval file
switch iFormat
    case 'mat'    
        intervalData = readIntervalFileMat(iFile);
    case 'narrowpeak'
        intervalData = readNarrowPeak(iFile);
    case 'bed'
        intervalData = readBed(iFile);
    case 'summit'
        intervalData = readSummit(iFile);
    case {'gff','gtf'}
        intervalData = readGff(iFile);
    otherwise
        error('ERROR: Unsupported file format %s for file %s' , iFormat , iFile);
end

if strcmpi(readMode,'minimal')    
    colNames = get( intervalData , 'VarNames');
    minimalSet = {'chr','start','stop','strand','summit'};
    retColNames = intersect( minimalSet , colNames );
    intervalData = intervalData( : , retColNames );
end

end