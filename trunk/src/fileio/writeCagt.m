function [ ] = writeCagt( cagtFileName , intervalData , signalData )        
% writes out signal data in CAGT format
% function [ ] = writeCagt( cagtFileName , intervalData , signalData )
% -------------
% INPUT ARGUMENTS
% -------------
% cagtFileName<string> : path to CAGT output file
%
% intervalData<dataset>
%   .chr[nominal] : chromosome names
%   .start[double]: start positions (1-based)
%   .stop[double]: stop positions (1-based)
%   .strand[nominal]: +/-
%   .summit[double]: absolute position of summit (1-based) [OPTIONAL]
%
% signalData: cell array of variable length single vectors OR matrix of singles

numIntervals = size(intervalData,1);
assert( (numIntervals == size(signalData,1)) , 'extractSignal:writeCagt' , 'ERROR: number of rows in intervalData dont match number of rows in signalData' );

%% Check file extension
[~ , ~ , fExt] = fileparts(cagtFileName);
if strcmpi( fExt , '.gz' )
    oFile = regexprep( cagtFileName , '\.gz$' , '');
    isGz = 1;
else
    oFile = cagtFileName;
    isGz = 0;
end

%% Write CAGT format
% chr\tstart\tstop\tsignalVals(separated by ,)
fp = fopen( oFile , 'w' );
assert ( (fp ~= -1) , 'extractSignal:writeCagt' , 'ERROR: cannot open file - ', cagtFileName );
if iscell(signalData)
    for irow = 1:numIntervals
        sigVal = sprintf( '%0.1f,' , signalData{irow} );
        sigVal = sigVal(1:end-1);
        fprintf( fp , '%s\t%d\t%d\t%s\n' , char(intervalData.chr(irow)) , intervalData.start(irow) , intervalData.stop(irow) , sigVal );
    end
else
    for irow = 1:numIntervals
        sigVal = sprintf( '%0.1f,' , signalData(irow,:) );
        sigVal = sigVal(1:end-1);
        fprintf( fp , '%s\t%d\t%d\t%s\n' , char(intervalData.chr(irow)) , intervalData.start(irow) , intervalData.stop(irow) , sigVal );
    end    
end
fclose(fp);

%% Gzip file if necessary
if isGz
    [status,~] = system( sprintf('gzip %s' , oFile) );
    assert( (status==0) , 'extractSignal:writeCagt' , 'ERROR: Unable to gzip file %s' , oFile );
end

end