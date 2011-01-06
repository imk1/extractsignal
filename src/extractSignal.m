function varargout = extractSignal(varargin)
% Extracts signal from a genome-wide signal file corresponding to intervals in an interval File
% [signal , adjustedIntervals , intervals ] = extractSignal( intervalFile , trackFile , varargin )
% --------------------
% MANDATORY ARGUMENTS
% --------------------
% intervalFile: path of interval file (first argument)
% trackFile: path of track file (second argument)
% --------------------
% OPTIONAL param/value pairs
% --------------------
% if : format of interval file
%                            : def: 'narrowpeak'
%                            :      'bed','mat','summit','gff','gtf'
% tf : format of track file
%                            : def: 'mat'
% tv : track variable prefix (only relevant for mat files) {def: 'signal_'}
% us : use summit information to readjust intervals by slopping
%                            : def: 'false' (slop intervals using the ends of the intervals)
%                            :      'midpoint' (use midpoint of interval as summit to slop)
%                            :      'peak' (use summit column in narrowpeak or summit format as summit to slop)
% sl : increment/decrement left end of all intervals
%                            : positive/negative values will increase/decrease interval length
%                            : If '-ss' is set, then left is with respect to strand (i.e. 5' end of interval)
%                            : def: 100 if summit information is being used ie. -us!=false OR if=summit (can only be set to > 0)
%                            : def: 0 if '-us=false' (can be set to positive or negative values)
% sr : increment/decrement right end of all intervals
%                            : positive/negative values will increase/decrease interval length
%                            : If '-ss' is set, then right is with respect to strand (i.e. 3' end of interval)
%                            : def: 100 if summit information is being used ie. -us!=false OR if=summit (can only be set to > 0)
%                            : def: 0 if '-us=false' (can be set to positive or negative values)
% ss : slop ends of intervals relative to strand info. If set to true left=5' and right=3'
%                            : def: false
%                            :      true
% fw : treat intervals as if they are fixed width intervals and do implicit adjustment (extend or contract to mode interval length AFTER slop adjustment)
%                            : def: false
%                            :      true
% o  : output file name
%                            : def: [intervalFileDir]/[trackFileRoot]_VS_[intervalFileRoot].[outputFileFormat]
% of : output file format
%                            : def: 'mat'
%                            : 'cagt'
% ov : output file variable name (valid for mat format)
%                            : def: genvarname(esig_[trackFileRoot]_VS_[intervalFileRoot])
% mf : metaFunction to be applied to each signal vector corresponding to each interval
%                            : def: signal : get signal vector
%                            :      mean: get mean of signal vector
%                            :      std : standard deviation
%                            :      var : variance
%                            :      geomean : geometric mean
%                            :      harmean : harmonic mean
%                            :      trimmean : trimmed mean (removing top and bottom x/2 %)
%                            :      sum: sum
%                            :      max: max and maxpos
%                            :      min: min and minpos
%                            :      iqr: interquartile range
%                            :      range: range
%                            :      zscore : normalize vector
%                            :      skewness : skewness of distribution
%                            :      kurtosis : kurtosis of distribution
%                            :      quantile : multiple parameters allowed
%                            :      samplerate : sample interval every 'x' bp
%                            :      samplewidth: sample interval uniformly to get 'x' bp intervals
%                            :      samplepos: sample specific positions
% mp : parameters for the metaFunctions
%                            : def: std,var: 1 {0,1} 0 means std is obtained by dividing by N-1
%                            : def: zscore: 1 {0,1} 0 means std is obtained by dividing by N-1
%                            : def: trimmean: 5 [0,100] x/2 % of top and bottom values will be removed before computing mean
%                            : def: samplerate: 10 [>=1]
%                            : def: samplewidth: 0 [>=0] 0 means full width
%                            : def: samplepos: 0.5 [>=0] decimal values means fractional position from start, values >=1 refers to actual coordinates (Can be a vector of values)
%                            : def: quantile: 0.5 [0,1] (Can be a vector of values)
%                            : def: otherwise: NaN
% ms : smoothing bandwidth for triweight kernel density smoothing
%                            : def: 1 (no smoothing) >=1

%% Arguments for deployed version
% extractSignal(vargin)
% --help, -h : help
%
% ARGUMENTS
% -i= / -intervalFile=       : path to interval file
% -if= / -intervalFileFormat=: format of interval file
%                            : def: narrowpeak
%                            :      bed,mat,summit,gff,gtf
% -t= / -trackFile=          : path to track file
% -tf= / -trackFileFormat=   : format of track file
%                            : def: mat
% -tv= / trackFileVarPrefix= : track variable prefix (only relevant for mat files) {def: signal_}
% -us= / useSummit=          : use summit information to readjust intervals by slopping
%                            : def: false (slop intervals using the ends of the intervals)
%                            :      midpoint (use midpoint of interval as summit to slop)
%                            :      peak (use summit column in narrowpeak or summit format as summit to slop)
% -sl= / -slopLeft=          : increment/decrement left end of all intervals
%                            : positive/negative values will increase/decrease interval length
%                            : If '-ss' is set, then left is with respect to strand (i.e. 5' end of interval)
%                            : def: 100 if summit information is being used ie. -us!=false OR if=summit (can only be set to > 0)
%                            : def: 0 if '-us=false' (can be set to positive or negative values)
% -sr= / -slopRight=         : increment/decrement right end of all intervals
%                            : positive/negative values will increase/decrease interval length
%                            : If '-ss' is set, then right is with respect to strand (i.e. 3' end of interval)
%                            : def: 100 if summit information is being used ie. -us!=false OR if=summit (can only be set to > 0)
%                            : def: 0 if '-us=false' (can be set to positive or negative values)
% -ss= / -slopStrand=        : slop ends of intervals relative to strand info. If set to true left=5' and right=3'
%                            : def: false
%                            :      true
% -fw= / -fixWidth=          : treat intervals as if they are fixed width intervals and do implicit adjustment (extend or contract to mode interval length AFTER slop adjustment)
%                            : def: false
%                            :      true
% -o= / -outputFileName=     : output file name
%                            : def: [intervalFileDir]/[trackFileRoot]_VS_[intervalFileRoot].[outputFileFormat]
% -of= / -outputFileFormat=  : output file format
%                            : def: mat
%                            : def: cagt
% -ov= / -outputFileVar=     : output file variable name (valid for mat format)
%                            : def: genvarname(esig_[trackFileRoot]_VS_[intervalFileRoot])
% -mf= / -metaFuncName=      : metaFunction to be applied to each signal vector corresponding to each interval
%                            : def: signal : get signal vector
%                            :      mean: get mean of signal vector
%                            :      std : standard deviation
%                            :      var : variance
%                            :      geomean : geometric mean
%                            :      harmean : harmonic mean
%                            :      trimmean : trimmed mean (removing top and bottom x/2 %)
%                            :      sum: sum
%                            :      max: max and maxpos
%                            :      min: min and minpos
%                            :      iqr: interquartile range
%                            :      range: range
%                            :      zscore : normalize vector
%                            :      skewness : skewness of distribution
%                            :      kurtosis : kurtosis of distribution
%                            :      quantile : multiple parameters allowed
%                            :      samplerate : sample interval every 'x' bp
%                            :      samplewidth: sample interval uniformly to get 'x' bp intervals
%                            :      samplepos: sample specific positions
% -mp= / -metaFuncParams=    : parameters for the metaFunctions
%                            : def: std,var: 1 {0,1} 0 means std is obtained by dividing by N-1
%                            : def: zscore: 1 {0,1} 0 means std is obtained by dividing by N-1
%                            : def: trimmean: 5 [0,100] x/2 % of top and bottom values will be removed before computing mean
%                            : def: samplerate: 10 [>=1]
%                            : def: samplewidth: 0 [>=0] 0 means full width
%                            : def: samplepos: 0.5 [>=0] multiple arguments allowed, decimal values means fractional position from start, values >=1 refers to actual coordinates,
%                            : def: quantile: 0.5 [0,1] multiple arguments allowed
%                            : def: otherwise: NaN
% -ms= / -metaFuncSmooth=    : smoothing bandwidth for triweight kernel density smoothing
%                            : def: 1 (no smoothing) >=1

% --------------

%% Parse input parameters
% iParams<struct>
%   intervalFile.name<string>
%               .format<string> {def:'narrowpeak','bed','mat','summit','gff','gtf'}
%   trackFile.name<string>
%            .format<string> {def:'mat'}
%            .varPrefix<string> def: signal_
%   slop.summit<string> {def:'false','midpoint','peak'}
%       .left<double> def: 0 if no summit / 100 if summit
%       .right<double> def: 0 if no summit / 100 if summit
%       .strand<string> {def:false,true}
%   fixWidth<string> {def:false,true}
%   outFile.name<string>
%          .format<string> {def:'mat','cagt'}
%          .varName<string>
%   metaFunc.name<string> {def:'signal','mean','std','var','geomean','harmean','trimmean','sum','quantile','max','min','iqr','range','samplerate','samplewidth','samplepos','zscore','skewness','kurtosis'}
%   metaFunc.params[double]
%   metaFunc.smooth[double]
% --------------
allArgs = varargin;

try
    if isdeployed
        iParams = processInputArgumentsDeployed( allArgs );
    else
        iParams = processInputArgumentsNative( allArgs{:} );
    end      
catch ME
    if strcmp(ME.identifier,'normalExit:normalExit')
        return;
    else
        rethrow(ME);
    end
end

clear allArgs;
iParams.tempFile = [tempname(),'.mat'];

% Check if output file exists
if exist(iParams.outFile.name,'file')
    replacedOutFile = [iParams.outFile.name,'_',datestr(now(),30),'.old'];
    fprintf( 2 , 'WARNING: Output file %s exists and will be replaced. Older version is renamed %s\n' , iParams.outFile.name, replacedOutFile );
    movefile(iParams.outFile.name , replacedOutFile );
end

% --------------
%% Set up metaFunc
% metaFunc<functionHandle>
% --------------
metaFunc = initializeMetaFunc( iParams.metaFunc.name , iParams.metaFunc.params , iParams.metaFunc.smooth);

% --------------
%% Read interval File
% intervalData<dataset>: chr[nominal],start[double],stop[double],strand[nominal], OPTIONAL:summit[double]
% --------------
intervalData = readIntervalFile(iParams.intervalFile.name , iParams.intervalFile.format , 'minimal');

% Check summit information
if ismember('summit' , intervalData.Properties.VarNames)
    if any( intervalData.summit < intervalData.start ) && strcmp(iParams.slop.summit , 'peak')
        error('ERROR: Atleast one interval has invalid summit information (-1) in narrowPeak file. Cannot use -us=peak option for invalid summit information');
    end
end

intervalData.Properties.UserData = iParams;
save(iParams.tempFile,'intervalData','iParams'); % Temporarily save original interval data

% --------------
%% Adjust/slop intervals
% --------------
adjIntervalData = slopIntervalData(intervalData , iParams.slop.left , iParams.slop.right , iParams.slop.strand , iParams.slop.summit);
clear intervalData;

% --------------
%% Get signal data
% --------------
if iParams.fixWidth
    [adjIntervalData , outStruct.(iParams.outFile.varName) ] = extractSignalFixedLenIntervals( adjIntervalData , metaFunc , iParams.trackFile.name , iParams.trackFile.varPrefix );
else
    [adjIntervalData , outStruct.(iParams.outFile.varName) ] = extractSignalVarLenIntervals( adjIntervalData , metaFunc , iParams.trackFile.name , iParams.trackFile.varPrefix );
end
adjIntervalData.Properties.UserData = iParams;
save(iParams.tempFile,'adjIntervalData','-append'); % Temporarily save adjusted interval data
clear adjIntervalData;

% --------------
%% Ouput/save signal data
% --------------
switch iParams.outFile.format
    case {'mat','cagt'}
        [success,msg,msgId] = copyfile(iParams.tempFile , iParams.outFile.name);
        assert(success,msgId,msg);
        save( iParams.outFile.name , '-append', '-struct', 'outStruct' );
end

if isdeployed
    varargout = {};
else
    if nargout > 0        
        varargout{1} = outStruct.(iParams.outFile.varName);
        load(iParams.tempFile,'adjIntervalData');
        varargout{2} = adjIntervalData;
        load(iParams.tempFile,'intervalData');
        varargout{3} = intervalData;        
    else
        varargout = {};
    end
end

% Delete temp file
delete(iParams.tempFile);

end


% #############################################################################################
% AUXILIARY FUNCTIONS
% #############################################################################################

% ===============================================================================================

function helpLine = getUsageHelp()
%GETUSAGEHELP generates help/usage information for extractSignal
% function helpLine = getUsageHelp()

% hlnum = hlnum + 1; helpLine{hlnum} = '\n';
hlnum = 0;
hlnum = hlnum + 1; helpLine{hlnum} = '-------------------------------------------------------------------------\n';
hlnum = hlnum + 1; helpLine{hlnum} = 'Program: extractSignal (extract signal data from genome-wide signal/coverage data using interval data and apply optional meta functions on the signal vectors)\n';
hlnum = hlnum + 1; helpLine{hlnum} = 'Version: 0.1 (MATLAB r2009b)\n';
hlnum = hlnum + 1; helpLine{hlnum} = 'Contact: Anshul Kundaje (akundaje@stanford.edu)\n';
hlnum = hlnum + 1; helpLine{hlnum} = '-------------------------------------------------------------------------\n';
hlnum = hlnum + 1; helpLine{hlnum} = '\n';

hlnum = hlnum + 1; helpLine{hlnum} = 'Usage: extractSignal <options> -i=<intervalFile> -t=<trackFile>\n';
hlnum = hlnum + 1; helpLine{hlnum} = '\n';

hlnum = hlnum + 1; helpLine{hlnum} = '-------------------------------------------------------------------------\n';
hlnum = hlnum + 1; helpLine{hlnum} = 'Mandatory arguments:\n';
hlnum = hlnum + 1; helpLine{hlnum} = '-------------------------------------------------------------------------\n';
hlnum = hlnum + 1; helpLine{hlnum} = '\n';

hlnum = hlnum + 1; helpLine{hlnum} = '-i= / -intervalFile=        : path to interval file\n';
hlnum = hlnum + 1; helpLine{hlnum} = '-t= / -trackFile=           : path to track file\n';
hlnum = hlnum + 1; helpLine{hlnum} = '\n';

hlnum = hlnum + 1; helpLine{hlnum} = '-------------------------------------------------------------------------\n';
hlnum = hlnum + 1; helpLine{hlnum} = 'Optional arguments:\n';
hlnum = hlnum + 1; helpLine{hlnum} = '-------------------------------------------------------------------------\n';
hlnum = hlnum + 1; helpLine{hlnum} = '\n';

hlnum = hlnum + 1; helpLine{hlnum} = '--help, -h                  : help\n';
hlnum = hlnum + 1; helpLine{hlnum} = '\n';

hlnum = hlnum + 1; helpLine{hlnum} = '-if= / -intervalFileFormat= : format of interval file\n';
hlnum = hlnum + 1; helpLine{hlnum} = '                            : def: narrowpeak\n';
hlnum = hlnum + 1; helpLine{hlnum} = '                            :      bed,mat,summit,gff,gtf\n';
hlnum = hlnum + 1; helpLine{hlnum} = '\n';

hlnum = hlnum + 1; helpLine{hlnum} = '-tf= / -trackFileFormat=    : format of track file\n';
hlnum = hlnum + 1; helpLine{hlnum} = '                            : def: mat\n';
hlnum = hlnum + 1; helpLine{hlnum} = '\n';

hlnum = hlnum + 1; helpLine{hlnum} = '-tv= / -trackFileVarPrefix= : track variable prefix (only relevant for mat files)\n';
hlnum = hlnum + 1; helpLine{hlnum} = '                            : def: signal_\n';
hlnum = hlnum + 1; helpLine{hlnum} = '\n';

hlnum = hlnum + 1; helpLine{hlnum} = '-us= / -useSummit=          : use summit information to readjust intervals by slopping\n';
hlnum = hlnum + 1; helpLine{hlnum} = '                            : def: false (slop intervals using the ends of the intervals)\n';
hlnum = hlnum + 1; helpLine{hlnum} = '                            :      midpoint (use midpoint of interval as summit to slop)\n';
hlnum = hlnum + 1; helpLine{hlnum} = '                            :      peak (use summit column in narrowpeak or summit format as summit to slop)\n';
hlnum = hlnum + 1; helpLine{hlnum} = '\n';

hlnum = hlnum + 1; helpLine{hlnum} = '-sl= / -slopLeft=           : increment/decrement left end of all intervals\n';
hlnum = hlnum + 1; helpLine{hlnum} = '                            : positive/negative values will increase/decrease interval length\n';
hlnum = hlnum + 1; helpLine{hlnum} = '                            : If -ss is set, then left is with respect to strand (i.e. 5 prime end of interval)\n';
hlnum = hlnum + 1; helpLine{hlnum} = '                            : def: 100 if summit information is being used ie. -us!=false OR if=summit (can only be set to > 0)\n';
hlnum = hlnum + 1; helpLine{hlnum} = '                            : def: 0 if -us=false (can be set to positive or negative values)\n';
hlnum = hlnum + 1; helpLine{hlnum} = '\n';

hlnum = hlnum + 1; helpLine{hlnum} = '-sr= / -slopRight=          : increment/decrement right end of all intervals\n';
hlnum = hlnum + 1; helpLine{hlnum} = '                            : positive/negative values will increase/decrease interval length\n';
hlnum = hlnum + 1; helpLine{hlnum} = '                            : If -ss is set, then right is with respect to strand (i.e. 3 prime end of interval)\n';
hlnum = hlnum + 1; helpLine{hlnum} = '                            : def: 100 if summit information is being used ie. -us!=false OR if=summit (can only be set to > 0)\n';
hlnum = hlnum + 1; helpLine{hlnum} = '                            : def: 0 if -us=false (can be set to positive or negative values)\n';
hlnum = hlnum + 1; helpLine{hlnum} = '\n';

hlnum = hlnum + 1; helpLine{hlnum} = '-ss= / -slopStrand=         : slop ends of intervals relative to strand info. If set to true left=5prime and right=3prime\n';
hlnum = hlnum + 1; helpLine{hlnum} = '                            : def: false\n';
hlnum = hlnum + 1; helpLine{hlnum} = '                            :      true\n';
hlnum = hlnum + 1; helpLine{hlnum} = '\n';

hlnum = hlnum + 1; helpLine{hlnum} = '-fw= / -fixWidth=           : treat intervals as if they are fixed width intervals and do implicit adjustment (extend or contract to mode interval length AFTER slop adjustment)\n';
hlnum = hlnum + 1; helpLine{hlnum} = '                            : def: false\n';
hlnum = hlnum + 1; helpLine{hlnum} = '                            :      true\n';
hlnum = hlnum + 1; helpLine{hlnum} = '\n';

hlnum = hlnum + 1; helpLine{hlnum} = '-o= / -outputFileName=      : output file name\n';
hlnum = hlnum + 1; helpLine{hlnum} = '                            : def: [intervalFileDir]/[trackFileRoot]_VS_[intervalFileRoot].[outputFileFormat]\n';
hlnum = hlnum + 1; helpLine{hlnum} = '\n';

hlnum = hlnum + 1; helpLine{hlnum} = '-of= / -outputFileFormat=   : output file format\n';
hlnum = hlnum + 1; helpLine{hlnum} = '                            : def: mat\n';
hlnum = hlnum + 1; helpLine{hlnum} = '                            : def: cagt\n';
hlnum = hlnum + 1; helpLine{hlnum} = '\n';

hlnum = hlnum + 1; helpLine{hlnum} = '-ov= / -outputFileVar=      : output file variable name (valid for mat format)\n';
hlnum = hlnum + 1; helpLine{hlnum} = '                            : def: genvarname(esig_[trackFileRoot]_VS_[intervalFileRoot])\n';
hlnum = hlnum + 1; helpLine{hlnum} = '\n';

hlnum = hlnum + 1; helpLine{hlnum} = '-mf= / -metaFuncName=       : metaFunction to be applied to each signal vector corresponding to each interval\n';
hlnum = hlnum + 1; helpLine{hlnum} = '                            : def: signal\n';
hlnum = hlnum + 1; helpLine{hlnum} = '                            :      mean: get mean of signal vector\n';
hlnum = hlnum + 1; helpLine{hlnum} = '                            :      std : standard deviation\n';
hlnum = hlnum + 1; helpLine{hlnum} = '                            :      var : variance\n';
hlnum = hlnum + 1; helpLine{hlnum} = '                            :      geomean : geometric mean\n';
hlnum = hlnum + 1; helpLine{hlnum} = '                            :      harmean : harmonic mean\n';
hlnum = hlnum + 1; helpLine{hlnum} = '                            :      trimmean : trimmed mean (removing top and bottom x/2 percent)\n';
hlnum = hlnum + 1; helpLine{hlnum} = '                            :      sum: sum\n';
hlnum = hlnum + 1; helpLine{hlnum} = '                            :      max: max and maxpos\n';
hlnum = hlnum + 1; helpLine{hlnum} = '                            :      min: min and minpos\n';
hlnum = hlnum + 1; helpLine{hlnum} = '                            :      iqr: interquartile range\n';
hlnum = hlnum + 1; helpLine{hlnum} = '                            :      range: range\n';
hlnum = hlnum + 1; helpLine{hlnum} = '                            :      zscore : normalize vector\n';
hlnum = hlnum + 1; helpLine{hlnum} = '                            :      skewness : skewness of distribution\n';
hlnum = hlnum + 1; helpLine{hlnum} = '                            :      kurtosis : kurtosis of distribution\n';
hlnum = hlnum + 1; helpLine{hlnum} = '                            :      mean,std,var,geomean,harmean,trimmean,sum,max,min,iqr,range,zscore,skewness,kurtosis\n';
hlnum = hlnum + 1; helpLine{hlnum} = '                            :      quantile : multiple parameters allowed\n';
hlnum = hlnum + 1; helpLine{hlnum} = '                            :      samplerate : sample interval every x bp\n';
hlnum = hlnum + 1; helpLine{hlnum} = '                            :      samplewidth: sample interval uniformly to get x bp intervals\n';
hlnum = hlnum + 1; helpLine{hlnum} = '                            :      samplepos: sample specific positions\n';
hlnum = hlnum + 1; helpLine{hlnum} = '\n';

hlnum = hlnum + 1; helpLine{hlnum} = '-mp= / -metaFuncParams=     : parameters for the metaFunctions\n';
hlnum = hlnum + 1; helpLine{hlnum} = '                            : def: std,var: 1 {0,1} 0 means std is obtained by dividing by N-1\n';
hlnum = hlnum + 1; helpLine{hlnum} = '                            : def: trimmean: 5 [0,100] (removing top and bottom x/2 percent)\n';
hlnum = hlnum + 1; helpLine{hlnum} = '                            : def: samplerate: 10 [>=1] sample every x bp\n';
hlnum = hlnum + 1; helpLine{hlnum} = '                            : def: samplewidth: 0 [>=0] sample to get specific width (0 means full width)\n';
hlnum = hlnum + 1; helpLine{hlnum} = '                            : def: quantile: 0.5 [0,1] multiple arguments allowed\n';
hlnum = hlnum + 1; helpLine{hlnum} = '                            : def: samplepos: 0.5 [>=0] multiple arguments allowed, decimal values means fractional position from start, values >=1 refers to actual coordinates\n';
hlnum = hlnum + 1; helpLine{hlnum} = '                            : def: otherwise: NaN\n';
hlnum = hlnum + 1; helpLine{hlnum} = '\n';

hlnum = hlnum + 1; helpLine{hlnum} = '-ms= / -metaFuncSmooth=     : smoothing bandwidth for triweight kernel density smoothing\n';
hlnum = hlnum + 1; helpLine{hlnum} = '                            : def: 1 (no smoothing) >=1\n';
hlnum = hlnum + 1; helpLine{hlnum} = '\n';
hlnum = hlnum + 1; helpLine{hlnum} = '-------------------------------------------------------------------------\n';

helpLine = cell2mat(helpLine);
end

% ===============================================================================================

function iParams = processInputArgumentsDeployed(allArgs)
% Parses input arguments
% iParams<struct>
%   intervalFile.name<string>
%               .format<string> {def:'narrowpeak','bed','mat','summit','gff','gtf'}
%   trackFile.name<string>
%            .format<string> {def:'mat'}
%            .varPrefix<string> def: signal_
%   slop.summit<string> {def:'false','midpoint','peak'}
%       .left<double> def: 0 if no summit / 100 if summit
%       .right<double> def: 0 if no summit / 100 if summit
%       .strand<string> {def:false,true}
%   fixWidth<string> {def:false,true}
%   outFile.name<string>
%          .format<string> {def:'mat','cagt'}
%          .varName<string>
%   metaFunc.name<string> {def:'signal','mean','std','var','geomean','harmean','trimmean','sum','quantile','max','min','iqr','range','samplerate','samplepos','samplewidth','zscore','skewness','kurtosis','quantile'}
%   metaFunc.params[double]
%   metaFunc.smooth<double>

% --------------
% Check number of parameters and --help, -h
% --------------
helpLine = getUsageHelp();
minArgin = 2; % minimum number of mandatory arguments
if ( numel(allArgs) < minArgin ) || any(ismember( {'--help','-h'} , allArgs ))
    fprintf( 2 , helpLine );
    error('normalExit:normalExit','This is a normal exit');
end

% --------------
% interval file '-i=' '-intervalFile='
% --------------
iParams.intervalFile.name = processStringArg( allArgs ,  {'-i=','-intervalFile='} , '' , 1 , {} , 'file' , 0 );
% Get file extension
[iDir , iFileRoot , iFileExt] = fileparts( regexprep(iParams.intervalFile.name ,'\.gz$' , '') );
validIntervalFormats = {'narrowpeak','bed','mat','summit','gff','gtf'}; % valid interval file formats
if ~ismember( iFileExt , lower(validIntervalFormats) )
    iFileExt = 'narrowpeak';
end

% --------------
% interval file format '-if=' '-intervalFileFormat='
% {'narrowpeak','bed','mat','summit','gff','gtf'}
% --------------
iParams.intervalFile.format = processStringArg( allArgs ,  {'-if=','-intervalFileFormat='} , iFileExt , 0 , validIntervalFormats , '' , 0 );

% --------------
% track File Name '-t=' '-trackFile='
% --------------
iParams.trackFile.name = processStringArg( allArgs ,  {'-t=','-trackFile='} , '' , 1 , {} , 'file' , 0 );
% Get file extension
[~ , tFileRoot , tFileExt] = fileparts( regexprep(iParams.trackFile.name ,'\.gz$' , '') );
validTrackFormats = {'mat'}; % valid track file formats
if ~ismember( tFileExt , lower(validTrackFormats) )
    tFileExt = 'mat';
end

% --------------
% track File Format '-tf=' '-trackFile.format='
% {'mat'}
% --------------
iParams.trackFile.format = processStringArg( allArgs ,  {'-tf=','-trackFileFormat='} , tFileExt , 0 , validTrackFormats , '' , 0 );

% --------------
% track variable prefix '-tv=' '-trackFile.varPrefix='
% --------------
iParams.trackFile.varPrefix = processStringArg( allArgs ,  {'-tv=','-trackFileVarPrefix='} , 'signal_' , 0 , {} , '' , 0 );

% --------------
% use Summit '-us=' '-useSummit='
% --------------
iParams.slop.summit = processStringArg( allArgs ,  {'-us=','-useSummit='} , 'false' , 0 , {'midpoint','peak','false'} , '' , 0 );

% --------------
% slop Left '-sl=' '-slopLeft='
% --------------
% Check if file format is 'summit' or if slop is necessary based on use of summit information
if ( ismember(iParams.intervalFile.format,{'summit'}) || ~strcmp(iParams.slop.summit,'false') )
    defaultSlopLength = 100;
    minSlop = 0;
    maxSlop = Inf;
else
    defaultSlopLength = 0;
    minSlop = -Inf;
    maxSlop = Inf;
end
iParams.slop.left = processNumericArg( allArgs ,  {'-sl=','-slopLeft='}, defaultSlopLength , 0 , [minSlop,maxSlop] , 0 );

% --------------
% slop Right '-sr=' '-slopRight='
% --------------
iParams.slop.right = processNumericArg( allArgs ,  {'-sr=','-slopRight='}, defaultSlopLength , 0 , [minSlop,maxSlop] , 0 );

% --------------
% slop based on strand '-ss=' '-slopStrand='
% --------------
iParams.slop.strand = processStringArg( allArgs ,  {'-ss=','-slopStrand='} , 'false' , 0 , {'true','false'} , '' , 0 );
iParams.slop.strand = strcmpi( iParams.slop.strand , 'true' );
% --------------
% fixed Width intervals '-fw=' '-fixWidth='
% --------------
iParams.fixWidth = processStringArg( allArgs ,  {'-fw=','-fixWidth='} , 'false' , 0 , {'true','false'} , '' , 0 );
iParams.fixWidth = strcmpi( iParams.fixWidth , 'true' );
% --------------
% output File Format '-of=' '-outputFileFormat='
% {'mat','cagt'}
% --------------
validOutFileFormats = {'mat','cagt'};
iParams.outFile.format = processStringArg( allArgs ,  {'-of=','-outputFileFormat='} , 'mat' , 0 , validOutFileFormats , '' , 0 );

% --------------
% output File Name '-o=' '-outputFileName='
% --------------
defaultOutFileName = fullfile( iDir, [ tFileRoot , '_VS_' , iFileRoot , '.' , iParams.outFile.format] );
iParams.outFile.name = processStringArg( allArgs ,  {'-o=','-outputFileName='} , defaultOutFileName , 0 , {} , '' , 0 );

% --------------
% output File Variable name '-ov=' '-outputFileVar='
% --------------
defaultOutFileVar = genvarname( [ ...
    'esig_' , ...
    regexprep( tFileRoot , {'\.','-'} , '' ) , ...
    '_VS_' , ...
    regexprep( iFileRoot , {'\.','-'} , '' ) ] );

iParams.outFile.varName = genvarname( ...
    regexprep( ...
    processStringArg( allArgs ,  {'-ov=','-outputFileVar='} , defaultOutFileVar , 0 , {} , '' , 0 ) , ...
    {'\.','-'} , '' ) );

% --------------
% meta Function '-mf=' '-metaFuncName='
% --------------
validMetaFuncNames = {'mean','std','var','geomean','harmean','trimmean','sum', ...
    'quantile','max','min','iqr','range', ...
    'samplerate','samplewidth','samplepos','zscore','skewness','kurtosis', ...
    'signal'};
iParams.metaFunc.name = processStringArg( allArgs ,  {'-mf=','-metaFuncName='} , 'signal' , 0 , validMetaFuncNames , '' , 0 );

% --------------
% meta Function Parameters '-mp=' '-metaFuncParams='
% --------------
switch iParams.metaFunc.name
    case {'std','var','zscore'}
        iParams.metaFunc.params = round( processNumericArg( allArgs ,  {'-mp=','-metaFuncParams='}, 1 , 0 , [0,1] , 0 ) );
    case 'trimmean'
        iParams.metaFunc.params = processNumericArg( allArgs ,  {'-mp=','-metaFuncParams='}, 5 , 0 , [0,100] , 0 );
    case 'samplerate'
        iParams.metaFunc.params = round( processNumericArg( allArgs ,  {'-mp=','-metaFuncParams='}, 10 , 0 , [1,Inf] , 0 ) );
    case 'samplewidth'
        iParams.metaFunc.params = round( processNumericArg( allArgs ,  {'-mp=','-metaFuncParams='}, 0 , 0 , [0,Inf] , 0 ) );
    case 'samplepos'
        iParams.metaFunc.params = processNumericArg( allArgs ,  {'-mp=','-metaFuncParams='}, 0.5 , 0 , [0,Inf] , 1 );
    case 'quantile'
        iParams.metaFunc.params = processNumericArg( allArgs ,  {'-mp=','-metaFuncParams='}, 0.5 , 0 , [0,1] , 1 );
    otherwise
        iParams.metaFunc.params = NaN;
end

% --------------
% meta Function Smoothing BandwidthParameters -ms= / -metaFuncSmooth=
% --------------
iParams.metaFunc.smooth = round( processNumericArg( allArgs ,  {'-ms=','-metaFuncSmooth='}, 1 , 0 , [1,Inf] , 0 ) );

displayInputParams(iParams);

end

% ===============================================================================================

function iParams = processInputArgumentsNative( i , t , varargin )
% Process input arguments for native (within matlab calls)
% iParams<struct>
%   intervalFile.name<string>
%               .format<string> {def:'narrowpeak','bed','mat','summit','gff','gtf'}
%   trackFile.name<string>
%            .format<string> {def:'mat'}
%            .varPrefix<string> def: signal_
%   slop.summit<string> {def:'false','midpoint','peak'}
%       .left<double> def: 0 if no summit / 100 if summit
%       .right<double> def: 0 if no summit / 100 if summit
%       .strand<string> {def:false,true}
%   fixWidth<string> {def:false,true}
%   outFile.name<string>
%          .format<string> {def:'mat','cagt'}
%          .varName<string>
%   metaFunc.name<string> {def:'signal','mean','std','var','geomean','harmean','trimmean','sum','quantile','max','min','iqr','range','samplerate','samplewidth','samplepos','zscore','skewness','kurtosis'}
%   metaFunc.params[double]
%   metaFunc.smooth[double]

if nargin < 2
    help('extractSignal');
    error('normalExit:normalExit','This is a normal exit');
end

%% Set up input parser
p = inputParser;
p.KeepUnmatched = false;

% Mandatory arguments
p.addRequired( 'i' , @(x) ischar(x) && exist(x,'file') );
p.addRequired( 't' , @(x) ischar(x) && exist(x,'file') );

% Optional param/value pairs
validParams = {'narrowpeak','bed','mat','summit','gff','gtf'};
p.addParamValue( 'if' , 'narrowpeak' , ...
    @(x) ~isempty( validatestring( x , validParams , 'processInputArgumentsNative' , 'if' ) ) );

validParams = {'mat'};
p.addParamValue( 'tf' , 'mat' , ...
    @(x) ~isempty( validatestring( x , validParams , 'processInputArgumentsNative' , 'tf' ) ) );

p.addParamValue( 'tv' , 'signal_' , ...
    @(x) validateattributes( x , {'char'} , {'nonempty'} , 'processInputArgumentsNative' , 'tv' ) );

validParams = {'false','midpoint','peak'};
p.addParamValue( 'us' , 'false' , ...
    @(x) ~isempty( validatestring( x , validParams , 'processInputArgumentsNative' , 'us' ) ) );

p.addParamValue( 'sl' , NaN , ...
    @(x) validateattributes( x , {'numeric'} , {'finite','integer','scalar'} , 'processInputArgumentsNative' , 'sl' ) ); % Check

p.addParamValue( 'sr' , NaN , ...
    @(x) validateattributes( x , {'numeric'} , {'finite','integer','scalar'} , 'processInputArgumentsNative' , 'sr' ) ); % Check

p.addParamValue( 'ss' , false , ...
    @(x) validateattributes( x , {'logical'} , {'scalar'} , 'processInputArgumentsNative' , 'ss' ) );

p.addParamValue( 'fw' , false , ...
    @(x) validateattributes( x , {'logical'} , {'scalar'} , 'processInputArgumentsNative' , 'fw' ) );

p.addParamValue( 'o' , '' , ...
    @(x) validateattributes( x , {'char'} , {'nonempty'} , 'processInputArgumentsNative' , 'o' ) ); % Check

validParams = {'mat','cagt'};
p.addParamValue( 'of' , 'mat' , ...
    @(x) ~isempty( validatestring( x , validParams , 'processInputArgumentsNative' , 'of' ) ) );

p.addParamValue( 'ov' , '' , ...
    @(x) validateattributes( x , {'char'} , {'nonempty'} , 'processInputArgumentsNative' , 'ov' ) ); % Check

validParams = {'signal','mean','std','var','geomean','harmean','trimmean','sum','quantile','max','min','iqr','range','samplerate','samplepos','samplewidth','zscore','skewness','kurtosis','quantile'};
p.addParamValue( 'mf' , 'signal' , ...
    @(x) ~isempty( validatestring( x , validParams , 'processInputArgumentsNative' , 'mf' ) ) );

p.addParamValue( 'mp' , []); % Check

p.addParamValue( 'ms' , 1 , ...
    @(x) validateattributes( x , {'numeric'} , {'finite','integer','scalar','>=', 1} , 'processInputArgumentsNative' , 'ms' ) );

% Parse the input arguments
p.parse( i , t , varargin{:} );

%% Set up iParams
% iParams<struct>
%   intervalFile.name<string>
%               .format<string> {def:'narrowpeak','bed','mat','summit','gff','gtf'}
%   trackFile.name<string>
%            .format<string> {def:'mat'}
%            .varPrefix<string> def: signal_
%   slop.summit<string> {def:'false','midpoint','peak'}
%       .left<double> def: 0 if no summit / 100 if summit
%       .right<double> def: 0 if no summit / 100 if summit
%       .strand<string> {def:false,true}
%   fixWidth<string> {def:false,true}
%   outFile.name<string>
%          .format<string> {def:'mat','cagt'}
%          .varName<string>
%   metaFunc.name<string> {def:'signal','mean','std','var','geomean','harmean','trimmean','sum','quantile','max','min','iqr','range','samplerate','samplewidth','samplepos','zscore','skewness','kurtosis'}
%   metaFunc.params[double]
%   metaFunc.smooth[double]

iParams.intervalFile.name = p.Results.i;
iParams.intervalFile.format = p.Results.if;
iParams.trackFile.name = p.Results.t;
iParams.trackFile.format = p.Results.tf;
iParams.trackFile.varPrefix = p.Results.tv;
iParams.slop.summit = p.Results.us;
iParams.slop.left = p.Results.sl; % check
iParams.slop.right = p.Results.sr; % check
iParams.slop.strand = p.Results.ss;
iParams.fixWidth = p.Results.fw;
iParams.outFile.name = p.Results.o; % check
iParams.outFile.format = p.Results.of;
iParams.outFile.varName = p.Results.ov; % check
iParams.metaFunc.name = p.Results.mf;
iParams.metaFunc.params = p.Results.mp; % check
iParams.metaFunc.smooth = p.Results.ms;

%% Extra checks on arguments

% Check if summit information is used
checkSummit = ( ismember(iParams.intervalFile.format,{'summit'}) || ~strcmp(iParams.slop.summit,'false') );

% slop.left
if ismember('sl',p.UsingDefaults)       
    iParams.slop.left = 100 * checkSummit + 0 * ~checkSummit;    
end

% slop.right
if ismember('sr',p.UsingDefaults)       
    iParams.slop.right = 100 * checkSummit + 0 * ~checkSummit;
end

% if checkSummit sl and sr must be >=0
if checkSummit
    assert( all( [iParams.slop.left,iParams.slop.right] >= 0 ) , ...
        'ERROR: If interval Format is summit OR summit information is used slop left and right must be >= 0' );
end

% outFile.name
[iDir , iFileRoot , ~] = fileparts( regexprep(iParams.intervalFile.name ,'\.gz$' , '') );
[~ , tFileRoot , ~] = fileparts( regexprep(iParams.trackFile.name ,'\.gz$' , '') );
if ismember('o',p.UsingDefaults)    
    iParams.outFile.name = fullfile( iDir, [ tFileRoot , '_VS_' , iFileRoot , '.' , iParams.outFile.format] );
end

% outFile.varName
if ismember('ov',p.UsingDefaults)       
    iParams.outFile.varName = genvarname( [ ...
        'esig_' , ...
        regexprep( tFileRoot , {'\.','-'} , '' ) , ...
        '_VS_' , ...
        regexprep( iFileRoot , {'\.','-'} , '' ) ] );
end

% metaFunc.params
if ismember('op',p.UsingDefaults)
    switch iParams.metaFunc.name
        case {'std','var','zscore'}
            iParams.metaFunc.params = 1;
        case 'trimmean'
            iParams.metaFunc.params = 5;
        case 'samplerate'
            iParams.metaFunc.params = 10;
        case 'samplewidth'
            iParams.metaFunc.params = 0;
        case 'samplepos'
            iParams.metaFunc.params = 0.5;
        case 'quantile'
            iParams.metaFunc.params = 0.5;
        otherwise
            iParams.metaFunc.params = NaN;
    end    
end

switch iParams.metaFunc.name
    case {'std','var','zscore'}
        validateattributes( iParams.metaFunc.params , {'numeric'} , {'binary','scalar'} , 'processInputArgumentsNative' , 'mp' );
    case 'trimmean'
        validateattributes( iParams.metaFunc.params , {'numeric'} , {'scalar','>=',0,'<=',100} , 'processInputArgumentsNative' , 'mp' );
    case 'samplerate'
        validateattributes( iParams.metaFunc.params , {'numeric'} , {'scalar','integer','>=',1} , 'processInputArgumentsNative' , 'mp' );
    case 'samplewidth'
        validateattributes( iParams.metaFunc.params , {'numeric'} , {'scalar','integer','>=',0} , 'processInputArgumentsNative' , 'mp' );
    case 'samplepos'
        validateattributes( iParams.metaFunc.params , {'numeric'} , {'>=',0} , 'processInputArgumentsNative' , 'mp' );
    case 'quantile'
        validateattributes( iParams.metaFunc.params , {'numeric'} , {'>=',0,'<=',1} , 'processInputArgumentsNative' , 'mp' );
end

displayInputParams(iParams);

end


% ===============================================================================================

function [value] = processNumericArg( allArgs , attributeFlags , defaultVal , mandatory , valRange , multiArgs )
% Will process numeric arguments
% ------------------------------------------------------------------------------------------
% INPUT ARGUMENTS
% -----------------
% attributeFlags{<string>}/<string>: flag preceeding argument value eg. {'--mem='}
% defaultVal<double>: default value to use
% mandatory[0/1]: set to 1 if the field is mandatory
% valRange[double,double]: [min max] legal values
% multiArgs[0/1]: if set to 1 more than 1 instance of the attribute value pair is allowed
% -----------------
% OUTPUT ARGUMENTS
% -----------------
% value[double]: value of the numeric argument
% ------------------------------------------------------------------------------------------

attributeRegExp = strcat( '^' , cellstr(attributeFlags) , '|' ); % Convert cell array of attribute Flags into a combined regular expression
attributeRegExp = regexprep( [ attributeRegExp{:} ] , '\|$' , '' ); % Remove extra trailing |

matchIdx = find( cellfun( @length, regexp( allArgs , attributeRegExp ) ) ); % find index of attribute in allArgs
matchIdx = sort(matchIdx); % sort so that indices are in left to right order
nMatchIdx = numel(matchIdx); % number of matches

if (nMatchIdx==0)
    
    if mandatory % if mandatory throw error else set to default value
        error( 'ERROR: Missing mandatory argument %s\n' , attributeFlags{1} );
    else
        value = defaultVal;
    end
    
elseif ( ~multiArgs && (nMatchIdx > 1) )
    
    error( 'ERROR: Too many arguments of type %s\n', attributeFlags{1} );
    
else
    
    value = regexprep( allArgs(matchIdx) , attributeRegExp , '' ); % remove attributeFlag part
    value = regexprep( value , '^"|"$' , '' ); % remove leading and lagging quotes
    value = str2double(value); % convert to double
    % throw error if illegal value
    assert( all( ~isnan(value) ) , 'ERROR: Illegal value for attribute %s: Value MUST be numeric\n', attributeFlags{1} );
    % Check for illegal range of value
    if ~isempty(valRange)
        assert( (all( value >= valRange(1) ) && all( value <= valRange(2) ) ), 'ERROR: Illegal value for attribute %s: Legal range is [%f %f]\n', attributeFlags{1} , valRange(1) , valRange(2) );
    end
    
end

end

% ===============================================================================================

function [value] = processStringArg( allArgs, attributeFlags , defaultVal , mandatory , valRange , existCond , multiArgs )
% Will process string arguments
% ------------------------------------------------------------------------------------------
% INPUT ARGUMENTS
% -----------------
% attributeFlags{<string>}/<string>: flag preceeding argument value e.g. {'--mem='}
% defaultVal<string>: default value to use
% mandatory[0/1]: set to 1 if the field is mandatory
% valRange{}: cell arrray of legal string values. Set to {} if not used
% existCond<string>: (OPTIONAL) use the exist() function with the flag existCond. Useful for
%            checking if a file or directory exists. Set to '' if not used
% multiArgs[0/1]: if set to 1 more than 1 instance of the attribute value pair is allowed
% -----------------
% OUTPUT ARGUMENTS
% -----------------
% value{string}/<string>: value of the string argument
% ------------------------------------------------------------------------------------------

attributeRegExp = strcat( '^' , cellstr(attributeFlags) , '|' ); % Convert cell array of attribute Flags into a combined regular expression
attributeRegExp = regexprep( [attributeRegExp{:}] , '\|$' , '' ); % Remove extra trailing |

matchIdx = find( cellfun( @length , regexp( allArgs , attributeRegExp ) ) ); % find index of attribute in allArgs
matchIdx = sort(matchIdx); % sort so that indices are in left to right order
nMatchIdx = numel(matchIdx); % number of matches

if (nMatchIdx==0)
    
    if mandatory % if mandatory throw error else set to default value
        error( 'ERROR: Missing mandatory argument %s\n', attributeFlags{1} );
    else
        value = defaultVal;
    end
    
elseif ( ~multiArgs && (nMatchIdx > 1) )
    
    error( 'ERROR: Too many arguments of type %s\n', attributeFlags{1} );
    
else
    
    value = regexprep( allArgs(matchIdx) , attributeRegExp , '' ); % remove attributeFlag part
    value = regexprep( value , '^"|"$' , '' ); % remove leading and lagging quotes
    % throw error if illegal value
    if ~isempty(valRange)
        assert( all( ismember(value,valRange) ) , 'ERROR: Illegal value for attribute %s: Check help for valid values\n', attributeFlags{1} );
    end
    % check any exist() conditions if applicable
    if ~isempty(existCond)
        assert( all ( logical ( cellfun( @(x) exist(x,existCond), value, 'UniformOutput', true ) ) ) , ...
            'ERROR: Illegal value for attribute %s: %s does not exist\n', attributeFlags{1} , existCond );
    end
    
end

% Convert single cell into a string
if ( iscell(value) && ( numel(value)==1 ) )
    value = [ value{:} ];
end

end

% ===============================================================================================

function displayInputParams(iParams)
% displays input parameters
% iParams<struct>
%   intervalFile.name<string>
%               .format<string> {def:'narrowpeak','bed','mat','summit','gff','gtf'}
%   trackFile.name<string>
%            .format<string> {def:'mat'}
%            .varPrefix<string> def: signal_
%   slop.summit<string> {def:'false','midpoint','peak'}
%       .left<double> def: 0 if no summit / 100 if summit
%       .right<double> def: 0 if no summit / 100 if summit
%       .strand<string> {def:false,true}
%   fixWidth<string> {def:false,true}
%   outFile.name<string>
%          .format<string> {def:'mat','cagt'}
%          .varName<string>
%   metaFunc.name<string> {def:'signal','mean','std','var','geomean','harmean','trimmean','sum','quantile','max','min','iqr','range','samplerate','samplewidth','samplepos','zscore','skewness','kurtosis'}
%   metaFunc.params[double]
%   metaFunc.smooth[double]

fprintf('---------------------------------\nINPUT PARAMETERS\n---------------------------------\n');
fprintf('Interval file        : %s\n' , iParams.intervalFile.name);
fprintf('Interval file format : %s\n' , iParams.intervalFile.format);
fprintf('Track file           : %s\n' , iParams.trackFile.name);
fprintf('Track file format    : %s\n' , iParams.trackFile.format);
fprintf('Track file varPrefix : %s\n' , iParams.trackFile.varPrefix);
fprintf('Use Summit           : %s\n' , iParams.slop.summit);
fprintf('Slop Left            : %d\n' , iParams.slop.left);
fprintf('Slop Right           : %d\n' , iParams.slop.right);
fprintf('Slop Strand          : %d\n' , iParams.slop.strand);
fprintf('Fix Width            : %d\n' , iParams.fixWidth);
fprintf('Output file          : %s\n' , iParams.outFile.name);
fprintf('Output format        : %s\n' , iParams.outFile.format);
fprintf('Output variable Name : %s\n' , iParams.outFile.varName);
fprintf('Meta Function        : %s\n' , iParams.metaFunc.name);
fprintf('Meta Function Params : %g\n' , iParams.metaFunc.params(:));
fprintf('Smoothing bandwidth  : %d\n' , iParams.metaFunc.smooth);
fprintf('---------------------------------\n');
end