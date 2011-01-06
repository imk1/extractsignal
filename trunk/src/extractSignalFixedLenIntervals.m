function [adjInterval,extSignal] = extractSignalFixedLenIntervals(interval,metaFunc,trackFileName,trackVarPrefix)
% Takes in a set of intervals, treats them as fixed length intervals,
% extracts signal data from trackFileName and apply metaFunc function on the
% signal data
%
% function [adjInterval,extSignal] = extractSignalFixedLenIntervals(interval,metaFunc,trackFileName,trackVarPrefix)
% ====================
% INPUT ARGUMENTS:
% ====================
% interval: is a dataset containing the following variables
%     .chr[nominal]: chromosome name
%     .start<double>: start position (1-based)
%     .stop<double>: stop position (1-based)
%     .strand<nominal>: OPTIONAL +/- (If - then signal vector is flipped) 
%
% metaFunc: Either an empty matrix or cell array (which means just return signal vector)
%           OR
%           function handle that will operate on the signal vector corresponding to the interval
%           metaFunc will take in a matrix of signal values (each row
%           corresponds to an interval and will return a matrix or cell
%           array)
%
% trackFileName: .mat file containing signal data in chunks
%
% trackVarPrefix: name of new variable to be generated using metaFunc
% ====================
% OUTPUT ARGUMENTS:
% ====================
% adjInterval: is a 'corrected' version of interval where the start/stops have been adjusted to make them valid
%            An addition column 'correction' is added
%            correction==1: Pad the left of the signal vector with NaNs to fit intervalLen
%            correction==2: Pad the right of the signal vector with NaNs to fit intervalLen
%            correction==-1: invalid interval (set signal vector to NaN)
% extSignal: extracted signal processed through metaFunc. Has same number of rows as interval and adjInterval
% ================================================
% NOTE: How anomalous intervals are handled
% (1) Intervals with start < 1
%     - start == 1
%     - Signal vector is left-padded with NaNs to make it equal to the predominant interval length
% (2) Intervals with stop > chrLen
%     - stop == chrLen
%     - Signal vector is right-padded with NaNs to make it equal to the predominant interval length
% (3) Intervals with lengths smaller than the predominant interval length
%     - Signal vector is expanded symmetrically on left and right to make it equal to the predominant interval length
% (4) Intervals with lengths larger than the predominant interval length
%     - Signal vector is truncated symmetrically on left and right
% ================================================

%% Validate metaFunc
assert( (isa(metaFunc,'function_handle') || isempty(metaFunc)) , ...
    'extractSignalFixedLenIntervals:metaFuncError','metaFunc must be a function handle OR an empty matrix/cell array');

if isempty(metaFunc)
    metaFunc= @(x)( x );
end

%% Validate interval data structure
interval = validateIntervalDataStruct( interval );

%% Preprocess dataset (TRACK)
% chrLen.chr[nominal] : names of chromosomes
% chrLen.len[double] : length of chromosome
% chrLen.maxChunkId[double] : maximum chunk id for that chromosome
%
% chunkLen.true<double>: length of chunks in .mat file
% chunkLen.multiplier<double>: virtual chunk length mutliplier to use to speed up loading
% chunkLen.use<double> = mutiplier * true: chunk length to actually use

[adjInterval , chrLen , chunkLen] = preprocessSignalTrack( interval , trackFileName);
nIntervals = size( adjInterval , 1 ); % Number of intervals

%% Check and correct anomalous intervals
[adjInterval , intervalLen] = correctIntervals( adjInterval , chrLen );

%% Get chunk ids and find unique start/stop chunk pairs
% Set start/stop chunkIds. These chunkIds correspond to the actual chunks in the track data variable. However, they are computed based on chunkLen.use

adjInterval.startChunkId = ( ceil( adjInterval.start / chunkLen.use ) - 1 ) * chunkLen.multiplier + 1; % left most sub-chunk of startChunkId
adjInterval.stopChunkId = ceil( adjInterval.stop / chunkLen.use ) * chunkLen.multiplier; % right most sub-chunk of stopChunkId
adjInterval.stopChunkId = min( double( join( ...
    adjInterval, chrLen, ...
    'Keys','chr', ...
    'LeftVars',{'stopChunkId'}, ...
    'RightVars',{'maxChunkId'} ) ), ...
    [] , 2 ); % truncate stopChunkId to maxChunkId

% Compute relative start
adjInterval.relStart = adjInterval.start - ( ( adjInterval.startChunkId - 1 ) * chunkLen.true );
% Get unique start/stop chunk ids
[uniquePairs, ~ , pairMembership] = unique( adjInterval( :, {'chr','startChunkId','stopChunkId'} ) );
nPairs = size(uniquePairs,1);
fprintf('Number of lookups = %d\n',nPairs);

%% For each start/stop chunk - load, concatenate, select, applyMetaFunc
for iPair = 1:nPairs
    if rem(iPair,10)==0
        fprintf('%d..',iPair);
    end
    currPair = uniquePairs(iPair,:);
    
    % Get all chunk names from startChunkId to stopChunkId
    chunkNames = strcat( ...
        {trackVarPrefix}, ...
        char(currPair.chr), ...
        '_',  ...
        num2str( ( currPair.startChunkId : currPair.stopChunkId )' , '%-d' ) );
    % Load chunks IN ORDER into a structure
    dataVec = load( trackFileName , chunkNames{:} );
    % Convert to ROW vector of doubles
    dataVec = transpose( cell2mat( struct2cell( dataVec ) ) );
    
    % Find members of current pair
    membersBinIdx = (pairMembership == iPair); % Use this to fill values in the adjInterval
    currDataVal = adjInterval(membersBinIdx,:); % adjInterval slice for current pair
    nMembers = size(currDataVal,1); % Number of members
    signalMatrix = nan(nMembers,intervalLen,'single'); % Initialize signalMatrix (all NaNs). Each signal is a ROW VECTOR
    
    % Get normal intervals from this set
    tempBinIdx = (currDataVal.correction == 0);  % normal intervals index into currDataVal/signalMatrix
    nNormPair = nnz(tempBinIdx);
    if nNormPair
        dataVecIdx = bsxfun(@plus,currDataVal.relStart(tempBinIdx),(0:(intervalLen-1))); % [relStart, ... ,relStart + intervalLen-1]
        signalMatrix(tempBinIdx,:) = dataVec(dataVecIdx); % Extract vectors for each member
        % signalMatrix(tempBinIdx,:) = arrayfun(@(x) dataVec(x : (x+intervalLen-1)) , currDataVal.relStart(tempBinIdx) , 'UniformOutput', true); % Extract vectors for each member
    end
    clear dataVecIdx;
    
    % Process anomalous short intervals
    % correction==1: Pad the left of the signal vector with NaNs to fit intervalLen
    % correction==2: Pad the right of the signal vector with NaNs to fit intervalLen
    if (nNormPair ~= nMembers)
        tempIdx = find(ismember(currDataVal.correction,[1,2])); % anomalous intervals index into currDataVal/signalMatrix
        nAnoPair = numel(tempIdx);
        for count = 1:nAnoPair
            anoIdx = tempIdx(count); % index into currDataVal/signalMatrix
            anoPair = currDataVal(anoIdx,:); % current interval
            anoLen = anoPair.stop - anoPair.start + 1; % interval length
            dataVecIdx = (anoPair.relStart : anoPair.relStart+anoLen-1); % index into dataVec
            padIdxStart = (anoPair.correction==1)*(intervalLen-anoLen+1) + (anoPair.correction==2)*1; % padding start index in signal vector
            padIdxStop = (anoPair.correction==1)*intervalLen + (anoPair.correction==2)*anoLen; % padding stop index in signal vector
            signalMatrix(anoIdx,(padIdxStart:padIdxStop)) = dataVec(dataVecIdx);
        end
    end
    
    % Reverse signal vectors corresponding to - strand intervals
    revBinIdx = (currDataVal.strand == '-');
    if nnz(revBinIdx)
        signalMatrix(revBinIdx,:) = signalMatrix( revBinIdx , (end:-1:1) );
    end
    clear revBinIdx;
    
    % Apply metaFunc on signalMatrix
    if (iPair~=1)
        extSignal(membersBinIdx,:) = metaFunc(signalMatrix);
    else
        funcOut = metaFunc(signalMatrix);
        [nrFuncOut,ncFuncOut] = size(funcOut); % size of funcOut
        assert(((nrFuncOut==nMembers) && (ncFuncOut >=1)), ...
            'extractSignalFixedLenIntervals:metaFuncError','size of output of metaFunc is not appropriate'); % check size of output
        if isa(funcOut , 'dataset')
            extSignal = datasetfun( @(x) repmat(x,nIntervals,1) , funcOut(1,:) , 'DatasetOutput' , true );
        else
            extSignal = repmat(funcOut(1,:),nIntervals,1);  % Initialize adjInterval.(dataVarName)
        end
        extSignal(membersBinIdx,:) = funcOut;
        clear funcOut;        
    end
end
fprintf('\n');

%% Clean up adjInterval
adjInterval(:,{'relStart','startChunkId','stopChunkId'}) = []; % Remove unwanted variables

end

% #############################################################################################
% AUXILIARY FUNCTIONS
% #############################################################################################
function [interval] = validateIntervalDataStruct( interval )
% Validates interval data structure

% ---------------------------
% interval must be a dataset
% ---------------------------
assert( isa( interval , 'dataset' ), ...
    'extractSignalFixedLenIntervals:argError', ...
    'Argument Error:interval: must be of class dataset with variables: chr,start,stop');

intervalProps = summary( interval );
% ---------------------------
% chr, start, stop and strand (OPTIONAL) must exist
% ---------------------------
assert( all( ismember( {'chr','start','stop'} , {intervalProps.Variables.Name} ) ), ...
    'extractSignalFixedLenIntervals:argError', ...
    'Argument Error:interval: missing one or more of the required variables: chr,start,stop,strand');
% If strand info does not exist set it to default '+'
if ~ismember( 'strand' , {intervalProps.Variables.Name} )
    interval.strand = nominal( repmat( '+' , length(interval.chr) , 1 ) );
end
end

% ===============================================================================================

function [adjInterval , chrLen , chunkLen] = preprocessSignalTrack( interval , trackFileName)
% Preprocess Signal Track dataset

% ---------------------------
% Load dataset parameters
% ---------------------------
try
    dataParams = load( trackFileName , 'globalParameters' , 'specificParameters' ); % load globalParameters, specificParameters
catch %#ok<CTCH>
    error('extractSignalFixedLenIntervals:dataParamError', ...
        'Unable to load globalParameters and specificParameters from %s', ...
        trackFileName);
end

% ---------------------------
% Get chunk length
% ---------------------------
chunkLen.true = dataParams.globalParameters.oChunkLen; % length of chunks in .mat file
chunkLen.multiplier = max( 1 , round( 5e7 / chunkLen.true ) ); % virtual chunk length mutliplier to use to speed up loading
chunkLen.use = chunkLen.multiplier * chunkLen.true; % chunk length to actually use

% ---------------------------
% Create chrLen table (variable names are chromosome names)
% chrLen.chr[nominal] : names of chromosomes
% chrLen.len[double] : length of chromosome
% chrLen.maxChunkId[double] : maximum chunk id for that chromosome
% ---------------------------
chrLen = dataset( ...
    { nominal( dataParams.specificParameters.chrNames' ) , 'chr' }, ...
    { dataParams.specificParameters.chrLen , 'len' } );
chrLen.maxChunkId = ceil( chrLen.len / chunkLen.true );

% ---------------------------
% Initialize adjInterval which is a copy of interval
% ---------------------------
adjInterval = interval( : , {'chr','start','stop','strand'} );

% ---------------------------
% Check that 'chr','start','stop' and 'strand' are the correct class
% ---------------------------
if ~isa( adjInterval.chr , 'nominal' )
    adjInterval.chr = nominal(adjInterval.chr);
end

if ~isa( adjInterval.start , 'double' )
    adjInterval.start = double(adjInterval.start);
end

if ~isa( adjInterval.stop , 'double' )
    adjInterval.stop = double(adjInterval.stop);
end

if ~isa( adjInterval.strand , 'nominal' )
    adjInterval.strand = nominal(adjInterval.strand);
end

% Check that chromosome names in interval are present in chrLen
assert( all( islevel( getlevels(adjInterval.chr) , getlevels(chrLen.chr) ) ), ...
    'extractSignalFixedLenIntervals:chrError', ...
    'Chromosome names dont match in interval argument and trackFileName');

end

% ===============================================================================================

function [adjInterval , intervalLen] = correctIntervals( adjInterval , chrLen )
% NOTE: How anomalous intervals are handled
% (1) Intervals with start < 1
%     - start == 1
%     - Signal vector is left-padded with NaNs to make it equal to the predominant interval length
% (2) Intervals with stop > chrLen
%     - stop == chrLen
%     - Signal vector is right-padded with NaNs to make it equal to the predominant interval length
% (3) Intervals with lengths smaller than the predominant interval length
%     - Signal vector is expanded symmetrically on left and right to make it equal to the predominant interval length
% (4) Intervals with lengths larger than the predominant interval length
%     - Signal vector is truncated symmetrically on left and right

nIntervals = size( adjInterval , 1 ); % Number of intervals

% ---------------------------
% Check that all stop positions are > start positions else swap positions
% ---------------------------

flipBinIdx = (adjInterval.stop < adjInterval.start);
adjInterval(flipBinIdx,{'start','stop'}) = adjInterval(flipBinIdx,{'stop','start'});
clear flipBinIdx;

intervalLen = mode(adjInterval.stop - adjInterval.start + 1); % predominant interval Length

% ---------------------------
% Check lengths of intervals and adjust short and long intervals
% ---------------------------
diffLen = ( (adjInterval.stop - adjInterval.start + 1) - intervalLen ) / 2; % 1/2 the difference in length from mode of intervalLength
invalidBinIdx = (diffLen ~= 0); % short and long intervals
if nnz(invalidBinIdx)
    adjInterval.start(invalidBinIdx) = adjInterval.start(invalidBinIdx) + round(diffLen(invalidBinIdx)); % shift start to the left (if short) or right (if long)
    adjInterval.stop(invalidBinIdx) = adjInterval.start(invalidBinIdx) + intervalLen - 1; % adjust stop accordingly
end
clear diffLen;

% ---------------------------
% Initialize Interval correction flag
% correction==1: Pad the left of the signal vector with NaNs to fit intervalLen
% correction==2: Pad the right of the signal vector with NaNs to fit intervalLen
% correction==-1: invalid interval (set signal vector to NaN)
% ---------------------------
adjInterval.correction = sparse(nIntervals,1);

% ---------------------------
% Check that start positions are >= 1 (correction flag = 1)
% ---------------------------
invalidBinIdx = (adjInterval.start < 1);
adjInterval.start(invalidBinIdx) = 1;
adjInterval.correction(invalidBinIdx) = 1; % set correction flag

% ---------------------------
% Check that stop positions are >= 1 (correction flag = -1)
% ---------------------------
invalidBinIdx = (adjInterval.stop < 1);
adjInterval.stop(invalidBinIdx) = 1;
adjInterval.correction(invalidBinIdx) = -1; % set correction flag

% ---------------------------
% Check that stop positions are <= the corresponding chromosome length (correction flag = 2)
% ---------------------------
tempLen = double( join( ...
    adjInterval , chrLen, ...
    'Keys' , 'chr', ...
    'LeftVars' , 'stop', ...
    'RightVars' , {'len'} ) ); % column 1 is stop, column 2 is chr_length
[adjInterval.stop , invalidBinIdx] = min( tempLen , [] , 2 );
invalidBinIdx = ( invalidBinIdx == 2 );
adjInterval.correction( invalidBinIdx ) = 2; % set correction flag

% ---------------------------
% Check that start positions are <= the corresponding chromosome length (correction flag = -1)
% ---------------------------
[ adjInterval.start , invalidBinIdx ] = min( [ adjInterval.start , tempLen(:,2) ], [], 2);
invalidBinIdx = (invalidBinIdx==2);
adjInterval.correction(invalidBinIdx) = -1; % set correction flag
clear invalidBinIdx tempLen;

end