function metaFuncHandle = initializeMetaFunc( metaFuncName, metaFuncParams, smoothBw )
% Generates function handle to meta function that will process each interval
% metaFuncName<string> {def:'none','mean','std','var','geomean','harmean','trimmean','sum','quantile','max','min','iqr','range','samplerate','samplewidth','zscore','skewness','kurtosis'}
% metaFuncParams[double] : parameters for metaFunc
% smoothBw[double]: smoothing bandwidth
% -mf= / -metaFuncName=      : metaFunction to be applied to each signal vector corresponding to each interval
%                            : def: none
%                            :      'mean','std','mode','var','geomean','harmean','trimmean','sum','max','min','iqr','range','zscore','skewness','kurtosis'
%                            :      'quantile' : multiple parameters allowed
%                            :      'samplerate' : sample interval every 'x' bp
%                            :      'samplewidth': sample interval uniformly to get 'x' bp intervals
%                            :      'samplepos': sample specific positions
% -mp= / -metaFuncParams=    : parameters for the metaFunctions
%                            : def: std,var: 1 {0,1} 0 means divide by N-1
%                            : def: zscore: 1 {0,1} 0 means std is obtained by dividing by N-1
%                            : def: trimmean: 5 [0,100]
%                            : def: samplerate: 10 [>=1]
%                            : def: samplewidth: 0 [>=0] 0 means full width
%                            : def: samplepos: 0.5 [>=0] multiple arguments allowed, decimal values means fractional length offset from start, values >=1 refers to actual coordinates,
%                            : def: quantile: 0.5 [0,1] multiple arguments allowed
%                            : def: otherwise: NaN
% NOTE: Each function MUST be able to take in a row vector or a matrix (X) where each row represents a signal vector
%       and perform an operation on each vector to output an array of data structures with the 
%       SAME NUMBER OF ROWS AS X and a FIXED NUMBER OF COLUMNS. The number of columns can be different from the 
%       number of cols of X. The data returned by the function can be an array of built in data types like
%       like ints/doubles etc. OR structs/cells or it can be a dataset array.
% NOTE: Make sure the functions can handle all NaN vectors

switch metaFuncName
    case 'mean'
        metaFuncHandle = @(x) nanmean( smoothSignal(x,smoothBw) , 2 );
    case 'mode'
        metaFuncHandle = @(x) mode( smoothSignal(x,smoothBw) , 2 );
    case 'std'
        metaFuncHandle = @(x) nanstd( smoothSignal(x,smoothBw) , metaFuncParams(1) , 2 );
    case 'var'
        metaFuncHandle = @(x) nanvar( smoothSignal(x,smoothBw) , metaFuncParams(1) , 2 );
    case 'geomean'
        metaFuncHandle = @(x) exp( nanmean( log(smoothSignal(x,smoothBw)) , 2 ) );
    case 'harmean'
        metaFuncHandle = @(x) 1 ./nanmean( 1./smoothSignal(x,smoothBw) , 2 ) ;
    case 'trimmean'
        metaFuncHandle = @(x) trimmean( smoothSignal(x,smoothBw) , metaFuncParams(1) , 'round' , 2 );
    case 'sum'
        metaFuncHandle = @(x) nansum( smoothSignal(x,smoothBw) , 2 );
    case 'max'
        metaFuncHandle = @(x) getMax( smoothSignal(x,smoothBw) );
    case 'min'
        metaFuncHandle = @(x) getMin( smoothSignal(x,smoothBw) );
    case 'iqr'
        metaFuncHandle = @(x) iqr( smoothSignal(x,smoothBw) , 2 );
    case 'range'
        metaFuncHandle = @(x) range( smoothSignal(x,smoothBw) , 2 );
    case 'zscore'
        metaFuncHandle = @(x) nanZscore( smoothSignal(x,smoothBw) , metaFuncParams(1) , 2 );
    case 'skewness'
        metaFuncHandle = @(x) skewness( smoothSignal(x,smoothBw) , 1 , 2 );
    case 'kurtosis'
        metaFuncHandle = @(x) kurtosis( smoothSignal(x,smoothBw) , 1 , 2 );
    case 'quantile'
        metaFuncHandle = @(x) quantile( smoothSignal(x,smoothBw) , metaFuncParams,2 );
    case 'samplerate'
        metaFuncHandle = @(x) sampleRate( smoothSignal(x,smoothBw) , metaFuncParams(1) );
    case 'samplewidth'
        metaFuncHandle = @(x) sampleWidth( smoothSignal(x,smoothBw) , metaFuncParams(1) );
    case 'samplepos'
        metaFuncHandle = @(x) samplePos( smoothSignal(x,smoothBw) , metaFuncParams );
    otherwise
        metaFuncHandle = @(x) smoothSignal(x,smoothBw) ;
end

end