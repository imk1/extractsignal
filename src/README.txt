Author: Anshul Kundaje
Email: akundaje _at_ stanford _dot_ edu
Jan 16, 2011

-------------------------------------------------------------------------------
SUPPORTED PLATFORMS AND REQUIRED SOFTWARE
-------------------------------------------------------------------------------
- linux x64

- The native matlab code runs on matlab2009b and above.
  Start matlab and type the following before using the code .
  
  cd src/
  addpath(genpath(pwd));
  help extractSignal
  
- The compiled matlab code in /bin runs on unix x64 machines. It requires the 
matlab2010b run time environment (MCR) to be installed. Email me if you need the
MCR.
 
NOTE: The MCR is free and you do not need a matlab license to install or run it.

- NOTE: the syntax for calling the compiled code vs native matlab code are different

-------------------------------------------------------------------------------
INTRODUCTION
-------------------------------------------------------------------------------

extractSignal extracts signal data corresponding to user defined intervals from
genome-wide signal tracks. Several input formats are supported for the interval
data as well as the signal data.

You can also specify meta functions that operate on the signal vectors corresponding 
to each interval. eg. get the maximum value of signal for each interval.

You can also perform specific manipulations of the intervals themselves such as 
expand or contract the intervals or automatic correction to force all intervals
to be equal to the predominant interval length in the dataset.

The code handles invalid intervals elegantly by returning NaNs as the output signal.

The order of operations is as follows:
- Read in interval data (automatic correction and identification of invalid intervals)
- Perform explicit interval manipulation as requested by the user
- Optional automatic interval adjustment to force all intervals to be equal to the 
  predominant size in the interval file.
- Signal is then extracted for the adjusted intervals
- Signal is optionally smoothed if requested by the user
- Optional meta function is applied to the signal vector from each interval

NOTE:
- All arguments are specified as [argument_name]=[argument_value] . No spaces are allowed.
- The output signal can contain NaN (not a number). This can happen if the interval is invalid OR
  if the signal itself was missing/not measured throughout the interval.
  
-------------------------------------------------------------------------------
COMPILED BINARY USAGE:
-------------------------------------------------------------------------------

extractSignal [arguments]

--help, -h : help

====================
MANDATORY ARGUMENTS
====================
-i= / -intervalFile=       : path to interval file
-t= / -trackFile=          : path to track file


====================
OPTIONAL ARGUMENTS
====================
------------
Input File formats
------------
-if= / -intervalFileFormat=: format of interval file
                           : def: narrowpeak
                           :      bed,mat,summit,gff,gtf
-tf= / -trackFileFormat=   : format of track file
                           : def: mat
-tv= / trackFileVarPrefix= : track variable prefix (only relevant for mat signal files) {def: signal_}

------------
Interval manipulation specifications
------------

-us= / useSummit=          : use summit information to readjust intervals by slopping
                           : def: false (slop intervals using the ends of the intervals)
                           :      midpoint (use midpoint of interval as summit to slop)
                           :      peak (use summit column in narrowpeak or summit format as summit to slop)

-sl= / -slopLeft=          : increment/decrement left end of all intervals
                           : positive/negative values will increase/decrease interval length
                           : If '-ss' is set, then left is with respect to strand (i.e. 5' end of interval)
                           : def: 100 if summit information is being used ie. -us!=false OR if=summit (can only be set to > 0)
                           : def: 0 if '-us=false' (can be set to positive or negative values)

-sr= / -slopRight=         : increment/decrement right end of all intervals
                           : positive/negative values will increase/decrease interval length
                           : If '-ss' is set, then right is with respect to strand (i.e. 3' end of interval)
                           : def: 100 if summit information is being used ie. -us!=false OR if=summit (can only be set to > 0)
                           : def: 0 if '-us=false' (can be set to positive or negative values)

-ss= / -slopStrand=        : slop ends of intervals relative to strand info. If set to true left=5' and right=3'
                           : def: false
                           :      true

-fw= / -fixWidth=          : treat intervals as if they are fixed width intervals and do implicit adjustment (extend or contract to mode interval length AFTER slop adjustment)
                           : def: false
                           :      true
                           
------------
Output file specifications
------------                           
-o= / -outputFileName=     : output file name
                           : def: [intervalFileDir]/[trackFileRoot]_VS_[intervalFileRoot].[outputFileFormat]

-of= / -outputFileFormat=  : output file format
                           : def: mat
                           :      cagt Format: chr\tstart\stop\tsignalVals (separated by ,)

-ov= / -outputFileVar=     : output file variable name (valid for mat format)
                           : def: genvarname(esig_[trackFileRoot]_VS_[intervalFileRoot])
                           
------------
Meta function specifications
------------
-mf= / -metaFuncName=      : metaFunction to be applied to each signal vector corresponding to each interval
                           : def: signal : get signal vector
                           :      mean: get mean of signal vector
                           :      std : standard deviation
                           :      var : variance
                           :      geomean : geometric mean
                           :      harmean : harmonic mean
                           :      trimmean : trimmed mean (removing top and bottom x/2 %)
                           :      sum: sum
                           :      max: max and maxpos
                           :      min: min and minpos
                           :      iqr: interquartile range
                           :      range: range
                           :      zscore : normalize vector
                           :      skewness : skewness of distribution
                           :      kurtosis : kurtosis of distribution
                           :      quantile : multiple parameters allowed
                           :      samplerate : sample interval every 'x' bp
                           :      samplewidth: sample interval uniformly to get 'x' bp intervals
                           :      samplepos: sample specific positions relative to interval start=1

-mp= / -metaFuncParams=    : parameters for the metaFunctions
                           : def: std,var: 1 {0,1} 0 means std is obtained by dividing by N-1
                           : def: zscore: 1 {0,1} 0 means std is obtained by dividing by N-1
                           : def: trimmean: 5 [0,100] x/2 of top and bottom values will be removed before computing mean
                           : def: samplerate: 10 [>=1]
                           : def: samplewidth: 0 [>=0] 0 means full width
                           : def: samplepos: 0.5 [>=0] multiple arguments allowed
                           :      decimal values <1 means fractional position from start
                           :      values >=1 refers to actual coordinates,
                           : def: quantile: 0.5 [0,1] multiple arguments allowed
                           : def: otherwise: NaN

-ms= / -metaFuncSmooth=    : smoothing bandwidth for triweight kernel density smoothing
                           : def: 1 (no smoothing) >=1

-------------------------------------------------------------------------------
Meta Function usage
-------------------------------------------------------------------------------
Meta function names are specified using the -mf argument.
eg: -mf=quantile

Some meta functions can also take in parameters. These parameters are specified 
using the -mp option.
eg. -mp=0.5

Not all meta functions take in parameters. Also, some meta functions can take in 
multiple parameters. You simply specify these using multiple -mp arguments. 
eg. -mp=0.5 -mp=0.99

NOTE: the order of multiple -mp arguments is very important.
eg. The above order will produce an output with the 0.5 quantile followed by the
0.99 quantile

Currently, you can specify only one metafunction for a single run.
