**Author**: Anshul Kundaje

**Email**: anshul _at_ kundaje _dot_ net

extractSignal extracts signal data corresponding to user defined intervals from
genome-wide signal tracks. Several input formats are supported for the interval
data as well as the signal data. A significant advantage of extractSignal is that
it is pretty fast especially when used with mat format genome-wide signal files
which allow quick random access. The plan is to support bigWig, hdf5 and some
text based formats such as bedGraph.

You can also specify meta functions that operate on the signal vectors corresponding
to each interval. eg. get the maximum value of signal for each interval.

You can also perform specific manipulations of the intervals themselves such as
expand or contract the intervals or automatic correction to force all intervals
to be equal to the predominant interval length in the dataset.

The code handles invalid intervals elegantly by returning NaNs as the output signal.

The order of operations is as follows:
  * Read in interval data (automatic correction and identification of invalid intervals)
  * Perform explicit interval manipulation as requested by the user
  * Optional automatic interval adjustment to force all intervals to be equal to the  predominant size in the interval file.
  * Signal is then extracted for the adjusted intervals
  * Signal is optionally smoothed if requested by the user
  * Optional meta function is applied to the signal vector from each interval



## SUPPORTED PLATFORMS AND REQUIRED SOFTWARE ##
  * linux x64

  * The native matlab code runs on matlab2009b and above.
Start matlab and type the following before using the code .
```
  cd src/
  addpath(genpath(pwd));
  help extractSignal
```

  * The compiled matlab code in /bin runs on unix x64 machines. It requires the matlab2010b run time environment (MCR) to be installed. MCR installation instructions
are given below.

**NOTE**: The MCR is free and you do not need a matlab license to install or run it.

**NOTE**: the syntax for calling the compiled code vs native matlab code are different

## MCR INSTALLTION INSTRUCTIONS ##

In order to run this code and/or any MATLAB compiled code, you will need the
MATLAB runtime library. Please use the correct version of the MCR.
This version of the executable was compiled using MCR V7.14 which is equivalent to the MATLAB R2010b release.

Download the installer for unix x64 from here http://www.broadinstitute.org/~anshul/softwareRepo/MCR2010b.bin

If you haven't installed the MCR, you MUST do that using this command
```
./MCRInstaller.bin -console
```
If you need to specify a specific temp directory then also use the option
```
-is:tempdir <tempdirname>
```
The installer will prompt you to select the directory (`<MCR_ROOT>`) you want to
install the MCR into. e.g. /home/akundaje/software/mcroot

**NOTE**: Make sure your installation directory has write permissions and has atleast 500 MB of disk space.

The installation should go smoothly with the above command.
However, if you are interested in other installation options you can consult
http://www.mathworks.com/access/helpdesk/help/toolbox/compiler/bru23df-1.html

**NOTE**: You need to install the MCR ONLY once on the machine/cluster you plan to
run MATLAB compiled code.

If you want to uninstall the MCR , follow this procedure:

Navigate to your MCR installation directory using the cd command.
cd into the _uninst directory
Run the uninstaller.bin program.
{{
./uninstaller.bin -console
}}_

## Setting paths ##
You need to set the following environment variables for the compiled MATLAB code to run correctly. These environment variables MUST be set before calling any matlab compiled code.

You can add the following lines to your .bashrc or .cshrc file if you want to avoid settings these variables every time you want to run the code

If you are using the bash shell or modifying .bashrc then use
```
MCRROOT=<MCR_ROOT>/v714
LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRROOT}/runtime/glnxa64
LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRROOT}/bin/glnxa64
LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRROOT}/sys/os/glnxa64
MCRJRE=${MCRROOT}/sys/java/jre/glnxa64/jre/lib/amd64
LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRJRE}/native_threads
LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRJRE}/server
LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRJRE}
XAPPLRESDIR=${MCRROOT}/X11/app-defaults
export LD_LIBRARY_PATH
export XAPPLRESDIR
```

If you are using the csh shell or modifying .cshrc then use
```
setenv MCRROOT <MCR_ROOT>/v714
setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:${MCRROOT}/runtime/glnxa64
setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:${MCRROOT}/bin/glnxa64
setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:${MCRROOT}/sys/os/glnxa64 
setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:${MCRROOT}/sys/java/jre/glnxa64/jre/lib/amd64/native_threads
setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:${MCRROOT}/sys/java/jre/glnxa64/jre/lib/amd64/server
setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:${MCRROOT}/sys/java/jre/glnxa64/jre/lib/amd64
setenv XAPPLRESDIR ${MCRROOT}/X11/app-defaults
```

You can add the /bin directory to your $PATH environment variable by adding the
following line to your .bashrc file
```
export PATH="[absolute_path_to_bin_directory]:${PATH}"
```
Now you can check if extractSignal works by typing
```
extractSignal -h
```

If you dont add the /bin directory to your $PATH. Then you need to go to the /bin
directory and call extract signal as
```
./extractSignal -h
```

## IMPORTANT NOTES ##
  * All arguments are specified as [argument\_name](argument_name.md)=[argument\_value](argument_value.md) . No spaces are allowed.
  * The output signal can contain NaN (not a number). This can happen if the interval is invalid OR if the signal itself was missing/not measured throughout the interval.

## COMPILED BINARY USAGE ##
```
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
```

## Meta Function usage ##

Meta function names are specified using the -mf argument.
```
eg: -mf=quantile
```
Some meta functions can also take in parameters. These parameters are specified
using the -mp option.
```
eg. -mp=0.5
```
Not all meta functions take in parameters. Also, some meta functions can take in multiple parameters. You simply specify these using multiple -mp arguments.
```
eg. -mp=0.5 -mp=0.99
```

**NOTE**: the order of multiple -mp arguments is very important.
eg. The above order will produce an output with the 0.5 quantile followed by the
0.99 quantile

Currently, you can specify only one metafunction for a single run.