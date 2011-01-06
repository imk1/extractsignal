This directory contains functions to extract data from genome-wide data tracks stored in .mat format

*******************************
Structure of the .mat file
*******************************
The .mat file contains genome wide values for some genomic signal track(s).

There is one value per track for every nucleotide. 
The value is always a double BUT it can represent characters, binary, discrete or continuous data.
Illegal values are stored as NaN.

The genome-wide signal track(s) are stored as chunked matrices where the rows represent the nucleotide position and the column represents one or more tracks. 
Typically, they are simply column vectors. All chunks in a .mat dataset are of a fixed length (#rows) except the last chunk in a chromosome (it can be smaller or larger).

Each chunk is named <prefix>_<chrname>_<chunkid>
<prefix>: defines what type of data it is e.g. signal, maxtags etc.

*******************************
Functions
*******************************
# Extract signal vector corresponding to a single interval
extractData.m

# Extract signal matrix corresponding to a several fixed length intervals
extractFixedLengthBatch.m

# Extract signal cell arrays corresponding to a several variable length intervals
extractVarLengthBatch.m
