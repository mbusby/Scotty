Scotty
======
These are the files that power the backend of the Scotty web app available at http://euler.bc.edu/marthlab/scotty.

These files were written using an older version of Matlab.  At some point Matlab depricated one of the functions used 
which will cause running errors on newer versions (Matlab 2013+) of Matlab.

To use these functions on newer versions of Matlab where it says: 

RandStream.setDefaultStream(s);   

you have to change it to 

RandStream.setGlobalStream(s);


These are located in the following lines:

getProbabilitiesFromFit Line 23 

getProbSequencedMonteCarloByReads Line 7
