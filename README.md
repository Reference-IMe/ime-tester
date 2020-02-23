# IMePACK

Inhibition Method PACKage

A collection of linear algebra routines.

##### Testers
Some code for testing is provided.
IMe routines are compared with those of LAPACK, ScaLAPACK and FTLA distributions.
These distributions are comprised in this repository, slightly modified in the makefile and/or in some code
- ScaLAPACK's makefile to accept variable for pointing to the local LAPACK lib in this repository
- ScaLAPACK's `pgemraux.c` in `REDIST` to fix the "xxmr2d:out of memory" error (see [report](https://software.intel.com/en-us/forums/intel-math-kernel-library/topic/286499) and [fix](http://icl.cs.utk.edu/lapack-forum/viewtopic.php?f=2&t=465&p=1537&hilit=xxmr2d#p1537))
- FTLA's `ftla_dcsum.c` to add external reference to a MPI communicator used in the caller code

###### Notes
`src` is the *current* source folder

`src-vX.Y` is the version *X.Y* source folder
