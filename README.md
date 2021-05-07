# IMePACK

Inhibition Method PACKage

A collection of linear algebra routines.

## Testers
Some code for testing is provided.
IMe routines are compared with those of LAPACK, ScaLAPACK, FTLA distributions and home-made Gaussian Elimination.
These distributions are comprised in this repository, slightly modified in the makefile and/or in some code.

### ScaLAPACK
- makefile to accept variable for pointing to the local LAPACK lib in this repository
- `pgemraux.c` in `REDIST` to fix the "xxmr2d:out of memory" error (see [report](https://software.intel.com/en-us/forums/intel-math-kernel-library/topic/286499) and [fix](http://icl.cs.utk.edu/lapack-forum/viewtopic.php?f=2&t=465&p=1537&hilit=xxmr2d#p1537))

### FTLA
- `ftla_dcsum.c` to add external reference to an MPI communicator used in the caller code

## Notes
Safe constraints (exceptions are found, but no clear rule) on input parameters for FTLA are:
 - number of processes must be perfect square
 - matrix rank must be multiple of the square root of the number of processes (i.e. the number of columns or rows per process is the quotient)
 - the number of columns or rows per process must be multiple of the blocking factor
 