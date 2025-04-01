# IMe tester
Test suite for the [Inhibition Method Library](https://github.com/Reference-IMe/ime-lib)

If you use this work, please cite:
> D. Loreti, M. Artioli and A. Ciampolini, "Rollback-free recovery for a high performance dense linear solver with reduced memory footprint," IEEE Transactions on Parallel and Distributed Systems, vol. 35, no. 7, pp. 1307-1319, July 2024, [doi: 10.1109/TPDS.2024.3400365](https://doi.org/10.1109/TPDS.2024.3400365)
> D. Loreti, M. Artioli and A. Ciampolini, "Solving Linear Systems on High Performance Hardware with Resilience to Multiple Hard Faults," 2020 International Symposium on Reliable Distributed Systems (SRDS), Shanghai, China, 2020, pp. 266-275, [doi: 10.1109/SRDS51746.2020.00034](https://dx.doi.org/10.1109/SRDS51746.2020.00034)
> M. Artioli, D. Loreti and A. Ciampolini, "Fault Tolerant High Performance Solver for Linear Equation Systems," 2019 38th Symposium on Reliable Distributed Systems (SRDS), Lyon, France, 2019, pp. 113-11309, [doi: 10.1109/SRDS47363.2019.00022](https://doi.org/10.1109/SRDS47363.2019.00022)

## The suite
It's a tool for comparing the performance of IMe with some other linear solvers.
Type:
- `tester --help`	: to show command line options.
- `tester --list` : to show the list of testable routines.

### Testers
IMe routines are compared with those of LAPACK, ScaLAPACK, FTLA distributions and home-made Gaussian Elimination.
These distributions are downloaded from their repository and slightly modified in the makefile and/or in some code to better fit the test suite.

#### ScaLAPACK
Modified:
- makefile to accept variable for pointing to the local LAPACK lib in this repository
- `pgemraux.c` in `REDIST` to fix the "xxmr2d:out of memory" error (see [report](https://software.intel.com/en-us/forums/intel-math-kernel-library/topic/286499) and [fix](http://icl.cs.utk.edu/lapack-forum/viewtopic.php?f=2&t=465&p=1537&hilit=xxmr2d#p1537))

Added own made checkpoint-and-restore versions to deal with fault tolerance.

#### FTLA
Modified:
- `ftla_dcsum.c` to add external reference to an MPI communicator used in the caller code

Notes: safe constraints (exceptions are found, but no clear rule) on input parameters for FTLA are:
 - number of processes must be perfect square
 - matrix rank must be multiple of the square root of the number of processes (i.e. the number of columns or rows per process is the quotient)
 - the number of columns or rows per process must be multiple of the blocking factor
 
