### Hints for achieving high performance on a distributed memory computer

source: [netlib](http://netlib.org/scalapack/slug/node106.html#SECTION04511000000000000000)

- nprocs = N\*M/1'000'000 for an N\*M matrix. This provides a local matrix of size approximately 1000 by 1000.  
- Do not try to solve a small problem on too many processors. 
- Do not exceed physical memory. 
- Use an efficient data distribution. 
- Block size (i.e., MB,NB) = 64.  
- Square processor grid.
- Use efficient machine-specific BLAS (not the Fortran 77 reference implementation BLAS) and BLACS (nondebug, BLACSDBGLVL=0 in Bmake.inc).

at Cineca: [example code](http://www.hpc.cineca.it/content/scalapack-solution-0) and [how to compile](https://hpc-forge.cineca.it/files/ScuolaCalcoloParallelo_WebDAV/public/anno-2016/12_Advanced_School/scalapack_short.pdf)

### distributing matrices
https://stackoverflow.com/questions/30167724/how-to-use-pdgemr2d-to-copy-distributed-matrix-in-total-to-all-processes/30358361#30358361
- problem with 32 bit integers when large matrices
- low performance of pdgemr2d