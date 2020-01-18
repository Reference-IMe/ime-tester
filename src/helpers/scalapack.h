/*
 * scalapack.h
 *
 *  Created on: Feb 3, 2019
 *      Author: marcello
 */

// reference: http://icl.cs.utk.edu/lapack-forum/viewtopic.php?t=271#
//            http://aragorn.pb.bialystok.pl/~mars/tutorial/scalapack/

#ifndef __SCALAPACK_H__
#define __SCALAPACK_H__

#define SCALAPACKNB 4 // blocking factor

extern void descinit_( int* desc, int* m, int* n, int* mb, int* nb, int* irsrc, int* icsrc, int* ictxt, int* lld, int* info);
extern int  numroc_( int* n, int* nb, int* iproc, int* isrcproc, int* nprocs);
extern void pdgemr2d_ (int* m , int* n , double* a , int* ia , int* ja , int* desca , double* b , int* ib , int* jb , int* descb , int* ictxt );
extern void pigemr2d_ (int* m , int* n , int* a , int* ia , int* ja , int* desca , int* b , int* ib , int* jb , int* descb , int* ictxt );
extern void pdgesv_( int* n, int* nrhs, double* A, int* ia, int* ja, int* desca, int* ipiv, double* B, int* ib, int* jb, int* descb, int* info);
extern void pdgetrf_ ( int *m, int *n, double *A, int *ia, int *ja, int *descA, int *ipiv, int *info );
extern void pdgetrs_( char *trans, int*n, int *nrhs, double *A, int *ia, int *ja, int *descA, int *ipiv, double *B, int *ib, int *jb, int *descB, int *info );
extern void pdtran_( int *m, int *n, double *ALPHA, double *A, int *ia, int *ja, int *descA, double *BETA, double *B, int *ib, int *jb, int *descB);

extern void pdgeqrf_ ( int *m, int *n, double *A, int *ia, int *ja, int *descA, double *tau, double *work, int *lwork, int *info );


// ScaLAPACK modified
extern void pdgetrf_cp_ ( int *m, int *n, double *A, int *ia, int *ja, int *descA, double *Acp, int *iacp, int *jacp, int *descAcp, int *ipiv, int *ipivcp, int *nipiv, int *cpfreq, int *jfault, int* ictxt, int *info );
extern void pdgeqrf_cp_ ( int *m, int *n, double *A, int *ia, int *ja, int *descA, double *Acp, int *iacp, int *jacp, int *descAcp, double *tau, double *taucp, int* ltau, double *work, double *workcp, int *lwork, int *cpfreq, int *jfault, int* ictxt, int *info );
//extern void pdgetf2_cp_ (int *m, int *n, double *A, int *ia, int *ja, int *descA, int *ipiv, int *info );

// FTLA
extern void pdmatgen_ (int *, char *, char *, int *, int *, int *, int *, double *, int *, int *, int *, int *, int *, int *, int *, int *, int *, int *, int *, int *);
extern void pdgeqrrv_ (int *m, int *n, double *A, int *ia, int *ja, int *descA, double*tau, double *work);

#endif /* __SCALAPACK_H__ */
