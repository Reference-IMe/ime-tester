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

extern void descinit_(int* desc, int* m, int* n, int* mb, int* nb, int* irsrc, int* icsrc, int* ictxt, int* lld, int* info);
extern int  numroc_  (int* n, int* nb, int* iproc, int* isrcproc, int* nprocs);
extern int  indxg2l_ (int* indxglob, int* nb, int* iproc, int* isrcproc, int* nprocs);
extern int  indxl2g_ (int* indxloc, int* nb, int* iproc, int* isrcproc, int* nprocs);
extern int  indxg2p_ (int* indxglob, int* nb, int* iproc, int* isrcproc, int* nprocs);

extern void pigemr2d_(int* m, int* n, int* A, int* ia, int* ja, int* descA, int* B, int* ib, int* jb, int* descB, int* ictxt);

extern void pdgemr2d_(int* m, int* n, double* A, int* ia, int* ja, int* descA, double* B, int* ib, int* jb, int* descB, int* ictxt);
extern void pdgesv_  (int* n, int* nrhs, double* A, int* ia, int* ja, int* descA, int* ipiv, double* B, int* ib, int* jb, int* descB, int* info);
extern void pdgetrf_ (int* m, int* n, double* A, int* ia, int* ja, int* descA, int* ipiv, int* info);
extern void pdgetrs_ (char* trans, int*n, int* nrhs, double* A, int* ia, int* ja, int* descA, int* ipiv, double* B, int* ib, int* jb, int* descB, int* info);
extern void pdtran_  (int* m, int* n, double* alpha, double* A, int* ia, int* ja, int* descA, double* beta, double* B, int* ib, int* jb, int* descB);
extern void pdgeqrf_ (int* m, int* n, double* A, int* ia, int* ja, int* descA, double* tau, double* work, int* lwork, int* info);
extern void pdlaset_ (char* uplo, int* m, int* n, double* alpha, double* beta, double* A, int* ia, int* ja, int* descA);
extern void pdgemm_  (char* transa, char* transb, int* m, int* n, int* k, double* alpha, double* A, int* ia, int* ja, int* descA, double* B, int* ib, int* jb, int* descB, double* beta, double* C, int* ic, int* jc, int* descC);
extern void pdorgqr_ (int* m, int* n, int* k, double* A, int* ia, int* ja, int* descA, double* tau, double* work, int* lwork, int* info);
extern void pdgesvd_ (char* jobu, char* jobvt, int* m, int* n, double* A, int* ia, int* ja, int* desAa, double* S, double* U, int* iu, int* ju, int* descU, double* VT, int* ivt, int* jvt, int* descVT, double* work, int* lwork, int* info);
extern void pdger_ ( int* m, int* n, double* alpha, double* X, int *ix, int* jx, int* descX, int* incx, double* Y, int* iy, int* jy, int* descY, int* incy, double* A, int* ia, int* ja, int* descA );
extern void pdlacpy_ ( char* uplo, int* m, int* n, double* A, int* ia, int* ja, int* descA, double* B, int* ib, int* jb, int* descB);
//extern void pdtrsm_ (char* side, char* uplo, char* transa, char* diag, int* m, int* n, double* alpha, double* A, int* ia, int* ja, int* descA, double* B, int* ib, int* jb, int* descB);
extern void  pdtrmm_ (char* side, char* uplo, char* trans, char* diag, int* m, int* n, double* alpha, double* A, int* ia, int* ja, int* descA, double* B, int* ib, int* jb, int* descB);

extern void psgemr2d_(int* m, int* n, float* A, int* ia, int* ja, int* descA, float* B, int* ib, int* jb, int* descB, int* ictxt);
extern void psgesv_  (int* n, int* nrhs, float* A, int* ia, int* ja, int* descA, int* ipiv, float* B, int* ib, int* jb, int* descB, int* info);
extern void psgetrf_ (int* m, int* n, float* A, int* ia, int* ja, int* descA, int* ipiv, int* info);
extern void psgetrs_ (char* trans, int*n, int* nrhs, float* A, int* ia, int* ja, int* descA, int* ipiv, float* B, int* ib, int* jb, int* descB, int* info);
extern void pstran_  (int* m, int* n, float* alpha, float* A, int* ia, int* ja, int* descA, float* beta, float* B, int* ib, int* jb, int* descB);
extern void psgeqrf_ (int* m, int* n, float* A, int* ia, int* ja, int* descA, float* tau, float* work, int* lwork, int* info);
extern void pslaset_ (char* uplo, int* m, int* n, float* alpha, float* beta, float* A, int* ia, int* ja, int* descA);
extern void psgemm_  (char* transa, char* transb, int* m, int* n, int* k, float* alpha, float* A, int* ia, int* ja, int* descA, float* B, int* ib, int* jb, int* descB, float* beta, float* C, int* ic, int* jc, int* descC);
extern void psorgqr_ (int* m, int* n, int* k, float* A, int* ia, int* ja, int* descA, float* tau, float* work, int* lwork, int* info);
extern void psgesvd_ (char* jobu, char* jobvt, int* m, int* n, float* A, int* ia, int* ja, int* desAa, float* S, float* U, int* iu, int* ju, int* descU, float* VT, int* ivt, int* jvt, int* descVT, float* work, int* lwork, int* info);
extern void psger_ ( int* m, int* n, float* alpha, float* X, int *ix, int* jx, int* descX, int* incx, float* Y, int* iy, int* jy, int* descY, int* incy, float* A, int* ia, int* ja, int* descA );
extern void pslacpy_ ( char* uplo, int* m, int* n, float* A, int* ia, int* ja, int* descA, float* B, int* ib, int* jb, int* descB);
//extern void pstrsm_ (char* side, char* uplo, char* transa, char* diag, int* m, int* n, float* alpha, float* A, int* ia, int* ja, int* descA, float* B, int* ib, int* jb, int* descB);
extern void  pstrmm_ (char* side, char* uplo, char* trans, char* diag, int* m, int* n, float* alpha, float* A, int* ia, int* ja, int* descA, float* B, int* ib, int* jb, int* descB);


// ScaLAPACK modified
extern void pdgetrf_cp_ (int* m, int* n, double* A, int* ia, int* ja, int* descA, double* Acp, int* iacp, int* jacp, int* descAcp, int* ipiv, int* ipivcp, int* nipiv, int* cpfreq, int* jfault, int* ictxt, int* info);
extern void pdgeqrf_cp_ (int* m, int* n, double* A, int* ia, int* ja, int* descA, double* Acp, int* iacp, int* jacp, int* descAcp, double* tau, double* taucp, int* ltau, double* work, double* workcp, int* lwork, int* cpfreq, int* jfault, int* ictxt, int* info);
extern void pdgetrf_cs_ (int* m, int* n, double* A, int* ia, int* ja, int* descA, double* Acp, int* iacp, int* jacp, int* descAcp, int* ipiv, int* ipivcp, int* nipiv, int* cpfreq, int* jfault, int* ictxt, int* info);
extern void pdgetrf_cpx_ (int* m, int* n, double* A, int* ia, int* ja, int* descA, double* Acp, int* iacp, int* jacp, int* descAcp, int* ipiv, int* ipivcp, int* nipiv, int* cpfreq, int* faultnum, int* faultlist, int* jfault, int* recovery, int* ictxt, int* cprocs, int* info);
extern void pdgesvnopiv_  (int* n, int* nrhs, double* A, int* ia, int* ja, int* descA, int* ipiv, double* B, int* ib, int* jb, int* descB, int* info);

extern void psgetrf_cp_ (int* m, int* n, float* A, int* ia, int* ja, int* descA, float* Acp, int* iacp, int* jacp, int* descAcp, int* ipiv, int* ipivcp, int* nipiv, int* cpfreq, int* jfault, int* ictxt, int* info);
extern void psgetrf_cpx_ (int* m, int* n, float* A, int* ia, int* ja, int* descA, float* Acp, int* iacp, int* jacp, int* descAcp, int* ipiv, int* ipivcp, int* nipiv, int* cpfreq, int* faultnum, int* faultlist, int* jfault, int* recovery, int* ictxt, int* cprocs, int* info);


// FTLA
extern void pdmatgen_ (int*, char*, char*, int*, int*, int*, int*, double*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*);
extern void pdgeqrrv_ (int* m, int* n, double* A, int* ia, int* ja, int* descA, double* tau, double* work);
extern int  pdtrsm_ (char* side, char* uplo, char* trans, char* diag, int* m, int* n, double* alpha, double* A, int* ia, int* ja, int* descA, double* B, int* ib, int* jb, int* descB);
extern void pdormqr_(char* side, char* trans, int* m, int* n, int* k, double* A, int* ia, int* ja, int* descA, double* tau, double* C, int* ic, int* jc, int* descC, double* work, int* lwork, int* info);

#endif /* __SCALAPACK_H__ */
