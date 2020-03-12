/*
 * lapack.h
 *
 *  Created on: Feb 1, 2020
 *      Author: marcello
 */


#ifndef __LAPACK_H__
#define __LAPACK_H__

extern void dgesv_( int* n, int* nrhs, double* a, int* lda, int* ipiv, double* b, int* ldb, int* info );
extern void dgeqrfp_(int* m, int* n, double* a, int* lda, double* tau, double* work, int* lwork, int* info);
extern void dorgqr_(int* m, int* n, int* k, double* a, int* lda, double* tau, double* work, int* lwork, int* info);
extern void dgemm_(char* transA, char* transB, int* m, int* n, int* k, double* alpha, double* A, int* lda, double* B, int* ldb, double* beta, double* C, int* ldc);
extern void dgesvd_( char* jobu, char* jobvt, int* m, int* n, double* a,
                int* lda, double* s, double* u, int* ldu, double* vt, int* ldvt,
                double* work, int* lwork, int* info );
#endif /* __LAPACK_H__ */
