/*
 * lapack.h
 *
 *  Created on: Feb 1, 2020
 *      Author: marcello
 */


#ifndef __LAPACK_H__
#define __LAPACK_H__

extern void dgesv_  (int* n, int* nrhs, double* A, int* lda, int* ipiv, double* B, int* ldb, int* info );
extern void dgeqrfp_(int* m, int* n, double* A, int* lda, double* tau, double* work, int* lwork, int* info);
extern void dorgqr_ (int* m, int* n, int* k, double* A, int* lda, double* tau, double* work, int* lwork, int* info);
extern void dgemm_  (char* transA, char* transB, int* m, int* n, int* k, double* alpha, double* A, int* lda, double* B, int* ldb, double* beta, double* C, int* ldc);
extern void dgesvd_ (char* jobu, char* jobvt, int* m, int* n, double* A, int* lda, double* S, double* U, int* ldu, double* Vt, int* ldvt, double* work, int* lwork, int* info);
extern void dgetrf_ (int* m, int *n, double* A, int* lda, int* ipiv, int* info);
extern void dgetri_ (int *n, double* A, int* lda, int* ipiv, double* work, int* lwork, int* info);

extern double dlange_(char *norm, int* m, int* n, double* A, int* lda, double* work);

#endif /* __LAPACK_H__ */
