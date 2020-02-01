/*
 * lapack.h
 *
 *  Created on: Feb 1, 2020
 *      Author: marcello
 */


#ifndef __LAPACK_H__
#define __LAPACK_H__

extern void dgesv( int* n, int* nrhs, double* a, int* lda, int* ipiv, double* b, int* ldb, int* info );

#endif /* __LAPACK_H__ */
