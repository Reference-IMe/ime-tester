#ifndef __INCLUDE_UTIL_MATRIX_H__
#define __INCLUDE_UTIL_MATRIX_H__
/*  
 * Copyright (c) 2011-2013 The University of Tennessee and The University                                                                          
 *                         of Tennessee Research Foundation.  All rights
 *                         reserved.
 *
 *
 * $COPYRIGHT$
 * 
 * Additional copyrights may follow
 * 
 * $HEADER$
 */
#include <complex.h> 

void sprintsmatrixf_ (int * m, int * n, float * A, int * ia, int * ja, int * descA, char * NAME, int * len);
void sprintsmatrix (int m, int n, float * A, int ia, int ja, int * descA, char * NAME, int len);
void sprintmatrixf_(float * A, int * descA, char * NAME, int * len);
void sprintmatrix (float * A, int * descA, char * NAME, int len);

void dprintsmatrixf_ (int * m, int * n, double * A, int * ia, int * ja, int * descA, char * NAME, int * len);
void dprintsmatrix (int m, int n, double * A, int ia, int ja, int * descA, char * NAME, int len);
void dprintmatrixf_(double * A, int * descA, char * NAME, int * len);
void dprintmatrix (double * A, int * descA, char * NAME, int len);

void cprintsmatrixf_ (int * m, int * n, float complex * A, int * ia, int * ja, int * descA, char * NAME, int * len);
void cprintsmatrix (int m, int n, float complex * A, int ia, int ja, int * descA, char * NAME, int len);
void cprintmatrixf_(float complex * A, int * descA, char * NAME, int * len);
void cprintmatrix (float complex * A, int * descA, char * NAME, int len);

void zprintsmatrixf_ (int * m, int * n, double complex * A, int * ia, int * ja, int * descA, char * NAME, int * len);
void zprintsmatrix (int m, int n, double complex * A, int ia, int ja, int * descA, char * NAME, int len);
void zprintmatrixf_(double complex * A, int * descA, char * NAME, int * len);
void zprintmatrix (double complex * A, int * descA, char * NAME, int len);

#endif /* __INCLUDE_UTIL_PRINTMATRIX_H__ */
