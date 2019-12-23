#ifndef __INCLUDE_UTIL_INJECT_H__
#define __INCLUDE_UTIL_INJECT_H__
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

#include "ftla_ftwork.h"
#include <complex.h>

int *create_error_list( int N, int NB, int f, char strategy );
int sinject_errors( int *errors, int k, float *A, int *descA, ftla_work_t *ftw );
int dinject_errors( int *errors, int k, double *A, int *descA, ftla_work_t *ftw );
int cinject_errors( int *errors, int k, float complex *A, int *descA, ftla_work_t *ftw );
int zinject_errors( int *errors, int k, double complex *A, int *descA, ftla_work_t *ftw );

#endif /* __INCLUDE_UTIL_INJECT_H__ */
