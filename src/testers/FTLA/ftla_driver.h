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

extern int ftla_psgetrf( int *m, int *n, float *a, int *ia, int *ja, int *descA, int *ipiv, int *info, int *ftwork );
extern int ftla_psgeqrf( int *m, int *n, float *a, int *ia, int *ja, int *descA, float *tau, float *work, int *lwork, int *info, int *ftwork );
extern int ftla_pdgetrf( int *m, int *n, double *a, int *ia, int *ja, int *descA, int *ipiv, int *info, int *ftwork );
extern int ftla_pdgeqrf( int *m, int *n, double *a, int *ia, int *ja, int *descA, double *tau, double *work, int *lwork, int *info, int *ftwork );
extern int ftla_pcgetrf( int *m, int *n, float complex *a, int *ia, int *ja, int *descA, int *ipiv, int *info, int *ftwork );
extern int ftla_pcgeqrf( int *m, int *n, float complex *a, int *ia, int *ja, int *descA, float complex *tau, float complex *work, int *lwork, int *info, int *ftwork );
extern int ftla_pzgetrf( int *m, int *n, double complex *a, int *ia, int *ja, int *descA, int *ipiv, int *info, int *ftwork );
extern int ftla_pzgeqrf( int *m, int *n, double complex *a, int *ia, int *ja, int *descA, double *tau, double *work, int *lwork, int *info, int *ftwork );
