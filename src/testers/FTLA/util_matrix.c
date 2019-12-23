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

/**********************************************
 **   print out matrix for check             **
 *********************************************/

#include "util_matrix.h"
#include <stdlib.h>
#include <string.h>
#include "slp.h"

//#define PRINTMATRIX

#ifdef PRINTMATRIX
static int i0=0;
#endif 

void sprintsmatrixf_ (int * m, int * n, float * A, int * ia, int * ja, int * descA, char * NAME, int * len) {
#ifdef PRINTMATRIX
    int out=0; /* 0=stderr, 6=stdout */
    int nb = descA[NB_];
    typeof(A) work = malloc(nb * sizeof(*A));
    pslaprnt_ (m, n, A, ia, ja, descA, &i0, &i0, NAME, &out, work, *len);
    free (work);
#endif
}
void sprintsmatrix (int m, int n, float * A, int ia, int ja, int * descA, char * NAME, int len) {
#ifdef PRINTMATRIX
    sprintsmatrixf_(&m, &n, A, &ia, &ja, descA, NAME, &len);
#endif
}
void sprintmatrixf_(float * A, int * descA, char * NAME, int * len) {
#ifdef PRINTMATRIX
    sprintmatrix (A, descA, NAME, *len);
#endif
}
void sprintmatrix (float * A, int * descA, char * NAME, int len) {
#ifdef PRINTMATRIX
    int m = descA[M_], n = descA[N_];
    sprintsmatrix (m, n, A, 1, 1, descA, NAME, len);
#endif
}

void dprintsmatrixf_ (int * m, int * n, double * A, int * ia, int * ja, int * descA, char * NAME, int * len) {
#ifdef PRINTMATRIX
    int out=0; /* 0=stderr, 6=stdout */
    int nb = descA[NB_];
    typeof(A) work = malloc(nb * sizeof(*A));
    pdlaprnt_ (m, n, A, ia, ja, descA, &i0, &i0, NAME, &out, work, *len);
    free (work);
#endif
}
void dprintsmatrix (int m, int n, double * A, int ia, int ja, int * descA, char * NAME, int len) {
#ifdef PRINTMATRIX
    dprintsmatrixf_(&m, &n, A, &ia, &ja, descA, NAME, &len);
#endif
}
void dprintmatrixf_(double * A, int * descA, char * NAME, int * len) {
#ifdef PRINTMATRIX
    dprintmatrix (A, descA, NAME, *len);
#endif
}
void dprintmatrix (double * A, int * descA, char * NAME, int len) {
#ifdef PRINTMATRIX
    int m = descA[M_], n = descA[N_];
    dprintsmatrix (m, n, A, 1, 1, descA, NAME, len);
#endif
}

void cprintsmatrixf_ (int * m, int * n, float complex * A, int * ia, int * ja, int * descA, char * NAME, int * len) {
#ifdef PRINTMATRIX
    int out=0; /* 0=stderr, 6=stdout */
    int nb = descA[NB_];
    typeof(A) work = malloc(nb * sizeof(*A));
    pclaprnt_ (m, n, A, ia, ja, descA, &i0, &i0, NAME, &out, work, *len);
    free (work);
#endif
}
void cprintsmatrix (int m, int n, float complex * A, int ia, int ja, int * descA, char * NAME, int len) {
#ifdef PRINTMATRIX
    cprintsmatrixf_(&m, &n, A, &ia, &ja, descA, NAME, &len);
#endif
}
void cprintmatrixf_(float complex * A, int * descA, char * NAME, int * len) {
#ifdef PRINTMATRIX
    cprintmatrix (A, descA, NAME, *len);
#endif
}
void cprintmatrix (float complex * A, int * descA, char * NAME, int len) {
#ifdef PRINTMATRIX
    int m = descA[M_], n = descA[N_];
    cprintsmatrix (m, n, A, 1, 1, descA, NAME, len);
#endif
}

void zprintsmatrixf_ (int * m, int * n, double complex * A, int * ia, int * ja, int * descA, char * NAME, int * len) {
#ifdef PRINTMATRIX
    int out=0; /* 0=stderr, 6=stdout */
    int nb = descA[NB_];
    typeof(A) work = malloc(nb * sizeof(*A));
    pzlaprnt_ (m, n, A, ia, ja, descA, &i0, &i0, NAME, &out, work, *len);
    free (work);
#endif
}
void zprintsmatrix (int m, int n, double complex * A, int ia, int ja, int * descA, char * NAME, int len) {
#ifdef PRINTMATRIX
    zprintsmatrixf_(&m, &n, A, &ia, &ja, descA, NAME, &len);
#endif
}
void zprintmatrixf_(double complex * A, int * descA, char * NAME, int * len) {
#ifdef PRINTMATRIX
    zprintmatrix (A, descA, NAME, *len);
#endif
}
void zprintmatrix (double complex * A, int * descA, char * NAME, int len) {
#ifdef PRINTMATRIX
    int m = descA[M_], n = descA[N_];
    zprintsmatrix (m, n, A, 1, 1, descA, NAME, len);
#endif
}

