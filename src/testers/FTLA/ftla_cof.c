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


#include "ftla_cof.h"
#include <stdio.h>
#include <unistd.h>
#include <sys/errno.h>
#include <limits.h>
#include "slp.h"

#define COF_FTLA_BASENAME "ftla_cof_"

static int i0=0;

static int cof_geqrw( int sdcz, void *A, int *descA, void *tau, void *work, int lwork, ftla_work_t *ftw ) {
    // grid parameters
    int ctxt = descA[CTXT_];
    int nprow, npcol, myrow, mycol;
    Cblacs_gridinfo( ctxt, &nprow, &npcol, &myrow, &mycol );
    int   nb = descA[NB_];
    int np_A = numroc_( &descA[M_], &nb, &myrow, &i0, &nprow );
    int nq_A = numroc_( &descA[N_], &nb, &mycol, &i0, &npcol );

    FILE * pFile;
    char filename[PATH_MAX];
    snprintf( filename, PATH_MAX, COF_FTLA_BASENAME "%d.%d.ckpt", myrow, mycol );
    pFile = fopen ( filename , "rb" );
    if( pFile == NULL ) {
        if( errno != ENOENT ) {
            perror( "ftla:cof_load:fopen" );
        }
        return -errno;
    }
    //printf ("(%d,%d): loading from %s\n", myrow, mycol, filename);
    // load A
    fread( A, sdcz, np_A*nq_A, pFile );
    // read local copy
    fread( ftw->pcopy.Pc, sdcz, ftw->pcopy.nrPc*ftw->pcopy.ncPc, pFile );
    // load tau 
    fread( tau, sdcz, nq_A, pFile );
    // load work 
    fread( work, sdcz, lwork, pFile );
    fclose( pFile );
    return 0;
}

static int cof_geqrr( int sdcz, void *A, int *descA, void *tau, void *work, int lwork, ftla_work_t *ftw ) {
    // grid parameters
    int ctxt = descA[CTXT_];
    int nprow, npcol, myrow, mycol;
    Cblacs_gridinfo( ctxt, &nprow, &npcol, &myrow, &mycol );
    int   nb = descA[NB_];
    int np_A = numroc_( &descA[M_], &nb, &myrow, &i0, &nprow );
    int nq_A = numroc_( &descA[N_], &nb, &mycol, &i0, &npcol );
    
    FILE * pFile;
    char filename[PATH_MAX];
    snprintf( filename, PATH_MAX, COF_FTLA_BASENAME "%d.%d.ckpt", myrow, mycol );
    //printf ("(%d,%d): dumping to %s\n", myrow, mycol, filename);
    pFile = fopen ( filename , "wb" );
    if( pFile == NULL ) {
        perror( "ftla:cof_save:fopen" );
        return -errno;
    }
    //printf ("(%d,%d): saving to %s\n", myrow, mycol, filename);
    // dump A
    fwrite( A, sdcz, np_A*nq_A, pFile );
    // dump local copy
    fwrite( ftw->pcopy.Pc, sdcz, ftw->pcopy.nrPc*ftw->pcopy.ncPc, pFile );
    // dump tau 
    fwrite( tau, sdcz, nq_A, pFile );
    // dump work 
    fwrite( work, sdcz, lwork, pFile );
    fclose( pFile );
    return 0;
}

int Cftla_cof_sgeqrw( float *A, int *descA, float *tau, float *work, int lwork, ftla_work_t *ftw ) {
    return cof_geqrw( sizeof(float), A, descA, tau, work, lwork, ftw );
}
int Cftla_cof_sgeqrr( float *A, int *descA, float *tau, float *work, int lwork, ftla_work_t *ftw ) {
    return cof_geqrr( sizeof(float), A, descA, tau, work, lwork, ftw );
}
int Cftla_cof_dgeqrw( double *A, int *descA, double *tau, double *work, int lwork, ftla_work_t *ftw ) {
    return cof_geqrw( sizeof(double), A, descA, tau, work, lwork, ftw );
}
int Cftla_cof_dgeqrr( double *A, int *descA, double *tau, double *work, int lwork, ftla_work_t *ftw ) {
    return cof_geqrr( sizeof(double), A, descA, tau, work, lwork, ftw );
}
int Cftla_cof_cgeqrw( float *A, int *descA, float *tau, float *work, int lwork, ftla_work_t *ftw ) {
    return cof_geqrw( 2*sizeof(float), A, descA, tau, work, lwork, ftw );
}
int Cftla_cof_cgeqrr( float *A, int *descA, float *tau, float *work, int lwork, ftla_work_t *ftw ) {
    return cof_geqrr( 2*sizeof(float), A, descA, tau, work, lwork, ftw );
}
int Cftla_cof_zgeqrw( double *A, int *descA, double *tau, double *work, int lwork, ftla_work_t *ftw ) {
    return cof_geqrw( 2*sizeof(double), A, descA, tau, work, lwork, ftw );
}
int Cftla_cof_zgeqrr( double *A, int *descA, double *tau, double *work, int lwork, ftla_work_t *ftw ) {
    return cof_geqrr( 2*sizeof(double), A, descA, tau, work, lwork, ftw );
}

int ftla_cof_sgeqrw_( float *A, int *descA, float *tau, float *work, int *lwork, int *ftw ) {
    return cof_geqrw( sizeof(float), A, descA, tau, work, *lwork, (ftla_work_t*)ftw );
}
int ftla_cof_sgeqrr_( float *A, int *descA, float *tau, float *work, int *lwork, int *ftw ) {
    return cof_geqrr( sizeof(float), A, descA, tau, work, *lwork, (ftla_work_t*)ftw );
}
int ftla_cof_dgeqrw_( double *A, int *descA, double *tau, double *work, int *lwork, int *ftw ) {
    return cof_geqrw( sizeof(double), A, descA, tau, work, *lwork, (ftla_work_t*)ftw );
}
int ftla_cof_dgeqrr_( double *A, int *descA, double *tau, double *work, int *lwork, int *ftw ) {
    return cof_geqrr( sizeof(double), A, descA, tau, work, *lwork, (ftla_work_t*)ftw );
}
int ftla_cof_cgeqrw_( float *A, int *descA, float *tau, float *work, int *lwork, int *ftw ) {
    return cof_geqrw( 2*sizeof(float), A, descA, tau, work, *lwork, (ftla_work_t*)ftw );
}
int ftla_cof_cgeqrr_( float *A, int *descA, float *tau, float *work, int *lwork, int *ftw ) {
    return cof_geqrr( 2*sizeof(float), A, descA, tau, work, *lwork, (ftla_work_t*)ftw );
}
int ftla_cof_zgeqrw_( double *A, int *descA, double *tau, double *work, int *lwork, int *ftw ) {
    return cof_geqrw( 2*sizeof(double), A, descA, tau, work, *lwork, (ftla_work_t*)ftw );
}
int ftla_cof_zgeqrr_( double *A, int *descA, double *tau, double *work, int *lwork, int *ftw ) {
    return cof_geqrr( 2*sizeof(double), A, descA, tau, work, *lwork, (ftla_work_t*)ftw );
}


void Cftla_cof_cleanup( ftla_work_t *ftw ) {
    int ctxt = ftw->descA[CTXT_];
    int nprow, npcol, myrow, mycol;
    char filename[PATH_MAX];

    Cblacs_gridinfo( ctxt, &nprow, &npcol, &myrow, &mycol );
    snprintf( filename, PATH_MAX, COF_FTLA_BASENAME "%d.%d.ckpt", myrow, mycol );
    unlink( filename );
}
void ftla_cof_cleanup_( int *ftw ) {
    Cftla_cof_cleanup( (ftla_work_t*)ftw );
}
