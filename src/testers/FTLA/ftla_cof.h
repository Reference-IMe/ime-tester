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

int Cftla_cof_sgeqrw( float *A, int *descA, float *tau, float *work, int lwork, ftla_work_t *ftwork );
int ftla_cof_sgeqrw_( float *A, int *descA, float *tau, float *work, int *lwork, int *ftwork );
int Cftla_cof_sgeqrr( float *A, int *descA, float *tau, float *work, int lwork, ftla_work_t *ftwork );
int ftla_cof_sgeqrr_( float *A, int *descA, float *tau, float *work, int *lwork, int *ftwork );

int Cftla_cof_dgeqrw( double *A, int *descA, double *tau, double *work, int lwork, ftla_work_t *ftwork );
int ftla_cof_dgeqrw_( double *A, int *descA, double *tau, double *work, int *lwork, int *ftwork );
int Cftla_cof_dgeqrr( double *A, int *descA, double *tau, double *work, int lwork, ftla_work_t *ftwork );
int ftla_cof_dgeqrr_( double *A, int *descA, double *tau, double *work, int *lwork, int *ftwork );

int Cftla_cof_cgeqrw( float *A, int *descA, float *tau, float *work, int lwork, ftla_work_t *ftwork );
int ftla_cof_cgeqrw_( float *A, int *descA, float *tau, float *work, int *lwork, int *ftwork );
int Cftla_cof_cgeqrr( float *A, int *descA, float *tau, float *work, int lwork, ftla_work_t *ftwork );
int ftla_cof_cgeqrr_( float *A, int *descA, float *tau, float *work, int *lwork, int *ftwork );

int Cftla_cof_zgeqrw( double *A, int *descA, double *tau, double *work, int lwork, ftla_work_t *ftwork );
int ftla_cof_zgeqrw_( double *A, int *descA, double *tau, double *work, int *lwork, int *ftwork );
int Cftla_cof_zgeqrr( double *A, int *descA, double *tau, double *work, int lwork, ftla_work_t *ftwork );
int ftla_cof_zgeqrr_( double *A, int *descA, double *tau, double *work, int *lwork, int *ftwork );

void Cftla_cof_cleanup( ftla_work_t *ftwork );
void ftla_cof_cleanup_( int *ftwork );
