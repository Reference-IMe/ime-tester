#ifndef __INCLUDE_FTLA_FTWORK_H__
#define __INCLUDE_FTLA_FTWORK_H__
/*  
 * Copyright (c) 2013-2013 The University of Tennessee and The University                                                                          
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

#include <mpi.h>

typedef struct {
    MPI_Comm rowcomm;
    MPI_Comm colcomm;
    void *S;
    void *R;  
} ftla_csum_comms_t;

typedef struct {
    int descPc[9];
    int nrPc;
    int ncPc;
    void *Pc;
} ftla_panel_lcopy_t;

typedef struct {
    int errstep; /* when 0, no error needs to be recovered */
    int errproc;
    int curop;
    int errop;
    int descA[9];
    int Mc;
    int Nc;
    ftla_panel_lcopy_t pcopy;
    ftla_csum_comms_t comms;
} ftla_work_t;

extern void Cftla_work_construct( int reentry, int *descA, int Mc, int Nc, ftla_work_t *ftwork );
extern void ftla_work_construct_( int * reentry, int *descA, int *Mc, int *Nc, int *ftwork );
extern void Cftla_work_destruct( ftla_work_t *ftwork );
extern void ftla_work_destruct_( int *ftwork );

extern int  Cftla_procfailed( ftla_work_t *ftwork, int je );
extern int  ftla_procfailed_( int *ftwork, int *je );
extern int  Cftla_replay( ftla_work_t *ftwork );
extern int  ftla_replay_( int *ftwork );
extern int  Cftla_repairstep( ftla_work_t *ftwork );
extern int  ftla_repairstep_( int *ftwork );

#endif
