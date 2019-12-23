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

#include "util_inject.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include "util_matrix.h"
#include "slp.h"

static int i0=0;

int *create_error_list( int N, int NB, int f, char strategy ) {
    int i;
    int *errors = malloc((3*f + 2) * sizeof(int));
    
    if( 0 == f ) errors[0] = 0;
    else errors[0] = 1;
    if( 'l' == strategy ) {
        /* inject an error for f consecutive iterations from the start 
           on processor 1,k+1 */
        for( i=0; i < f; i++ ) {
            int ii = i*3 + 1;
            errors[ii]   = NB*i+1;
            errors[ii+1] = 1;
            errors[ii+2] = i+1;
        }
        errors[f*3+1] = -1;
    }
    else if( 'r' == strategy ) {
        /* inject an error for f consecutive iterations from the N-k-1'th
           on processor 1,k+1 */
        int kmin = N/NB - f - 1;
        assert( kmin >= 0 );
        for( i=0; i < f; i++ ) {
            int ii = i*3 + 1;
            errors[ii]   = NB*(kmin+i)+1;
            errors[ii+1] = 1;
            errors[ii+2] = i+1;
        }
        errors[f*3+1] = -1;
    }
    else if( 'e' == strategy ) {
        int kspace = N/NB/(f+1);
        //printf("%d %d %d -> %d",N,NB,f+1,kspace);
        assert( kspace >= 1 );
        for( i=0; i < f; i++ ) {
            int ii = i*3 + 1;
            errors[ii]   = NB*(i+1)*kspace+1;
            errors[ii+1] = 1;
            errors[ii+2] = i+1;
        }
        errors[f*3+1] = -1;
    }
    else assert( strategy != strategy );
    
    //for( i=0; i<f; i++ ) printf("## %d (%d,%d)\t", errors[i*3+1], errors[i*3+2], errors[i*3+3]);
    return errors;
}

int sinject_errors( int *errors, int k, float *A, int *descA, ftla_work_t *ftw ) {
    int e = errors[0];
    int Perr = errors[e+1];
    int Qerr = errors[e+2];    
    int ctxt = descA[CTXT_];
    int nprow, npcol, myrow, mycol;
    Cblacs_gridinfo( ctxt, &nprow, &npcol, &myrow, &mycol );
    
    if( ftw->curop>ftw->errop && k==errors[e] ) {
        Cblacs_barrier( ctxt, "a" );
        
        //ftw->errstep = k;
        ftw->errproc = Cblacs_pnum( ctxt, Perr, Qerr );
        errors[0]+=3;
        sprintmatrix(A, descA, "Af_orig", 7); 
        if( myrow==Perr && mycol==Qerr ) {
            int np_A = numroc_( &descA[M_], &descA[NB_], &myrow, &i0, &nprow );
            int nq_A = numroc_( &descA[N_], &descA[NB_], &mycol, &i0, &npcol );                
        //    fprintf(stderr, "### Injecting a failure at iteration %d process [%d, %d] (op %d)\n", k, Perr, Qerr, ftw->curop );
            memset( A, 0, np_A*nq_A*sizeof(*A) );
            memset( ftw->pcopy.Pc, 0, (ftw->pcopy.nrPc)*(ftw->pcopy.ncPc)*sizeof(*A) );
        }
        Cblacs_barrier( ctxt, "a" );
        return -1;
    }
    return 0;
}

int dinject_errors( int *errors, int k, double *A, int *descA, ftla_work_t *ftw ) {
    int e = errors[0];
    int Perr = errors[e+1];
    int Qerr = errors[e+2];    
    int ctxt = descA[CTXT_];
    int nprow, npcol, myrow, mycol;
    Cblacs_gridinfo( ctxt, &nprow, &npcol, &myrow, &mycol );
    
    if( ftw->curop>ftw->errop && k==errors[e] ) {
        Cblacs_barrier( ctxt, "a" );
        
        //ftw->errstep = k;
        ftw->errproc = Cblacs_pnum( ctxt, Perr, Qerr );
        errors[0]+=3;
        dprintmatrix(A, descA, "Af_orig", 7); 
        if( myrow==Perr && mycol==Qerr ) {
            int np_A = numroc_( &descA[M_], &descA[NB_], &myrow, &i0, &nprow );
            int nq_A = numroc_( &descA[N_], &descA[NB_], &mycol, &i0, &npcol );                
            //printf("### Injecting a failure at iteration %d process [%d, %d] (op %d)  -ftw->curop=%d ftw->errop=%d k=%d ,errors[e]=%d \n", k, Perr, Qerr, ftw->curop, ftw->curop, ftw->errop, k,errors[e] );
            printf("### Injecting a failure at iteration %d process [%d, %d] (op %d)\n", k, Perr, Qerr, ftw->curop );
            memset( A, 0, np_A*nq_A*sizeof(*A) );
            memset( ftw->pcopy.Pc, 0, (ftw->pcopy.nrPc)*(ftw->pcopy.ncPc)*sizeof(*A) );
        }
        Cblacs_barrier( ctxt, "a" );
        return -1;
    }
    return 0;
}

int cinject_errors( int *errors, int k, float complex *A, int *descA, ftla_work_t *ftw ) {
    int e = errors[0];
    int Perr = errors[e+1];
    int Qerr = errors[e+2];    
    int ctxt = descA[CTXT_];
    int nprow, npcol, myrow, mycol;
    Cblacs_gridinfo( ctxt, &nprow, &npcol, &myrow, &mycol );
    
    if( ftw->curop>ftw->errop && k==errors[e] ) {
        Cblacs_barrier( ctxt, "a" );
        
        //ftw->errstep = k;
        ftw->errproc = Cblacs_pnum( ctxt, Perr, Qerr );
        errors[0]+=3;
        cprintmatrix(A, descA, "Af_orig", 7); 
        if( myrow==Perr && mycol==Qerr ) {
            int np_A = numroc_( &descA[M_], &descA[NB_], &myrow, &i0, &nprow );
            int nq_A = numroc_( &descA[N_], &descA[NB_], &mycol, &i0, &npcol );                
        //    fprintf(stderr, "### Injecting a failure at iteration %d process [%d, %d] (op %d)\n", k, Perr, Qerr, ftw->curop );
            memset( A, 0, np_A*nq_A*sizeof(*A) );
            memset( ftw->pcopy.Pc, 0, (ftw->pcopy.nrPc)*(ftw->pcopy.ncPc)*sizeof(*A) );
        }
        Cblacs_barrier( ctxt, "a" );
        return -1;
    }
    return 0;
}

int zinject_errors( int *errors, int k, double complex *A, int *descA, ftla_work_t *ftw ) {
    int e = errors[0];
    int Perr = errors[e+1];
    int Qerr = errors[e+2];    
    int ctxt = descA[CTXT_];
    int nprow, npcol, myrow, mycol;
    Cblacs_gridinfo( ctxt, &nprow, &npcol, &myrow, &mycol );
    
    if( ftw->curop>ftw->errop && k==errors[e] ) {
        Cblacs_barrier( ctxt, "a" );
        
        //ftw->errstep = k;
        ftw->errproc = Cblacs_pnum( ctxt, Perr, Qerr );
        errors[0]+=3;
        zprintmatrix(A, descA, "Af_orig", 7); 
        if( myrow==Perr && mycol==Qerr ) {
            int np_A = numroc_( &descA[M_], &descA[NB_], &myrow, &i0, &nprow );
            int nq_A = numroc_( &descA[N_], &descA[NB_], &mycol, &i0, &npcol );                
        //    fprintf(stderr, "### Injecting a failure at iteration %d process [%d, %d] (op %d)\n", k, Perr, Qerr, ftw->curop );
            memset( A, 0, np_A*nq_A*sizeof(*A) );
            memset( ftw->pcopy.Pc, 0, (ftw->pcopy.nrPc)*(ftw->pcopy.ncPc)*sizeof(*A) );
        }
        Cblacs_barrier( ctxt, "a" );
        return -1;
    }
    return 0;
}
