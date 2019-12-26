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

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "../../helpers/Cblacs.h"
#include "../../helpers/scalapack.h"
#include "../FTLA/util_matrix.h"
#include "../FTLA/util_inject.h"
#include "../FTLA/ftla_ftwork.h"
#include "../FTLA/ftla_cof.h"
#include "../FTLA/ftla_driver.h"

static int i0=0, i1=1;
static double p0=0e0, p1=1e0, m1=-1e0;
static double eps;
static int rank, np;

#define DTYPE_      0
#define CTXT_       1
#define M_          2
#define N_          3
#define MB_         4
#define NB_         5
#define RSRC_       6
#define CSRC_       7
#define LLD_        8
#define TAG_L    1021
#define TAG_LL    904


static void create_matrix (int ctxt, int seed, double **A, int *descA, int M, int N, int nb, int *np_A, int *nq_A);

int *errors;

//int test_FTLA_ftdqr_new(int P, int Q, int Fmin, int Fmax, int Finc)
int calc_FTLA_ftdqr_new(int rows, double* A_global, int P, int Q, int Fmin, int Fmax, int Finc)
{
    int ictxt, ictxt_global, info;
    int PxQ = P * Q;
    int nprow, npcol, myrow, mycol;
    ftla_work_t ftwork;

    int NB=SCALAPACKNB;
    int M, N, Nc, Ne, S=1;
    int Fstrat='e', F; // Fmin=0, Fmax=0, Finc=1;
    double *A0=NULL, *A=NULL;    
    int descA[9], descA_global[9];
    double *B0=NULL, *X0=NULL, *Xo=NULL, *Xf=NULL;
    int descX[9];

    double t1, t2, To, Tf;
    double Rom, Rfm, Rob, Rfb, Rox, Rfx, KA;
    int err = 0;
    int i, start, end, step;
    start = end = step = rows;

    {/* init BLACS */
        Cblacs_pinfo( &rank, &np );
        //print0( ("Setting up MPI and BLACS...") );
        Cblacs_get( -1, 0, &ictxt );
        Cblacs_gridinit( &ictxt, "Row", P, Q );
        Cblacs_gridinfo( ictxt, &nprow, &npcol, &myrow, &mycol );
        //print0( ("done\tPxQ=%dx%d, NB=%d\n", nprow, npcol, NB) );
        //eps = pdlamch_( &ictxt, "Epsilon" );
        Cblacs_get( -1, 0, &ictxt_global );
        Cblacs_gridinit( &ictxt_global, "Row", i1, i1 );
    }

    //print0( ("M\tNe\t : time(s)\tgflops/proc\t(realflops)\t: |A-A~|\t|B-AX~|\t\t|X-X~|\n") );

    for( i = start; i <= end; i += step) {
        {/* allocate matrices */
            /* determine checksum size, generate A matrix */
            N = M = i;
            Nc = numroc_( &N, &NB, &mycol, &i0, &npcol ); //LOCc(N_A) 
            MPI_Allreduce( MPI_IN_PLACE, &Nc, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD );
#ifndef NO_EXTRAFLOPS
            Ne = N + Nc*2;// + ((Nc/NB)%Q==0)*NB;
#else
            Ne = N;
#endif
            //print0( ("%-7d\t%-7d\t", M, Ne) );

            int n=N;
            int nb=NB;
            int zero=i0;
            int one=i1;
            int lld_global;
            int nr;
            nr = numroc_( &N, &NB, &myrow, &i0, &nprow ); //LOCr(N_A)
            int nce = numroc_( &Ne, &NB, &mycol, &i0, &npcol );
            int lld=MAX( 1 , nr );
            int nc = Nc;

            //A = malloc(nr*nce*sizeof(double));

        	// Descriptors (local)
        	//descinit_( descA, &n, &nce, &nb, &nb, &zero, &zero, &ictxt, &lld, &info );

        	if (rank==0)
        	{
        		// Descriptors (global)
        		lld_global = n;
        		descinit_( descA_global, &n, &n, &one, &one, &zero, &zero, &ictxt_global, &lld_global, &info );
        	}
        	else
        	{
        		// Descriptors (global, for non-root nodes)
        		for (i=0; i<9; i++)
        		{
        			descA_global[i]=0;
        		}
        		descA_global[1]=-1;
        	}

            //create_matrix( ictxt, 800, &A0, descA, M, Ne, NB, NULL, NULL );
            //create_matrix( ictxt, 0,   &A,  NULL, M, Ne, NB, NULL, NULL );
            create_matrix( ictxt, 0,   &A,  descA, M, Ne, NB, NULL, NULL );
            //pdlacpy_( "All", &M, &N, A0, &i1, &i1, descA, A, &i1, &i1, descA );
            //printf("\n fin qua \n");
        	// spread matrices
			pdgemr2d_(&n, &n, A_global, &one, &one, descA_global, A, &one, &one, descA, &ictxt);

            /* allocate local buffer for the Q-wide local panel copy */
            create_matrix( ictxt, 0, (typeof(&A))&(ftwork.pcopy.Pc), ftwork.pcopy.descPc, M, (Q+2)*NB, NB, &(ftwork.pcopy.nrPc), &(ftwork.pcopy.ncPc) );
        }

        {/* call resilient QR */
            int err=0;
            int lwork=-1;
            double lazywork;
            pdgeqrf_( &M, &Ne, NULL, &i1, &i1, descA, NULL, &lazywork, &lwork, &info );
            lwork = (int)lazywork;
            double *work = (double*)malloc( lwork*sizeof(double) );
            double *tau  = (double*)malloc( Nc*sizeof(double) );
            //printf("\n fin qui \n");
#ifdef INJECT        
          for( F = Fmin; F<=Fmax; F+=Finc ) {
            errors = create_error_list( M, NB, F, Fstrat );
#endif

            Cftla_work_construct( 0, descA, 0, Ne-N, &ftwork );

            do { // call ftpdgeqrf until we complete w/o a failure
#ifdef USE_CoF
              if(err) Cftla_cof_dgeqrr( A, descA, tau, work, lwork, &ftwork );
#endif
              //printf("\n fin quo \n");
              err = ftla_pdgeqrf( &M, &Ne, A, &i1, &i1, descA, tau, work, &lwork, &info, (int*)&ftwork );

#ifdef USE_CoF
              if(err) Cftla_cof_dgeqrr( A, descA, tau, work, lwork, &ftwork );
#endif
            } while(err);

            //checkerror( info, 0 );

            pdgemr2d_ (&N, &N, A, &i1, &i1, descA, A_global, &i1, &i1, descA_global, &ictxt);

            Cftla_cof_cleanup( &ftwork );
            Cftla_work_destruct( &ftwork );
            free( errors );

#ifdef INJECT
          }
#endif
            free( work );
            free( tau );
        }
        {/* Cleanup */
            if( NULL != A0 ) free( A0 );
            if( NULL != A  ) free( A );
            A0 = A = NULL;
            if( NULL != B0 ) free( B0 );
            if( NULL != X0 ) free( X0 );
            if( NULL != Xo ) free( Xo );
            if( NULL != Xf ) free( Xf );
            B0 = X0 = Xo = Xf = NULL;

            if( NULL != ftwork.pcopy.Pc ) free( ftwork.pcopy.Pc);
            ftwork.pcopy.Pc = NULL;
        }
    }

    fflush( stdout );
    //Cblacs_gridexit( ictxt );
    //MPI_Finalize();
    return 0;
}

//#define MIN(a,b) ((a>b)?b:a)
//#define MAX(a,b) ((a>b)?a:b)

/*
 * produce distributed matrix,  
 */
void create_matrix( int ctxt, int seed, double **A, int *descA, int M, int N, int nb, int *np_A, int *nq_A)
{
    int info;
    int nprow, npcol, myrow, mycol;
    Cblacs_gridinfo( ctxt, &nprow, &npcol, &myrow, &mycol );

    // allocate the generator matrix and check matrix
    int np_iA = numroc_( &M, &nb, &myrow, &i0, &nprow );
    int nq_iA = numroc_( &N, &nb, &mycol, &i0, &npcol );

    if (np_iA*nq_iA!=0)
    {
        *A = malloc(np_iA*nq_iA*sizeof(**A)) ;
        if (*A == NULL) Cblacs_abort( ctxt, 10 );
        memset (*A, 0, np_iA*nq_iA*sizeof(**A));
    }
    else *A = NULL;

    if (descA != NULL)
    {
        int itemp = MAX( 1, np_iA );
        descinit_( descA, &M, &N, &nb, &nb, &i0, &i0, &ctxt, &itemp, &info );
        if (info != 0) Cblacs_abort( ctxt, 12 );
    }

    if (seed)
    {
        // fill in random numbers
        pdmatgen_ (&ctxt, "N", "N", &M, &N, &nb, &nb, *A,
                descA+8, descA+6, descA+7, 
                &seed, &i0, &np_iA, &i0, &nq_iA, 
                &myrow, &mycol, &nprow, &npcol);
    }

    /* set np and nq */
    if (np_A != NULL)
        *np_A = np_iA;
    if (nq_A != NULL)
        *nq_A = nq_iA;
}

