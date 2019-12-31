/*
 * FTLA_pDGETRF.h
 *
 *  Created on: Dec 27, 2019
 *      Author: marcello
 */

// code modified from ftdtr_main.c (create by Daniela in FTLA)

#include <mpi.h>
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
#include "create_matrix.h"
#include "errors.h"

extern void create_matrix (int ctxt, int seed, double **A, int *descA, int M, int N, int nb, int *np_A, int *nq_A);

extern int *errors;

int FTLA_ftdtr_calc(int rows, double* A_global, int rank, int cprocs, int sprocs)
{
	int i0=0, i1=1;
	int i;
    int NB=SCALAPACKNB;
    int M, N, Nc, Nr, Ne, S=1;
    ftla_work_t ftwork;

	// MPI
	int ndims = 2, dims[2] = {0,0};
	int P, Q;
	MPI_Dims_create(cprocs, ndims, dims);
	P = dims[0];
	Q = dims[1];
    // BLACS
    int ictxt, ictxt_global, info;
    int myrow, mycol, lld;
    // faults
    int Fstrat='e', F; // Fmin=0, Fmax=0, Finc=1;
    int Fmin, Fmax;
    int Finc = 1;
    Fmin= Fmax = sprocs;
    int err = 0;
    // matrices
    double* A=NULL;
    int descA[9], descA_global[9];

    {/* init BLACS */
        Cblacs_get( -1, 0, &ictxt );
        Cblacs_gridinit( &ictxt, "Row", P, Q );
        Cblacs_gridinfo( ictxt, &P, &Q, &myrow, &mycol );
        Cblacs_get( -1, 0, &ictxt_global );
        Cblacs_gridinit( &ictxt_global, "Row", i1, i1 );
    }

	{/* allocate matrices */
		/* determine checksum size, generate A matrix */
		N = M = rows;
		Nc = numroc_( &N, &NB, &mycol, &i0, &Q ); //LOCc(N_A)
		Nr = numroc_( &N, &NB, &myrow, &i0, &P ); //LOCr(N_A)
		lld = MAX( 1 , Nr );
		MPI_Allreduce( MPI_IN_PLACE, &Nc, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD );

#ifndef NO_EXTRAFLOPS
		Ne = N + Nc*2;// + ((Nc/NB)%Q==0)*NB;
#else
		Ne = N;
#endif

		// Descriptors (local)
		descinit_( descA, &N, &Ne, &NB, &NB, &i0, &i0, &ictxt, &lld, &info );

		if (rank==0)
		{
			// Descriptors (global)
			descinit_( descA_global, &N, &N, &i1, &i1, &i0, &i0, &ictxt_global, &N, &info );
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
		create_matrix( ictxt, 0,   &A,  descA, M, Ne, NB, NULL, NULL );

		/* allocate local buffer for the Q-wide local panel copy */
		create_matrix( ictxt, 0, (typeof(&A))&(ftwork.pcopy.Pc), ftwork.pcopy.descPc, M, (Q+2)*NB, NB, &(ftwork.pcopy.nrPc), &(ftwork.pcopy.ncPc) );

		MPI_Barrier(MPI_COMM_WORLD);

		// spread matrices
		pdgemr2d_(&N, &N, A_global, &i1, &i1, descA_global, A, &i1, &i1, descA, &ictxt);
	}

	{/* call resilient QR */
		int err=0;
		int* ipiv = (int*)malloc(Ne*sizeof(int) );

#ifdef INJECT
		for( F = Fmin; F<=Fmax; F+=Finc )
		{
			errors = create_error_list( M, NB, F, Fstrat );
#endif

			Cftla_work_construct( 0, descA, 0, Ne-N, &ftwork );

			do { // call ftpdgetrf until we complete w/o a failure

#ifdef USE_CoF
				if(err) Cftla_cof_dgeqrr( A, descA, tau, work, lwork, &ftwork );
#endif

				err = ftla_pdgetrf( &M, &Ne, A, &i1, &i1, descA, ipiv, &info, (int*)&ftwork );

#ifdef USE_CoF
				if(err) Cftla_cof_dgeqrr( A, descA, tau, work, lwork, &ftwork );
#endif

			} while(err);

			// collect matrices
			pdgemr2d_ (&N, &N, A, &i1, &i1, descA, A_global, &i1, &i1, descA_global, &ictxt);

			// cleanup
			Cftla_cof_cleanup( &ftwork );
			Cftla_work_destruct( &ftwork );
			free( errors );

#ifdef INJECT
		}
#endif

		free(ipiv);
	}

	{/* Cleanup */
		if( NULL != A  ) free( A );
		A = NULL;
		if( NULL != ftwork.pcopy.Pc ) free( ftwork.pcopy.Pc);
		ftwork.pcopy.Pc = NULL;
	}

    fflush( stdout );

	MPI_Barrier(MPI_COMM_WORLD);
    return 0;
}

