/*
 * FTLA_pDGEQRF.h
 *
 *  Created on: Dec 27, 2019
 *      Author: marcello
 */

// code modified from ftdqr_main.c of FTLA

#include <mpi.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include "../../helpers/Cblacs.h"
#include "../../helpers/scalapack.h"
#include "../tester_structures.h"
#include "commons.h"
#include "create_matrix.h"
#include "ftla-rSC13.mod/ftla_cof.h"
#include "ftla-rSC13.mod/ftla_driver.h"
#include "ftla-rSC13.mod/ftla_ftwork.h"
#include "ftla-rSC13.mod/util_inject.h"
#include "ftla-rSC13.mod/util_matrix.h"

extern int *errors;

extern MPI_Comm ftla_current_comm;

test_output FTLA_ftdqr(	int rows, double* A_global, int NB,			\
						int mpi_rank, int cprocs, int sprocs,		\
						int nprow, int npcol, int myrow, int mycol,	\
						int ctxt, int ctxt_root)
{
	test_output result = EMPTY_OUTPUT;

	result.total_start_time = time(NULL);

	int i;
	int i0=0;
	int i1=1;
    int M, N, Nc, Ne;
	int info;
    ftla_work_t ftwork;

	if (mpi_rank>=cprocs)
	{
		MPI_Comm_split(MPI_COMM_WORLD, MPI_UNDEFINED, MPI_UNDEFINED, &ftla_current_comm);
	}
	else
	{
		MPI_Comm_split(MPI_COMM_WORLD, 1, mpi_rank, &ftla_current_comm);
	}

    // faults
    int Fstrat='e', F; // Fmin=0, Fmax=0, Finc=1;
    int Fmin, Fmax;
    int Finc = 1;
    Fmin= Fmax = sprocs;

    // matrices
    double* A=NULL;
    int descA[9];
    int descA_global[9];

	{/* allocate matrices */
		/* determine checksum size, generate A matrix */
		N = M = rows;
		Nc = numroc_( &N, &NB, &mycol, &i0, &npcol ); //LOCc(N_A)
		//Nr = numroc_( &N, &NB, &myrow, &i0, &P ); //LOCr(N_A)
		//lld = MAX( 1 , Nr );
		MPI_Allreduce( MPI_IN_PLACE, &Nc, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD );

#ifndef NO_EXTRAFLOPS
		Ne = N + Nc*2;// + ((Nc/NB)%npcol==0)*NB;
#else
		Ne = N;
#endif

		if (mpi_rank < cprocs)	// only calc nodes have a local copy of submatrix A
		{
			// Descriptors (local)
			descinit_( descA, &N, &Ne, &NB, &NB, &i0, &i0, &ctxt, &Nc, &info );
		}
		else
		{
			for (i=0; i<9; i++)
			{
				descA[i]=0;
			}
			descA[1]=-1;
		}

		if (mpi_rank==0)
		{
			// Descriptors (global)
			descinit_( descA_global, &N, &N, &i1, &i1, &i0, &i0, &ctxt_root, &N, &info );
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

		if (mpi_rank < cprocs)	// only calc nodes initialize matrices
		{
			create_matrix( ctxt, 0,   &A,  descA, M, Ne, NB, NULL, NULL );

			/* allocate local buffer for the npcol-wide local panel copy */
			create_matrix( ctxt, 0, (typeof(&A))&(ftwork.pcopy.Pc), ftwork.pcopy.descPc, M, (npcol+2)*NB, NB, &(ftwork.pcopy.nrPc), &(ftwork.pcopy.ncPc) );

			// spread matrices
			pdgemr2d_(&N, &N, A_global, &i1, &i1, descA_global, A, &i1, &i1, descA, &ctxt);
		}
	}

	if (mpi_rank < cprocs)
	{/* call resilient QR */
		int err=0;
		int lwork=-1;
		double lazywork;
		pdgeqrf_( &M, &Ne, NULL, &i1, &i1, descA, NULL, &lazywork, &lwork, &info );
		lwork = (int)lazywork;
		double *work = (double*)malloc( lwork*sizeof(double) );
		double *tau  = (double*)malloc( Nc*sizeof(double) );

		result.core_start_time = time(NULL);

#ifdef INJECT        
		for( F = Fmin; F<=Fmax; F+=Finc )
		{
			errors = create_error_list( M, NB, F, Fstrat );
#endif

			Cftla_work_construct( 0, descA, 0, Ne-N, &ftwork );

			do { // call ftpdgeqrf until we complete w/o a failure

#ifdef USE_CoF
				if(err) Cftla_cof_dgeqrr( A, descA, tau, work, lwork, &ftwork );
#endif

				err = ftla_pdgeqrf( &M, &Ne, A, &i1, &i1, descA, tau, work, &lwork, &info, (int*)&ftwork );

#ifdef USE_CoF
				if(err) Cftla_cof_dgeqrr( A, descA, tau, work, lwork, &ftwork );
#endif

			} while(err);

			result.core_end_time = time(NULL);
			result.exit_code = info;

			// collect matrices
			pdgemr2d_ (&N, &N, A, &i1, &i1, descA, A_global, &i1, &i1, descA_global, &ctxt);

			// cleanup
			Cftla_cof_cleanup( &ftwork );
			Cftla_work_destruct( &ftwork );
			free( errors );

#ifdef INJECT
		}
#endif

		free( work );
		free( tau );
	}
	else
	{
		result.core_start_time = time(NULL);
		Cftla_work_construct( 0, descA, 0, Ne-N, &ftwork );
		result.core_end_time = time(NULL);
	}

	if (mpi_rank < cprocs)
	{/* Cleanup */
		if( NULL != A  ) free( A );
		A = NULL;
		if( NULL != ftwork.pcopy.Pc ) free( ftwork.pcopy.Pc);
		ftwork.pcopy.Pc = NULL;
	}

    fflush( stdout );

	MPI_Barrier(MPI_COMM_WORLD);

	result.total_end_time = time(NULL);

	return result;
}

