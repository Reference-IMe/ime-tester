/*
 * FTLA_pDGEQRF.h
 *
 *  Created on: May 15, 2020
 *      Author: marcello
 */

// code modified from ftdqr_main.c of FTLA

#include <mpi.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include "../../helpers/macros.h"
#include "../../helpers/Cblacs.h"
#include "../../helpers/scalapack.h"
#include "../tester_structures.h"
#include "commons.h"
#include "create_matrix.h"
#include "ftla-rSC13/ftla_cof.h"
#include "ftla-rSC13/ftla_driver.h"
#include "ftla-rSC13/ftla_ftwork.h"
#include "ftla-rSC13/util_inject.h"
#include "ftla-rSC13/util_matrix.h"

extern int *errors;

extern MPI_Comm ftla_current_comm;

test_output FTLA_ftdqr(	int n, double* A_global, double* B_global, int nb,
						int mpi_rank, int cprocs, int sprocs,
						int nprow, int npcol, int myrow, int mycol,
						int ctxt, int ctxt_root)
{
	test_output result = EMPTY_OUTPUT;

	result.total_start_time = time(NULL);

	int i;
	int i0=0;
	int i1=1;
	double d0 = 0.0;
	double d1 = 1.0;
	int info;
	ftla_work_t ftwork;
	int lwork;
	double lazywork;
	double* work;
	double* tau;

	if (mpi_rank>=cprocs)
	{
		MPI_Comm_split(MPI_COMM_WORLD, MPI_UNDEFINED, MPI_UNDEFINED, &ftla_current_comm);
	}
	else
	{
		MPI_Comm_split(MPI_COMM_WORLD, 1, mpi_rank, &ftla_current_comm);
	}

	// faults
	int Fstrat = 'e', F;
	/*
	int Fmin, Fmax, Finc;

	Fmin = 0;
	Fmax = sprocs;
	Finc = 1;
	*/

	// matrices
	int nc, nr, ne;
	double* A;
	double* At;
	int m = 1; // B is a vector that will hold the solution vector for checking purposes
	int ncrhs, nrrhs;
	double *B;

	int descA_global[9];
	int descB_global[9];
	int descA[9];
	int descAt[9];
	int descB[9];

	/* allocate matrices */

		/* A */
		nc = numroc_( &n, &nb, &mycol, &i0, &npcol );
		nr = numroc_( &n, &nb, &myrow, &i0, &nprow );
		//MPI_Allreduce( MPI_IN_PLACE, &nc, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD );
		//MPI_Allreduce( MPI_IN_PLACE, &nr, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD );

		/* B */
		ncrhs = numroc_( &m, &nb, &mycol, &i0, &npcol ); // one column vector
		nrrhs = numroc_( &n, &nb, &myrow, &i0, &nprow );

#ifndef NO_EXTRAFLOPS
		ne = n + nc*2;// + ((nc/nb)%npcol==0)*nb;
#else
		ne = n;
#endif
		if (mpi_rank < cprocs)	// only calc nodes have a local copy of submatrices A and B
		{
			// Descriptors (local)
			descinit_( descB, &n, &m, &nb, &nb, &i0, &i0, &ctxt, &nr, &info );
			// descA inited below
		}
		else
		{
			for (i=0; i<9; i++)
			{
				descA[i]=0;
				descAt[i]=0;
			}
			descA[1]=-1;
			descAt[1]=-1;
		}

		if (mpi_rank==0)
		{
			// Descriptors (global)
			descinit_( descA_global, &n, &n, &i1, &i1, &i0, &i0, &ctxt_root, &n, &info );
			descinit_( descB_global, &n, &m, &i1, &i1, &i0, &i0, &ctxt_root, &n, &info );
		}
		else
		{
			// Descriptors (global, for non-root nodes)
			for (i=0; i<9; i++)
			{
				descA_global[i]=0;
				descB_global[i]=0;
			}
			descA_global[1]=-1;
			descB_global[1]=-1;
		}

		if (mpi_rank < cprocs)	// only calc nodes initialize matrices
		{
			B = malloc(nrrhs*ncrhs*sizeof(double));
			create_matrix( ctxt, 0, &A,  descA,  n, ne, nb, NULL, NULL );
			create_matrix( ctxt, 0, &At, descAt, n, ne, nb, NULL, NULL );

			/* allocate local buffer for the npcol-wide local panel copy */
			create_matrix( ctxt, 0, (typeof(&A))&(ftwork.pcopy.Pc), ftwork.pcopy.descPc, n, (npcol+2)*nb, nb, &(ftwork.pcopy.nrPc), &(ftwork.pcopy.ncPc) );

			// spread matrices
			pdgemr2d_(&n, &n, A_global, &i1, &i1, descA_global, A, &i1, &i1, descA, &ctxt);

			// transpose system matrix
			pdtran_(&n, &n, &d1, A, &i1, &i1, descA, &d0, At, &i1, &i1, descAt);
		}
		else
		{
			B  = NULL;
			A  = NULL;
			At = NULL;
			ftwork.pcopy.Pc = NULL;
		}

	/* call resilient QR */
		if (mpi_rank < cprocs)
		{
			int err=0;
			lwork=-1;
			pdgeqrf_( &n, &ne, NULL, &i1, &i1, descAt, NULL, &lazywork, &lwork, &info );
			lwork = (int)lazywork;

			work = (double*)malloc( lwork*sizeof(double) );
			tau  = (double*)malloc( nc*sizeof(double) );

			result.core_start_time = time(NULL);

#ifdef 		INJECT
			//for( F = Fmin; F<=Fmax; F+=Finc )
			{
				F=sprocs;
				errors = create_error_list( n, nb, F, Fstrat );
#endif

				Cftla_work_construct( 0, descAt, 0, ne-n, &ftwork );

				do { // call ftpdgeqrf until we complete w/o a failure

#ifdef USE_CoF
					if(err) Cftla_cof_dgeqrr( At, descAt, tau, work, lwork, &ftwork );
#endif

					err = ftla_pdgeqrf( &n, &ne, At, &i1, &i1, descAt, tau, work, &lwork, &info, (int*)&ftwork );

#ifdef USE_CoF
					if(err) Cftla_cof_dgeqrr( At, descAt, tau, work, lwork, &ftwork );
#endif

				} while(err);

				// cleanup
				Cftla_cof_cleanup( &ftwork );
				Cftla_work_destruct( &ftwork );
				free( errors );

#ifdef INJECT
			}
#endif

			result.core_end_time = time(NULL);
			result.exit_code = info;

			// transpose back
			pdtran_(&n, &n, &d1, At, &i1, &i1, descAt, &d0, A, &i1, &i1, descA);

			// collect matrices
			pdgemr2d_ (&n, &n, A, &i1, &i1, descA, A_global, &i1, &i1, descA_global, &ctxt);

			free( work );
			//free( tau ); // do NOT clear, as tau is used in n.r.e. estimation
		}
		else
		{
			work = NULL;
			tau  = NULL;

			result.core_start_time = time(NULL);
			//Cftla_work_construct( 0, descAt, 0, ne-n, &ftwork );
			result.core_end_time = time(NULL);
		}

		MPI_Barrier(MPI_COMM_WORLD);

		result.total_end_time = time(NULL);

	/*
	 * estimate n.r.e
	 */
		if (mpi_rank < cprocs)
		{
			/*
			 * calc solution vector for n.r.e
			 *
			 * https://en.wikipedia.org/wiki/QR_decomposition#Using_for_solution_to_linear_inverse_problems
			 *
			 * A=Q.R
			 * Q.R.x=B
			 * x=R^-1.(Q'.B)
			 *
			 * QR -> A		already calculated by pdgeqrf_
			 * Q'.B -> B	calc with pdormqr_
			 * R^-1.B) -> B	calc with pdtrsm_
			 *
			 */

			pdgemr2d_(&n, &m, B_global, &i1, &i1, descB_global, B, &i1, &i1, descB, &ctxt);

			lwork = -1;
			pdormqr_( "L", "T", &n, &m, &n, At, &i1, &i1, descAt, tau, B, &i1, &i1, descB, &lazywork, &lwork, &info );

			lwork = (int) lazywork;
			work = (double*) malloc( lwork * sizeof(double) );

			pdormqr_( "L", "T", &n, &m, &n, At, &i1, &i1, descAt, tau, B, &i1, &i1, descB, work, &lwork, &info);
			pdtrsm_("L", "U", "N", "N", &n, &m, &d1, At, &i1, &i1, descAt, B, &i1, &i1, descB);

			// collect result
			pdgemr2d_(&n, &m, B, &i1, &i1, descB, B_global, &i1, &i1, descB_global, &ctxt);
		}

	/* Cleanup */
		NULLFREE(A);
		NULLFREE(At);
		NULLFREE(B);
		NULLFREE(work);
		NULLFREE(tau);
		NULLFREE(ftwork.pcopy.Pc);

	fflush(stdout);
	MPI_Barrier(MPI_COMM_WORLD);
	return result;
}
