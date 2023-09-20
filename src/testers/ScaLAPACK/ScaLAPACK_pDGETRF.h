/*
 * ScaLAPACK_pDGETRF.h
 *
 *  Created on: May 15, 2020
 *      Author: marcello
 */

#include <mpi.h>
#include <time.h>
#include "../../helpers/macros.h"
#include "../../helpers/vector.h"
#include "../../helpers/Cblacs.h"
#include "../../helpers/matrix_basic.h"
#include "../../helpers/scalapack.h"
#include "../tester_structures.h"


test_output ScaLAPACK_pDGETRF(int n, double* A_global, double* B_global, int nb,
								int mpi_rank, int cprocs,
								int nprow, int npcol, int myrow, int mycol,
								int context, int context_global)
{
	test_output result = EMPTY_OUTPUT;

	result.total_start_time = time(NULL);

	// general
	int i;
	int i0 = 0;
	int i1 = 1;
	int m = 1;
	double d0 = 0.0;
	double d1 = 1.0;
	int info;
	int *ipiv;

	// matrix
	int nr, nc;
	double *A;
	double *At;
	int ncrhs, nrrhs;
	double *B;
	int ncrhst, nrrhst;
	double *Bt;

	int descA_global[9];
	int descB_global[9];
	int descA[9];
	int descAt[9];
	int descB[9];
	int descBt[9];

	int lld, lldt;


	if (mpi_rank < cprocs)
	{
		// Computation of local matrix size
		nc = numroc_( &n, &nb, &mycol, &i0, &npcol );
		nr = numroc_( &n, &nb, &myrow, &i0, &nprow );
		lld = MAX( 1 , nr );
		A = malloc(nr*nc*sizeof(double));
		At = malloc(nr*nc*sizeof(double));

		ncrhs = numroc_( &m, &nb, &mycol, &i0, &npcol );
		nrrhs = numroc_( &n, &nb, &myrow, &i0, &nprow );
		B = malloc(nrrhs*ncrhs*sizeof(double));

		ncrhst = numroc_( &n, &nb, &mycol, &i0, &npcol );
		nrrhst = numroc_( &m, &nb, &myrow, &i0, &nprow );
		lldt = MAX( 1 , nrrhst );
		Bt = malloc(nrrhst*ncrhst*sizeof(double));

		ipiv = malloc((lld+nb)*sizeof(int));

		// Descriptors (local)
		descinit_( descA, &n, &n, &nb, &nb, &i0, &i0, &context, &lld, &info );
		descinit_( descAt, &n, &n, &nb, &nb, &i0, &i0, &context, &lld, &info );
		descinit_( descB, &n, &m, &nb, &nb, &i0, &i0, &context, &lld, &info );
		descinit_( descBt, &m, &n, &nb, &nb, &i0, &i0, &context, &lldt, &info );

		if (mpi_rank==0)
		{
			// Descriptors (global)
			descinit_( descA_global, &n, &n, &i1, &i1, &i0, &i0, &context_global, &n, &info );
			descinit_( descB_global, &n, &m, &i1, &i1, &i0, &i0, &context_global, &n, &info );
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

		// spread A
		pdgemr2d_(&n, &n, A_global, &i1, &i1, descA_global, A, &i1, &i1, descA, &context);

		// transpose A
		pdtran_(&n, &n, &d1, A, &i1, &i1, descA, &d0, At, &i1, &i1, descAt);

		// LU factorization
		result.core_start_time = time(NULL);
		pdgetrf_( &n, &n, At, &i1, &i1, descAt, ipiv, &info );
		result.core_end_time = time(NULL);
		result.exit_code = info;

		// transpose back A
		pdtran_(&n, &n, &d1, At, &i1, &i1, descAt, &d0, A, &i1, &i1, descA);

		// get back A
		pdgemr2d_ (&n, &n, A, &i1, &i1, descA, A_global, &i1, &i1, descA_global, &context);

		result.total_end_time = time(NULL);

		/*
		 * evaluate n.r.e.
		 */
		// spread B
		pdgemr2d_(&n, &m, B_global, &i1, &i1, descB_global, B, &i1, &i1, descB, &context);

		pdgetrs_("N",  &n, &m, At, &i1, &i1, descAt, ipiv, B, &i1, &i1, descB, &info  );
			/*
			 * no need of transposing back B, being an 1D vector
			 */
			/*
			// re-transpose result
			pdtran_(&m, &n, &d1, B, &i1, &i1, descB, &d0, Bt, &i1, &i1, descBt);

			// collect result
			if (mpi_rank==0)
			{
				// Adapt descriptor for B to accept transposed matrix
				descinit_( descB_global, &m, &n, &i1, &i1, &i0, &i0, &context_global, &m, &info );
			}

			// get back result
			pdgemr2d_(&m, &n, Bt, &i1, &i1, descBt, B_global, &i1, &i1, descB_global, &context);
			 */
		pdgemr2d_(&n, &m, B, &i1, &i1, descB, B_global, &i1, &i1, descB_global, &context);
	}
	else
	{
		A  = NULL;
		At = NULL;
		B  = NULL;
		Bt = NULL;
		ipiv = NULL;

		result.core_start_time = time(NULL);
		result.core_end_time = time(NULL);
		result.exit_code = 0;
	}

	// cleanup
	NULLFREE(A);
	NULLFREE(At);
	NULLFREE(B);
	NULLFREE(Bt);
	NULLFREE(ipiv);

	MPI_Barrier(MPI_COMM_WORLD);
	return result;
}
