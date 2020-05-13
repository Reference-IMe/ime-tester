/*
 * ScaLAPACK_pDGEQRF.h
 *
 *  Created on: Dec 27, 2019
 *      Author: marcello
 */

#include <mpi.h>
#include <time.h>
#include "../../helpers/macros.h"
#include "../../helpers/matrix.h"
#include "../../helpers/vector.h"
#include "../../helpers/Cblacs.h"
#include "../../helpers/scalapack.h"
#include "../tester_structures.h"

/*
 * QR decomposition by ScaLAPACK
 *
 * WARNING: WITHOUT transposition of the input matrix A
 *
 */

test_output ScaLAPACK_pDGEQRF(	int n, double* A_global, int nb,			\
								int mpi_rank, int cprocs,					\
								int nprow, int npcol, int myrow, int mycol,	\
								int context, int context_global)
{
	test_output result = EMPTY_OUTPUT;

	result.total_start_time = time(NULL);

	/*
	 * n = system rank (A_global n x n)
	 */

	// general
	int i;
	int i0 = 0;
	int i1 = 1;
	int info;
	double* work;
	double* tau;

	// matrix
	int nr;
	int nc;
	int lld;
	int lld_global;
	double* A;
	int descA_global[9];
	int descA[9];


	if (mpi_rank < cprocs)
	{
		// Computation of local matrix size
		nr = numroc_( &n, &nb, &myrow, &i0, &nprow );
		nc = numroc_( &n, &nb, &mycol, &i0, &npcol );
		lld = MAX( 1 , nr );
		A = malloc(nr*nc*sizeof(double));
		tau = malloc( nc*sizeof(double) );

		// Descriptors (local)
		descinit_( descA, &n, &n, &nb, &nb, &i0, &i0, &context, &lld, &info );

		if (mpi_rank==0)
		{
			// Descriptors (global, for root node)
			lld_global = n;
			descinit_( descA_global, &n, &n, &i1, &i1, &i0, &i0, &context_global, &lld_global, &info );
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

		// spread matrices
		pdgemr2d_(&n, &n, A_global, &i1, &i1, descA_global, A, &i1, &i1, descA, &context);

		// init work space
		int lwork=-1;
		double lazywork;
		pdgeqrf_(  &n, &n, A, &i1, &i1, descA, NULL, &lazywork, &lwork, &info );
		lwork = (int)lazywork;
		work = malloc( lwork*sizeof(double) );

		// QR factorization
		result.core_start_time = time(NULL);
		pdgeqrf_(  &n, &n, A, &i1, &i1, descA, tau, work, &lwork, &info );
		result.core_end_time = time(NULL);
		result.exit_code = info;

		pdgemr2d_ (&n, &n, A, &i1, &i1, descA, A_global, &i1, &i1, descA_global, &context);

		// cleanup
		free(A);
		free(work);
		free(tau);
	}

	MPI_Barrier(MPI_COMM_WORLD);

	result.total_end_time = time(NULL);

	return result;
}
