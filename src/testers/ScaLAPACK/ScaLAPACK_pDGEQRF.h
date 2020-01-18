/*
 * ScaLAPACK_pDGEQRF.h
 *
 *  Created on: Dec 27, 2019
 *      Author: marcello
 */

#include <mpi.h>
#include "../../helpers/macros.h"
#include "../../helpers/matrix.h"
#include "../../helpers/vector.h"
#include "../../helpers/Cblacs.h"
#include "../../helpers/scalapack.h"


void ScaLAPACK_pDGEQRF_calc(int n, double* A_global, int nb, int mpi_rank, int cprocs)
{
	/*
	 * n = system rank (A_global n x n)
	 */

	// general
	int i;						//iterators
	int zero = 0, one = 1;			//numbers
	// MPI
	int ndims = 2, dims[2] = {0,0};
	// BLACS/SCALAPACK
	int nprow, npcol, info, ic = -1, context, context_global, myrow, mycol;
	int descA_global[9], descA[9];
	char order = 'R';
	// MATRIX
	int nr, nc, lld, lld_global;
	double* A;
	double* work;
	double* tau;

	// Initialize a default BLACS context and the processes grid
	MPI_Dims_create(cprocs, ndims, dims);
	nprow = dims[0];
	npcol = dims[1];
	Cblacs_get( ic, zero, &context );
	Cblacs_gridinit( &context, &order, nprow, npcol );
	Cblacs_get( ic, zero, &context_global );
	Cblacs_gridinit( &context_global, &order, one, one );
	Cblacs_gridinfo( context, &nprow, &npcol, &myrow, &mycol );

	if (mpi_rank < cprocs)
	{
		// Computation of local matrix size
		//nb = SCALAPACKNB;
		nr = numroc_( &n, &nb, &myrow, &zero, &nprow );
		nc = numroc_( &n, &nb, &mycol, &zero, &npcol );
		lld = MAX( 1 , nr );
		A = malloc(nr*nc*sizeof(double));
		tau = malloc( nc*sizeof(double) );

		// Descriptors (local)
		descinit_( descA, &n, &n, &nb, &nb, &zero, &zero, &context, &lld, &info );

		if (mpi_rank==0)
		{
			// Descriptors (global)
			lld_global = n;
			descinit_( descA_global, &n, &n, &one, &one, &zero, &zero, &context_global, &lld_global, &info );
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
		pdgemr2d_(&n, &n, A_global, &one, &one, descA_global, A, &one, &one, descA, &context);

		// init work space
		int lwork=-1;
		double lazywork;
		pdgeqrf_(  &n, &n, A, &one, &one, descA, NULL, &lazywork, &lwork, &info );
		lwork = (int)lazywork;
		work = malloc( lwork*sizeof(double) );

		// QR factorization
		pdgeqrf_(  &n, &n, A, &one, &one, descA, tau, work, &lwork, &info );

		pdgemr2d_ (&n, &n, A, &one, &one, descA, A_global, &one, &one, descA_global, &context);

		// cleanup
		free(A);
		free(work);
		free(tau);
	}

	MPI_Barrier(MPI_COMM_WORLD);
}
