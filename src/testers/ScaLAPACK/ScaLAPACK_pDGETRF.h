/*
 * ScaLAPACK_pDGETRF.h
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


void ScaLAPACK_pDGETRF_calc(int n, double* A_global, int nb, int mpi_rank, int cprocs)
{
	/*
	 * n = system rank (A_global n x n)
	 */

	// general
	int i;				//iterators
	int zero = 0, one = 1;	//numbers
	//int nprocs = cprocs + sprocs;
	// MPI
	int ndims = 2, dims[2] = {0,0};
	// BLACS/SCALAPACK
	int nprow, npcol, info, ic = -1, context, context_global, myrow, mycol;
	int descA_global[9], descA[9];
	char order = 'R';
	// MATRIX
	int nr, nc, lld, lld_global;
	double *A;
	int *ipiv;

	// Initialize a default BLACS context and the processes grid
	MPI_Dims_create(cprocs, ndims, dims);
	nprow = dims[0];
	npcol = dims[1];
	Cblacs_get( ic, zero, &context );
	Cblacs_gridinit( &context, &order, nprow, npcol );
	Cblacs_get( ic, zero, &context_global );
	Cblacs_gridinit( &context_global, &order, one, one );
	//Cblacs_get( ic, zero, &context_all );
	//Cblacs_gridinit( &context_all, &order, one, nprocs );
	Cblacs_gridinfo( context, &nprow, &npcol, &myrow, &mycol );

	if (mpi_rank < cprocs)
	{
		// Computation of local matrix size
		//nb = SCALAPACKNB;
		nr = numroc_( &n, &nb, &myrow, &zero, &nprow );
		nc = numroc_( &n, &nb, &mycol, &zero, &npcol );
		lld = MAX( 1 , nr );
		A = malloc(nr*nc*sizeof(double));
		ipiv = malloc((lld+nb)*sizeof(int));

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

		// LU factorization
		pdgetrf_(      &n, &n, A, &one, &one, descA, ipiv, &info );
		pdgemr2d_ (&n, &n, A, &one, &one, descA, A_global, &one, &one, descA_global, &context);

		// cleanup
		free(A);
		free(ipiv);
	}

		//Close BLACS environment
		//Cblacs_barrier( context, "All");	// not working! why?
		//MPI_Barrier(MPI_COMM_WORLD);		// working but not needed
		//Cblacs_gridexit( context );		// not needed if calling blacs_exit
		//Cblacs_gridexit( context_global );// not needed if calling blacs_exit
		//Cblacs_exit( one );				// argument not 0: it is assumed the user will continue using the machine after the BLACS are done
											// error, if main function called more tha once, why?
	MPI_Barrier(MPI_COMM_WORLD);
}
