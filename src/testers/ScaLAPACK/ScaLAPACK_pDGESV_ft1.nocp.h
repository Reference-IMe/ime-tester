/*
 * ScaLAPACK_pDGESV_ft1.nocp.h
 *
 *  Created on: Dec 28, 2019
 *      Author: marcello
 */

#include <mpi.h>
#include <time.h>
#include "../../helpers/macros.h"
#include "../IMe/lib/src/helpers/vector_basic.h"
#include "../../helpers/Cblacs.h"
#include "../IMe/lib/src/helpers/matrix_basic.h"
#include "../../helpers/scalapack.h"
#include "../tester_structures.h"

// solve with separate calls to factorization and solution
// to prepare for a checkpointing version
// checkpointing NOT implemented => old, unmaintained version

exit_status ScaLAPACK_pDGESV_ft1(int n, double* A_global, int m, double* B_global, int nb, int mpi_rank, int cprocs, int sprocs, int failing_level, int checkpoint_freq)
{
	exit_status result = EMPTY_OUTPUT;

	result.total_start_time = time(NULL);

	/*
	 * n = system rank (A_global n x n)
	 * m = num. of r.h.s (B_global n x m)
	 */

	// to avoid warnings, until checkpointing is implemented
	failing_level = failing_level+0;

	// general
	int i;				//iterators
	int zero = 0, one = 1;	//numbers
	double dzero = 0.0, done = 1.0;
	int nprocs = cprocs + sprocs;
	// MPI
	int ndims = 2, dims[2] = {0,0};
	// BLACS/SCALAPACK
	int nprow, npcol, info, ic = -1, context, context_global, context_all, myrow, mycol;
	int descA_global[9], descB_global[9], descA[9], descAt[9], descB[9], descBt[9];
	char order = 'R';
	// MATRIX
	int nr, nc, ncrhs, nrrhs, lld, lld_global, lldt;
	int ncrhst, nrrhst;
	double *A, *At, *B, *Bt;
	int *ipiv;
	
	// Initialize a default BLACS context and the processes grid
	MPI_Dims_create(cprocs, ndims, dims);
	nprow = dims[0];
	npcol = dims[1];
	Cblacs_get( ic, zero, &context );
	Cblacs_gridinit( &context, &order, nprow, npcol );
	Cblacs_get( ic, zero, &context_global );
	Cblacs_gridinit( &context_global, &order, one, one );
	Cblacs_get( ic, zero, &context_all );
	Cblacs_gridinit( &context_all, &order, one, nprocs );
	Cblacs_gridinfo( context, &nprow, &npcol, &myrow, &mycol );

	if (mpi_rank < cprocs)
	{
		//printf("\nI %d have to",myrank);

		// Computation of local matrix size
		//nb = SCALAPACKNB;
		nc = numroc_( &n, &nb, &mycol, &zero, &npcol );
		nr = numroc_( &n, &nb, &myrow, &zero, &nprow );
		lld = MAX( 1 , nr );
		A = malloc(nr*nc*sizeof(double));
		At = malloc(nr*nc*sizeof(double));

		ncrhs = numroc_( &m, &nb, &mycol, &zero, &npcol );
		nrrhs = numroc_( &n, &nb, &myrow, &zero, &nprow );
		B = malloc(nrrhs*ncrhs*sizeof(double));

		ncrhst = numroc_( &n, &nb, &mycol, &zero, &npcol );
		nrrhst = numroc_( &m, &nb, &myrow, &zero, &nprow );
		lldt = MAX( 1 , nrrhst );
		Bt = malloc(nrrhst*ncrhst*sizeof(double));

		ipiv = malloc((lld+nb)*sizeof(int));

		// Descriptors (local)
		descinit_( descA, &n, &n, &nb, &nb, &zero, &zero, &context, &lld, &info );
		descinit_( descAt, &n, &n, &nb, &nb, &zero, &zero, &context, &lld, &info );
		descinit_( descB, &n, &m, &nb, &nb, &zero, &zero, &context, &lld, &info );
		descinit_( descBt, &m, &n, &nb, &nb, &zero, &zero, &context, &lldt, &info );

		if (mpi_rank==0)
		{
			// Descriptors (global)
			lld_global = n;
			descinit_( descA_global, &n, &n, &one, &one, &zero, &zero, &context_global, &lld_global, &info );
			descinit_( descB_global, &n, &m, &one, &one, &zero, &zero, &context_global, &lld_global, &info );
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

		// spread matrices
		pdgemr2d_(&n, &n, A_global, &one, &one, descA_global, A, &one, &one, descA, &context);
		pdgemr2d_(&n, &m, B_global, &one, &one, descB_global, B, &one, &one, descB, &context);

		// transpose system matrix
		pdtran_(&n, &n, &done, A, &one, &one, descA, &dzero, At, &one, &one, descAt);

		// Linear system equations solver
		// pdgesv_(  &n, &m, A, &one, &one, descA, ipiv, B, &one, &one, descB, &info );
		// split in LU factorization + solve (pdgetrf + pdgetrs) to introduce checkpointing
		result.core_start_time = time(NULL);
		pdgetrf_(      &n, &n, At, &one, &one, descAt, ipiv, &info );
		pdgetrs_("N",  &n, &m, At, &one, &one, descAt, ipiv, B, &one, &one, descB, &info  );
	    result.core_end_time = time(NULL);
		result.exit_code = info;

		// re-transpose result
		pdtran_(&m, &n, &done, B, &one, &one, descB, &dzero, Bt, &one, &one, descBt);

		// collect result
		if (mpi_rank==0)
		{
			// Descriptors (global)
			descinit_( descB_global, &m, &n, &one, &one, &zero, &zero, &context_global, &m, &info );
		}
		pdgemr2d_(&m, &n, Bt, &one, &one, descBt, B_global, &one, &one, descB_global, &context);

		free(A);
		free(B);
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

	result.total_end_time = time(NULL);

	return result;
}
