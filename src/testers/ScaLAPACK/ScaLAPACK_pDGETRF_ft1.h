/*
 * ScaLAPACK_pDGETRF_ft1.h
 *
 *  Created on: Dec 28, 2019
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
 * LU decomposition by ScaLAPACK
 *
 * WARNING: WITHOUT transposition of the input matrix A
 *
 */

test_output ScaLAPACK_pDGETRF_ft1(int n, double* A_source, int nb, \
									int mpi_rank, int cprocs, int sprocs, int failing_level, int checkpoint_freq, \
									int nprow, int npcol, int myrow, int mycol, \
									int context_distributed, int context_source, int context_all, int context_cp)
{
	test_output result = EMPTY_OUTPUT;

	result.total_start_time = time(NULL);

	/*
	 * n = system rank (A_global n x n)
	 */

	// general
	int i;					//iterators
	int zero = 0, one = 1;	//numbers
	int nprocs = cprocs + sprocs;

	// MPI
	//int ndims = 2, dims[2] = {0,0};
	// BLACS/SCALAPACK
	int info;
	//int ic = -1, context_distributed, context_source, context_cp, context_all;
	//int nprow, npcol, myrow, mycol;
	//int tmprow,tmpcol,tmpmyrow, tmpmycol;
	int descA_source[9], descA[9], descA_cp[9];
	//int descB_global[9], descB[9];
	//char order = 'R';
	// MATRIX
	int nr, nc, lldA;
	int nIPIV;
	double *A;
	double *A_cp;
	int *IPIV;
	int *IPIV_cp;

	/*
	// Initialize a default BLACS context and the processes grid
	MPI_Dims_create(cprocs, ndims, dims);
	nprow = dims[0];
	npcol = dims[1];

	// context for all processes
	Cblacs_get( ic, zero, &context_all );
	Cblacs_gridinit( &context_all, &order, one, nprocs );

	// context for distributed matrix A
	Cblacs_get( ic, zero, &context_distributed );
	Cblacs_gridinit( &context_distributed, &order, nprow, npcol );
	Cblacs_gridinfo( context_distributed, &nprow, &npcol, &myrow, &mycol );

	// context for source global matrix A (A_source)
	Cblacs_get( ic, zero, &context_source );
	Cblacs_gridinit( &context_source, &order, one, one );
	Cblacs_gridinfo( context_source, &tmprow, &tmpcol, &tmpmyrow, &tmpmycol );

	// context for checkpointing global matrix (A_cp)
	Cblacs_get( ic, zero, &context_cp );
	int map_cp[1];
	map_cp[0]=cprocs;
	Cblacs_gridmap( &context_cp, map_cp, one, one, one);
	Cblacs_gridinfo( context_cp, &tmprow, &tmpcol, &tmpmyrow, &tmpmycol );
	*/

	// Computation of local matrix size
	nr = numroc_( &n, &nb, &myrow, &zero, &nprow );
	nc = numroc_( &n, &nb, &mycol, &zero, &npcol );
	lldA = MAX( 1 , nr );
	MPI_Allreduce( MPI_IN_PLACE, &lldA, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD );

	nIPIV = (lldA+nb);

	// global matrices
	if (mpi_rank==0)
	{
		// Descriptors (global, for root node)
		descinit_( descA_source, &n, &n, &one, &one, &zero, &zero, &context_source, &n, &info );
	}
	else
	{
		// Allocation not needed
		A_source=NULL;

		// Descriptors (global, for non-root nodes)
		for (i=0; i<9; i++)
		{
			descA_source[i]=0;
		}
		descA_source[1]=-1;
	}

	if (mpi_rank==cprocs) // spare (checkpointing) node
	{
		// Descriptors
		A_cp = malloc(n*n*sizeof(double));
		descinit_( descA_cp, &n, &n, &one, &one, &zero, &zero, &context_cp, &n, &info );

		//IPIV_cp=malloc(nIPIV*cprocs*sizeof(int)); // with cprocs is not good because MPI_GATHER wants a buffer for everyone
		IPIV_cp=malloc(nIPIV*nprocs*sizeof(int));
	}
	else
	{
		// Allocation not needed
		A_cp=NULL;
		IPIV_cp=NULL;

		// Descriptors (global, for non-spare nodes)
		for (i=0; i<9; i++)
		{
			descA_cp[i]=0;
		}
		descA_cp[1]=-1;
		descA_cp[4]=nb; // *** important!
		descA_cp[5]=nb; // ***
	}

	if (mpi_rank < cprocs) // non-spare nodes
	{
		// Allocation
		A = malloc(nr*nc*sizeof(double));
		IPIV = malloc(nIPIV*sizeof(int));

		// Descriptors (local)
		// all processes have to know descA (not only non-spare nodes)
		descinit_( descA, &n, &n, &nb, &nb, &zero, &zero, &context_distributed, &lldA, &info );
	}
	else
	{
		// Allocation
		A=NULL;
		IPIV = malloc(nIPIV*sizeof(int)); // also allocated on the spare proc because MPI_GATHER wants a buffer for everyone

		// Descriptors
		for (i=0; i<9; i++)
		{
			descA[i]=0;
		}
		// all processes have to know something about descA (not only non-spare nodes)
		// can't use descinint, due to illegal values of spare process not belonging to the right context
		descA[1]=-1;
		descA[4]=nb;
		descA[5]=nb;
	}

	// spread matrices
	pdgemr2d_(&n, &n, A_source, &one, &one, descA_source, A, &one, &one, descA, &context_all);

	if (failing_level>=0)
	{
		failing_level=n-failing_level;
	}

	// LU factorization
	result.core_start_time = time(NULL);
	pdgetrf_cp_  (&n, &n, A, &one, &one, descA, A_cp, &one, &one, descA_cp, IPIV, IPIV_cp, &nIPIV, &checkpoint_freq, &failing_level, &context_all, &info );
	result.core_end_time = time(NULL);
	result.exit_code = info;

	pdgemr2d_ (&n, &n, A, &one, &one, descA, A_source, &one, &one, descA_source, &context_all);

	if (mpi_rank < cprocs)
	{
		// cleanup
		free(A);
		free(IPIV);
	}
	else
	{
		free(A_cp);
		free(IPIV_cp);
		free(IPIV);
	}

	MPI_Barrier(MPI_COMM_WORLD);

	result.total_end_time = time(NULL);

	return result;
}
