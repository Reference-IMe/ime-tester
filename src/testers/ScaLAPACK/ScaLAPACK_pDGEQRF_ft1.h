/*
 * ScaLAPACK_pDGEQRF_ft1.h
 *
 *  Created on: Jan 8, 2020
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

test_output ScaLAPACK_pDGEQRF_ft1(int n, double* A_source, int nb, \
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
	int i;				//iterators
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

	//char order = 'R';
	// MATRIX
	int nr, nc, lldA;
	double* A;
	double* A_cp;
	double* work;
	double* work_cp;
	double* tau;
	double* tau_cp;
	int lwork=-1, ltau;
	double lazywork;

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

	ltau = lldA;

	// global matrices
	if (mpi_rank==0)
	{
		// Descriptors (global, for root node)
		descinit_( descA_source, &n, &n, &one, &one, &zero, &zero, &context_source, &n, &info );
	}
	else
	{
		// Descriptors (global, for non-root nodes)
		A_source=NULL;
		for (i=0; i<9; i++)
		{
			descA_source[i]=0;
		}
		descA_source[1]=-1;
	}

	// locally distributed matrices
	if (mpi_rank < cprocs)				// non-spare nodes
	{
		// Allocation
		A = malloc(nr*nc*sizeof(double));

		// Descriptors (local)
		// all processes have to know descA (not only non-spare nodes)
		descinit_( descA, &n, &n, &nb, &nb, &zero, &zero, &context_distributed, &lldA, &info );

		// init work space
		pdgeqrf_(  &n, &n, A, &one, &one, descA, NULL, &lazywork, &lwork, &info );
		lwork = (int)lazywork;
	}
	else								// spare node
	{
		// Allocation not needed
		A=NULL;

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

	// workspace
	MPI_Allreduce( MPI_IN_PLACE, &lwork, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD );
	work = malloc( lwork*sizeof(double) ); // also allocated on the spare proc because MPI_GATHER wants a buffer for everyone
	tau = malloc( ltau*sizeof(double) );

	// checkpointing matrices
	if (mpi_rank==cprocs)	// spare (checkpointing) node
	{
		// Allocation
		A_cp = malloc(n*n*sizeof(double));
		work_cp=malloc(lwork*nprocs*sizeof(double)); // with cprocs instead of nprocs is not good because MPI_GATHER wants a buffer for everyone!
		tau_cp=malloc(ltau*nprocs*sizeof(double));

		// Descriptors
		descinit_( descA_cp, &n, &n, &one, &one, &zero, &zero, &context_cp, &n, &info );
	}
	else					// non-spare nodes
	{
		// Allocation not needed
		A_cp=NULL;
		work_cp=NULL;
		tau_cp=NULL;

		// Descriptors (global, for non-spare nodes)
		for (i=0; i<9; i++)
		{
			descA_cp[i]=0;
		}
		descA_cp[1]=-1;
		descA_cp[4]=nb; // *** important!
		descA_cp[5]=nb; // ***
	}

	// spread matrices
	pdgemr2d_(&n, &n, A_source, &one, &one, descA_source, A, &one, &one, descA, &context_all);

	if (failing_level>=0)
	{
		failing_level=n-failing_level;
	}

	// QR factorization
	result.core_start_time = time(NULL);
	pdgeqrf_cp_  (&n, &n, A, &one, &one, descA, A_cp, &one, &one, descA_cp, tau, tau_cp, &ltau, work, work_cp, &lwork, &checkpoint_freq, &failing_level, &context_all, &info );
	result.core_end_time = time(NULL);
	result.exit_code = info;

	// get matrix back
	pdgemr2d_ (&n, &n, A, &one, &one, descA, A_source, &one, &one, descA_source, &context_all);

	if (mpi_rank < cprocs)
	{
		// cleanup
		free(A);
		free(work);
		free(tau);
	}
	else
	{
		free(A_cp);
		free(work_cp);
		free(tau_cp);
	}

	MPI_Barrier(MPI_COMM_WORLD);

	result.total_end_time = time(NULL);

	return result;
}
