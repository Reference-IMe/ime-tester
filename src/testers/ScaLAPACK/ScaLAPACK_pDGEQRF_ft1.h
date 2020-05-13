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

test_output ScaLAPACK_pDGEQRF_ft1(	int n, double* A_global, int nb,												\
									int mpi_rank, int cprocs, int sprocs, int failing_level, int checkpoint_freq,	\
									int nprow, int npcol, int myrow, int mycol,										\
									int ctxt, int ctxt_root, int ctxt_onerow, int ctxt_cp)
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
	int nprocs = cprocs + sprocs;
	int info;
	double* work;
	double* work_cp;
	double* tau;
	double* tau_cp;
	int lwork=-1;
	int ltau;
	double lazywork;

	// matrix
	int nr;
	int nc;
	int lldA;
	double* A;
	double* A_cp;
	int descA_global[9];
	int descA[9];
	int descA_cp[9];

	// Computation of local matrix size
	nr = numroc_( &n, &nb, &myrow, &i0, &nprow );
	nc = numroc_( &n, &nb, &mycol, &i0, &npcol );

	lldA = MAX( 1 , nr );
	MPI_Allreduce( MPI_IN_PLACE, &lldA, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD );

	ltau = lldA;

	// global matrices
	if (mpi_rank==0)
	{
		// Descriptors (global, for root node)
		descinit_( descA_global, &n, &n, &i1, &i1, &i0, &i0, &ctxt_root, &n, &info );
	}
	else
	{
		// Descriptors (global, for non-root nodes)
		A_global=NULL;
		for (i=0; i<9; i++)
		{
			descA_global[i]=0;
		}
		descA_global[1]=-1;
	}

	// locally distributed matrices
	if (mpi_rank < cprocs)				// non-spare nodes
	{
		// Allocation
		A = malloc(nr*nc*sizeof(double));

		// Descriptors (local)
		// all processes have to know descA (not only non-spare nodes)
		descinit_( descA, &n, &n, &nb, &nb, &i0, &i0, &ctxt, &lldA, &info );

		// init work space
		pdgeqrf_(  &n, &n, A, &i1, &i1, descA, NULL, &lazywork, &lwork, &info );
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
		descinit_( descA_cp, &n, &n, &i1, &i1, &i0, &i0, &ctxt_cp, &n, &info );
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
	pdgemr2d_(&n, &n, A_global, &i1, &i1, descA_global, A, &i1, &i1, descA, &ctxt_onerow);

	if (failing_level>=0)
	{
		failing_level=n-failing_level;
	}

	// QR factorization
	result.core_start_time = time(NULL);
	pdgeqrf_cp_  (&n, &n, A, &i1, &i1, descA, A_cp, &i1, &i1, descA_cp, tau, tau_cp, &ltau, work, work_cp, &lwork, &checkpoint_freq, &failing_level, &ctxt_onerow, &info );
	result.core_end_time = time(NULL);
	result.exit_code = info;

	// get matrix back
	pdgemr2d_ (&n, &n, A, &i1, &i1, descA, A_global, &i1, &i1, descA_global, &ctxt_onerow);

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
