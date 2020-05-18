/*
 * ScaLAPACK_pDGESV_ft1.h
 *
 *  Created on: May 17, 2020
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


test_output ScaLAPACK_pDGEQRF_ft1(	int n, double* A_global, double* B_global, int nb,
									int mpi_rank, int cprocs, int sprocs, int failing_level, int checkpoint_freq,
									int nprow, int npcol, int myrow, int mycol,
									int context, int context_global, int context_all, int context_cp)
{
	test_output result = EMPTY_OUTPUT;

	result.total_start_time = time(NULL);

	// general
	int i;
	int i0 = 0;
	int i1 = 1;
	double d0 = 0.0;
	double d1 = 1.0;
	int nprocs = cprocs + sprocs;
	int info;
	
	double* work;
	double* work_cp;
	double* tau;
	double* tau_cp;
	int lwork;
	int ltau;
	double lazywork;

	// matrix
	int nr, nc;
	double *A;
	double *At;
	double *A_cp;
	int m = 1; // B is a vector that will hold the solution vector for checking purposes
	int ncrhs;
	int nrrhs;
	//int ncrhst;
	int nrrhst;
	double *B;

	int descA_global[9];
	int descB_global[9];
	int descA[9];
	int descAt[9];
	int descA_cp[9];
	int descB[9];

	int lld;
	int lldt;

	// Computation of local matrix size
	nc = numroc_( &n, &nb, &mycol, &i0, &npcol );
	nr = numroc_( &n, &nb, &myrow, &i0, &nprow );
	MPI_Allreduce( MPI_IN_PLACE, &nc, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD );
	MPI_Allreduce( MPI_IN_PLACE, &nr, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD );
	lld = MAX( 1 , nr );
	ltau = MAX( 1 , nc );

	ncrhs = numroc_( &m, &nb, &mycol, &i0, &npcol );
	nrrhs = numroc_( &n, &nb, &myrow, &i0, &nprow );

	//ncrhst = numroc_( &n, &nb, &mycol, &i0, &npcol );
	nrrhst = numroc_( &m, &nb, &myrow, &i0, &nprow );
	lldt = MAX( 1 , nrrhst );
	MPI_Allreduce( MPI_IN_PLACE, &lldt, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD );

	// global matrices
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

	// locally distributed matrices
	if (mpi_rank < cprocs)				// non-spare nodes
	{
		// Allocation
		A  = malloc(nr*nc*sizeof(double));
		At = malloc(nr*nc*sizeof(double));
		B  = malloc(nrrhs*ncrhs*sizeof(double));

		// Descriptors (local)
		descinit_( descA, &n, &n, &nb, &nb, &i0, &i0, &context, &lld, &info );
		descinit_( descAt, &n, &n, &nb, &nb, &i0, &i0, &context, &lld, &info );
		descinit_( descB, &n, &m, &nb, &nb, &i0, &i0, &context, &lld, &info );

		// init work space
		lwork=-1;
		pdgeqrf_( &n, &n, A, &i1, &i1, descA, NULL, &lazywork, &lwork, &info );
		lwork = (int)lazywork;
	}
	else								// spare node
	{
		// Allocation not needed
		A  = NULL;
		At = NULL;
		B  = NULL;

		// Descriptors
		for (i=0; i<9; i++)
		{
			descA[i]=0;
			descB[i]=0;
			descAt[i]=0;
		}
							// all processes have to know something about locally distributed matrices (not only non-spare nodes)
							// can't use descinint, due to illegal values of spare process not belonging to the right context
		descA[1]=-1;
		descA[4]=nb;
		descA[5]=nb;

		descAt[1]=-1;
		descAt[4]=nb;
		descAt[5]=nb;

		descB[1]=-1;
		descB[4]=nb;
		descB[5]=nb;

		lwork=-1;
	}

	// workspace
	MPI_Allreduce( MPI_IN_PLACE, &lwork, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD );
	MPI_Allreduce( MPI_IN_PLACE, &ltau, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD );
	work = malloc( lwork*sizeof(double) ); // also allocated on the spare proc because MPI_GATHER wants a buffer for everyone
	tau  = malloc( ltau*sizeof(double) );

	// checkpointing matrices
	if (mpi_rank==cprocs)	// spare (checkpointing) node
	{
		// Allocation
		A_cp    = malloc(n*n*sizeof(double));
		work_cp = malloc(lwork*nprocs*sizeof(double)); // with cprocs instead of nprocs is not good because MPI_GATHER wants a buffer for everyone!
		tau_cp  = malloc(ltau*nprocs*sizeof(double));

		// Descriptors
		descinit_( descA_cp, &n, &n, &i1, &i1, &i0, &i0, &context_cp, &n, &info );
	}
	else					// non-spare nodes
	{
		// Allocation not needed
		A_cp    = NULL;
		work_cp = NULL;
		tau_cp  = NULL;

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
	pdgemr2d_(&n, &n, A_global, &i1, &i1, descA_global, A, &i1, &i1, descA, &context_all);
	pdgemr2d_(&n, &m, B_global, &i1, &i1, descB_global, B, &i1, &i1, descB, &context_all);

	// transpose system matrix
	if (mpi_rank<cprocs) pdtran_(&n, &n, &d1, A, &i1, &i1, descA, &d0, At, &i1, &i1, descAt);

	result.core_start_time = time(NULL);

	if (failing_level>=0)
	{
		failing_level=n-failing_level;
	}

	// QR factorization
	result.core_start_time = time(NULL);
	pdgeqrf_cp_(&n, &n, At, &i1, &i1, descAt, A_cp, &i1, &i1, descA_cp, tau, tau_cp, &ltau, work, work_cp, &lwork, &checkpoint_freq, &failing_level, &context_all, &info );
	result.core_end_time = time(NULL);
	result.exit_code = info;

	// transpose back
	if (mpi_rank<cprocs) pdtran_(&n, &n, &d1, At, &i1, &i1, descAt, &d0, A, &i1, &i1, descA);

	// get matrix back
	pdgemr2d_(&n, &n, A, &i1, &i1, descA, A_global, &i1, &i1, descA_global, &context_all);

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

		pdgemr2d_(&n, &m, B_global, &i1, &i1, descB_global, B, &i1, &i1, descB, &context);

		lwork = -1;
		pdormqr_( "L", "T", &n, &m, &n, At, &i1, &i1, descAt, tau, B, &i1, &i1, descB, &lazywork, &lwork, &info );
		lwork = (int) lazywork;
		work = (double*) malloc( lwork * sizeof(double) );

		pdormqr_( "L", "T", &n, &m, &n, At, &i1, &i1, descAt, tau, B, &i1, &i1, descB, work, &lwork, &info);

		pdtrsm_("L", "U", "N", "N", &n, &m, &d1, At, &i1, &i1, descAt, B, &i1, &i1, descB);

		// collect result
		pdgemr2d_(&n, &m, B, &i1, &i1, descB, B_global, &i1, &i1, descB_global, &context);
	}

	// cleanup
	NULLFREE(A);
	NULLFREE(At);
	NULLFREE(B);
	NULLFREE(work);
	NULLFREE(tau);
	NULLFREE(work_cp);
	NULLFREE(tau_cp);
	NULLFREE(A_cp);

	MPI_Barrier(MPI_COMM_WORLD);
	return result;
}
