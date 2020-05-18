/*
 * ScaLAPACK_pDGETRF_ft1.h
 *
 *  Created on: May 15, 2020
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


test_output ScaLAPACK_pDGETRF_ft1(	int n, double* A_global, double* B_global, int nb,
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
	int m = 1;
	double d0 = 0.0;
	double d1 = 1.0;
	int nprocs = cprocs + sprocs;
	int info;
	int *ipiv;
	int *ipiv_cp;
	int nipiv;
	
	// matrix
	int nr, nc;
	double *A;
	double *At;
	int ncrhs, nrrhs;
	double *B;
	int ncrhst, nrrhst;
	double *Bt;
	double *A_cp;

	int descA_global[9];
	int descB_global[9];
	int descA[9];
	int descAt[9];
	int descB[9];
	int descBt[9];
	int descA_cp[9];

	int lld, lldt;

	// Computation of local matrix size
	nc = numroc_( &n, &nb, &mycol, &i0, &npcol );
	nr = numroc_( &n, &nb, &myrow, &i0, &nprow );
	lld = MAX( 1 , nr );
	MPI_Allreduce( MPI_IN_PLACE, &lld, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD );

	nipiv = (lld+nb);
	ipiv = malloc(nipiv*sizeof(int));// also allocated on the spare proc because MPI_GATHER wants a buffer for everyone

	ncrhs = numroc_( &m, &nb, &mycol, &i0, &npcol );
	nrrhs = numroc_( &n, &nb, &myrow, &i0, &nprow );

	ncrhst = numroc_( &n, &nb, &mycol, &i0, &npcol );
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

	// checkpointing matrices
	if (mpi_rank==cprocs) // spare (checkpointing) node
	{
		// Allocation
		A_cp    = malloc(n*n*sizeof(double));
		ipiv_cp = malloc(nipiv*nprocs*sizeof(int)); //ipiv_cp=malloc(nIPIV*cprocs*sizeof(int)); // with cprocs is not good because MPI_GATHER wants a buffer for everyone

		// Descriptors
		descinit_( descA_cp, &n, &n, &i1, &i1, &i0, &i0, &context_cp, &n, &info );
	}
	else				// non-spare nodes
	{
		// Allocation not needed
		A_cp    = NULL;
		ipiv_cp = NULL;

		// Descriptors
		for (i=0; i<9; i++)
		{
			descA_cp[i]=0;
		}
		descA_cp[1]=-1;
		descA_cp[4]=nb; // *** important!
		descA_cp[5]=nb; // ***
	}

	// locally distributed matrices
	if (mpi_rank < cprocs)				// non-spare nodes
	{
		// Allocation
		A  = malloc(nr*nc*sizeof(double));
		At = malloc(nr*nc*sizeof(double));
		B  = malloc(nrrhs*ncrhs*sizeof(double));
		Bt = malloc(nrrhst*ncrhst*sizeof(double));

		// Descriptors (local)
		descinit_( descA, &n, &n, &nb, &nb, &i0, &i0, &context, &lld, &info );
		descinit_( descAt, &n, &n, &nb, &nb, &i0, &i0, &context, &lld, &info );
		descinit_( descB, &n, &m, &nb, &nb, &i0, &i0, &context, &lld, &info );
		descinit_( descBt, &m, &n, &nb, &nb, &i0, &i0, &context, &lldt, &info );
	}
	else								// spare node
	{
		// Allocation not needed
		A  = NULL;
		At = NULL;
		B  = NULL;
		Bt = NULL;

		// Descriptors
		for (i=0; i<9; i++)
		{
			descA[i]=0;
			descB[i]=0;
			descAt[i]=0;
			descBt[i]=0;
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

		descBt[1]=-1;
		descBt[4]=nb;
		descBt[5]=nb;
	}

	// spread matrices
	pdgemr2d_(&n, &n, A_global, &i1, &i1, descA_global, A, &i1, &i1, descA, &context_all);
	pdgemr2d_(&n, &m, B_global, &i1, &i1, descB_global, B, &i1, &i1, descB, &context_all);

	// transpose system matrix
	if (mpi_rank<cprocs) pdtran_(&n, &n, &d1, A, &i1, &i1, descA, &d0, At, &i1, &i1, descAt);

	result.core_start_time = time(NULL);

	// Linear system equations solver

	if (failing_level>=0)
	{
		failing_level=n-failing_level;
	}

	// checkpointed factorization called by everyone
	pdgetrf_cp_(&n, &n, At, &i1, &i1, descAt, A_cp, &i1, &i1, descA_cp, ipiv, ipiv_cp, &nipiv, &checkpoint_freq, &failing_level, &context_all, &info );

	result.core_end_time = time(NULL);
	result.exit_code = info;

	if (mpi_rank < cprocs)
	{
		// transpose back A
		pdtran_(&n, &n, &d1, At, &i1, &i1, descAt, &d0, A, &i1, &i1, descA);

		// get back A
		pdgemr2d_ (&n, &n, A, &i1, &i1, descA, A_global, &i1, &i1, descA_global, &context);

		result.total_end_time = time(NULL);

		/*
		 * evaluate n.r.e.
		 */
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
		result.total_end_time = time(NULL);
	}

	// cleanup
	NULLFREE(A);
	NULLFREE(At);
	NULLFREE(B);
	NULLFREE(Bt);
	NULLFREE(ipiv);
	NULLFREE(ipiv_cp);
	NULLFREE(A_cp);

	MPI_Barrier(MPI_COMM_WORLD);
	return result;
}
