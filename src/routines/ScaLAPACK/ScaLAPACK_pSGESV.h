/*
 * ScaLAPACK_pDGESV.h
 *
 *  Created on: Dec 27, 2019
 *      Author: marcello
 */

#include <mpi.h>
#include <time.h>
#include "../../helpers/macros.h"
#include "../IMe/lib/src/helpers/vector_basic.h"
#include "../../helpers/Cblacs.h"
#include "../IMe/lib/src/helpers/matrix_basic.h"
#include "../../helpers/scalapack.h"
#include "../../tester_structures.h"


exit_status ScaLAPACK_pSGESV(int n, float* A_global, int m, float* B_global, int nb,
								int mpi_rank, int cprocs,
								int nprow, int npcol, int myrow, int mycol,
								int context, int context_global)
{
	exit_status result = EMPTY_OUTPUT;

	result.total_start_time = time(NULL);

	/*
	 * n = system rank (A_global n x n)
	 * m = num. of r.h.s (B_global n x m)
	 */

	// general
	int i;
	int i0 = 0;
	int i1 = 1;
	float d0 = 0.0;
	float d1 = 1.0;
	int info;
	int *ipiv;

	// matrix
	int nr, nc;
	float *A;
	float *At;
	int ncrhs, nrrhs;
	float *B;
	int ncrhst, nrrhst;
	float *Bt;

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
		A  = malloc(nr*nc*sizeof(float));
		At = malloc(nr*nc*sizeof(float));

		ncrhs = numroc_( &m, &nb, &mycol, &i0, &npcol );
		nrrhs = numroc_( &n, &nb, &myrow, &i0, &nprow );
		B = malloc(nrrhs*ncrhs*sizeof(float));

		ncrhst = numroc_( &n, &nb, &mycol, &i0, &npcol );
		nrrhst = numroc_( &m, &nb, &myrow, &i0, &nprow );
		lldt = MAX( 1 , nrrhst );
		Bt = malloc(nrrhst*ncrhst*sizeof(float));

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

		// spread matrices
		psgemr2d_(&n, &n, A_global, &i1, &i1, descA_global, A, &i1, &i1, descA, &context);
		psgemr2d_(&n, &m, B_global, &i1, &i1, descB_global, B, &i1, &i1, descB, &context);

		// transpose system matrix
		pstran_(&n, &n, &d1, A, &i1, &i1, descA, &d0, At, &i1, &i1, descAt);

		// linear system equations solver
		result.core_start_time = time(NULL);
		psgesv_( &n, &m, At, &i1, &i1, descAt, ipiv, B, &i1, &i1, descB, &info );
		result.core_end_time = time(NULL);
		result.exit_code = info;

		// re-transpose result
		pstran_(&m, &n, &d1, B, &i1, &i1, descB, &d0, Bt, &i1, &i1, descBt);

		// collect result
		if (mpi_rank==0)
		{
			// Descriptors (global)
			descinit_( descB_global, &m, &n, &i1, &i1, &i0, &i0, &context_global, &m, &info );
		}
		psgemr2d_(&m, &n, Bt, &i1, &i1, descBt, B_global, &i1, &i1, descB_global, &context);
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

	result.total_end_time = time(NULL);

	return result;
}
