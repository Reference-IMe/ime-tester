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
 * WARNING: WITH transposition of the input (and output) matrix A to accomodate Fortran/C translation
 *
 */

test_output ScaLAPACK_pDGEQRF(	int n, double* A_global, double* B_global, int nb,	\
								int mpi_rank, int cprocs,							\
								int nprow, int npcol, int myrow, int mycol,			\
								int context, int context_global)
{
	test_output result = EMPTY_OUTPUT;

	result.total_start_time = time(NULL);

	// general
	int i;
	int i0 = 0;
	int i1 = 1;
	double d0 = 0.0;
	double d1 = 1.0;
	int info;
	double* work;
	double* tau;
	int lwork;
	double lazywork;

	// matrix
	int nr;
	int nc;
	double* A;
	double* At;
	int m=1; // B is a vector that will hold the solution vector for checking purposes
	int ncrhs, nrrhs;
	double *B;

	int descA_global[9];
	int descB_global[9];
	int descA[9];
	int descAt[9];
	int descB[9];
	int lld;
	int lld_global;

	if (mpi_rank < cprocs)
	{
		// Computation of local matrix size
		nr = numroc_( &n, &nb, &myrow, &i0, &nprow );
		nc = numroc_( &n, &nb, &mycol, &i0, &npcol );
		lld = MAX( 1 , nr );
		A = malloc(nr*nc*sizeof(double));
		At = malloc(nr*nc*sizeof(double));
		tau = malloc( nc*sizeof(double) );

		ncrhs = numroc_( &i1, &nb, &mycol, &i0, &npcol ); // one column vector
		nrrhs = numroc_( &n, &nb, &myrow, &i0, &nprow );
		B = malloc(nrrhs*ncrhs*sizeof(double));

		// Descriptors (local)
		descinit_( descA, &n, &n, &nb, &nb, &i0, &i0, &context, &lld, &info );
		descinit_( descAt, &n, &n, &nb, &nb, &i0, &i0, &context, &lld, &info );
		descinit_( descB, &n, &m, &nb, &nb, &i0, &i0, &context, &lld, &info );

		if (mpi_rank==0)
		{
			// Descriptors (global, for root node)
			lld_global = n;
			descinit_( descA_global, &n, &n, &i1, &i1, &i0, &i0, &context_global, &lld_global, &info );
			descinit_( descB_global, &n, &m, &i1, &i1, &i0, &i0, &context_global, &lld_global, &info );
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
		pdgemr2d_(&n, &n, A_global, &i1, &i1, descA_global, A, &i1, &i1, descA, &context);

		// transpose system matrix
		pdtran_(&n, &n, &d1, A, &i1, &i1, descA, &d0, At, &i1, &i1, descAt);

		// init work space
		lwork=-1;
		pdgeqrf_(  &n, &n, At, &i1, &i1, descAt, NULL, &lazywork, &lwork, &info );
		lwork = (int)lazywork;
		work = malloc( lwork*sizeof(double) );

		// QR factorization
		result.core_start_time = time(NULL);
		pdgeqrf_(  &n, &n, At, &i1, &i1, descAt, tau, work, &lwork, &info );
		result.core_end_time = time(NULL);
		result.exit_code = info;

		// transpose back
		pdtran_(&n, &n, &d1, At, &i1, &i1, descAt, &d0, A, &i1, &i1, descA);

		// get back
		pdgemr2d_ (&n, &n, A, &i1, &i1, descA, A_global, &i1, &i1, descA_global, &context);

		free(work);
	}
	else
	{
		A  = NULL;
		At = NULL;
		B  = NULL;
		tau  = NULL;
		work = NULL;
	}

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

		/*
		free(A);
		free(At);
		free(B);
		free(work);
		free(tau);
		*/
	}

	// cleanup
	NULLFREE(A);
	NULLFREE(At);
	NULLFREE(B);
	NULLFREE(work);
	NULLFREE(tau);

	MPI_Barrier(MPI_COMM_WORLD);
	return result;
}
