/*
 * test_ScaLAPACK_pDGETRF.h
 *
 *  Created on: Dec 27, 2019
 *      Author: marcello
 */

#include <mpi.h>
#include <time.h>
#include "../helpers/info.h"
#include "../helpers/macros.h"
#include "../helpers/matrix.h"
#include "ScaLAPACK/ScaLAPACK_pDGETRF.h"

duration_t test_ScaLAPACK_pDGETRF(const char* label, int verbosity, int rows, int cols, int nb, int rank, int cprocs)
{
	duration_t timing, timing_max;
	result_info info;

	double* A;
	//double* bb;

	if (rank==0)
	{
		A=AllocateMatrix1D(rows, cols);
		//bb=AllocateMatrix1D(rows, nrhs);
		FillMatrixT1D(A, rows, cols);
		//OneMatrix1D(bb, rows, nrhs);

		if (verbosity>2)
		{
			printf("\n\n Matrix A:\n");
			PrintMatrix1D(A, rows, cols);
			//printf("\n Vector b:\n");
			//PrintMatrix1D(bb, rows, nrhs);
		}
	}

	// Scalapack_pDGETRF(rows, A, nrhs, bb, rank, cprocs, sprocs); // for consistency checking
	info = ScaLAPACK_pDGETRF_calc(rows, A, nb, rank, cprocs);

	if (rank==0)
	{
		if (verbosity>1)
		{
			printf("\nThe %s factorization is:\n",label);
			PrintMatrix1D(A, rows, cols);
		}
		DeallocateMatrix1D(A);
		//DeallocateMatrix1D(bb);
	}

	TEST_END(info, timing, timing_max);
}
