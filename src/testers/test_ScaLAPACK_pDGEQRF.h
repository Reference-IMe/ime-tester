/*
 * test_ScaLAPACK_pDGEQRF.h
 *
 *  Created on: Dec 5, 2019
 *      Author: marcello
 */

#include <mpi.h>
#include <time.h>
#include "../helpers/info.h"
#include "../helpers/macros.h"
#include "../helpers/matrix.h"
#include "ScaLAPACK/ScaLAPACK_pDGEQRF.h"

duration_t test_ScaLAPACK_pDGEQRF(const char* label, int verbosity, int rows, int cols, int nb, int rank, int cprocs)
{
	duration_t timing, timing_max;
	test_output info;

	double* A;

	if (rank==0)
	{
		A=AllocateMatrix1D(rows, cols);
		FillMatrixT1D(A, rows, cols);

		if (verbosity>2)
		{
			printf("\n\n Matrix A:\n");
			PrintMatrix1D(A, rows, cols);
		}
	}

	info = ScaLAPACK_pDGEQRF_calc(rows, A, nb, rank, cprocs);

	if (rank==0)
	{
		if (verbosity>1)
		{
			printf("\nThe %s factorization is:\n",label);
			PrintMatrix1D(A, rows, cols);
		}
		DeallocateMatrix1D(A);
	}

	TEST_END(info, timing, timing_max);
}
