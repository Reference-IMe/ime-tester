/*
 * test_ScaLAPACK_pDGESV.h
 *
 *  Created on: Dec 27, 2019
 *      Author: marcello
 */

#include <mpi.h>
#include <time.h>
#include "../helpers/info.h"
#include "../helpers/macros.h"
#include "../helpers/matrix.h"
#include "../helpers/matrix_advanced.h"
#include "ScaLAPACK/ScaLAPACK_pDGESV.h"

duration_t test_ScaLAPACK_pDGESV(const char* label, int verbosity, int rows, int cols, int nrhs, int nb, int rank, int cprocs)
{
	duration_t timing, timing_max;
	result_info info;

	double* A;
	double* bb;

	if (rank==0)
	{
		A=AllocateMatrix1D(rows, cols);
		bb=AllocateMatrix1D(rows, nrhs);
		//FillMatrix1D(A, rows, cols);
		//ReferenceMatrix1D(A, rows, cols);
		RandomMatrix1D(A, rows, cols, 1);
		OneMatrix1D(bb, rows, nrhs);

		if (verbosity>2)
		{
			printf("\n\n Matrix A:\n");
			PrintMatrix1D(A, rows, cols);
			printf("\n Vector b:\n");
			PrintMatrix1D(bb, rows, nrhs);
		}

		//OrthogonalizeMatrix1D(A, rows, cols);
		RandomSquareMatrix1D_cnd(A, rows, 1, 10);
		printf("\n\n GOOD-Matrix A:\n");
		PrintMatrix1D(A, rows, cols);
		printf("\n\n CND A: %g\n",ConditionNumber1D( A, rows, cols));
	}

	info = ScaLAPACK_pDGESV_calc(rows, A, nrhs, bb, nb, rank, cprocs);

	if (rank==0)
	{
		if (verbosity>1)
		{
			printf("\nThe %s solution is:\n",label);
			PrintMatrix1D(bb, rows, nrhs);
		}
		DeallocateMatrix1D(A);
		DeallocateMatrix1D(bb);
	}

	TEST_END(info, timing, timing_max);
}
