/*
 * test_ScaLAPACK_pDGESV_cp_ft1.h
 *
 *  Created on: Dec 5, 2019
 *      Author: marcello
 */

#include <mpi.h>
#include <time.h>
#include "../helpers/info.h"
#include "../helpers/macros.h"
#include "../helpers/matrix.h"
#include "ScaLAPACK/ScaLAPACK_pDGESV_cp_ft1_sim.h"

duration_t test_ScaLAPACK_pDGESV_cp_ft1_sim(const char* label, int verbosity, int rows, int cols, int nb, int nrhs, int rank, int cprocs, int sprocs, int failing_level)
{
	duration_t timing, timing_max;
	result_info info;

	double* A;
	double* bb;

	if (rank==0)
	{
		A=AllocateMatrix1D(rows, cols);
		bb=AllocateMatrix1D(rows, nrhs);
		FillMatrixT1D(A, rows, cols);
		OneMatrix1D(bb, rows, nrhs);

		if (verbosity>2)
		{
			printf("\n\n Matrix A:\n");
			PrintMatrix1D(A, rows, cols);
			printf("\n Vector b:\n");
			PrintMatrix1D(bb, rows, nrhs);
		}
	}

	info = ScaLAPACK_pDGESV_cp_ft1_sim(rows, A, nrhs, bb, nb, rank, cprocs, sprocs, failing_level);

	if (rank==0)
	{
		if (verbosity>1)
		{
			printf("\nThe %s solution is:\n",label);
			//PrintMatrix1D(bb, nrhs, rows);
			PrintMatrix1D(bb, nrhs, rows);
		}
		DeallocateMatrix1D(A);
		DeallocateMatrix1D(bb);
	}

	TEST_END(info, timing, timing_max);
}
