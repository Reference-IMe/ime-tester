/*
 * test_IMe_pviDGESV_cs.h
 *
 *  Created on: Dec 5, 2019
 *      Author: marcello
 */

#include <mpi.h>
#include <time.h>
#include "../helpers/info.h"
#include "../helpers/macros.h"
#include "../helpers/matrix.h"
#include "../ft-simulated/pviDGESV_WO_cs.h"

duration_t test_IMe_pviDGESV_cs(const char* label, int verbosity, int rows, int cols, int nrhs, int rank, int sprocs)
{
	duration_t timing, timing_max;
	result_info info;

	double** A2;
	double** bb;
	double** xx;

	xx=AllocateMatrix2D(rows, nrhs, CONTIGUOUS);
	bb=AllocateMatrix2D(rows, nrhs, CONTIGUOUS);

	if (rank==0)
	{
		A2=AllocateMatrix2D(rows,cols,CONTIGUOUS);

		FillMatrix2D(A2, rows, cols);

		OneMatrix2D(bb, rows, nrhs);

		if (verbosity>2)
		{
			printf("\n\n Matrix A:\n");
			PrintMatrix2D(A2, rows, cols);
			printf("\n Vector b:\n");
			PrintMatrix2D(bb, rows, nrhs);
		}
	}
	else
	{
		A2=AllocateMatrix2D(0, 0, CONTIGUOUS);
	}

	info = pviDGESV_WO_cs(rows, A2, nrhs, bb, xx, MPI_COMM_WORLD, sprocs);

	if (rank==0 && verbosity>1)
	{
		printf("\nThe %s solution is:\n",label);
		PrintMatrix2D(xx, rows, nrhs);
	}

	DeallocateMatrix2D(xx,rows,CONTIGUOUS);
	DeallocateMatrix2D(bb,rows,CONTIGUOUS);
	if (rank==0)
	{
		DeallocateMatrix2D(A2,rows,CONTIGUOUS);
	}
	else
	{
		DeallocateMatrix2D(A2, 0, CONTIGUOUS);
	}

	TEST_END(info, timing, timing_max);
}
