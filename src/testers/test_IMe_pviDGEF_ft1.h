/*
 * test_IMe_pviDGESV_ft1.h
 *
 *  Created on: Nov 20, 2019
 *      Author: marcello
 */

#include <mpi.h>
#include <time.h>
#include "../helpers/macros.h"
#include "../helpers/matrix.h"
#include "../ft-simulated/pviDGEF_WO_ft1.h"
#include "../helpers/tester_structures.h"

duration_t test_IMe_pviDGEF_ft1_sim(const char* label, int verbosity, int rows, int cols, int rank, int sprocs, int failing_rank, int failing_level)
{
	duration_t timing, timing_max;
	test_output info;

	int r,c;
	double** A2;
	double** K;

	if (rank==0)
	{
		A2=AllocateMatrix2D(rows,cols,CONTIGUOUS);
		FillMatrix2D(A2, rows, cols);

		K=AllocateMatrix2D(rows, cols, CONTIGUOUS);

		if (verbosity>2)
		{
			printf("\n\n Matrix A:\n");
			PrintMatrix2D(A2, rows, cols);
		}
	}
	else
	{
		A2=AllocateMatrix2D(0, 0, CONTIGUOUS);
		K=AllocateMatrix2D(0, 0, CONTIGUOUS);
	}

	info = pviDGEF_WO_ft1_sim(rows, A2, K, MPI_COMM_WORLD, sprocs, failing_rank, failing_level);

	// clean K for output
	if (rank==0)
	{
		for (r=0;r<rows;r++)
		{
			for (c=r;c<cols;c++)
			{
				K[r][c]=0;
			}
		}
	}

	if (rank==0 && verbosity>1)
	{
		printf("\nThe %s X factor is:\n",label);
		PrintMatrix2D(A2, rows, cols);
		printf("\nThe %s K factor is:\n",label);
		PrintMatrix2D(K, rows, cols);
	}

	if (rank==0)
	{
		DeallocateMatrix2D(A2, rows, CONTIGUOUS);
		DeallocateMatrix2D(K, rows, CONTIGUOUS);
	}
	else
	{
		DeallocateMatrix2D(A2, 0, CONTIGUOUS);
		DeallocateMatrix2D(K, 0, CONTIGUOUS);
	}

	TEST_END(info, timing, timing_max);
}
