/*
 * test_ScaLAPACK_pDGETRF_cp_ft1.h
 *
 *  Created on: Dec 27, 2019
 *      Author: marcello
 */

#include <mpi.h>
#include <time.h>
#include "../helpers/info.h"
#include "../helpers/macros.h"
#include "../helpers/matrix.h"
#include "ScaLAPACK/ScaLAPACK_pDGETRF_cp_ft1_sim.h"

run_info test_ScaLAPACK_pDGETRF_cp_ft1_sim(const char* label, int verbosity, int n, double* A_ref, int nb, int rank, int cprocs, int sprocs, int failing_level, int checkpoint_freq)
{
	run_info process_info = {0, 0, 0, 0};
	run_info team_info = {0, 0, 0, 0};
	result_info info = {0, 0, 0, 0, 0, 0};

	double* A;
	//double* bb;

	if (rank==0)
	{
		A=AllocateMatrix1D(n, n);
		//bb=AllocateMatrix1D(rows, nrhs);
		//FillMatrixT1D(A, rows, cols);
		//OneMatrix1D(bb, rows, nrhs);
		CopyMatrix1D(A_ref, A, n, n);

		if (verbosity>2)
		{
			printf("\n\n Matrix A:\n");
			PrintMatrix1D(A, n, n);
			//printf("\n Vector b:\n");
			//PrintMatrix1D(bb, rows, nrhs);
		}
	}

	info = ScaLAPACK_pDGETRF_cp_ft1_sim(n, A, nb, rank, cprocs, sprocs, failing_level, checkpoint_freq);

	if (rank==0)
	{
		if (verbosity>1)
		{
			printf("\nThe %s factorization is:\n",label);
			PrintMatrix1D(A, n, n);
		}
		DeallocateMatrix1D(A);
		//DeallocateMatrix1D(bb);
	}

	TEST_END(info, process_info, team_info);
}
