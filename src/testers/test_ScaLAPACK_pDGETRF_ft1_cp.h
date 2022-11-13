/*
 * test_ScaLAPACK_pDGESV_ft1.h
 *
 *  Created on: Dec 5, 2019
 *      Author: marcello
 */

#include <mpi.h>
#include <time.h>
#include "../helpers/macros.h"
#include "../helpers/matrix.h"
#include "ScaLAPACK/ScaLAPACK_pDGETRF_ft1.h"
#include "tester_structures.h"

test_result test_ScaLAPACK_pDGETRF_ft1(const char check, const char* label, int verbosity, parallel_env env, test_input input, int fault_tolerance, int faulty_procs, int failing_level, int checkpoint_freq)
{
	test_result rank_result = TEST_NOT_RUN;
	test_result team_result = TEST_NOT_RUN;
	test_output output      = EMPTY_OUTPUT;

	int i;

	double* A;
	double* bb;

	int rank_calc_procs;

	rank_calc_procs=sqrt(env.calc_procs);

	if (check)
	{
		#include "test_ScaLAPACK_pre-check_ft.inc"
	}
	else
	{
		if (env.mpi_rank==0)
		{
			A=AllocateMatrix1D(input.n, input.n);
			bb=AllocateMatrix1D(input.n, 1);

			CopyMatrix1D(input.A_ref, A, input.n, input.n);

			for (i=0;i<input.n;i++)
			{
				bb[i] = input.b_ref[i];
			}
			if (verbosity>2)
			{
				printf("\n\n Matrix A:\n");
				PrintMatrix1D(A, input.n, input.n);
			}
		}
		else
		{
			A  = NULL;
			bb = NULL;
		}

		output = ScaLAPACK_pDGETRF_ft1(input.n, A, bb, input.scalapack_bf, env.mpi_rank, env.calc_procs, env.spare_procs,
										failing_level, checkpoint_freq,
										env.blacs_nprow, env.blacs_npcol, env.blacs_row, env.blacs_col,
										env.blacs_ctxt_grid, env.blacs_ctxt_root, env.blacs_ctxt_onerow, env.blacs_ctxt_spare[0]);

		if (env.mpi_rank==0)
		{
			#include "test_ScaLAPACK_post-check-factorization.inc"
			#include "test_ScaLAPACK_show-factorization.inc"
		}

		NULLFREE(A);
		NULLFREE(bb);
	}
	TEST_END(output, rank_result, team_result);
}
