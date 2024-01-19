/*
 * test_ScaLAPACK_pDGESV.h
 *
 *  Created on: Dec 27, 2019
 *      Author: marcello
 */

#include <mpi.h>
#include <time.h>
#include "../helpers/macros.h"
#include "../helpers/matrix_advanced.h"
#include "../helpers/matrix_basic.h"
#include "ScaLAPACK/ScaLAPACK_pDGETRF.h"
#include "tester_structures.h"

test_result test_ScaLAPACK_pDGETRF (	const char check,
										const char* tag,
										int verbosity,
										parallel_env env,
										test_input input	)
{
	test_result rank_result = TEST_NOT_RUN;
	test_result team_result = TEST_NOT_RUN;
	test_output output      = EMPTY_OUTPUT;

	sds label=sdsempty();
	TAG2LABEL(tag,label);

	int i;

	double* A;
	double* bb;

	if (check)
	{
		#include "test_ScaLAPACK_pre-check.inc"
	}
	else
	{
		if (env.mpi_rank==0)
		{
			A=AllocateMatrix1D_double(input.n, input.n);
			bb=AllocateMatrix1D_double(input.n, 1);

			CopyMatrix1D_double(input.A_ref_d, A, input.n, input.n);

			for (i=0;i<input.n;i++)
			{
				bb[i] = input.b_ref_d[i];
			}
			if (verbosity>2)
			{
				printf("\n\n Matrix A:\n");
				PrintMatrix1D_double(A, input.n, input.n);
			}
		}
		else
		{
			A  = NULL;
			bb = NULL;
		}

		output = ScaLAPACK_pDGETRF(input.n, A, bb, input.scalapack_bf, env.mpi_rank, env.calc_procs,
									env.blacs_nprow, env.blacs_npcol, env.blacs_row, env.blacs_col,
									env.blacs_ctxt_grid, env.blacs_ctxt_root);

		if (env.mpi_rank==0)
		{
			#include "test_ScaLAPACK_post-check-factorization.inc"
			#include "test_ScaLAPACK_show-factorization.inc"
		}

		NULLFREE(A);
		NULLFREE(bb);
	}
	sdsfree(label);
	TEST_END(output, rank_result, team_result);
}
