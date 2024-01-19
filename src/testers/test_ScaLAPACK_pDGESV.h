/*
 * test_ScaLAPACK_pDGESV.h
 *
 *  Created on: Dec 27, 2019
 *      Author: marcello
 */

#include <mpi.h>
#include <time.h>
#include "../constants.h"
#include "../helpers/macros.h"
#include "../helpers/matrix_advanced.h"
#include "../helpers/matrix_basic.h"
#include "../helpers/simple_dynamic_strings/sds.h"
#include "tester_structures.h"
#include "ScaLAPACK/ScaLAPACK_pDGESV.h"

test_result test_ScaLAPACK_pDGESV (	const char check,
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

	int i,j;

	double* A;
	double* bb;
	double* xx_ref;

	if (check)
	{
		#include "test_ScaLAPACK_pre-check.inc"
	}
	else
	{
		if (env.mpi_rank==0)
		{
			A=AllocateMatrix1D_double(input.n, input.n);
			bb=AllocateMatrix1D_double(input.n, input.nrhs);
			xx_ref=AllocateMatrix1D_double(input.n, input.nrhs);

			CopyMatrix1D_double(input.A_ref_d, A, input.n, input.n);
			for (i=0;i<input.n;i++)
			{
				for (j=0;j<input.nrhs;j++)
				{
					bb[i*input.nrhs+j] = input.b_ref_d[i];
					xx_ref[i*input.nrhs+j] = input.x_ref_d[i];
				}
			}

			if (verbosity>2)
			{
				printf("\n\n Matrix A:\n");
				PrintMatrix1D_double(A, input.n, input.n);
				printf("\n Vector b:\n");
				PrintMatrix1D_double(bb, input.n, input.nrhs);
			}
		}
		else
		{
			A      = NULL;
			bb     = NULL;
			xx_ref = NULL;
		}

		output = ScaLAPACK_pDGESV(input.n, A, input.nrhs, bb, input.scalapack_bf, env.mpi_rank, env.calc_procs,
									env.blacs_nprow, env.blacs_npcol, env.blacs_row, env.blacs_col,
									env.blacs_ctxt_grid, env.blacs_ctxt_root);

		if (env.mpi_rank==0)
		{
			#define TYPE REAL_DOUBLE
			#include "test_ScaLAPACK_post-check-solution.inc"
			#include "test_ScaLAPACK_show-solution.inc"
			#undef TYPE
		}

		NULLFREE(A);
		NULLFREE(bb);
		NULLFREE(xx_ref);
	}
	sdsfree(label);
	TEST_END(output, rank_result, team_result);
}
