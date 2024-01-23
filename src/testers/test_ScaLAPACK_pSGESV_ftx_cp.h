/*
 * test_ScaLAPACK_pDGESV_ftx_cp.h
 *
 *  Created on: Dec 5, 2019
 *      Author: marcello
 */

#include <mpi.h>
#include <time.h>
#include "IMe/constants.h"
#include "../helpers/macros.h"
#include "../helpers/matrix_basic.h"
#include "../helpers/simple_dynamic_strings/sds.h"
#include "ScaLAPACK/ScaLAPACK_pDGESV_ftx_cp.h"
#include "tester_structures.h"

test_result test_ScaLAPACK_pDGESV_ftx_cp (	const char check,
											const char* tag,
											int verbosity,
											parallel_env env,
											test_input input,
											int fault_tolerance,
											int faulty_procs,
											int failing_rank,
											int failing_level,
											int checkpoint_freq	)
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

	int sqrt_calc_procs;
	sqrt_calc_procs=sqrt(env.calc_procs);

	int* failing_rank_list;

	int rank_calc_procs;

	rank_calc_procs=sqrt(env.calc_procs);

	if (check)
	{
		#include "test_ScaLAPACK_pre-check_ft.inc"
	}
	else
	{
		// failure list
		if (sqrt_calc_procs - (failing_rank % sqrt_calc_procs) < faulty_procs)
		{
			failing_rank=(floor(failing_rank / sqrt_calc_procs) +1 )*sqrt_calc_procs - faulty_procs;
		}
		failing_rank_list = malloc(faulty_procs*sizeof(int));

		for (i=failing_rank; i < failing_rank + faulty_procs; i++)
		{
			failing_rank_list[i-failing_rank]=i;
		}

		if (failing_level == 0)
		{
			failing_level = 1;
		}
		else if (failing_level >= (input.n - 1) )
		{
			failing_level = input.n - 2;
		}

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

		if (fault_tolerance < 1)
		{
			// if fault tolerance is disabled, force disable checkpoint and recovery
			output = ScaLAPACK_pDGESV_ftx_cp(input.n, A, input.nrhs, bb, input.scalapack_bf, env.mpi_rank, env.calc_procs, env.spare_procs,
													-1, faulty_procs, failing_rank_list, failing_level, 0,
													env.blacs_nprow, env.blacs_npcol, env.blacs_row, env.blacs_col,
													env.blacs_ctxt_grid, env.blacs_ctxt_root, env.blacs_ctxt_onerow, env.blacs_ctxt_spare);
		}
		else
		{
			// if fault tolerance is enabled, pass recovery option as set by user
			output = ScaLAPACK_pDGESV_ftx_cp(input.n, A, input.nrhs, bb, input.scalapack_bf, env.mpi_rank, env.calc_procs, env.spare_procs,
													checkpoint_freq, faulty_procs, failing_rank_list, failing_level, 1,
													env.blacs_nprow, env.blacs_npcol, env.blacs_row, env.blacs_col,
													env.blacs_ctxt_grid, env.blacs_ctxt_root, env.blacs_ctxt_onerow, env.blacs_ctxt_spare);
		}


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
