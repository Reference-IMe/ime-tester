/*
 * test_ScaLAPACK_pDGEQRF_cp_ft1.h
 *
 *  Created on: Jan 8, 2020
 *      Author: marcello
 */

#include <mpi.h>
#include <time.h>
#include "../helpers/macros.h"
#include "../helpers/matrix.h"
#include "ScaLAPACK/ScaLAPACK_pDGEQRF_ft1.h"
#include "tester_structures.h"

test_result test_ScaLAPACK_pDGEQRF_ft1(const char* label, int verbosity, parallel_env env, test_input input, int failing_level, int checkpoint_freq)
{
	test_result rank_result = TEST_NOT_RUN;
	test_result team_result = TEST_NOT_RUN;
	test_output output = EMPTY_OUTPUT;

	double* A;
	//double* bb;

	if (env.mpi_rank==0)
	{
		A=AllocateMatrix1D(input.n, input.n);
		//bb=AllocateMatrix1D(rows, nrhs);
		//FillMatrixT1D(A, rows, cols);
		//OneMatrix1D(bb, rows, nrhs);
		CopyMatrix1D(input.A_ref, A, input.n, input.n);

		if (verbosity>2)
		{
			printf("\n\n Matrix A:\n");
			PrintMatrix1D(A, input.n, input.n);
			//printf("\n Vector b:\n");
			//PrintMatrix1D(bb, rows, nrhs);
		}
	}

	output = ScaLAPACK_pDGEQRF_ft1(input.n, A, input.scalapack_bf, env.mpi_rank, input.calc_procs, input.spare_procs, \
									failing_level, checkpoint_freq, \
									env.blacs_nprow, env.blacs_npcol, env.blacs_row, env.blacs_col, \
									env.blacs_ctxt_grid, env.blacs_ctxt_root, env.blacs_ctxt_onerow, env.blacs_ctxt_spare);

	if (env.mpi_rank==0)
	{
		// check exit condition
		if (output.exit_code!=0)
		{
			printf("\n** Dangerous exit code.. (%d)**\n",output.exit_code);
		}
		// TODO calc error
		output.norm_rel_err = -99;

		if (verbosity>1)
		{
			printf("\nThe %s factorization is:\n",label);
			PrintMatrix1D(A, input.n, input.n);
			printf("\n with exit code     %d\n",output.exit_code);
			printf("      norm.rel.err. %f\n",output.norm_rel_err);
		}
		DeallocateMatrix1D(A);
	}

	TEST_END(output, rank_result, team_result);
}
