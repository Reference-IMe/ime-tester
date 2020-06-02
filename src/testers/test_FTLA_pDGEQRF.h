/*
 * test_FTLA_pDGEQRF.h
 *
 *  Created on: Dec 28, 2019
 *      Author: marcello
 */

#include <mpi.h>
#include <time.h>
#include "../helpers/macros.h"
#include "../helpers/matrix.h"
#include "tester_structures.h"
#include "FTLA/FTLA_pDGEQRF.h"

test_result test_FTLA_pDGEQRF(const char* label, int verbosity, parallel_env env, test_input input, int faults)
{
	test_result rank_result = TEST_NOT_RUN;
	test_result team_result = TEST_NOT_RUN;
	test_output output      = EMPTY_OUTPUT;

	double* A;
	double* bb;
	int i;

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

	output = FTLA_ftdqr(input.n, A, bb, input.scalapack_bf, env.mpi_rank, input.calc_procs,
							faults,
							env.blacs_nprow, env.blacs_npcol, env.blacs_row, env.blacs_col,
							env.blacs_ctxt_grid, env.blacs_ctxt_root);

	if (env.mpi_rank==0)
	{
		// check exit condition
		if (output.exit_code!=0)
		{
			printf("\n** Dangerous exit code.. (%d)**\n",output.exit_code);
		}
		// calc error
		if (input.calc_nre) rank_result.norm_rel_err = NormwiseRelativeError1D(bb, input.x_ref, input.n, 1);

		if (verbosity>1)
		{
			printf("\nThe %s factorization is:\n",label);
			PrintMatrix1D(A, input.n, input.n);
			printf("\n with exit code     %d\n",output.exit_code);
			printf("      norm.rel.err. %f\n",rank_result.norm_rel_err);
		}
	}

	NULLFREE(A);
	NULLFREE(bb);

	TEST_END(output, rank_result, team_result);
}
