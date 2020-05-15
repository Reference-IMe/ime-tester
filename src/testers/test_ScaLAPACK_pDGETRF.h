/*
 * test_ScaLAPACK_pDGESV.h
 *
 *  Created on: Dec 27, 2019
 *      Author: marcello
 */

#include <mpi.h>
#include <time.h>
#include "../helpers/macros.h"
#include "../helpers/matrix.h"
#include "../helpers/matrix_advanced.h"
#include "ScaLAPACK/ScaLAPACK_pDGETRF.h"
#include "tester_structures.h"

test_result test_ScaLAPACK_pDGETRF(const char* label, int verbosity, parallel_env env, test_input input)
{
	test_result rank_result = TEST_NOT_RUN;
	test_result team_result = TEST_NOT_RUN;
	test_output output = EMPTY_OUTPUT;

	int i,j;

	double* A;
	double* bb;
	double* xx_ref;

	if (env.mpi_rank==0)
	{
		A=AllocateMatrix1D(input.n, input.n);
		bb=AllocateMatrix1D(input.n, input.nrhs);
		xx_ref=AllocateMatrix1D(input.n, input.nrhs);

		CopyMatrix1D(input.A_ref, A, input.n, input.n);
		for (i=0;i<input.n;i++)
		{
			for (j=0;j<input.nrhs;j++)
			{
				bb[i*input.nrhs+j] = input.b_ref[i];
				xx_ref[i*input.nrhs+j] = input.x_ref[i];
			}
		}


		if (verbosity>2)
		{
			printf("\n\n Matrix A:\n");
			PrintMatrix1D(A, input.n, input.n);
			printf("\n Vector b:\n");
			PrintMatrix1D(bb, input.n, input.nrhs);
		}
	}

	output = ScaLAPACK_pDGETRF(input.n, A, bb, input.scalapack_bf, env.mpi_rank, input.calc_procs, \
								env.blacs_nprow, env.blacs_npcol, env.blacs_row, env.blacs_col, \
								env.blacs_ctxt_grid, env.blacs_ctxt_root);

	if (env.mpi_rank==0)
	{
		// check exit condition
		if (output.exit_code!=0)
		{
			printf("\n** Dangerous exit code.. (%d)**\n",output.exit_code);
		}
		// calc error
		output.norm_rel_err = NormwiseRelativeError1D(bb, xx_ref, input.n, input.nrhs);

		if (verbosity>1)
		{
			printf("\nThe %s factorization is:\n",label);
			PrintMatrix1D(A, input.n, input.n);
			printf("\n with exit code     %d\n",output.exit_code);
			printf("      norm.rel.err. %f\n",output.norm_rel_err);
		}
		DeallocateMatrix1D(A);
		DeallocateMatrix1D(bb);
		DeallocateMatrix1D(xx_ref);
	}

	TEST_END(output, rank_result, team_result);
}
