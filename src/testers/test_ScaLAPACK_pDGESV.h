/*
 * test_ScaLAPACK_pDGESV.h
 *
 *  Created on: Dec 27, 2019
 *      Author: marcello
 */

#include <mpi.h>
#include <time.h>
#include "../helpers/info.h"
#include "../helpers/macros.h"
#include "../helpers/matrix.h"
#include "../helpers/matrix_advanced.h"
#include "ScaLAPACK/ScaLAPACK_pDGESV.h"

test_result test_ScaLAPACK_pDGESV(const char* label, int verbosity, test_input input, int rank)
{
	test_result rank_result = {0, 0, 0, 0};
	test_result team_result = {0, 0, 0, 0};
	test_output output = {0, 0, 0, 0, 0, 0};

	int i,j;

	double* A;
	double* bb;
	double* xx_ref;

	if (rank==0)
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

	output = ScaLAPACK_pDGESV_calc(input.n, A, input.nrhs, bb, input.scalapack_bf, rank, input.calc_procs);

	if (rank==0)
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
			printf("\nThe %s solution is:\n",label);
			PrintMatrix1D(bb, input.n, input.nrhs);
			printf("\n with exit code     %d\n",output.exit_code);
			printf("      norm.rel.err. %f\n",output.norm_rel_err);
		}
		DeallocateMatrix1D(A);
		DeallocateMatrix1D(bb);
		DeallocateMatrix1D(xx_ref);
	}

	TEST_END(output, rank_result, team_result);
}
