/*
 * test_FTLA_pDGESV_TRF.h
 *
 *  Created on: Aug 9, 2021
 *      Author: marcello
 */

#include <mpi.h>
#include <time.h>
#include "../../helpers/macros.h"
#include "../../tester_structures.h"
#include "../IMe/lib/src/helpers/matrix_basic.h"
#include "FTLA_pDGESV_TRF.h"

test_result test_FTLA_pDGESV_TRF (	const char check,
									const char* tag,
									int verbosity,
									parallel_env env,
									test_input input,
									int faults			)
{
	test_result rank_result = TEST_NOT_RUN;
	test_result team_result = TEST_NOT_RUN;
	exit_status output      = EMPTY_OUTPUT;

	sds label=sdsempty();
	TAG2LABEL(tag,label);

	double* A;
	double* bb;
	int i;

	int sqrt_calc_procs;

	if (check)
	{
		if (env.mpi_rank==0)
		{
			if (input.scalapack_bf < 1)
			{
				DISPLAY_ERR(label,"the blocking factor has to be greater than 0");
			}
			else
			{
				if (input.scalapack_bf < 64)
				{
					DISPLAY_WRN(label,"blocking factor < 64")
				}
				if (IS_SQUARE(env.calc_procs))
				{
					sqrt_calc_procs=sqrt(env.calc_procs);

					if (IS_MULT(input.n, sqrt_calc_procs))
					{
						if (input.n / sqrt_calc_procs > 0)
						{
							if (env.spare_procs > 0)
							{
								DISPLAY_WRN(label,"can run also with FT enabled or spare processes allocated, but calc. processes will differ from total processes")
							}
							if (IS_MULT(input.n / sqrt_calc_procs, input.scalapack_bf))
							{
								DISPLAY_MSG(label,"OK");
								output.exit_code = 0;
							}
							else
							{
								DISPLAY_ERR(label,"the number of columns (rows) per calc. process has to be a multiple of the blocking factor");
							}
						}
						else DISPLAY_ERR(label,"the number of columns (rows) per calc. process has to be greater than 0");
					}
					else DISPLAY_ERR(label,"the number of columns (rows) has to be a multiple of the square root of calc. processes");
				}
				else DISPLAY_ERR(label,"the number of calc. process has to be square");
			}
		}
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
				printf("\n Vector b:\n");
				PrintMatrix1D_double(bb, input.n, input.nrhs);
			}
		}
		else
		{
			A  = NULL;
			bb = NULL;
		}

		output = FTLA_ftdtr_sv(input.n, A, bb, input.scalapack_bf, env.mpi_rank, env.calc_procs,
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
			if (input.calc_nre) rank_result.norm_rel_err = NormwiseRelativeError1D_double(bb, input.x_ref_d, input.n, 1);

			if (verbosity>1)
			{
				printf("\nThe %s solution is:\n",tag);
				PrintMatrix1D_double(bb, input.n, input.nrhs);
				printf("\n with exit code     %d\n",output.exit_code);
				printf("      norm.rel.err. %.17f\n",rank_result.norm_rel_err);
			}
		}

		NULLFREE(A);
		NULLFREE(bb);
	}
	sdsfree(label);
	TEST_END(output, rank_result, team_result);
}
