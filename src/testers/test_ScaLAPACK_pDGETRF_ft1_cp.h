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
		if (env.mpi_rank==0)
		{
			if (faulty_procs > fault_tolerance)
			{
				DISPLAY_ERR(label,"requested fault occurrences greater than fault tolerance (single fault routine)");
			}
			else if (fault_tolerance < 1)
			{
				DISPLAY_ERR(label,"fault tolerance disabled ('-ft 0')");
			}
			else if (fault_tolerance > 1)
			{
				DISPLAY_ERR(label,"requested fault tolerance too high (single fault routine)");
			}
			else
			{
				if (input.scalapack_bf < 64)
				{
					DISPLAY_WRN(label,"blocking factor < 64");
				}
				if (env.spare_procs > 1)
				{
					DISPLAY_ERR(label,"too many spare processes allocated (single fault routine)");
				}
				else if (env.spare_procs > 0)
				{
					if (IS_SQUARE(env.calc_procs))
					{
						if (IS_MULT(input.n, rank_calc_procs))
						{
							if (IS_MULT(input.n / rank_calc_procs, input.scalapack_bf))
							{
								DISPLAY_MSG(label,"OK");
								output.exit_code = 0;
							}
							else
							{
								DISPLAY_ERR(label,"the number of columns (rows) per calc. process has to be a multiple of the blocking factor");
							}
						}
						else
						{
							DISPLAY_ERR(label,"the matrix size has to be a multiple of the calc. processes per rank");
						}
					}
					else
					{
						DISPLAY_ERR(label,"the number of the calc. processes must be square");
					}
				}
				else
				{
					DISPLAY_ERR(label,"fault tolerance enabled, but no requested spare processes");
				}
			}
		}
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
			// check exit condition
			if (output.exit_code!=0)
			{
				printf("\n** Dangerous exit code.. (%d)**\n",output.exit_code);
			}
			// calc error
			if (input.calc_nre) rank_result.norm_rel_err = NormwiseRelativeError1D(bb, input.x_ref, input.n, input.nrhs);

			if (verbosity>1)
			{
				printf("\nThe %s factorization is:\n",label);
				PrintMatrix1D(A, input.n, input.n);
				printf("\n with exit code     %d\n",output.exit_code);
				printf("      norm.rel.err. %.17f\n",rank_result.norm_rel_err);
			}
		}

		NULLFREE(A);
		NULLFREE(bb);
	}
	TEST_END(output, rank_result, team_result);
}
