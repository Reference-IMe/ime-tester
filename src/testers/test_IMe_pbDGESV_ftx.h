#include <mpi.h>
#include <time.h>
#include <unistd.h>
#include "tester_structures.h"
#include "../helpers/macros.h"
#include "../helpers/matrix.h"
#include "../helpers/matrix_advanced.h"

#include "../pbDGESV_CO.bf1.ftx.h"


test_result test_IMe_pbDGESV_ftx(const char check, const char* label, const char* variant, int verbosity, parallel_env env, test_input input)
{
	test_result rank_result = TEST_NOT_RUN;
	test_result team_result = TEST_NOT_RUN;
	test_output output      = EMPTY_OUTPUT;

	int i,j;

	double** A2;
	double*  A2_1D;
	double** bb;
	double** xx;
	double*  xx_ref;

	if (check)
	{
		//TODO: ckecking rules
		/*
		if (rank==0)
		{
			if (input.ime_bf < 1)
			{
				DISPLAY_ERR(label,"the blocking factor has to be greater than 0");
			}
			else
			{
				if (IS_MULT(input.n, input.calc_procs))
				{
					if (input.n / input.calc_procs > 0)
					{
						if (input.spare_procs > 0)
						{
							DISPLAY_WRN(label,"can run also with FT enabled, but calc. processes differ from total processes")
						}
						if (IS_MULT(input.n / input.calc_procs, input.ime_bf))
						{
							DISPLAY_MSG(label,"OK");
							output.exit_code = 0;
						}
						else
						{
							DISPLAY_ERR(label,"the number of columns per calc. process has to be a multiple of the blocking factor");
						}
					}
					else DISPLAY_ERR(label,"the number of columns per calc. process has to be greater than 0");
				}
				else DISPLAY_ERR(label,"the number of columns has to be a multiple of the calc. processes");
			}
		}
		*/ output.exit_code = 0;
	}
	else
	{
		if (env.mpi_rank < input.calc_procs)
		{
			xx=AllocateMatrix2D(input.n, input.nrhs, CONTIGUOUS);
			bb=AllocateMatrix2D(input.n, input.nrhs, CONTIGUOUS);

			if (env.mpi_rank==0)
			{
				xx_ref=AllocateMatrix1D(input.n, input.nrhs);
				A2=AllocateMatrix2D(input.n, input.n, CONTIGUOUS);
				A2_1D=&A2[0][0];
				CopyMatrix1D(input.A_ref, A2_1D, input.n, input.n);

				for (i=0;i<input.n;i++)
				{
					for (j=0;j<input.nrhs;j++)
					{
						bb[i][j] = input.b_ref[i];
						xx_ref[i*input.nrhs+j] = input.x_ref[i];
					}
				}

				if (verbosity>2)
				{
					printf("\n\n Matrix A:\n");
					PrintMatrix2D(A2, input.n, input.n);
					printf("\n Vector b:\n");
					PrintMatrix2D(bb, input.n, input.nrhs);
				}
			}
			else
			{
				A2=AllocateMatrix2D(1, 1, CONTIGUOUS); // to avoid segmentation fault in mpi collectives with 2D arrays
				xx_ref=NULL;
			}
		}
		else
		{
			A2=AllocateMatrix2D(1, 1, CONTIGUOUS); // to avoid segmentation fault in mpi collectives with 2D arrays
			xx_ref=NULL;
			xx=NULL;
			bb=NULL;
		}


		if ( strcmp( variant, "PB-CO-bf1-ftx/0") == 0) output = pbDGESV_CO_bf1_ftx (A2, bb, xx, input, env);
		else
		{
			DISPLAY_ERR(label,"not yet implemented! UNDEFINED BEHAVIOUR!");
		}

		if (env.mpi_rank==0)
		{
			// check exit condition
			if (output.exit_code!=0)
			{
				printf("\n** Dangerous exit code.. (%d)**\n",output.exit_code);
			}
			// calc error
			if (input.calc_nre) rank_result.norm_rel_err = NormwiseRelativeError1D(&xx[0][0], xx_ref, input.n, input.nrhs);

			if (verbosity>1)
			{
				printf("\nThe %s solution is:\n",label);
				PrintMatrix2D(xx, input.n, input.nrhs);
				printf("\n with exit code     %d\n",output.exit_code);
				printf("      norm.rel.err. %f\n",rank_result.norm_rel_err);
			}
		}

		//cleanup
		if (env.mpi_rank < input.calc_procs)
		{
			if (env.mpi_rank==0)
			{
				DeallocateMatrix2D(A2, input.n, CONTIGUOUS);
				DeallocateMatrix1D(xx_ref);
			}
			else
			{
				DeallocateMatrix2D(A2, 1, CONTIGUOUS);
			}
			DeallocateMatrix2D(xx, input.n, CONTIGUOUS);
			DeallocateMatrix2D(bb, input.n, CONTIGUOUS);
		}
	}
	TEST_END(output, rank_result, team_result);
}
