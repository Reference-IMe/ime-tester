#include <mpi.h>
#include <time.h>
#include <unistd.h>
#include "tester_structures.h"
#include "../helpers/macros.h"
#include "../helpers/matrix.h"
#include "../helpers/matrix_advanced.h"

#include "../pbDGESV_CO.bf1.ftx.h"
#include "../pbDGESV_CO.dev.h"

test_result test_IMe_pbDGESV_ftx(const char check, const char* label, const char* variant, int verbosity, parallel_env env, test_input input, int fault_protection, int fault_number, int failing_rank, int failing_level, int recovery)
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

	int sqrt_calc_procs;
	sqrt_calc_procs=sqrt(input.calc_procs);

	int* failing_rank_list;

	if (check)
	{
		if (env.mpi_rank==0)
		{
			if (fault_protection*sqrt_calc_procs < input.spare_procs)
			{
				DISPLAY_ERR(label,"not enough spare processes for the requested fault tolerance level");
			}
			if (input.ime_bf < 1)
			{
				DISPLAY_ERR(label,"the blocking factor has to be greater than 0");
			}
			else
			{
				if (IS_SQUARE(input.calc_procs))
				{
					if (IS_MULT(input.n, sqrt_calc_procs))
					{
						if (input.n / sqrt_calc_procs > 0)
						{
							if (sqrt_calc_procs - (failing_rank % sqrt_calc_procs) < fault_protection)
							{
								DISPLAY_WRN(label,"has first faulty rank too high: lowering..");
								failing_rank=(floor(failing_rank/sqrt_calc_procs)+1)*sqrt_calc_procs-fault_protection;
							}
							if (input.spare_procs > 0)
							{
								printf("     Faulty ranks:");
								for (i=failing_rank; i<failing_rank+fault_protection; i++)
								{
									printf(" %d",i);
								}
								printf("\n");
							}
							if (IS_MULT(input.n / sqrt_calc_procs, input.ime_bf))
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

		// failure list
		if (sqrt_calc_procs - (failing_rank % sqrt_calc_procs) < fault_protection)
		{
			failing_rank=(floor(failing_rank/sqrt_calc_procs)+1)*sqrt_calc_procs-fault_protection;
		}
		failing_rank_list = malloc(fault_protection*sizeof(int));

		for (i=failing_rank; i<failing_rank+fault_protection; i++)
		{
			failing_rank_list[i-failing_rank]=i;
		}

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

			 if ( strcmp( variant, "dev"            ) == 0) output = pbDGESV_CO_dev (A2, bb, xx, input, env, fault_protection, fault_number, failing_rank_list, failing_level, recovery);
		else if ( strcmp( variant, "PB-CO-bf1-ftx/0") == 0) output = pbDGESV_CO_bf1_ftx (A2, bb, xx, input, env);
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
