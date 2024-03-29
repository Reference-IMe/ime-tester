#include <mpi.h>
#include <time.h>
#include <unistd.h>
#include "lib/src/psgesv-co-ft.h"
#include "lib/src/helpers/matrix_basic.h"
#include "../../helpers/matrix_advanced.h"
#include "../../helpers/simple_dynamic_strings/sds.h"
#include "../../helpers/macros.h"
#include "../../tester_structures.h"


test_result test_IMe_pSGESV_CO_FT (	const char check,
									const char* tag,
									const char* variant,
									int verbosity,
									parallel_env env,
									test_input input,
									int fault_tolerance,
									int faulty_procs,
									int failing_rank,
									int failing_level,
									int recovery		)
{
	test_result rank_result = TEST_NOT_RUN;
	test_result team_result = TEST_NOT_RUN;
	exit_status output      = EMPTY_OUTPUT;

	sds label=sdsempty();
	TAG2LABEL(tag,label);

	int i,j;

	float** A2;
	float*  A2_1D;
	float** bb;
	float** xx;
	float*  xx_ref;

	int sqrt_calc_procs;
	sqrt_calc_procs=sqrt(env.calc_procs);

	int* failing_rank_list;

	if (check)
	{
		if (env.mpi_rank==0)
		{
			if (fault_tolerance < 1)
			{
				DISPLAY_WRN(label,"fault tolerance disabled: never recovering");
			}
			else if (faulty_procs > fault_tolerance)
			{
				DISPLAY_ERR(label,"fault tolerance set but requested fault occurrences greater than fault tolerance level");
			}

			if (faulty_procs == 0)
			{
				DISPLAY_WRN(label,"no fault will be actually injected: never failing");
			}

			if (fault_tolerance*sqrt_calc_procs > env.spare_procs)
			{
				DISPLAY_ERR(label,"not enough spare processes for the requested fault tolerance level");
			}
			else if (fault_tolerance*sqrt_calc_procs < env.spare_procs)
			{
				DISPLAY_ERR(label,"in this experimental version spare processes have to match the requested fault tolerance level");
			}
			else if (input.ime_bf < 1)
			{
				DISPLAY_ERR(label,"the blocking factor has to be greater than 0");
			}
			else
			{
				if (IS_SQUARE(env.calc_procs))
				{
					if (IS_MULT(input.n, sqrt_calc_procs))
					{
						if (input.n / sqrt_calc_procs > 0)
						{
							if (sqrt_calc_procs - (failing_rank % sqrt_calc_procs) < faulty_procs)
							{
								DISPLAY_WRN(label,"has first faulty rank too high: lowering..");
								failing_rank=(floor(failing_rank/sqrt_calc_procs)+1)*sqrt_calc_procs - faulty_procs;
							}
							if (failing_level == 0)
							{
								DISPLAY_WRN(label,"has failing level at last one: not allowed, correcting to last-but-one..");
								failing_level = 1;
							}
							else if (failing_level >= (input.n - 1) )
							{
								DISPLAY_WRN(label,"has failing level at first one: not allowed, correcting to second one..");
								failing_level = input.n - 2;
							}
							if (faulty_procs > 0)
							{
								printf("     %-30s","Faulty ranks:");
								for (i=failing_rank; i<failing_rank+faulty_procs; i++)
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

		if (env.mpi_rank < env.calc_procs)
		{
			xx=AllocateMatrix2D_float(input.n, input.nrhs, CONTIGUOUS);
			bb=AllocateMatrix2D_float(input.n, input.nrhs, CONTIGUOUS);

			if (env.mpi_rank==0)
			{
				xx_ref=AllocateMatrix1D_float(input.n, input.nrhs);
				A2=AllocateMatrix2D_float(input.n, input.n, CONTIGUOUS);
				A2_1D=&A2[0][0];
				CopyMatrix1D_float(input.A_ref_s, A2_1D, input.n, input.n);

				for (i=0;i<input.n;i++)
				{
					for (j=0;j<input.nrhs;j++)
					{
						bb[i][j] = input.b_ref_s[i];
						xx_ref[i*input.nrhs+j] = input.x_ref_s[i];
					}
				}

				if (verbosity>2)
				{
					printf("\n\n Matrix A:\n");
					PrintMatrix2D_float(A2, input.n, input.n);
					printf("\n Vector b:\n");
					PrintMatrix2D_float(bb, input.n, input.nrhs);
				}
			}
			else
			{
				A2=AllocateMatrix2D_float(1, 1, CONTIGUOUS); // to avoid segmentation fault in mpi collectives with 2D arrays
				xx_ref=NULL;
			}
		}
		else
		{
			A2=AllocateMatrix2D_float(1, 1, CONTIGUOUS); // to avoid segmentation fault in mpi collectives with 2D arrays
			xx_ref=NULL;
			xx=NULL;
			bb=NULL;
		}

		if ( strcmp( variant, "PB-CO-BF1-FT") == 0)
		{
			if (fault_tolerance < 1)
			{
				// if fault tolerance is disabled, force disable recovery
				output = psgesv_co_ft (A2, bb, xx, input, env, faulty_procs, failing_rank_list, failing_level, 0);
			}
			else
			{
				// if fault tolerance is enabled, pass recovery option as set by user
				output = psgesv_co_ft (A2, bb, xx, input, env, faulty_procs, failing_rank_list, failing_level, recovery);
			}
		}
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
			if (input.calc_nre) rank_result.norm_rel_err = NormwiseRelativeError1D_float(&xx[0][0], xx_ref, input.n, input.nrhs);

			if (verbosity>1)
			{
				printf("\nThe %s solution is:\n",tag);
				PrintMatrix2D_float(xx, input.n, input.nrhs);
				printf("\n with exit code     %d\n",output.exit_code);
				printf("      norm.rel.err. %.17f\n",rank_result.norm_rel_err);
			}
		}

		//cleanup
		if (env.mpi_rank < env.calc_procs)
		{
			if (env.mpi_rank==0)
			{
				DeallocateMatrix2D_float(A2, input.n, CONTIGUOUS);
				DeallocateMatrix1D_float(xx_ref);
			}
			else
			{
				DeallocateMatrix2D_float(A2, 1, CONTIGUOUS);
			}
			DeallocateMatrix2D_float(xx, input.n, CONTIGUOUS);
			DeallocateMatrix2D_float(bb, input.n, CONTIGUOUS);
		}
		NULLFREE(failing_rank_list);
	}
	sdsfree(label);
	TEST_END(output, rank_result, team_result);
}
