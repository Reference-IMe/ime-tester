#include <mpi.h>
#include <time.h>
#include <unistd.h>
#include "../../helpers/macros.h"
#include "../../helpers/matrix_advanced.h"
#include "lib/src/helpers/matrix_basic.h"
#include "../../helpers/simple_dynamic_strings/sds.h"
#include "../../tester_structures.h"
#include "lib/src/pdgesv-co.h"

test_result test_IMe_pDGESV_CO (	const char check,
									const char* tag,
									const char* variant,
									int verbosity,
									parallel_env env,
									test_input input,
									int fault_tolerance	)
{
	test_result rank_result = TEST_NOT_RUN;
	test_result team_result = TEST_NOT_RUN;
	exit_status output      = EMPTY_OUTPUT;

	sds label=sdsempty();
	TAG2LABEL(tag,label);

	int i,j;

	double** A2;
	double*  A2_1D;
	double** bb;
	double** xx;
	double*  xx_ref;

	int sqrt_calc_procs;

	MPI_Comm comm_calc;

	int i_calc; // participating in ime calc = 1, checksumming = 0

	if (check)
	{
		if (env.mpi_rank==0)
		{
			if (input.ime_bf < 1)
			{
				DISPLAY_ERR(label,"the blocking factor has to be greater than 0");
			}
			else
			{
				if (IS_SQUARE(env.calc_procs))
				{
					sqrt_calc_procs=sqrt(env.calc_procs);

					if (IS_MULT(input.n, sqrt_calc_procs))
					{
						if (input.n / sqrt_calc_procs > 0)
						{
							if (env.spare_procs > 0 || fault_tolerance > 0)
							{
								if (verbosity>0) DISPLAY_WRN(label,"can run also with fault tolerance set or spare processes allocated, but calc. processes will differ from total processes and no faults will be injected");
							}
							if (IS_MULT(input.n / sqrt_calc_procs, input.ime_bf))
							{
								if (verbosity>0) DISPLAY_MSG(label,"OK");
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
		if (env.mpi_rank >= env.calc_procs)
		{
			i_calc=0;
			MPI_Comm_split(MPI_COMM_WORLD, MPI_UNDEFINED, MPI_UNDEFINED, &comm_calc); // checksumming procs don't belong to calc communicator
		}
		else
		{
			i_calc=1;
			MPI_Comm_split(MPI_COMM_WORLD, i_calc, env.mpi_rank, &comm_calc); // calc procs belong to calc communicator
		}

		if (i_calc)
		{
			xx=AllocateMatrix2D_double(input.n, input.nrhs, CONTIGUOUS);
			bb=AllocateMatrix2D_double(input.n, input.nrhs, CONTIGUOUS);

			if (env.mpi_rank==0)
			{
				xx_ref=AllocateMatrix1D_double(input.n, input.nrhs);
				A2=AllocateMatrix2D_double(input.n, input.n, CONTIGUOUS);
				A2_1D=&A2[0][0];
				CopyMatrix1D_double(input.A_ref_d, A2_1D, input.n, input.n);

				for (i=0;i<input.n;i++)
				{
					for (j=0;j<input.nrhs;j++)
					{
						bb[i][j] = input.b_ref_d[i];
						xx_ref[i*input.nrhs+j] = input.x_ref_d[i];
					}
				}

				if (verbosity>2)
				{
					printf("\n\n Matrix A:\n");
					PrintMatrix2D_double(A2, input.n, input.n);
					printf("\n Vector b:\n");
					PrintMatrix2D_double(bb, input.n, input.nrhs);
				}
			}
			else
			{
				A2=AllocateMatrix2D_double(1, 1, CONTIGUOUS); // to avoid segmentation fault in mpi collectives with 2D arrays
				xx_ref=NULL;
			}

			if ( strcmp( variant, "PB-CO-BF1") == 0) output = pdgesv_co (input.ime_bf, input.n, A2, input.nrhs, bb, xx, comm_calc);
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
				if (input.calc_nre) rank_result.norm_rel_err = NormwiseRelativeError1D_double(&xx[0][0], xx_ref, input.n, input.nrhs);

				if (verbosity>1)
				{
					printf("\nThe %s solution is:\n",tag);
					PrintMatrix2D_double(xx, input.n, input.nrhs);
					printf("\n with exit code     %d\n",output.exit_code);
					printf("      norm.rel.err. %.17f\n",rank_result.norm_rel_err);
				}
			}
			//cleanup
			if (env.mpi_rank==0)
			{
				DeallocateMatrix2D_double(A2, input.n, CONTIGUOUS);
			}
			else
			{
				DeallocateMatrix2D_double(A2, 1, CONTIGUOUS);
			}
			DeallocateMatrix1D_double(xx_ref);
			DeallocateMatrix2D_double(xx, input.n, CONTIGUOUS);
			DeallocateMatrix2D_double(bb, input.n, CONTIGUOUS);

		}
		else
		{
			xx=NULL;
			bb=NULL;

			rank_result.total_time=0;
			rank_result.core_time=0;
		}

		if (env.spare_procs>0)
		{
			if (comm_calc != MPI_COMM_NULL)
			{
				MPI_Comm_free(&comm_calc);
			}
		}
	}
	sdsfree(label);
	TEST_END(output, rank_result, team_result);
}
