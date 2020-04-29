#include <mpi.h>
#include <time.h>
#include "../helpers/macros.h"
#include "../helpers/matrix.h"
#include "../helpers/matrix_advanced.h"
#include "tester_structures.h"
#include "../pviDGESV_WO.h"



test_result test_IMe_pviDGESV(const char* label, int verbosity, test_input input, int rank)
{
	test_result rank_result = TEST_NOT_RUN;
	test_result team_result = TEST_NOT_RUN;
	test_output output = EMPTY_OUTPUT;

	int i,j;

	double** A2;
	double* A2_1D;
	double** bb;
	double** xx;
	double* xx_ref;

	MPI_Comm comm_calc;

	int i_calc; // participating in ime calc = 1, checksumming = 0

	if (rank >= input.calc_procs)
	{
		i_calc=0;
		MPI_Comm_split(MPI_COMM_WORLD, MPI_UNDEFINED, MPI_UNDEFINED, &comm_calc); // checksumming procs don't belong to calc communicator
	}
	else
	{
		i_calc=1;
		MPI_Comm_split(MPI_COMM_WORLD, i_calc, rank, &comm_calc); // calc procs belong to calc communicator
	}

	if (i_calc)
	{
		xx=AllocateMatrix2D(input.n, input.nrhs, CONTIGUOUS);
		bb=AllocateMatrix2D(input.n, input.nrhs, CONTIGUOUS);

		if (rank==0)
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
			A2=AllocateMatrix2D(0, 0, CONTIGUOUS);
		}

		output = pviDGESV_WO(input.ime_bf, input.n, A2, input.nrhs, bb, xx, comm_calc);

		if (rank==0)
		{
			// check exit condition
			if (output.exit_code!=0)
			{
				printf("\n** Dangerous exit code.. (%d)**\n",output.exit_code);
			}
			// calc error
			output.norm_rel_err = NormwiseRelativeError1D(&xx[0][0], xx_ref, input.n, input.nrhs);

			if (verbosity>1)
			{
				printf("\nThe %s solution is:\n",label);
				PrintMatrix2D(xx, input.n, input.nrhs);
				printf("\n with exit code     %d\n",output.exit_code);
				printf("      norm.rel.err. %f\n",output.norm_rel_err);
			}
		}

		DeallocateMatrix2D(xx, input.n, CONTIGUOUS);
		DeallocateMatrix2D(bb, input.n, CONTIGUOUS);
		if (rank==0)
		{
			DeallocateMatrix2D(A2, input.n, CONTIGUOUS);
		}
		else
		{
			DeallocateMatrix2D(A2, 0, CONTIGUOUS);
		}
	}
	else
	{
		rank_result.total_time=0;
		rank_result.core_time=0;
	}

	if (input.spare_procs>0)
	{
		if (comm_calc != MPI_COMM_NULL)
		{
			MPI_Comm_free(&comm_calc);
		}
	}

	TEST_END(output, rank_result, team_result);
}
