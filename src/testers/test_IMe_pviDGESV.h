#include <mpi.h>
#include <time.h>
#include "../helpers/info.h"
#include "../helpers/macros.h"
#include "../helpers/matrix.h"
#include "../helpers/matrix_advanced.h"
#include "../pviDGESV_WO.h"


run_info test_IMe_pviDGESV(const char* label, int verbosity, int n, double* A_ref, double* x_ref, double* b_ref, int nrhs, int rank, int cprocs, int sprocs)
{
	run_info process_info, team_info;
	result_info info;

	int i,j;

	double** A2;
	double* A2_1D;
	double** bb;
	double** xx;
	double* xx_ref;

	MPI_Comm comm_calc;

	int i_calc; // participating in ime calc = 1, checksumming = 0

	if (rank>=cprocs)
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
		xx=AllocateMatrix2D(n, nrhs, CONTIGUOUS);
		bb=AllocateMatrix2D(n, nrhs, CONTIGUOUS);

		if (rank==0)
		{
			xx_ref=AllocateMatrix1D(n, nrhs);
			A2=AllocateMatrix2D(n, n, CONTIGUOUS);
			A2_1D=&A2[0][0];
			CopyMatrix1D(A_ref, A2_1D, n, n);

			for (i=0;i<n;i++)
			{
				for (j=0;j<nrhs;j++)
				{
					bb[i][j] = b_ref[i];
					xx_ref[i*nrhs+j] = x_ref[i];
				}
			}

			if (verbosity>2)
			{
				printf("\n\n Matrix A:\n");
				PrintMatrix2D(A2, n, n);
				printf("\n Vector b:\n");
				PrintMatrix2D(bb, n, nrhs);
			}
		}
		else
		{
			A2=AllocateMatrix2D(0, 0, CONTIGUOUS);
		}

		info = pviDGESV_WO(n, A2, nrhs, bb, xx, comm_calc);

		if (rank==0)
		{
			// check exit condition
			if (info.exit_code!=0)
			{
				printf("\n** Dangerous exit code.. (%d)**\n",info.exit_code);
			}
			// calc error
			info.norm_rel_err = NormwiseRelativeError1D(&xx[0][0], xx_ref, n, nrhs);

			if (verbosity>1)
			{
				printf("\nThe %s solution is:\n",label);
				PrintMatrix2D(xx, n, nrhs);
				printf("\n with exit code     %d\n",info.exit_code);
				printf("      norm.rel.err. %f\n",info.norm_rel_err);
			}
		}

		DeallocateMatrix2D(xx, n, CONTIGUOUS);
		DeallocateMatrix2D(bb, n, CONTIGUOUS);
		if (rank==0)
		{
			DeallocateMatrix2D(A2, n, CONTIGUOUS);
		}
		else
		{
			DeallocateMatrix2D(A2, 0, CONTIGUOUS);
		}
	}
	else
	{
		process_info.total_time=0;
		process_info.core_time=0;
	}

	if (sprocs>0)
	{
		if (comm_calc != MPI_COMM_NULL)
		{
			MPI_Comm_free(&comm_calc);
		}
	}

	TEST_END(info, process_info, team_info);
}
