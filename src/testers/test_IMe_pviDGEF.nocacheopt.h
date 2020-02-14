#include <mpi.h>
#include <time.h>
#include "../helpers/info.h"
#include "../helpers/macros.h"
#include "../helpers/matrix.h"
#include "../pviDGEF_WO.nocacheopt.h"


duration_t test_IMe_pviDGEF_nocacheopt(const char* label, int verbosity, int rows, int cols, int rank, int cprocs, int sprocs)
{
	duration_t timing, timing_max;
	result_info info;

	int r,c;
	double** A2;
	double** K;

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
		if (rank==0)
		{
			A2=AllocateMatrix2D(rows, cols, CONTIGUOUS);
			FillMatrix2D(A2, rows, cols);

			K=AllocateMatrix2D(rows, cols, CONTIGUOUS);

			if (verbosity>2)
			{
				printf("\n\n Matrix A:\n");
				PrintMatrix2D(A2, rows, cols);
			}
		}
		else
		{
			A2=AllocateMatrix2D(0, 0, CONTIGUOUS);
			K=AllocateMatrix2D(0, 0, CONTIGUOUS);
		}

		info = pviDGEF_WO_nocacheopt(rows, A2, K, comm_calc);

		// clean K for output
		if (rank==0)
		{
			for (r=0;r<rows;r++)
			{
				for (c=r;c<cols;c++)
				{
					K[r][c]=0;
				}
			}
		}

		if (rank==0 && verbosity>1)
		{
			printf("\nThe %s X factor is:\n",label);
			PrintMatrix2D(A2, rows, cols);
			printf("\nThe %s K factor is:\n",label);
			PrintMatrix2D(K, rows, cols);
		}

		if (rank==0)
		{
			DeallocateMatrix2D(A2, rows, CONTIGUOUS);
			DeallocateMatrix2D(K, rows, CONTIGUOUS);
		}
		else
		{
			DeallocateMatrix2D(A2, 0, CONTIGUOUS);
			DeallocateMatrix2D(K, 0, CONTIGUOUS);
		}
	}
	else
	{
		timing.total=0;
		timing.core=0;
	}

	if (sprocs>0)
	{
		if (comm_calc != MPI_COMM_NULL)
		{
			MPI_Comm_free(&comm_calc);
		}
	}

	TEST_END(info, timing, timing_max);
}
