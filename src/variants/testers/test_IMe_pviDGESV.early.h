#include <mpi.h>
#include <time.h>
#include "../helpers/info.h"
#include "../helpers/macros.h"
#include "../helpers/matrix.h"
#include "../pviDGESV_WO.early.h"


duration_t test_IMe_pviDGESV_early(const char* label, int verbosity, int rows, int cols, int nrhs, int rank, int cprocs, int sprocs)
{
	duration_t timing, timing_max;
	test_output info;

	double** A2;
	double** bb;
	double** xx;

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
		xx=AllocateMatrix2D(rows, nrhs, CONTIGUOUS);
		bb=AllocateMatrix2D(rows, nrhs, CONTIGUOUS);

		if (rank==0)
		{
			A2=AllocateMatrix2D(rows, cols, CONTIGUOUS);
			FillMatrix2D(A2, rows, cols);

			OneMatrix2D(bb, rows, nrhs);

			if (verbosity>2)
			{
				printf("\n\n Matrix A:\n");
				PrintMatrix2D(A2, rows, cols);
				printf("\n Vector b:\n");
				PrintMatrix2D(bb, rows, nrhs);
			}
		}
		else
		{
			A2=AllocateMatrix2D(0, 0, CONTIGUOUS);
		}

		info = pviDGESV_WO_early(rows, A2, nrhs, bb, xx, comm_calc);

		if (rank==0 && verbosity>1)
		{
			printf("\nThe %s solution is:\n",label);
			PrintMatrix2D(xx, rows, nrhs);
		}

		DeallocateMatrix2D(xx, rows, CONTIGUOUS);
		DeallocateMatrix2D(bb, rows, CONTIGUOUS);
		if (rank==0)
		{
			DeallocateMatrix2D(A2, rows, CONTIGUOUS);
		}
		else
		{
			DeallocateMatrix2D(A2, 0, CONTIGUOUS);
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
