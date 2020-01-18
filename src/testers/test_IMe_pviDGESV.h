#include <mpi.h>
#include <time.h>
#include "../helpers/matrix.h"
#include "../pviDGESV_WO.h"


double test_IMe_pviDGESV(const char* label, int verbosity, int rows, int cols, int nrhs, int rank, int cprocs, int sprocs)
{
	clock_t start, stop;
	double span, maxspan;
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

		//MPI_Barrier(MPI_COMM_WORLD);
		start = clock();

		pviDGESV_WO(rows, A2, nrhs, bb, xx, comm_calc);

		//MPI_Barrier(MPI_COMM_WORLD);
		stop = clock();

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
		start=stop=0;
	}

	if (sprocs>0)
	{
		if (comm_calc != MPI_COMM_NULL)
		{
			MPI_Comm_free(&comm_calc);
		}
	}

	span=(double)(stop - start);
    MPI_Reduce( &span, &maxspan, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD );

    MPI_Barrier(MPI_COMM_WORLD);
	return maxspan;
}
