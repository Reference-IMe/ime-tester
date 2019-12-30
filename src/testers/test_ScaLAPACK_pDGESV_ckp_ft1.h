#include "../helpers/matrix.h"
#include "ScaLAPACK.mod/ScaLAPACK_pDGESV_ckp_ft1_sim.h"

double test_Scalapack_pDGESV_ckp_ft1_sim(const char* label, int verbosity, int rows, int cols, int nrhs, int rank, int cprocs, int sprocs, int failing_rank, int failing_level)
{
	clock_t start, stop;
	double* A;
	double* bb;

	if (rank==0)
	{
		A=AllocateMatrix1D(rows, cols);
		bb=AllocateMatrix1D(rows, nrhs);
		FillMatrixT1D(A, rows, cols);
		OneMatrix1D(bb, rows, nrhs);

		if (verbosity>2)
		{
			printf("\n\n Matrix A:\n");
			PrintMatrix1D(A, rows, cols);
			printf("\n Vector b:\n");
			PrintMatrix1D(bb, rows, nrhs);
		}
	}

	MPI_Barrier(MPI_COMM_WORLD);
	start=clock();

	Scalapack_pDGESV_ckp_ft1_sim(rows, A, nrhs, bb, rank, cprocs, sprocs, failing_rank, failing_level);

	MPI_Barrier(MPI_COMM_WORLD);
	stop=clock();

	if (rank==0)
	{
		if (verbosity>1)
		{
			printf("\nThe %s solution is:\n",label);
			//PrintMatrix1D(bb, nrhs, rows);
			PrintMatrix1D(bb, nrhs, rows);
		}
		DeallocateMatrix1D(A);
		DeallocateMatrix1D(bb);
	}

	return (double)(stop - start);
}
