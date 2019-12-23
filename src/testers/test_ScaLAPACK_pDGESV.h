#include "../helpers/matrix.h"
#include "ScaLAPACK/Scalapack_pDGESV.h"

double test_Scalapack_pDGESV(const char* label, int verbosity, int rows, int cols, int nrhs, int rank, int cprocs)
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

	start=clock();

	Scalapack_pDGESV_calc(rows, A, nrhs, bb, rank, cprocs);

	stop=clock();

	if (rank==0)
	{
		if (verbosity>1)
		{
			printf("\nThe %s solution is:\n",label);
			PrintMatrix1D(bb, nrhs, rows);
		}
		DeallocateMatrix1D(A);
		DeallocateMatrix1D(bb);
	}

	return (double)(stop - start);
}
