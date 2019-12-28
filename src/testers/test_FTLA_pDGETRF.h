#include "../helpers/matrix.h"
#include "FTLA.mod/FTLA_pDGETRF.h"

double test_FTLA_pDGETRF(const char* label, int verbosity, int rows, int cols, int rank, int cprocs, int sprocs)
{
	clock_t start, stop;
	double* A;

	if (rank==0)
	{
		A=AllocateMatrix1D(rows, cols);
		FillMatrixT1D(A, rows, cols);

		if (verbosity>2)
		{
			printf("\n\n Matrix A:\n");
			PrintMatrix1D(A, rows, cols);
		}
	}

	start=clock();

	calc_FTLA_ftdtr(rows, A, rank, cprocs, sprocs);

	stop=clock();

	if (rank==0)
	{
		if (verbosity>1)
		{
			printf("\nThe %s factorization is:\n",label);
			PrintMatrix1D(A, rows, cols);
		}
		DeallocateMatrix1D(A);
	}

	return (double)(stop - start);
}
