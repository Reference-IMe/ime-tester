#include "../helpers/matrix.h"
#include "ScaLAPACK/Scalapack_pDGEQRF.h"

double test_Scalapack_pDGEQRF(const char* label, int verbosity, int rows, int cols, int rank, int cprocs, int sprocs)
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

	Scalapack_pDGEQRF_calc(rows, A, rank, cprocs, sprocs);

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
