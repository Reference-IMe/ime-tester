#include "../helpers/matrix.h"
#include "../1D/DGESV_WO_1D.h"

double test_IMe1D_DGESV(const char* label, int verbosity, int rows, int cols, int nrhs)
{
	clock_t start, stop;
	double* A;
	double* bb;
	double* xx;
	A=AllocateMatrix1D(rows, cols);
	bb=AllocateMatrix1D(rows, nrhs);
	xx=AllocateMatrix1D(rows, nrhs);

	FillMatrix1D(A, rows, cols);
	OneMatrix1D(bb, rows, nrhs);

	if (verbosity>2)
	{
		printf("\n\n Matrix A:\n");
		PrintMatrix1D(A, rows, cols);
		printf("\n Vector b:\n");
		PrintMatrix1D(bb, rows, nrhs);
	}

	start = clock();

	DGESV_WO_1D(rows, A, nrhs, bb, xx);

	stop = clock();

	if (verbosity>1)
	{
		printf("\nThe %s solution is:\n",label);
		PrintMatrix1D(xx, rows, nrhs);
	}

	DeallocateMatrix1D(A);
	DeallocateMatrix1D(bb);
	DeallocateVector(xx);

	return (double)(stop - start);
}
