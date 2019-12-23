#include "../helpers/matrix.h"
#include "LAPACK/Lapack_DGESV.h"

double test_Lapack_DGESV(const char* label, int verbosity, int rows, int cols, int nrhs)
{
	clock_t start, stop;
	double* A;
	double* bb;
	A=AllocateMatrix1D(rows, cols);
	bb=AllocateMatrix1D(rows, nrhs);

	FillMatrixT1D(A, rows, cols);	// input has to be transposed
	//TODO: fill normally and transpose inside the timed section

	OneMatrix1D(bb, rows, nrhs);

	if (verbosity>2)
	{
		printf("\n\n Matrix A:\n");
		PrintMatrix1D(A, rows, cols);
		printf("\n Vector b:\n");
		PrintMatrix1D(bb, rows, nrhs);
	}

	start=clock();

	Lapack_DGESV(rows, A, nrhs, bb);

	stop=clock();

	if (verbosity>1)
	{
		printf("\nThe %s solution is:\n",label);
		PrintMatrix1D(bb, nrhs, rows);	// result is transposed
		//TODO: transpose the result inside the timed section
	}

	DeallocateMatrix1D(A);
	DeallocateMatrix1D(bb);

	return (double)(stop - start);
}
