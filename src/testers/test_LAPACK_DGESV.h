#include <time.h>

#include "../helpers/matrix_basic.h"
#include "LAPACK/LAPACK_DGESV.h"

double test_Lapack_DGESV(const char* label, int verbosity, int rows, int cols, int nrhs)
{
	clock_t start, stop;
	double* A;
	double* bb;
	A=AllocateMatrix1D_double(rows, cols);
	bb=AllocateMatrix1D_double(rows, nrhs);

	FillMatrixT1D_double(A, rows, cols);	// input has to be transposed
	//TODO: fill normally and transpose inside the timed section

	OneMatrix1D_double(bb, rows, nrhs);

	if (verbosity>2)
	{
		printf("\n\n Matrix A:\n");
		PrintMatrix1D_double(A, rows, cols);
		printf("\n Vector b:\n");
		PrintMatrix1D_double(bb, rows, nrhs);
	}

	start=time(NULL);

	Lapack_DGESV(rows, A, nrhs, bb);

	stop=time(NULL);

	if (verbosity>1)
	{
		printf("\nThe %s solution is:\n",label);
		PrintMatrix1D_double(bb, nrhs, rows);	// result is transposed
		//TODO: transpose the result inside the timed section
	}

	DeallocateMatrix1D_double(A);
	DeallocateMatrix1D_double(bb);

	return (double)(stop - start);
}
