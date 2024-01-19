#include <time.h>

#include "../DGESV-WO.h"
#include "../helpers/matrix_basic.h"

double test_IMe_DGESV (	const char* tag,
						int verbosity,
						int rows,
						int cols,
						int nrhs		)
{
	clock_t start, stop;

	double** A;
	double** bb;
	double** xx;
	A=AllocateMatrix2D_double(rows, cols, CONTIGUOUS);
	bb=AllocateMatrix2D_double(rows, nrhs, CONTIGUOUS);
	xx=AllocateMatrix2D_double(rows, nrhs, CONTIGUOUS);

	FillMatrix2D_double(A, rows, cols);
	OneMatrix2D_double(bb, rows, nrhs);

	if (verbosity>2)
	{
		printf("\n\n Matrix A:\n");
		PrintMatrix2D_double(A, rows, cols);
		printf("\n Vector b:\n");
		PrintMatrix2D_double(bb, rows, nrhs);
	}

	start = time(NULL);

	DGESV_WO(rows, A, nrhs, bb, xx);

	stop = time(NULL);

	if (verbosity>1)
	{
		printf("\nThe %s solution is:\n",tag);
		PrintMatrix2D_double(xx, rows, nrhs);
	}

	DeallocateMatrix2D_double(A,rows,CONTIGUOUS);
	DeallocateMatrix2D_double(bb,rows,CONTIGUOUS);
	DeallocateMatrix2D_double(xx,rows,CONTIGUOUS);

	return (double)(stop - start);
}
