#include <time.h>
#include "../helpers/matrix.h"
#include "GaussianElimination/GE.h"
#include "GaussianElimination/BackSubst.h"

double test_GaussianElimination(const char* label, int verbosity, int rows, int cols, int nrhs)
{
	clock_t start, stop;
	double** A;
	double** bb;
	double** xx;

	A=AllocateMatrix2D(rows,cols,CONTIGUOUS);
	bb=AllocateMatrix2D(rows,nrhs,CONTIGUOUS);
	xx=AllocateMatrix2D(rows,nrhs,CONTIGUOUS);

	FillMatrix2D(A, rows, cols);
	OneMatrix2D(bb, rows, nrhs);

	if (verbosity>2)
	{
		printf("\n\n Matrix A:\n");
		PrintMatrix2D(A, rows, cols);
		printf("\n Vector b:\n");
		PrintMatrix2D(bb, rows, nrhs);
	}

	start =time(NULL);

	GaussianElimination(rows, A, nrhs, bb);
	BackSubstitution(rows, A, nrhs, bb, xx);

	stop = time(NULL);

	if (verbosity>1)
	{
		printf("\nThe %s solution is:\n",label);
		PrintMatrix2D(xx, rows, nrhs);
	}

	DeallocateMatrix2D(A,rows,CONTIGUOUS);
	DeallocateMatrix2D(bb,rows,CONTIGUOUS);
	DeallocateMatrix2D(xx,rows,CONTIGUOUS);

	return (double)(stop - start);
}
