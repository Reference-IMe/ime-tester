#include "../helpers/matrix.h"
#include "../ft-simulated/pviDGESV_WO_ft1.h"

double test_IMe_pviDGESV_ft1(const char* label, int verbosity, int rows, int cols, int nrhs, int rank, int cprocs, int sprocs)
{
	clock_t start, stop;
	double** A2;
	double** bb;
	double** xx;

	xx=AllocateMatrix2D(rows, nrhs, CONTIGUOUS);
	bb=AllocateMatrix2D(rows, nrhs, CONTIGUOUS);

	if (rank==0)
	{
		A2=AllocateMatrix2D(rows,cols,CONTIGUOUS);

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

	start = clock();

	pviDGESV_WO_ft1(rows, A2, nrhs, bb, xx, rank, cprocs, sprocs);


	stop = clock();

	if (rank==0 && verbosity>1)
	{
		printf("\nThe %s solution is:\n",label);
		PrintMatrix2D(xx, rows, nrhs);
	}

	DeallocateMatrix2D(xx,rows,CONTIGUOUS);
	DeallocateMatrix2D(bb,rows,CONTIGUOUS);
	if (rank==0)
	{
		DeallocateMatrix2D(A2,rows,CONTIGUOUS);
	}
	else
	{
		DeallocateMatrix2D(A2, 0, CONTIGUOUS);
	}

	return (double)(stop - start);
}
