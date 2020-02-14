#include <time.h>
#include "../helpers/matrix.h"
#include "GaussianElimination/pGE.h"
#include "GaussianElimination/BackSubst.h"

double test_pGaussianElimination(const char* label, int verbosity, int rows, int cols, int nrhs, int rank, int nprocs)
{
	clock_t start, stop;
	double** A2;
	double** bb;
	double** xx;

	if (rank==0)
	{
		A2=AllocateMatrix2D(rows,cols,CONTIGUOUS);
		xx=AllocateMatrix2D(rows,nrhs,CONTIGUOUS);
		bb=AllocateMatrix2D(rows,nrhs,CONTIGUOUS);

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
	else // avoid problem (bug?) in mpi when calling collectives with undefined pointers in buffers even if unused!
	{
		A2=AllocateMatrix2D(0,0,CONTIGUOUS);
		bb=AllocateMatrix2D(0,0,CONTIGUOUS);
		xx=AllocateMatrix2D(0,0,CONTIGUOUS);
	}

	MPI_Barrier(MPI_COMM_WORLD);
	start=time(NULL);

	pGaussianElimination(rows, A2, nrhs, bb, rank, nprocs);

	if (rank==0)
	{
		BackSubstitution(rows, A2, nrhs, bb, xx);

		MPI_Barrier(MPI_COMM_WORLD);
		stop=time(NULL);

		if (verbosity>1)
		{
			printf("\nThe %s solution is:\n",label);
			PrintMatrix2D(xx, rows, nrhs);
		}
		DeallocateMatrix2D(bb,rows,CONTIGUOUS);
		DeallocateMatrix2D(xx,rows,CONTIGUOUS);
		DeallocateMatrix2D(A2,rows,CONTIGUOUS);
	}
	else
	{
		MPI_Barrier(MPI_COMM_WORLD);
		stop=time(NULL);

		DeallocateMatrix2D(bb,0,CONTIGUOUS);
		DeallocateMatrix2D(xx,0,CONTIGUOUS);
		DeallocateMatrix2D(A2,0,CONTIGUOUS); // due to mpi problem
	}

	MPI_Barrier(MPI_COMM_WORLD);
	return (double)(stop - start);
}
