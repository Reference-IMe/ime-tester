/*
 * test_ScaLAPACK_pDGESV.h
 *
 *  Created on: Dec 27, 2019
 *      Author: marcello
 */

#include <mpi.h>
#include <time.h>
#include "../helpers/matrix.h"
#include "ScaLAPACK/ScaLAPACK_pDGESV.h"

double test_ScaLAPACK_pDGESV(const char* label, int verbosity, int rows, int cols, int nrhs, int nb, int rank, int cprocs)
{
	clock_t start, stop;
	double span, maxspan;
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

	//MPI_Barrier(MPI_COMM_WORLD);
	start=clock();

	ScaLAPACK_pDGESV_calc(rows, A, nrhs, bb, nb, rank, cprocs);

	//MPI_Barrier(MPI_COMM_WORLD);
	stop=clock();

	span=(double)(stop - start);
    MPI_Reduce( &span, &maxspan, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD );

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

	MPI_Barrier(MPI_COMM_WORLD);
	return maxspan;
}
