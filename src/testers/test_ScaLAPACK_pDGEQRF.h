/*
 * test_ScaLAPACK_pDGEQRF.h
 *
 *  Created on: Dec 5, 2019
 *      Author: marcello
 */

#include <mpi.h>
#include <time.h>
#include "../helpers/matrix.h"
#include "ScaLAPACK/ScaLAPACK_pDGEQRF.h"

double test_ScaLAPACK_pDGEQRF(const char* label, int verbosity, int rows, int cols, int nb, int rank, int cprocs)
{
	clock_t start, stop;
	double span, maxspan;
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

	//MPI_Barrier(MPI_COMM_WORLD);
	start=clock();

	ScaLAPACK_pDGEQRF_calc(rows, A, nb, rank, cprocs);

	//MPI_Barrier(MPI_COMM_WORLD);
	stop=clock();

	span=(double)(stop - start);
    MPI_Reduce( &span, &maxspan, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD );

	if (rank==0)
	{
		if (verbosity>1)
		{
			printf("\nThe %s factorization is:\n",label);
			PrintMatrix1D(A, rows, cols);
		}
		DeallocateMatrix1D(A);
	}

	MPI_Barrier(MPI_COMM_WORLD);
	return maxspan;
}
