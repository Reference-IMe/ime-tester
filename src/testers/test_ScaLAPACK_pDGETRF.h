/*
 * test_ScaLAPACK_pDGETRF.h
 *
 *  Created on: Dec 27, 2019
 *      Author: marcello
 */

#include "../helpers/matrix.h"
#include "ScaLAPACK/ScaLAPACK_pDGETRF.h"

double test_ScaLAPACK_pDGETRF(const char* label, int verbosity, int rows, int cols, int rank, int cprocs, int sprocs)
{
	clock_t start, stop;
	double* A;
	//double* bb;

	if (rank==0)
	{
		A=AllocateMatrix1D(rows, cols);
		//bb=AllocateMatrix1D(rows, nrhs);
		FillMatrixT1D(A, rows, cols);
		//OneMatrix1D(bb, rows, nrhs);

		if (verbosity>2)
		{
			printf("\n\n Matrix A:\n");
			PrintMatrix1D(A, rows, cols);
			//printf("\n Vector b:\n");
			//PrintMatrix1D(bb, rows, nrhs);
		}
	}

	MPI_Barrier(MPI_COMM_WORLD);
	start=clock();

	// Scalapack_pDGETRF(rows, A, nrhs, bb, rank, cprocs, sprocs); // for consistency checking
	ScaLAPACK_pDGETRF_calc(rows, A, rank, cprocs, sprocs);

	MPI_Barrier(MPI_COMM_WORLD);
	stop=clock();

	if (rank==0)
	{
		if (verbosity>1)
		{
			printf("\nThe %s factorization is:\n",label);
			PrintMatrix1D(A, rows, cols);
		}
		DeallocateMatrix1D(A);
		//DeallocateMatrix1D(bb);
	}

	return (double)(stop - start);
}
