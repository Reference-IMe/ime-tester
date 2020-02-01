/*
 * test_ScaLAPACK_pDGEQRF_cp_ft1.h
 *
 *  Created on: Jan 8, 2020
 *      Author: marcello
 */

#include <mpi.h>
#include <time.h>
#include "../helpers/matrix.h"
#include "ScaLAPACK/ScaLAPACK_pDGEQRF_cp_ft1_sim.h"

double test_ScaLAPACK_pDGEQRF_cp_ft1_sim(const char* label, int verbosity, int rows, int cols, int nb, int rank, int cprocs, int sprocs, int failing_level, int checkpoint_freq)
{
	clock_t start, stop;
	double span, maxspan;
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

	//MPI_Barrier(MPI_COMM_WORLD);
	start=clock();

	ScaLAPACK_pDGEQRF_cp_ft1_sim(rows, A, nb, rank, cprocs, sprocs, failing_level, checkpoint_freq);

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
		//DeallocateMatrix1D(bb);
	}

	MPI_Barrier(MPI_COMM_WORLD);
	return maxspan;
}
