#include <mpi.h>
#include <time.h>
#include "../helpers/macros.h"
#include "../helpers/matrix.h"
#include "../helpers/matrix_advanced.h"
#include "tester_structures.h"


void test_dummy(const char* label, int verbosity, test_input input, int rank, MPI_Comm comm_world)
{
	double** xx;
	double x;


	if (rank==0)
	{
		if (verbosity>1)
		{
			printf("\n\n     Opening MPI communication channels..\n");
		}
	}

	xx=AllocateMatrix2D(1, input.calc_procs+input.spare_procs, CONTIGUOUS);
	MPI_Scatter (&xx[0][0], 1, MPI_DOUBLE, &x, 1, MPI_DOUBLE, 0, comm_world);
	DeallocateMatrix2D(xx, 1, CONTIGUOUS);
	if (rank==0)
	{
		if (verbosity>1)
		{
			printf("     ..comm_world\n");
		}
	}

	MPI_Comm comm_calc;
	int i_calc; // participating in ime calc = 1, checksumming = 0

	if (rank >= input.calc_procs)
	{
		i_calc=0;
		MPI_Comm_split(comm_world, MPI_UNDEFINED, MPI_UNDEFINED, &comm_calc); // checksumming procs don't belong to calc communicator
	}
	else
	{
		i_calc=1;
		MPI_Comm_split(comm_world, i_calc, rank, &comm_calc); // calc procs belong to calc communicator
	}

	if (i_calc)
	{
		xx=AllocateMatrix2D(input.n, input.nrhs, CONTIGUOUS);
		MPI_Scatter (&xx[0][0], 1, MPI_DOUBLE, &x, 1, MPI_DOUBLE, 0, comm_calc);
		DeallocateMatrix2D(xx, input.n, CONTIGUOUS);
		if (rank==0)
		{
			if (verbosity>1)
			{
				printf("     ..comm_calc\n");
			}
		}
	}
	else
	{
		xx=NULL;
	}

	if (input.spare_procs>0)
	{
		if (comm_calc != MPI_COMM_NULL)
		{
			MPI_Comm_free(&comm_calc);
		}
	}

	if (rank==0)
	{
		if (verbosity>1)
		{
			printf("     ..done\n");
		}
	}

}
