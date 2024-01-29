#include <mpi.h>
#include <time.h>
#include "../helpers/macros.h"
#include "../helpers/matrix_advanced.h"
#include "../testers/IMe/lib/src/helpers/matrix_basic.h"
#include "tester_structures.h"


void test_init(const char* label, int verbosity, parallel_env env, test_input input)
{
	double** xx;
	double x;


	if (env.mpi_rank==0)
	{
		if (verbosity>1)
		{
			printf("     Opening MPI communication channels..\n");
		}
	}

	xx=AllocateMatrix2D_double(1, env.calc_procs+env.spare_procs, CONTIGUOUS);
	xx[0][0]=0;
	MPI_Scatter (&xx[0][0], 1, MPI_DOUBLE, &x, 1, MPI_DOUBLE, 0, env.mpi_comm);
	DeallocateMatrix2D_double(xx, 1, CONTIGUOUS);
	if (env.mpi_rank==0)
	{
		if (verbosity>1)
		{
			printf("     ..comm_world\n");
		}
	}

	MPI_Comm comm_calc;
	int i_calc; // participating in ime calc = 1, checksumming = 0

	if (env.mpi_rank >= env.calc_procs)
	{
		i_calc=0;
		MPI_Comm_split(env.mpi_comm, MPI_UNDEFINED, MPI_UNDEFINED, &comm_calc); // checksumming procs don't belong to calc communicator
	}
	else
	{
		i_calc=1;
		MPI_Comm_split(env.mpi_comm, i_calc, env.mpi_rank, &comm_calc); // calc procs belong to calc communicator
	}

	if (i_calc)
	{
		xx=AllocateMatrix2D_double(input.n, input.nrhs, CONTIGUOUS);
		MPI_Scatter (&xx[0][0], 1, MPI_DOUBLE, &x, 1, MPI_DOUBLE, 0, comm_calc);
		DeallocateMatrix2D_double(xx, input.n, CONTIGUOUS);
		if (env.mpi_rank==0)
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

	if (env.spare_procs>0)
	{
		if (comm_calc != MPI_COMM_NULL)
		{
			MPI_Comm_free(&comm_calc);
		}
	}

	if (env.mpi_rank==0)
	{
		if (verbosity>1)
		{
			printf("     ..done\n");
		}
	}

}
