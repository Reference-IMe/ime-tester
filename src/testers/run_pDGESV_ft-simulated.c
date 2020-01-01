#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "../helpers/matrix.h"
#include "../helpers/vector.h"

#include "test_IMe_pviDGESV.h"
#include "test_IMe_pviDGESV_ft.h"

// globals
#define FRONTEND_RANK 0

int failing_rank;
int failed_rank;
int failing_level;
int faulty = 0;
int faulty_location=0;

MPI_Comm ftla_current_comm, current_world_comm, ccomm;

// needed for FT
#include <signal.h>
void common_error_handler(int s)
{
	int i,j ;
	int rank, world_size;

	// avoid race conditions (https://digilander.libero.it/uzappi/C/librerie/funzioni/signal.html)
	signal(s, SIG_IGN);

	MPI_Comm_rank(ftla_current_comm, &rank);
	MPI_Comm_size(ftla_current_comm, &world_size);

	//Localize fault
	MPI_Group failed_group;
	MPIX_Comm_failure_ack(ftla_current_comm);
	MPIX_Comm_failure_get_acked(ftla_current_comm, &failed_group);
	//MPI_Comm_set_errhandler(current_comm, MPI_ERRORS_RETURN);
	//signal(SIGUSR1, common_handler_error);

	MPI_Group current_group;
	MPI_Comm_group(ftla_current_comm, &current_group);
	int *ranks = (int *)malloc( world_size * sizeof(int) );
	int *ranks_out = (int *)malloc( 1 * sizeof(int) );
	for(i=0; i<world_size;i++)
		ranks[i]=i;
	ranks_out[0]=0;

	MPI_Group_translate_ranks(failed_group, 1, ranks, current_group, ranks_out );
	failed_rank = ranks_out[0];

	printf("\n%d knows that %d failed at %d\n",rank,failed_rank,faulty_location);
	if (failed_rank==FRONTEND_RANK)
	{
		//frontend_rank++;
	}

	MPI_Comm surviving_comm;
	MPI_Group tmp_group, surviving_group;
	MPI_Comm_group(ftla_current_comm, &tmp_group);
	MPI_Group_excl(tmp_group, 1, &failed_rank, &surviving_group);
	MPI_Comm_create_group(ftla_current_comm, surviving_group, 0, &surviving_comm);

	free(ranks);
	free(ranks_out);
	MPI_Group_free(&failed_group);
	MPI_Group_free(&tmp_group);
	MPI_Group_free(&surviving_group);

	MPI_Comm_free(&ftla_current_comm);
	ftla_current_comm=surviving_comm;

	signal(SIGUSR1, common_error_handler);
}



int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);

    int rank, totprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);		/* get current process id */
    MPI_Comm_size(MPI_COMM_WORLD, &totprocs);	/* get number of processes */

    // init
    int i,rep;

    double versionrun[20][100];
    const char* versionname[20];
    double versiontot[20];

    int n=atoi(argv[1]);
    int rows=n;
    int cols=n;

    int sprocs=atoi(argv[2]);		// number of processes to allocate for summing (0 = no fault tolerance)
    int cprocs=totprocs-sprocs;		// number of processes for real IMe calc

	current_world_comm=MPI_COMM_WORLD;
    if (sprocs==0)
    {

    	ccomm=MPI_COMM_WORLD;
    }
    else
    {
        // prepare for FT
    	// rank management: http://mpi.deino.net/mpi_functions/MPI_Group_range_excl.html

    	if (rank < cprocs)
    	{
    		MPI_Comm_split(MPI_COMM_WORLD, rank < cprocs, rank, &ccomm );
    		MPI_Comm_split(MPI_COMM_WORLD, MPI_UNDEFINED, rank, &scomm );
    	}
    	else
    	{
    		MPI_Comm_split(MPI_COMM_WORLD, MPI_UNDEFINED, rank, &ccomm );
    		MPI_Comm_split(MPI_COMM_WORLD, rank >= cprocs, rank, &scomm );
    	}
    	/*
        MPI_Group world_group;
    	MPI_Comm_group(MPI_COMM_WORLD, &world_group);
    	MPI_Comm_create(MPI_COMM_WORLD, world_group, &current_comm);
        MPI_Comm_rank(current_comm, &rank);		// get current process id

        MPI_Comm_size(current_comm, &totprocs);	// get number of processes
       	MPI_Comm_set_errhandler(current_comm, MPI_ERRORS_RETURN);
       	signal(SIGUSR1, common_error_handler);
       	*/
    }


    int crank;
    if (rank < cprocs)
    {
    	MPI_Comm_rank(ccomm, &crank);		// get current process id in calc communicator
    }
    else
    {
    	crank=-1;
    }
    //
    MPI_Barrier(MPI_COMM_WORLD);
    printf("world %d = calc %d\n",rank,crank);
    MPI_Barrier(MPI_COMM_WORLD);

    int repetitions=atoi(argv[3]);

    int verbose=atoi(argv[4]);

    failing_rank=atoi(argv[5]);
    failing_level=atoi(argv[6]);

    int nRHS=10;

	versionname[0]="IMe     1 ";
	versionname[1]="IMe-ft  1 ";

	int versions = 2;

	for (i=0; i<versions; i++)
	{
		versiontot[i] = 0;
	}

	if (rank==FRONTEND_RANK && verbose>0)
	{
		printf("\nMatrix size: %dx%d",n,n);
		printf("\nChecksum(s): %d",sprocs);
	}

    for (rep=0; rep<repetitions; rep++)
    {
    	if (rank==FRONTEND_RANK && verbose>0) {printf("\n Run #%d",rep+1);}

    	//////////////////////////////////////////////////////////////////////////////////

    	if (rank<cprocs)
    	{
    		//versionrun[0][rep]=test_IMe_pviDGESV(versionname[0], verbose, rows, cols, 1, ccomm);    		// IMe latest optimization with 1 rhs
    	}
    	//versionrun[1][rep]=test_IMe_pviDGESV_ft(versionname[1], verbose, rows, cols, 1, ccomm, sprocs);		// IMe latest optimization with 10 rhs

		//////////////////////////////////////////////////////////////////////////////////

    	if (rank==FRONTEND_RANK)
		{
			for (i=0; i<versions; i++)
			{
				versiontot[i] += versionrun[i][rep];
				if (verbose>0)
				{
					printf("\n%s    call    run time: %f clk", versionname[i], versionrun[i][rep]);
				}
			}
		}
    }

	if (rank==FRONTEND_RANK)
	{
		printf("\n\n Summary:");
		for (i=0; i<versions; i++)
		{
			printf("\n%s    Total   run time: %f clk\t\t%f s", versionname[i], versiontot[i], versiontot[i] / CLOCKS_PER_SEC);
		}
		printf("\n");
		for (i=0; i<versions; i++)
			{
				printf("\n%s    Average run time: %f clk\t\t%f s", versionname[i], versiontot[i]/repetitions, versiontot[i]/repetitions / CLOCKS_PER_SEC);
			}
		printf("\n");
	}

	MPI_Barrier(current_world_comm);
	printf("\nrank %d finished\n",rank);

	MPI_Finalize();
    return(0);
}
