#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "../helpers/matrix.h"
#include "../helpers/vector.h"

#include "test_IMe_pviDGESV.h"
#include "test_IMe_pviDGESV_swaploop.h"
#include "test_IMe_pviDGESV_unopt.h"
#include "test_IMe_pviDGESV_cs.h"
#include "test_IMe_pviDGESV_blind.h"
#include "test_IMe_pviDGESV_parinit.h"

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    int rank, totprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);		/* get current process id */
    MPI_Comm_size(MPI_COMM_WORLD, &totprocs);	/* get number of processes */

    int i,rep;

    double versionrun[20][100];
    const char* versionname[20];
    double versiontot[20];

    int n=atoi(argv[1]);

    int rows=n;
    int cols=n;
    int sprocs=atoi(argv[2]);		// number of processes to allocate for summing (0 = no fault tolerance)
    int cprocs=totprocs-sprocs;		// number of processes for real IMe calc
    int repetitions=atoi(argv[3]);
    int verbose=atoi(argv[4]);

    int nRHS=10;

	versionname[0]="IMe-uo  1 ";
	versionname[1]="IMe-uo  10";
	versionname[2]="IMe-sl  1 ";
	versionname[3]="IMe-sl  10";
	versionname[4]="IMe-pin 1 ";
	versionname[5]="IMe-pin 10";
	versionname[6]="IMe-cs  1 ";
	versionname[7]="IMe-cs  10";
	versionname[8]="IMe-bld 1 ";
	versionname[9]="IMe-bld 10";
	versionname[10]="IMe     1 ";
	versionname[11]="IMe     10";

	int versions = 12;

	for (i=0; i<versions; i++)
	{
		versiontot[i] = 0;
	}

	if (rank==0 && verbose>0)
	{
		printf("\nMatrix size: %dx%d",n,n);
		printf("\nCheckpoint : ");
		if(sprocs==0)
		{
			printf("no\n");
		}
		else
		{
			printf("yes\n");
		}
	}

    for (rep=0; rep<repetitions; rep++)
    {
    	if (rank==0 && verbose>0) {printf("\n Run #%d",rep+1);}

    	//////////////////////////////////////////////////////////////////////////////////

     	versionrun[0][rep]=test_IMe_pviDGESV_unopt(versionname[0], verbose, rows, cols, 1, rank, cprocs);  			// IMe unoptimized with 1 rhs
    	versionrun[1][rep]=test_IMe_pviDGESV_unopt(versionname[1], verbose, rows, cols, nRHS, rank, cprocs);		// IMe unoptimized with 10 rhs
    	versionrun[2][rep]=test_IMe_pviDGESV_swaploop(versionname[2], verbose, rows, cols, 1, rank, cprocs);    	// IMe swap loop with 1 rhs
     	versionrun[3][rep]=test_IMe_pviDGESV_swaploop(versionname[3], verbose, rows, cols, nRHS, rank, cprocs);		// IMe swap loop with 10 rh
    	versionrun[4][rep]=test_IMe_pviDGESV_parinit(versionname[4], verbose, rows, cols, 1, rank, cprocs);    				// IMe parallel init with 1 rhs
    	versionrun[5][rep]=test_IMe_pviDGESV_parinit(versionname[5], verbose, rows, cols, nRHS, rank, cprocs);				// IMe parallel init with 10 rhs
    	versionrun[6][rep]=test_IMe_pviDGESV_cs(versionname[6], verbose, rows, cols, 1, rank, cprocs, sprocs);    	// IMe checksumming with 1 rhs
    	versionrun[7][rep]=test_IMe_pviDGESV_cs(versionname[7], verbose, rows, cols, nRHS, rank, cprocs, sprocs);	// IMe checksumming with 10 rhs
    	versionrun[8][rep]=test_IMe_pviDGESV_blind(versionname[8], verbose, rows, cols, 1, rank, cprocs);    		// IMe blind with 1 rhs
    	versionrun[9][rep]=test_IMe_pviDGESV_blind(versionname[9], verbose, rows, cols, nRHS, rank, cprocs);		// IMe blind with 10 rhs
    	versionrun[10][rep]=test_IMe_pviDGESV(versionname[10], verbose, rows, cols, 1, rank, cprocs);    		// IMe latest optimization with 1 rhs
    	versionrun[11][rep]=test_IMe_pviDGESV(versionname[11], verbose, rows, cols, nRHS, rank, cprocs);		// IMe latest optimization with 10 rhs

		//////////////////////////////////////////////////////////////////////////////////

    	if (rank==0)
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

	if (rank==0)
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

	MPI_Finalize();
    return(0);
}
