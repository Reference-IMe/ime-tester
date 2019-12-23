#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "../helpers/matrix.h"
#include "../helpers/vector.h"

#include "test_IMe_pviDGESV.h"
#include "test_IMe_pviDGESV_cs.h"
#include "test_IMe_pviDGESV_ft1.h"
#include "test_ScaLAPACK_pDGESV_ckp_ft1.h"
#include "test_ScaLAPACK_pDGETRF_ckp_ft1.h"


int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    int main_rank, totprocs; //
    MPI_Comm_rank(MPI_COMM_WORLD, &main_rank);	//get current process id
    MPI_Comm_size(MPI_COMM_WORLD, &totprocs);	// get number of processes

    int i,rep;

    double versionrun[20][100];
    const char* versionname[20];
    double versiontot[20];

    int n;
    int rows;
    int cols;

    int sprocs;		// number of processes to allocate for summing (0 = no fault tolerance)
    int cprocs;		// number of processes for real IMe calc
    int repetitions;
    int verbose;
    int failing_rank;
    int failing_level;

    if ((argc != 5) && (argc != 7))
    {
    	if (main_rank==0)
    	{
			printf("ERR: Wrong number of parameters\n");
			printf("     %s matrix_size fault_tolerance repetitions verbosity\n",argv[0]);
			printf("     %s matrix_size fault_tolerance repetitions verbosity failing_rank failing_level\n",argv[0]);
    	}
    	MPI_Finalize();
        return(1);
    }

    n=atoi(argv[1]);
    rows=n;
    cols=n;
    sprocs=atoi(argv[2]);		// number of processes to allocate for summing (0 = no fault tolerance)
    cprocs=totprocs-sprocs;		// number of processes for real IMe calc
    repetitions=atoi(argv[3]);
    verbose=atoi(argv[4]);

    if (argc == 7)
    {
        failing_rank=atoi(argv[5]);
        failing_level=atoi(argv[6]);

		if (failing_rank < 0 || failing_rank >= totprocs)
		{
			if (main_rank==0)
			{
				printf("ERR: The failing rank is out of range\n");
			}
			MPI_Finalize();
			return(1);
		}
		if (failing_level < -1 || failing_level >= n) // if failing_level is -1, fault tolerance is enabled, but no fault will occur
		{
			if (main_rank==0)
			{
				printf("ERR: The failing level is out of range\n");
			}
			MPI_Finalize();
			return(1);
		}
    }

    if ((n % cprocs) != 0)
    {
    	if (main_rank==0)
    	{
    		printf("ERR: The size of the matrix has to be a multiple of the number of calc. nodes\n");
    	}
    	MPI_Finalize();
        return(1);
    }
    int nRHS=10;

	versionname[0]="IMe-cs     1 ";
	versionname[1]="IMe-cs     10";
	versionname[2]="IMe-ft1    1 ";
	versionname[3]="IMe-ft1    10";
	versionname[4]="SPK-SV-ft1 1 ";
	versionname[5]="SPK-SV-ft1 10";
	versionname[6]="SPK-LU-ft1   ";

	int versions = 7;

	for (i=0; i<versions; i++)
	{
		versiontot[i] = 0;
	}

	if (main_rank==0 && verbose>0)
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
    	if (main_rank==0 && verbose>0) {printf("\n Run #%d",rep+1);}

    	//////////////////////////////////////////////////////////////////////////////////

     	versionrun[0][rep]=test_IMe_pviDGESV_cs(versionname[0], verbose, rows, cols, 1, main_rank, cprocs, sprocs);    	// IMe checksumming with 1 rhs
    	versionrun[1][rep]=test_IMe_pviDGESV_cs(versionname[1], verbose, rows, cols, nRHS, main_rank, cprocs, sprocs);	// IMe checksumming with 10 rhs
    	versionrun[2][rep]=test_IMe_pviDGESV_ft1_sim(versionname[2], verbose, rows, cols, 1, main_rank, cprocs, sprocs, failing_rank, failing_level);	// IMe single FT with 1 rhs
    	versionrun[3][rep]=test_IMe_pviDGESV_ft1_sim(versionname[3], verbose, rows, cols, nRHS, main_rank, cprocs, sprocs, failing_rank, failing_level);// IMe single FT with 10 rhs
    	versionrun[4][rep]=test_Scalapack_pDGESV_ckp_ft1_sim(versionname[4], verbose, rows, cols, 1, main_rank, cprocs, sprocs, failing_rank, failing_level);
    	versionrun[5][rep]=test_Scalapack_pDGESV_ckp_ft1_sim(versionname[5], verbose, rows, cols, nRHS, main_rank, cprocs, sprocs, failing_rank, failing_level);
    	versionrun[6][rep]=test_Scalapack_pDGETRF_ckp_ft1_sim(versionname[6], verbose, rows, cols, 1, main_rank, cprocs, sprocs, failing_rank, failing_level);
    	//versionrun[7][rep]=test_Scalapack_pDGETRF_ckp_ft1_sim(versionname[7], verbose, rows, cols, nRHS, main_rank, cprocs, sprocs, failing_rank, failing_level);
    	//versionrun[2][rep]=test_IMe_pviDGESV(versionname[2], verbose, rows, cols, 1, main_rank, cprocs);    		// IMe latest optimization with 1 rhs
    	//versionrun[3][rep]=test_IMe_pviDGESV(versionname[3], verbose, rows, cols, nRHS, main_rank, cprocs);		// IMe latest optimization with 10 rhs

    	//////////////////////////////////////////////////////////////////////////////////

    	if (main_rank==0)
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

	if (main_rank==0)
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

	// slow down exit
	//sleep(3);
	//MPI_Barrier(MPI_COMM_WORLD);
	//printf("Done %d.\n",main_rank);

	MPI_Finalize();
    return(0);
}
