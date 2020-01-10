#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

#include "../helpers/matrix.h"
#include "../helpers/vector.h"
#include "test_IMe_pviDGESV.h"
#include "test_IMe_pviDGESV_cs.h"
#include "test_IMe_pviDGESV_ft1.h"
#include "test_ScaLAPACK_pDGESV.h"
#include "test_ScaLAPACK_pDGETRF.h"
#include "test_ScaLAPACK_pDGEQRF.h"
#include "test_FTLA_pDGEQRF.h"
#include "test_FTLA_pDGETRF.h"
#include "test_ScaLAPACK_pDGESV_cp_ft1.h"
#include "test_ScaLAPACK_pDGETRF_cp_ft1.h"
#include "test_ScaLAPACK_pDGEQRF_cp_ft1.h"


int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    int main_rank, totprocs; //
    MPI_Comm_rank(MPI_COMM_WORLD, &main_rank);	//get current process id
    MPI_Comm_size(MPI_COMM_WORLD, &totprocs);	// get number of processes

    int i,rep;

#define MAX_VERSIONS 30
#define MAX_RUNS 10

    double versionrun[MAX_VERSIONS][MAX_RUNS];
    const char* versionname[MAX_VERSIONS];
    double versiontot[MAX_VERSIONS];

    int n;		// matrix size
    int cnd;	// condition number
    int seed;	// seed for random matrix generation
    int file_name_len;
    char* file_name;
    int rows;
    int cols;
    int cleanup_interval;

    int sprocs;		// number of processes to allocate for summing (0 = no fault tolerance)
    int cprocs;		// number of processes for real IMe calc
    int repetitions;
    int verbose;
    int failing_rank;
    int failing_level;

    /*
     * default values
     */
    n=8;
    cnd=1;
    seed=1;
    verbose=1;			// minimum verbosity
    repetitions=1;
    sprocs=0;			// no fault tolerance enabled
    file_name_len=0;	// no output to file
    failing_rank=2;		// process 2 will fail
    failing_level=2;

    cleanup_interval=3;
    /*
     * read command line parameters
     */
	for( i = 1; i < argc; i++ )
	{
		if( strcmp( argv[i], "-n" ) == 0 ) {
			n = atoi(argv[i+1]);
			i++;
		}
		if( strcmp( argv[i], "-cnd" ) == 0 ) {
			cnd = atoi(argv[i+1]);
			i++;
		}
		if( strcmp( argv[i], "-seed" ) == 0 ) {
			seed = atoi(argv[i+1]);
			i++;
		}
		if( strcmp( argv[i], "-v" ) == 0 ) {
			verbose = atoi(argv[i+1]);
			i++;
		}
		if( strcmp( argv[i], "-r" ) == 0 ) {
			repetitions = atoi(argv[i+1]);
			i++;
		}
		if( strcmp( argv[i], "-ft" ) == 0 ) {
			sprocs = atoi(argv[i+1]);
			i++;
		}
		if( strcmp( argv[i], "-o" ) == 0 ) {
			file_name_len = strlen(argv[i+1]);
			file_name = malloc(file_name_len+1);
			strcpy(file_name, argv[i+1]);
			i++;
		}
		if( strcmp( argv[i], "-fr" ) == 0 ) {
			failing_rank = atoi(argv[i+1]);
			i++;
		}
		if( strcmp( argv[i], "-fl" ) == 0 ) {
			failing_level = atoi(argv[i+1]);
			i++;
		}
	}

	rows=n;
    cols=n;
    cprocs=totprocs-sprocs;		// number of processes for real IMe calc

    if (main_rank==0 && verbose>0)
    {
		printf("     Matrix size:                   %dx%d\n",rows,cols);
		printf("     Matrix condition number:       %d\n",cnd);
		printf("     Matrix random generation seed: %d\n",seed);

		printf("     Fault tolerance:               ");
		if (sprocs>0)
		{
			printf("enabled = %d\n",sprocs);
			printf("     Failing rank:                  %d\n",failing_rank);
			printf("     Failing level:                 %d\n",failing_level);
		}
		else
		{
			printf("disabled\n");
		}

		printf("     Run repetitions:               %d\n",repetitions);

		if (file_name_len>0)
		{
			printf("     Output file:                   %s\n",file_name);
		}
		else
		{
			printf("WRN: No output to file\n");
		}
    }

    if ((n % cprocs) != 0)
    {
    	if (main_rank==0)
    	{
    		printf("ERR: The size of the matrix has to be a multiple of the number (%d) of calc. nodes\n",cprocs);
    	}
    	MPI_Finalize();
        return(1);
    }
    int nRHS=10;
	versionname[ 0]= "IMe-SV        1 ";
	versionname[ 1]= "IMe-SV        10";
	versionname[ 2]= "IMe-SV-cs     1 ";
	versionname[ 3]= "IMe-SV-cs     10";
	versionname[ 4]= "IMe-SV-ft1/0  1 ";
	versionname[ 5]= "IMe-SV-ft1/0  10";
	versionname[ 6]= "IMe-SV-ft1/1  1 ";
	versionname[ 7]= "IMe-SV-ft1/1  10";
	versionname[ 8]= "SPK-SV        1 ";
	versionname[ 9]= "SPK-SV        10";
	versionname[10]= "SPK-SV-ft1    1 ";
	versionname[11]= "SPK-SV-ft1    10";

	versionname[12]= "IMe-XK-         ";
	versionname[13]= "IMe-XK-ft1/0    ";
	versionname[14]= "IMe-XK-ft1/1    ";

	versionname[15]= "SPK-LU          ";
	versionname[16]= "SPK-LU-ft1/0    ";
	versionname[17]= "SPK-LU-ft1/1    ";

	versionname[18]= "FTLA-LU-ft1/0   ";
	versionname[19]= "FTLA-LU-ft1/1   ";

	versionname[20]= "SPK-QR          ";
	versionname[21]= "SPK-QR-ft1/0    ";
	versionname[22]= "SPK-QR-ft1/1    ";

	versionname[23]= "FTLA-QR-ft1/0   ";
	versionname[24]= "FTLA-QR-ft1/1   ";


	int versions = 25;

	for (i=0; i<versions; i++)
	{
		versiontot[i] = 0;
	}


    for (rep=0; rep<repetitions; rep++)
    {
    	if (main_rank==0 && verbose>0) {printf("\n Run #%d",rep+1);}

    	//////////////////////////////////////////////////////////////////////////////////

     	versionrun[ 0][rep]=test_IMe_pviDGESV(versionname[0], verbose, rows, cols, 1, main_rank, cprocs, sprocs);		// vanilla IMe solve with 1 rhs
    	versionrun[ 1][rep]=test_IMe_pviDGESV(versionname[1], verbose, rows, cols, nRHS, main_rank, cprocs, sprocs);	// vanilla IMe solve with 10 rhs
     	versionrun[ 2][rep]=test_IMe_pviDGESV_cs(versionname[2], verbose, rows, cols, 1, main_rank, cprocs, sprocs);	// checksumming IMe solve with 1 rhs
    	versionrun[ 3][rep]=test_IMe_pviDGESV_cs(versionname[3], verbose, rows, cols, nRHS, main_rank, cprocs, sprocs);	// checksumming IMe solve with 10 rhs
    	versionrun[ 4][rep]=test_IMe_pviDGESV_ft1_sim(versionname[4], verbose, rows, cols, 1, main_rank, cprocs, sprocs, -1, -1);	// IMe single FT solve with 1 rhs and 0 faults
    	versionrun[ 5][rep]=test_IMe_pviDGESV_ft1_sim(versionname[5], verbose, rows, cols, nRHS, main_rank, cprocs, sprocs, -1, -1);// IMe single FT solve with 10 rhs and 0 faults

    	if (sprocs>0)
    	{
    		versionrun[ 6][rep]=test_IMe_pviDGESV_ft1_sim(versionname[6], verbose, rows, cols, 1, main_rank, cprocs, sprocs, failing_rank, failing_level);   // IMe single FT solve with 1 rhs
    		versionrun[ 7][rep]=test_IMe_pviDGESV_ft1_sim(versionname[7], verbose, rows, cols, nRHS, main_rank, cprocs, sprocs, failing_rank, failing_level);// IMe single FT solve with 10 rhs
    	}
    	else
    	{
    		versionrun[ 6][rep]=-1; // don't run IMe single FT solve with 1 rhs
    		versionrun[ 7][rep]=-1; // don't run IMe single FT solve with 10 rhs
    	}

    	versionrun[ 8][rep]=test_ScaLAPACK_pDGESV(versionname[8], verbose, rows, cols, 1, main_rank, cprocs, sprocs);		// SPK solve with 1 rhs
    	versionrun[ 9][rep]=test_ScaLAPACK_pDGESV(versionname[9], verbose, rows, cols, nRHS, main_rank, cprocs, sprocs);	// SPK solve with 10 rhs

    	if (sprocs>0)
    	{
    		versionrun[10][rep]=-99; // not yet SPKmod single FT solve with 1 rhs
    		versionrun[11][rep]=-99; // not yet SPKmod single FT solve with 10 rhs
    	}
    	else
    	{
    		versionrun[10][rep]=-1; // don't run SPKmod single FT solve with 1 rhs
    		versionrun[11][rep]=-1; // don't run SPKmod single FT solve with 10 rhs
    	}

     	versionrun[12][rep]=-99; // not yet IMe factorization only
    	versionrun[13][rep]=-99; // not yet IMe factorization only and 0 faults
    	versionrun[14][rep]=-99; // not yet IMe factorization only and 1 fault

    	versionrun[15][rep]=test_ScaLAPACK_pDGETRF(versionname[15], verbose, rows, cols, main_rank, cprocs, sprocs); // SPK LU factorization

		if (sprocs>0)
    	{
    		versionrun[16][rep]=test_ScaLAPACK_pDGETRF_cp_ft1_sim(versionname[16], verbose, rows, cols, main_rank, cprocs, sprocs, failing_rank, -n);	// SPKmod LU factorization single FT and 0 faults
    		versionrun[17][rep]=test_ScaLAPACK_pDGETRF_cp_ft1_sim(versionname[17], verbose, rows, cols, main_rank, cprocs, sprocs, failing_rank, failing_level);	// SPKmod LU factorization single FT and 1 fault
        	versionrun[18][rep]=test_FTLA_pDGETRF(versionname[18], verbose, rows, cols, main_rank, cprocs, 0);	// FTLA LU with 0 faults
        	versionrun[19][rep]=test_FTLA_pDGETRF(versionname[19], verbose, rows, cols, main_rank, cprocs, 1);	// FTLA LU with 1 fault
    	}
    	else
    	{
    		versionrun[16][rep]=-1;	// don't run SPKmod LU factorization single FT and 0 faults
    		versionrun[17][rep]=-1;	// don't run SPKmod LU factorization single FT and 1 fault
    		versionrun[18][rep]=-1;	// don't run FTLA LU with 0 faults
        	versionrun[19][rep]=-1;	// don't run FTLA LU with 1 fault
    	}

    	versionrun[20][rep]=test_ScaLAPACK_pDGEQRF(versionname[20], verbose, rows, cols, main_rank, cprocs, sprocs); // SPK LU factorization

    	if (sprocs>0)
    	{
    		versionrun[21][rep]=test_ScaLAPACK_pDGEQRF_cp_ft1_sim(versionname[21], verbose, rows, cols, main_rank, cprocs, sprocs, failing_rank, -n);	// SPKmod QR factorization single FT and 0 faults
    		versionrun[22][rep]=test_ScaLAPACK_pDGEQRF_cp_ft1_sim(versionname[22], verbose, rows, cols, main_rank, cprocs, sprocs, failing_rank, failing_level);	// SPKmod QR factorization single FT and 1 fault
        	versionrun[23][rep]=test_FTLA_pDGEQRF(versionname[23], verbose, rows, cols, main_rank, cprocs, 0);	// FTLA QR with 0 faults
        	versionrun[24][rep]=test_FTLA_pDGEQRF(versionname[24], verbose, rows, cols, main_rank, cprocs, 1);	// FTLA QR with 1 fault
    	}
    	else
    	{
    		versionrun[21][rep]=-1;	// don't run SPKmod QR factorization single FT and 0 faults
    		versionrun[22][rep]=-1;	// don't run SPKmod QR factorization single FT and 1 fault
    		versionrun[23][rep]=-1;	// don't run FTLA QR with 0 faults
    		versionrun[24][rep]=-1;	// don't run FTLA QR with 1 fault
    	}

    	//////////////////////////////////////////////////////////////////////////////////

    	if (main_rank==0)
		{
			for (i=0; i<versions; i++)
			{
				versiontot[i] += versionrun[i][rep];
				if (verbose>0)
				{
					printf("\n%s    call    run time: %10.0f clk", versionname[i], versionrun[i][rep]);
				}
			}
		}
    }

	if (main_rank==0)
	{
		printf("\n\n Summary:");
		for (i=0; i<versions; i++)
		{
			printf("\n%s    Total   run time: %10.0f clk", versionname[i], versiontot[i]); // in sec. versiontot[i] / CLOCKS_PER_SEC
		}
		printf("\n");
		for (i=0; i<versions; i++)
		{
			printf("\n%s    Average run time: %10.0f clk", versionname[i], (versiontot[i]/repetitions));
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
