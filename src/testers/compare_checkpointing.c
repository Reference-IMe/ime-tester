#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

#include "../helpers/matrix.h"
#include "../helpers/vector.h"

#include "test_ScaLAPACK_pDGETRF.h"
#include "test_ScaLAPACK_pDGEQRF.h"
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
    int checkpoint_skip_interval;

    int scalapack_iter;

	FILE* fp;

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
    checkpoint_skip_interval=-1; // -1=never, otherwise do at every (checkpoint_skip_interval+1) iteration (after the first)
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
		if( strcmp( argv[i], "-cp" ) == 0 ) {
			checkpoint_skip_interval = atoi(argv[i+1]);
			i++;
		}
	}

	rows=n;
    cols=n;
    cprocs=totprocs-sprocs;		// number of processes for real IMe calc
    scalapack_iter=ceil(rows/SCALAPACKNB);

    if (main_rank==0 && verbose>0)
    {
		//printf("     Matrix condition number:       %d\n",cnd);
		//printf("     Matrix random generation seed: %d\n",seed);
		printf("     Matrix size:                   %dx%d\n",rows,cols);
		printf("     IMe iterations:                %d\n",rows);
		printf("     SPK-like iterations:           %d\n",scalapack_iter);

		printf("     Fault tolerance:               ");
		if (sprocs>0)
		{
			printf("enabled = %d\n",sprocs);
			printf("     IMe failing rank:              %d\n",failing_rank);
			printf("     IMe failing level:             %d\n",failing_level);
			printf("     SPK-like failing level:        %d\n",n-failing_level);
			printf("     SPK-like failing iteration:    %d\n",(int)ceil((n-failing_level)/SCALAPACKNB));
			printf("     Checkpoint skip interval:      %d\n",checkpoint_skip_interval);

			printf("     Checkpoint freq.:              ");
			if (checkpoint_skip_interval<0)
			{
				printf("never\n");
			}
			else
			{
				if (checkpoint_skip_interval==0)
				{
					printf("always\n");
				}
				else
				{
					printf("every %d iterations\n",checkpoint_skip_interval+1);
				}
			}
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

	versionname[0]= "SPK-LU          ";
	versionname[1]= "SPK-LU-ft1/0    ";
	versionname[2]= "SPK-LU-ft1/1    ";

	versionname[3]= "SPK-QR          ";
	versionname[4]= "SPK-QR-ft1/0    ";
	versionname[5]= "SPK-QR-ft1/1    ";


	int versions = 6;

	for (i=0; i<versions; i++)
	{
		versiontot[i] = 0;
	}

#define fpinfo(string_label,integer_info) if (main_rank==0) {fprintf(fp,"info,%s,%d\n",string_label,integer_info);}
#define fpdata(track_num) if (main_rank==0) {fprintf(fp,"data,%s,%d,%.0f\n",versionname[track_num],rep+1,versionrun[ track_num][rep]);}

	time_t rawtime;
	struct tm *readtime;
	time ( &rawtime );
	readtime = localtime ( &rawtime );

	if (file_name_len>0 && main_rank==0)
	{
		fp=fopen(file_name,"w");
		fpinfo("year",readtime->tm_year+1900);
		fpinfo("month",readtime->tm_mon+1);
		fpinfo("day",readtime->tm_mday);
		fpinfo("hour",readtime->tm_hour);
		fpinfo("minute",readtime->tm_min);
		fpinfo("second",readtime->tm_sec);
		fpinfo("number of processes",totprocs);
		fpinfo("fault tolerance",sprocs);
		fpinfo("failing rank",failing_rank);
		fpinfo("failing level",failing_level);
		fpinfo("checkpoint skip interval",checkpoint_skip_interval);
		fpinfo("matrix size",n);
		fpinfo("repetitions",repetitions);

		fprintf(fp,"head,code name,run num. (0=avg),run time\n");
	}


    for (rep=0; rep<repetitions; rep++)
    {
    	if (main_rank==0 && verbose>0) {printf("\n Run #%d",rep+1);}

    	//////////////////////////////////////////////////////////////////////////////////

    	versionrun[0][rep]=test_ScaLAPACK_pDGETRF(versionname[0], verbose, rows, cols, main_rank, cprocs, sprocs); // SPK LU factorization
 		fpdata(0);

		versionrun[1][rep]=test_ScaLAPACK_pDGETRF_cp_ft1_sim(versionname[1], verbose, rows, cols, main_rank, cprocs, sprocs, failing_rank, -1, checkpoint_skip_interval);	// SPKmod LU factorization single FT and 0 faults
		fpdata(1);

		versionrun[2][rep]=test_ScaLAPACK_pDGETRF_cp_ft1_sim(versionname[2], verbose, rows, cols, main_rank, cprocs, sprocs, failing_rank, failing_level, checkpoint_skip_interval);	// SPKmod LU factorization single FT and 1 fault
		fpdata(2);

    	versionrun[3][rep]=test_ScaLAPACK_pDGEQRF(versionname[3], verbose, rows, cols, main_rank, cprocs, sprocs); // SPK LU factorization
    	fpdata(3);

		versionrun[4][rep]=test_ScaLAPACK_pDGEQRF_cp_ft1_sim(versionname[4], verbose, rows, cols, main_rank, cprocs, sprocs, failing_rank, -1, checkpoint_skip_interval);	// SPKmod QR factorization single FT and 0 faults
		fpdata(4);

		versionrun[5][rep]=test_ScaLAPACK_pDGEQRF_cp_ft1_sim(versionname[5], verbose, rows, cols, main_rank, cprocs, sprocs, failing_rank, failing_level, checkpoint_skip_interval);	// SPKmod QR factorization single FT and 1 fault
		fpdata(5);

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

			if (file_name_len>0)
			{
				fprintf(fp,"data,%s,%d,%.0f\n",versionname[i],0,(versiontot[i]/repetitions));
			}
		}
		printf("\n");
	}


	// slow down exit
	//sleep(3);
	//MPI_Barrier(MPI_COMM_WORLD);
	//printf("Done %d.\n",main_rank);

	if (file_name_len>0 && main_rank==0)
	{
		fclose(fp);
	}

	MPI_Finalize();
    return(0);
}
