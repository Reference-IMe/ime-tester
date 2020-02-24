/*
 * tester_head_p.c
 *
 *  Created on: Jan 17, 2020
 *      Author: marcello
 */

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

#include "../helpers/info.h"
#include "../helpers/matrix.h"
#include "../helpers/vector.h"
#include "../helpers/scalapack.h"

#define MAX_VERSIONS 30
#define MAX_RUNS 10

/*
 * header code for verbatim inclusion to create a code tester for some parallel versions
 */


void dsort(duration_t* duration, int n)
{
    /* Sort the given array number, of length n */
    int j, i;
    duration_t temp;

    for (i = 1; i < n; i++)
    {
        for (j = 0; j < n - i; j++)
        {
            if (duration[j].total > duration[j + 1].total )
            {
                temp = duration[j];
                duration[j] = duration[j + 1];
                duration[j + 1] = temp;
            }
        }
    }
}

duration_t dmedian(duration_t* duration, int n)
{
	dsort(duration, n);
	return duration[n/2];
}

duration_t not_run={-1, -1};
duration_t not_implemented={-99,-99};

int main(int argc, char **argv)
{
    int main_rank, totprocs, mpisupport; //

    //MPI_Init(&argc, &argv);
    MPI_Init_thread(&argc ,&argv, MPI_THREAD_FUNNELED, &mpisupport);

    MPI_Comm_rank(MPI_COMM_WORLD, &main_rank);	//get current process id
    MPI_Comm_size(MPI_COMM_WORLD, &totprocs);	// get number of processes

    int i,rep;

	int versions;
    const char* versionname[MAX_VERSIONS];
    duration_t versionrun[MAX_VERSIONS][MAX_RUNS];
    duration_t versiontot[MAX_VERSIONS];

    int n;		// matrix size
    int file_name_len;
    char* file_name;
    int rows;
    int cols;
    int nb;		// (scalapack) blocking factor

    //TODO
    //int cnd;	// condition number
    //int seed;	// seed for random matrix generation

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
    verbose=1;			// minimum verbosity
    repetitions=1;
    sprocs=0;			// no fault tolerance enabled
    file_name_len=0;	// no output to file
    failing_rank=2;		// process 2 will fail
    checkpoint_skip_interval=-1; // -1=never, otherwise do at every (checkpoint_skip_interval+1) iteration
    //cleanup_interval=3;
    nb=SCALAPACKNB;
    //TODO
    //cnd=1;
    //seed=1;

    /*
     * read command line parameters
     */
	for( i = 1; i < argc; i++ )
	{
		if( strcmp( argv[i], "-n" ) == 0 ) {
			n = atoi(argv[i+1]);
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
		if( strcmp( argv[i], "-nb" ) == 0 ) {
			nb = atoi(argv[i+1]);
			i++;
		}
		//TODO
		/*
			if( strcmp( argv[i], "-cnd" ) == 0 ) {
			cnd = atoi(argv[i+1]);
			i++;
		}
		if( strcmp( argv[i], "-seed" ) == 0 ) {
			seed = atoi(argv[i+1]);
			i++;
		}
		*/
	}

	rows=n;
    cols=n;
    cprocs=totprocs-sprocs;		// number of processes for real IMe calc
    scalapack_iter=(int)ceil(rows/nb);
    failing_level=n/2;

    if (main_rank==0 && verbose>0)
    {
		//TODO
    	//printf("     Matrix condition number:       %d\n",cnd);
		//printf("     Matrix random generation seed: %d\n",seed);
		printf("     Matrix size:                   %dx%d\n",rows,cols);
		printf("     IMe iterations:                %d\n",rows);
		printf("     SPK-like iterations:           %d\n",scalapack_iter);
		printf("     SPK-like blocking factor:      %d\n",nb);

		printf("     Fault tolerance:               ");
		if (sprocs>0)
		{
			printf("enabled = %d\n",sprocs);
			printf("     IMe failing rank:              %d\n",failing_rank);
			printf("     IMe failing level:             %d\n",failing_level);
			printf("     SPK-like failing level:        %d\n",n-failing_level);
			printf("     SPK-like failing iteration:    %d\n",(int)ceil((n-failing_level)/nb));
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

		if (n/cprocs < nb )
		{
			printf("WRN: Blocking factor probably too small\n");;
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

#define fpinfo(string_label,integer_info) if (fp!=NULL && main_rank==0) {fprintf(fp,"info,%s,%d\n",string_label,integer_info);}
#define fpdata(track_num) if (fp!=NULL && main_rank==0) { \
	fprintf(fp,"data,%s,%d,%.0f\n",versionname[track_num],rep+1,versionrun[ track_num][rep].total); \
	fprintf(fp,"data,%s%s,%d,%.0f\n",versionname[track_num],"(core)",rep+1,versionrun[ track_num][rep].core); \
	}


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
		//TODO
		//fpinfo("seed",seed);
		//fpinfo("condition number",cnd);
		fpinfo("matrix size",n);
		fpinfo("blocking factor",nb);
		fpinfo("repetitions",repetitions);

		fprintf(fp,"head,code name,run num. (0=avg),run time\n");
	}
	else
	{
		fp=NULL;
	}
