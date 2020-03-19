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
#include "../helpers/matrix_advanced.h"
#include "../helpers/vector.h"
#include "../helpers/lapack.h"
#include "../helpers/scalapack.h"

#define MAX_VERSIONS 30
#define MAX_RUNS 10

/*
 * header code for verbatim inclusion to create a code tester for some parallel versions
 */

char* faketrim(char* str)
{
	// credits to https://codeforwin.org/2016/04/c-program-to-trim-trailing-white-space-characters-in-string.html

	int index, i;

	/* Set default index */
	index = -1;

	/* Find last index of non-white space character */
	i = 0;
	while(str[i] != '\0')
	{
		if(str[i] != ' ' && str[i] != '\t' && str[i] != '\n')
		{
			index= i;
		}

		i++;
	}

	/* Mark next character to last non-white space character as NULL */
	str[index + 1] = '\0';
	return str;
}

void dsort(run_info* run, int n)
{
    /* Sort the given array number, of length n */
    int j, i;
    run_info temp;

    for (i = 1; i < n; i++)
    {
        for (j = 0; j < n - i; j++)
        {
            if (run[j].total_time > run[j + 1].total_time )
            {
                temp = run[j];
                run[j] = run[j + 1];
                run[j + 1] = temp;
            }
        }
    }
}

run_info dmedian(run_info* run, int n)
{
	dsort(run, n);
	return run[n/2];
}

run_info not_run={-1, -1, -1, -1};
run_info not_implemented={-99,-99, -1, -1};

int main(int argc, char **argv)
{
    int main_rank, totprocs, mpisupport; //

    //MPI_Init(&argc, &argv);
    MPI_Init_thread(&argc ,&argv, MPI_THREAD_FUNNELED, &mpisupport);

    MPI_Comm_rank(MPI_COMM_WORLD, &main_rank);	//get current process id
    MPI_Comm_size(MPI_COMM_WORLD, &totprocs);	// get number of processes

    int i,rep;

	int versions;
    char* versionname[MAX_VERSIONS];
    run_info versionrun[MAX_VERSIONS][MAX_RUNS];
    run_info versiontot[MAX_VERSIONS];

    int n;		// matrix size
    int file_name_len;
    char* file_name;
    int rows;
    int cols;
    int scalapack_nb;		// (scalapack) blocking factor
    int ime_nb;

    int cnd;	// condition number
    int seed;	// seed for random matrix generation

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
    verbose=1;					// minimum output verbosity (0 = none)
    repetitions=1;				// how many calls for each routine
    sprocs=0;					// no fault tolerance enabled
    file_name_len=0;			// no output to file
    failing_rank=2;				// process 2 will fail
    checkpoint_skip_interval=-1;// -1=never, otherwise do at every (checkpoint_skip_interval+1) iteration
    scalapack_nb=SCALAPACKNB;	// scalapack blocking factor, default defined in header
    ime_nb=1;					// ime blocking factor
    cnd=1;						// condition number for randomly generated matrices
    seed=1;						// seed for random generation

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
		if( strcmp( argv[i], "-spk-nb" ) == 0 ) {
			scalapack_nb = atoi(argv[i+1]);
			i++;
		}
		if( strcmp( argv[i], "-ime-nb" ) == 0 ) {
			ime_nb = atoi(argv[i+1]);
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
	}

	/*
	 * other default values, depending on inputs
	 */
	rows=n;
    cols=n;
    cprocs=totprocs-sprocs;		// number of processes for real IMe calc
    scalapack_iter=(int)ceil(rows/scalapack_nb);
    failing_level=n/2;

    /*
     * print summary to video
     */
    if (main_rank==0 && verbose>0)
    {
    	printf("     Matrix condition number:       %d\n",cnd);
		printf("     Matrix random generation seed: %d\n",seed);
		printf("     Matrix size:                   %dx%d\n",rows,cols);
		printf("     IMe iterations:                %d\n",rows);
		printf("     IMe blocking factor:           %d\n",ime_nb);
		printf("     SPK-like iterations:           %d\n",scalapack_iter);
		printf("     SPK-like blocking factor:      %d\n",scalapack_nb);

		printf("     Fault tolerance:               ");
		if (sprocs>0)
		{
			printf("enabled = %d\n",sprocs);
			printf("     IMe failing rank:              %d\n",failing_rank);
			printf("     IMe failing level:             %d\n",failing_level);
			printf("     SPK-like failing level:        %d\n",n-failing_level);
			printf("     SPK-like failing iteration:    %d\n",(int)ceil((n-failing_level)/scalapack_nb));
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

		if (n/cprocs < scalapack_nb )
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

	/*
	 * prepare input matrices
	 */
	int read_cnd;
	double* A_ref;
	double* x_ref;
	double* b_ref;
	char transA = 'T', transx = 'N';
	double one = 1.0, zero = 0.0;
	int m=1;
	if (main_rank==0)
	{
		A_ref = AllocateMatrix1D(n, n);
		x_ref = AllocateVector(n);
		b_ref = AllocateVector(n);
		RandomSquareMatrix1D_cnd(A_ref, n, seed, cnd);
		read_cnd = round(ConditionNumber1D(A_ref, n, n));
		if (read_cnd!=cnd && verbose>0)
		{
			printf("WRN: Condition number (%d) differs from read back (%d)\n",cnd,read_cnd);
		}
		FillVector(x_ref, n, 1);
		dgemm_(&transA, &transx, &n, &m, &n, &one, A_ref, &n, x_ref, &n, &zero, b_ref, &n);
	}

	/*
	 * print summary to file
	 */
	#define fpinfo(string_label,integer_info) if (fp!=NULL && main_rank==0) {fprintf(fp,"info,%s,%d\n",string_label,integer_info); \
	}
	#define fpdata(track_num) if (fp!=NULL && main_rank==0) { \
		fprintf(fp,"data,%s,%d,%d,%.0f,%f\n",versionname[track_num], rep+1,	versionrun[track_num][rep].exit_code,     \
																			versionrun[track_num][rep].total_time,    \
																			versionrun[track_num][rep].norm_rel_err); \
		fprintf(fp,"data,%s%s,%d,%d,%.0f,%f\n",versionname[track_num],"(core)", rep+1,	versionrun[track_num][rep].exit_code,     \
																						versionrun[track_num][rep].core_time,    \
																						versionrun[track_num][rep].norm_rel_err); \
	}

	time_t rawtime;
	struct tm *readtime;
	time ( &rawtime );
	readtime = localtime ( &rawtime );

	if (file_name_len>0 && main_rank==0)
	{
		fp=fopen(file_name,"w");

		fprintf(fp,"info,command,");
		for (i=0;i<argc;i++)
		{
			fprintf(fp,"%s ",argv[i]);
		}
		fprintf(fp,"\n");

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
		fpinfo("seed",seed);
		fpinfo("condition number set",cnd);
		fpinfo("condition number got",read_cnd);
		fpinfo("matrix size",n);
		fpinfo("scalapack blocking factor",scalapack_nb);
		fpinfo("ime blocking factor",ime_nb);
		fpinfo("repetitions",repetitions);

		fprintf(fp,"head,code name,run num. (0=avg,-1=mdn),exit code (0=ok),run time,rel.err.\n");
	}
	else
	{
		fp=NULL;
	}
	;





