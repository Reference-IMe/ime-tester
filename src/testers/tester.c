/*
 * tester.c
 *
 *  Created on: Apr 25, 2020
 *      Author: marcello
 */

#include <mpi.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

#include "../helpers/macros.h"
#include "../helpers/matrix.h"
#include "../helpers/matrix_advanced.h"
#include "../helpers/vector.h"
#include "../helpers/lapack.h"
#include "../helpers/scalapack.h"
#include "tester_labels.h"
#include "tester_routine.h"
#include "tester_structures.h"

#define MAX_VERSIONS 30
#define MAX_RUNS 10

int versionnumber_in(int n_all, char** all, char* selected)
{
	int i=n_all-1;
	while (i>=0)
	{
		if( strcmp( all[i], selected ) == 0 )
		{
			break;
		}
		i--;
	}
	return i;
}

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

void dsort(test_result* run, int n)
{
    /* Sort the given array number, of length n */
    int j, i;
    test_result temp;

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

test_result dmedian(test_result* run, int n)
{
	dsort(run, n);
	return run[n/2];
}

test_result not_run={-1, -1, -1, -1};
test_result not_implemented={-99,-99, -1, -1};

int main(int argc, char **argv)
{
	/*
	 * parallel environment setup
	 */
	int rank;			// mpi rank
	int totprocs;		// total num. of mpi ranks
	int thread_support; // level of provided thread support

	//MPI_Init(&argc, &argv);												// NOT working with OpenMP
	MPI_Init_thread(&argc ,&argv, MPI_THREAD_FUNNELED, &thread_support);	// OK for OpenMP
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);									// get current process id
	MPI_Comm_size(MPI_COMM_WORLD, &totprocs);								// get number of processes

    int ont;			// number of OpenMP threads set with OMP_NUM_THREADS
	int np;				// number of total processes

	if (getenv("OMP_NUM_THREADS")==NULL)
	{
		ont=1;
	}
	else
	{
		ont=atoi(getenv("OMP_NUM_THREADS"));
	}
	np=totprocs*ont;

	/*
	 * application variables
	 */
    int i,j,rep;

	int			versions_all;
	int			versions_selected;
    char*		versionname_all[MAX_VERSIONS];
    char*		versionname_selected[MAX_VERSIONS];
    int			versionnumber_selected[MAX_VERSIONS];
    test_result	versionrun[MAX_VERSIONS][MAX_RUNS];
    test_result	versiontot[MAX_VERSIONS];

    int n;					// matrix size
    int nrhs;				// num. of r.h.s
    int cnd;				// condition number
    int seed;				// seed for random matrix generation
    int rows;
    int cols;
    int scalapack_iter;		// scalapack total iterations
    int scalapack_nb;		// scalapack blocking factor
    int ime_nb;				// ime blocking factor

    int sprocs;				// number of processes to allocate for summing (0 = no fault tolerance)
    int cprocs;				// number of processes for real IMe calc
    int repetitions;
    int verbose;
    int failing_rank;
    int failing_level;
    int failing_level_override;
    int checkpoint_skip_interval;

    char* file_name;
    int   file_name_len;
	FILE* fp;

    /*
     * default values
     */
    n=8;
	nrhs=1;
    verbose=1;					// minimum output verbosity (0 = none)
    repetitions=1;				// how many calls for each routine
    sprocs=0;					// no fault tolerance enabled
    file_name_len=0;			// no output to file
    failing_rank=2;				// process 2 will fail
    failing_level_override=-1;
    checkpoint_skip_interval=-1;// -1=never, otherwise do at every (checkpoint_skip_interval+1) iteration
    scalapack_nb=SCALAPACKNB;	// scalapack blocking factor, default defined in header
    ime_nb=1;					// ime blocking factor
    cnd=1;						// condition number for randomly generated matrices
    seed=1;						// seed for random generation

    // list of testable routines (see tester_labels.h)
    versions_all = 0;
	versionname_all[versions_all++] = IME_SV;
	//versionname_all[versions_all++] = IME_SV_CHECKSUMMED;
	versionname_all[versions_all++] = IME_SV_FAULT_0_TOLERANT_1;
	versionname_all[versions_all++] = IME_SV_FAULT_1_TOLERANT_1;
	/*
	versionname_all[versions_all++] = IME_XK;
	versionname_all[versions_all++] = IME_XK_FAULT_0_TOLERANT_1;
	versionname_all[versions_all++] = IME_XK_FAULT_1_TOLERANT_1;
	*/
	versionname_all[versions_all++] = SPK_SV;
	versionname_all[versions_all++] = SPK_SV_FAULT_0_TOLERANT_1;
	versionname_all[versions_all++] = SPK_SV_FAULT_1_TOLERANT_1;

	versionname_all[versions_all++] = SPK_LU;
	versionname_all[versions_all++] = SPK_LU_FAULT_0_TOLERANT_1;
	versionname_all[versions_all++] = SPK_LU_FAULT_1_TOLERANT_1;

	versionname_all[versions_all++] = SPK_QR;
	versionname_all[versions_all++] = SPK_QR_FAULT_0_TOLERANT_1;
	versionname_all[versions_all++] = SPK_QR_FAULT_1_TOLERANT_1;

	versionname_all[versions_all++] = FTLA_LU_FAULT_0_TOLERANT_1;
	versionname_all[versions_all++] = FTLA_LU_FAULT_1_TOLERANT_1;

	versionname_all[versions_all++] = FTLA_QR_FAULT_0_TOLERANT_1;
	versionname_all[versions_all++] = FTLA_QR_FAULT_1_TOLERANT_1;

    /*
     * read command line parameters
     */
	for( i = 1; i < argc; i++ )
	{
		if( strcmp( argv[i], "-n" ) == 0 ) {
			n = atoi(argv[i+1]);
			i++;
		}
		if( strcmp( argv[i], "-nrhs" ) == 0 ) {
			nrhs = atoi(argv[i+1]);
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
			failing_level_override = 1;
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
		if( strcmp( argv[i], "--help" ) == 0 ) {
			if (rank==0)
			{
				printf("HELP\n");
			}
			MPI_Finalize();
			return(1);
		}
		if( strcmp( argv[i], "--run" ) == 0 ) {
			i++;
			// read functions to be tested
			for( j = i; j < argc; j++ )
			{
				versionname_selected[j-i]=argv[j];
			}
			versions_selected=j-i;
			i=j-1;
		}
		if( strcmp( argv[i], "--list" ) == 0 ) {
			i++;
			// print testable functions
			if (rank==0)
			{
				printf("Testable routines:\n");
				for( j = 0; j < versions_all; j++ )
				{
					printf("     %2d %s\n", j+1 ,versionname_all[j]);
				}
			}
			MPI_Finalize();
			return(1);
		}
	}

	/*
	 * other default values, depending on inputs
	 */
	rows=n;
    cols=n;
    cprocs=totprocs-sprocs;		// number of processes for real IMe calc
    scalapack_iter=(int)ceil(rows/scalapack_nb);

    if (failing_level_override<0) // if faulty level NOT set on command line
    {
        failing_level=n/2;			// faulty level/iteration, -1=none
    }

	/*
	 * print initial summary to video
	 */
	if (rank==0 && verbose>0)
	{
		printf("     Total processes:               %d\n",np);
		printf("     OMP threads:                   %d\n",ont);
		printf("     MPI ranks:                     %d\n",totprocs);
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

		printf("     Testing routines:\n");
		for( j = 0; j < versions_selected; j++ )
		{
			printf("      %2d %s\n", j+1 ,versionname_selected[j]);
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
			printf("WRN: ScaLAPACK blocking factor probably too small\n");;
		}
    }

	// check matrix size
    if ((n % cprocs) != 0)
    {
    	if (rank==0)
    	{
    		printf("ERR: The size of the matrix has to be a multiple of the number (%d) of calc. nodes\n",cprocs);
    	}
    	MPI_Finalize();
        return 1;
    }

    // check list of routines
	for (i=0; i<versions_selected; i++)
	{
		versionnumber_selected[i]=versionnumber_in(versions_all, versionname_all, versionname_selected[i]);
		if (versionnumber_selected[i]<0)
		{
	    	if (rank==0)
	    	{
	    		printf("ERR: Routine '%s' is unknown\n",versionname_selected[i]);
	    	}
	    	MPI_Finalize();
	        return 1;
		}
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
	if (rank==0)
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
	else
	{
		A_ref = NULL;
		x_ref = NULL;
		b_ref = NULL;
	}

	/*
	 * get time now
	 */
	time_t rawtime;
	struct tm *readtime;
	time ( &rawtime );
	readtime = localtime ( &rawtime );

	/*
	 * print initial summary to file
	 */
	#define fpinfo(string_label,integer_info) if (fp!=NULL && rank==0) {fprintf(fp,"info,%s,%d\n",string_label,integer_info); \
	}
	#define fpdata(track_num) if (fp!=NULL && rank==0) {																			\
		fprintf(fp,"data,%s,%d,%d,%.0f,%f\n",versionname_all[track_num], rep+1,	versionrun[track_num][rep].exit_code,				\
																			versionrun[track_num][rep].total_time,					\
																			versionrun[track_num][rep].norm_rel_err);				\
		fprintf(fp,"data,%s%s,%d,%d,%.0f,%f\n",versionname_all[track_num],"(core)", rep+1,	versionrun[track_num][rep].exit_code,	\
																						versionrun[track_num][rep].core_time,		\
																						versionrun[track_num][rep].norm_rel_err);	\
	}

	if (file_name_len>0 && rank==0)
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
		fpinfo("number of MPI ranks",cprocs);
		fpinfo("number of OMP threads",ont);
		fpinfo("number of processes",np);
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

	// init totals to 0 before accumulation
	for (i=0; i<versions_selected; i++)
	{
		versiontot[i].total_time   = 0;
		versiontot[i].core_time    = 0;
		versiontot[i].norm_rel_err = 0;
		versiontot[i].exit_code    = -1;
	}

	test_input routine_input = {
			n,
			A_ref,
			x_ref,
			b_ref,
			nrhs,
			cprocs,
			sprocs,
			ime_nb,
			scalapack_nb
	};

	/*
	 * main loop
	 */
	// runs are repeated
	for (rep=0; rep<repetitions; rep++)
	{
		if (rank==0 && verbose>0) {printf("\n Run #%d:\n",rep+1);}

		// every run calls some selected routines
		for (i=0; i<versions_selected; i++)
		{
			versionrun[i][rep]=tester_routine(versionname_selected[i], verbose, routine_input, rank, failing_rank, failing_level, checkpoint_skip_interval);
		}

		if (rank==0)
		{
			// accumulation of totals
			for (i=0; i<versions_selected; i++)
			{
				versiontot[i].total_time	+= versionrun[i][rep].total_time;
				versiontot[i].core_time		+= versionrun[i][rep].core_time;
				versiontot[i].norm_rel_err	+= versionrun[i][rep].norm_rel_err;
				if (verbose>0)
				{
					printf("%-20s    call    run time: %10.0f (%.0f)\ts\t nre: %f\n",	versionname_selected[i],		\
																					versionrun[i][rep].total_time,	\
																					versionrun[i][rep].core_time,	\
																					versionrun[i][rep].norm_rel_err	\
					);
				}
			}
		}
	}

	/*
	 * print final summary to video and to file
	 */
	if (rank==0)
	{
		printf("\n Summary:\n");

		// total
		for (i=0; i<versions_selected; i++)
		{
			printf("%-20s    Total   run time: %10.0f (%.0f)\ts\n",	versionname_selected[i],	\
																	versiontot[i].total_time,	\
																	versiontot[i].core_time		\
			);
		}
		printf("\n");

		// average
		for (i=0; i<versions_selected; i++)
		{
			printf("%-20s    Average run time: %10.0f (%.0f)\ts\t nre: %f\n",	versionname_selected[i],				\
																				versiontot[i].total_time/repetitions,	\
																				versiontot[i].core_time/repetitions,	\
																				versiontot[i].norm_rel_err/repetitions	\
			);
			if (file_name_len>0)
			{
				fprintf(fp,"data,%s,%d,%d,%.0f,%f\n",	versionname_selected[i], 0,				\
														versionrun[i][0].exit_code,				\
														versiontot[i].total_time/repetitions,	\
														versiontot[i].norm_rel_err/repetitions	\
				);
				fprintf(fp,"data,%s%s,%d,%d,%.0f,%f\n",	versionname_selected[i], "(core)", 0,	\
														versionrun[i][0].exit_code,				\
														versiontot[i].core_time/repetitions,	\
														versiontot[i].norm_rel_err/repetitions	\
				);
			}
		}
		printf("\n");

		// median
		for (i=0; i<versions_selected; i++)
		{
			dmedian( versionrun[i], repetitions);
			printf("%-20s    Median  run time: %10.0f (%.0f)\ts\t nre: %f\n",	versionname_selected[i],					\
																				versionrun[i][repetitions/2].total_time,	\
																				versionrun[i][repetitions/2].core_time,		\
																				versionrun[i][repetitions/2].norm_rel_err	\
			);
			if (file_name_len>0)
			{
				fprintf(fp,"data,%s,%d,%d,%.0f,%f\n",	versionname_selected[i], -1,				\
														versionrun[i][0].exit_code,					\
														versionrun[i][repetitions/2].total_time,	\
														versionrun[i][repetitions/2].norm_rel_err	\
				);
				fprintf(fp,"data,%s%s,%d,%d,%.0f,%f\n",	versionname_selected[i], "(core)",-1,		\
														versionrun[i][0].exit_code,					\
														versionrun[i][repetitions/2].core_time,		\
														versionrun[i][repetitions/2].norm_rel_err	\
				);
			}
		}
		printf("\n");
	}

	/*
	 * slow down exit
	 */
	//sleep(3);
	//MPI_Barrier(MPI_COMM_WORLD);
	//printf("Done %d.\n",rank);

	/*
	 * cleanup
	 */
	if (rank==0)
	{
		if (file_name_len>0)
		{
			fclose(fp);
		}
		DeallocateMatrix1D(A_ref);
		DeallocateVector(b_ref);
		DeallocateVector(x_ref);
	}

	MPI_Finalize();
	return 0;
}

