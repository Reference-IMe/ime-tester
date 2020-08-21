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
#include "../helpers/Cblacs.h"
#include "../helpers/simple_dynamic_strings/sds.h"

#include "tester_labels.h"
#include "tester_routine.h"
#include "tester_structures.h"
#include "test_dummy.h"

#define MAX_VERSIONS 50
#define MAX_RUNS 10

// TODO: rename auxiliary functions

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



int main(int argc, char **argv)
{

	/*
	 * *********************
	 * application variables
	 * *********************
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
    char verbose;
    int failing_rank;
    int failing_level;
    int failing_level_override;
    int checkpoint_skip_interval;

    sds matrix_output_base_name;
    sds matrix_output_file_name;
    sds matrix_input_base_name;
    sds matrix_input_file_name;
    sds test_output_file_name;

    char calc_nre;
    int output_to_file;
    int input_from_file;

	FILE* fp;

    char* command;

    /*
     * ************
     * input values
     * ************
     */

		/*
		 * FIXED defaults
		 */
		matrix_output_base_name = NULL;
		matrix_output_file_name = NULL;
		matrix_input_base_name  = NULL;
		matrix_input_file_name  = NULL;
		test_output_file_name   = NULL;
		fp = NULL;
		output_to_file = 0;			// no output to file
		input_from_file = 0;		// no input from file
		command="null";
		n=8;
		nrhs=1;
		verbose=1;					// minimum output verbosity (0 = none)
		repetitions=1;				// how many calls for each routine
		sprocs=0;					// no fault tolerance enabled
		failing_rank=2;				// process 2 will fail
		failing_level_override=-1;
		checkpoint_skip_interval=-1;// -1=never, otherwise do at every (checkpoint_skip_interval+1) iteration
		scalapack_nb=SCALAPACKNB;	// scalapack blocking factor, default defined in header
		ime_nb=1;					// ime blocking factor
		cnd=1;						// condition number for randomly generated matrices
		seed=1;						// seed for random generation
		calc_nre=1;					// calc (1) normwise relative error or not (0)

		/*
		 * list of testable routines (see tester_labels.h)
		 */

		// TODO: add a description for every routine

		versions_all = 0;
		versionname_all[versions_all++] = "dummy";

		versionname_all[versions_all++] =  IME_SV_WO_OAE;
		versionname_all[versions_all++] =  IME_SV_WO_OA;
		versionname_all[versions_all++] =  IME_SV_WO_OGE;
		versionname_all[versions_all++] =  IME_SV_WO_OG;

		versionname_all[versions_all++] =  IME_SV_WO_U1AE;
		versionname_all[versions_all++] =  IME_SV_WO_U1A;
		versionname_all[versions_all++] =  IME_SV_WO_U1GE;
		versionname_all[versions_all++] =  IME_SV_WO_U1G;

		versionname_all[versions_all++] =  IME_SV_WO_U2AE;
		versionname_all[versions_all++] =  IME_SV_WO_U2A;
		versionname_all[versions_all++] =  IME_SV_WO_U2GE;
		versionname_all[versions_all++] =  IME_SV_WO_U2G;

		versionname_all[versions_all++] =  IME_SV_WO_U3AE;
		versionname_all[versions_all++] =  IME_SV_WO_U3A;
		versionname_all[versions_all++] =  IME_SV_WO_U3GE;
		versionname_all[versions_all++] =  IME_SV_WO_U3G;

		//versionname_all[versions_all++] = IME_SV_WO_CHECKSUMMED;
		versionname_all[versions_all++] = IME_SV_WO_FAULT_0_TOLERANT_1;
		versionname_all[versions_all++] = IME_SV_WO_FAULT_1_TOLERANT_1;
		/*
		versionname_all[versions_all++] = IME_XK_WO;
		versionname_all[versions_all++] = IME_XK_WO_FAULT_0_TOLERANT_1;
		versionname_all[versions_all++] = IME_XK_WO_FAULT_1_TOLERANT_1;
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

		versions_selected = 0;

		/*
		 * read command line parameters
		 */
		// TODO: checking for valid arguments
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
				output_to_file = 1;
				test_output_file_name = sdsnew(argv[i+1]);
				i++;
			}
			if( strcmp( argv[i], "-i" ) == 0 ) {
				input_from_file = 1;
				matrix_input_base_name = sdsnew(argv[i+1]);
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
			if( strcmp( argv[i], "-no-nre" ) == 0 ) {
				calc_nre = 0;
				i++;
			}
			if( strcmp( argv[i], "--run" ) == 0 ) {
				command=argv[i];
				i++;
				// read functions to be tested
				for( j = i; j < argc; j++ )
				{
					versionname_selected[j-i]=argv[j];
				}
				versions_selected=j-i;
				//i=j-1;
				break;
			}
			if( strcmp( argv[i], "--list" ) == 0 ) {
				command=argv[i];
				break;
			}
			if( strcmp( argv[i], "--help" ) == 0 ) {
				command=argv[i];
				break;
			}
			if( strcmp( argv[i], "--save" ) == 0 ) {
				command=argv[i];
				i++;
				if (i < argc)
				{
					matrix_output_base_name = sdsnew(argv[i]);
				}
				else
				{
					matrix_output_base_name = NULL;
				}
				break;
			}
		}

		/*
		 * CALCULATED defaults (depending on command line inputs)
		 */
		rows=n;
		cols=n;
		scalapack_iter=(int)ceil(rows/scalapack_nb);

		if (failing_level_override<0) 	// if faulty level NOT set on command line
		{
			failing_level=n/2;			// faulty level/iteration, -1=none
		}

	/*
	 * **************************
	 * parallel environment setup
	 * **************************
	 */
		/*
		 * MPI + openMP
		 */
		int mpi_rank;				// mpi rank
		int mpi_procs;				// total num. of mpi ranks
		int mpi_thread_support; 	// level of provided thread support
		//MPI_Init(&argc, &argv);													// NOT working with OpenMP
		MPI_Init_thread(&argc ,&argv, MPI_THREAD_FUNNELED, &mpi_thread_support);	// OK for OpenMP
		MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);									// get current process id
		MPI_Comm_size(MPI_COMM_WORLD, &mpi_procs);									// get number of processes

		cprocs=mpi_procs-sprocs;	// number of MPI processes for real calc

		// OpenMP
		int omp_threads;			// number of OpenMP threads set with OMP_NUM_THREADS
		int np;						// number of total processes

		if (getenv("OMP_NUM_THREADS")==NULL)
		{
			omp_threads=1;
		}
		else
		{
			omp_threads=atoi(getenv("OMP_NUM_THREADS"));
		}

		np=mpi_procs*omp_threads;

		/*
		 * BLACS
		 */
		int ndims = 2, dims[2] = {0,0};
		MPI_Dims_create(cprocs, ndims, dims);

		int blacs_nprow, blacs_npcol;
		blacs_nprow = dims[0];
		blacs_npcol = dims[1];

		int i0 = 0;
		int i1 = 1;
		int ic = -1;
		char order = 'R';

		// general grid context
		int blacs_ctxt;
		Cblacs_get( ic, i0, &blacs_ctxt );
		Cblacs_gridinit( &blacs_ctxt, &order, blacs_nprow, blacs_npcol );

		// general single row context (1 x mpi_procs)
		int blacs_ctxt_onerow;
		Cblacs_get( ic, i0, &blacs_ctxt_onerow );
		Cblacs_gridinit( &blacs_ctxt_onerow, &order, i1, mpi_procs );

		// root node context (for holding global matrices) at {1,1}
		int blacs_ctxt_root;
		Cblacs_get( ic, i0, &blacs_ctxt_root );
		Cblacs_gridinit( &blacs_ctxt_root, &order, i1, i1 );

		int blacs_ctxt_cp;
		int blacs_row;
		int blacs_col;

		if (sprocs>0) // fault tolerance enabled
		{
			// context for the checkpointing node
			Cblacs_get( ic, i0, &blacs_ctxt_cp );
			int map_cp[1];
			map_cp[0]=cprocs;
			Cblacs_gridmap( &blacs_ctxt_cp, map_cp, i1, i1, i1);
		}
		else
		{
			blacs_ctxt_cp = -1;
		}

		// get coords in general grid context
		Cblacs_gridinfo( blacs_ctxt, &blacs_nprow, &blacs_npcol, &blacs_row, &blacs_col );

	/*
	 * ****************************
	 * common input data structures
	 * ****************************
	 */
		parallel_env routine_env = {
			mpi_rank,
			blacs_nprow,
			blacs_npcol,
			blacs_row,
			blacs_col,
			blacs_ctxt_onerow,
			blacs_ctxt,
			blacs_ctxt_root,
			blacs_ctxt_cp
		};

		test_input routine_input = {
				n,
				NULL,
				NULL,
				NULL,
				nrhs,
				cprocs,
				sprocs,
				ime_nb,
				scalapack_nb,
				calc_nre
		};

	/*
	 * ******************************
	 * print initial summary to video
	 * ******************************
	 */
	if (mpi_rank==0 && verbose>0)
	{
		printf("     Total processes:               %d\n",np);
		printf("     OMP threads:                   %d\n",omp_threads);
		printf("     MPI ranks:                     %d\n",mpi_procs);
		printf("     Matrix condition number:       %d\n",cnd);
		printf("     Matrix random generation seed: %d\n",seed);
		printf("     Matrix size:                   %dx%d\n",rows,cols);
		printf("     Calculate n.r.e.:              ");
			if (calc_nre)	{printf("yes\n");}
			else			{printf("no\n");}
		printf("     IMe iterations:                %d\n",rows);
		printf("     IMe blocking factor:           %d\n",ime_nb);
		printf("     SPK-like iterations:           %d\n",scalapack_iter);
		printf("     SPK-like blocking factor:      %d\n",scalapack_nb);

		printf("     Fault tolerance:               ");
		if (sprocs>0)
		{
			printf("enabled = %d\n",sprocs);
			printf("       Calc. processes:             %d\n",cprocs);
			printf("       Spare processes:             %d\n",sprocs);
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

		if (output_to_file)
		{
			printf("     Output file:                   %s\n",test_output_file_name);
		}
		else
		{
			printf("WRN: No output to file\n");
		}
	}

	/*
	 * ********
	 * checking
	 * ********
	 */

	// check list of selected routines
	for (i=0; i<versions_selected; i++)
	{
		versionnumber_selected[i]=versionnumber_in(versions_all, versionname_all, versionname_selected[i]);
		if (versionnumber_selected[i]<0)
		{
	    	if (mpi_rank==0)
	    	{
	    		printf("ERR: Routine '%s' is unknown\n",versionname_selected[i]);
	    	}
	    	MPI_Finalize();
	        return 3;
		}
	}

	// check specific condition for every selected routine
	j=0; // error accumulation
	for (i=0; i<versions_selected; i++)
	{
		j = j + (tester_routine(1, versionname_selected[i], verbose, routine_env, routine_input, mpi_rank, failing_rank, failing_level, checkpoint_skip_interval)).exit_code;
	}
	MPI_Bcast(&j,1,MPI_INT,0,MPI_COMM_WORLD);
	if (j != 0)
	{
    	MPI_Finalize();
        return 4;
	}
	// continue if no errors

	/*
	 * **********************
	 * prepare input matrices
	 * **********************
	 */
	int read_cnd = -99;
	double* A_ref;
	double* x_ref;
	double* b_ref;
	char transA = 'T', transx = 'N';
	double d1 = 1.0;
	double d0 = 0.0;
	int m=1;
	if (mpi_rank==0)
	{
		A_ref = AllocateMatrix1D(n, n);
		x_ref = AllocateVector(n);
		b_ref = AllocateVector(n);
	}
	else
	{
		A_ref = NULL;
		x_ref = NULL;
		b_ref = NULL;
	}
		/*
		 * matrices from file
		 */
		if (input_from_file)
		{
			if ( strcmp(command, "--save" ) == 0 )
			{
				if (mpi_rank == 0) printf("ERR: Cannot take input matrices with '--save' command.\n");
				MPI_Finalize();
				return 1;
			}
			if (mpi_rank==0)
			{
				if (verbose > 0) printf("     Loading reference matrices from binary files..\n");

				matrix_input_file_name = sdsdup(matrix_input_base_name);
				matrix_input_file_name = sdscat(matrix_input_file_name, ".X");
				fp=fopen(matrix_input_file_name,"rb");
				fread(x_ref,sizeof(double),n,fp);
				fclose(fp);
				sdsfree(matrix_input_file_name);
				if (verbose > 0) printf("     ..X\n");

				matrix_input_file_name = sdsdup(matrix_input_base_name);
				matrix_input_file_name = sdscat(matrix_input_file_name, ".B");
				fp=fopen(matrix_input_file_name,"rb");
				fread(b_ref,sizeof(double),n,fp);
				fclose(fp);
				sdsfree(matrix_input_file_name);
				if (verbose > 0) printf("     ..B\n");

				matrix_input_file_name = sdsdup(matrix_input_base_name);
				matrix_input_file_name = sdscat(matrix_input_file_name, ".A");
				fp=fopen(matrix_input_file_name,"rb");
				fread(A_ref,sizeof(double),n*n,fp);
				fclose(fp);
				fp=NULL;
				sdsfree(matrix_input_file_name);
				if (verbose > 0) printf("     ..A\n");
			}
		}
		/*
		 * generated matrices
		 */
		else
		{
			if (mpi_rank==0)
			{
				RandomSquareMatrix1D_cnd(A_ref, n, seed, cnd);
				read_cnd = round(ConditionNumber1D(A_ref, n, n));
				if (read_cnd!=cnd && verbose>0)
				{
					printf("WRN: Condition number (%d) differs from read back (%d)\n",cnd,read_cnd);
				}
				FillVector(x_ref, n, 1);
				dgemm_(&transA, &transx, &n, &m, &n, &d1, A_ref, &n, x_ref, &n, &d0, b_ref, &n);
			}
		}
		sdsfree(matrix_input_base_name);

		routine_input.A_ref = A_ref;
		routine_input.x_ref = x_ref;
		routine_input.b_ref = b_ref;
		/*
		{
				n,
				A_ref,
				x_ref,
				b_ref,
				nrhs,
				cprocs,
				sprocs,
				ime_nb,
				scalapack_nb,
				calc_nre
		};
		*/

	/*
	 * ********
	 * commands
	 * ********
	 */
	if ( strcmp(command, "--help" ) == 0 )			// print help summary
	{
		if (mpi_rank==0)
		{
			printf("Usage: tester [OPTION] command [TEST ROUTINE(S)]\n\n");
			printf("Commands are:\n");
			printf("  --help \t\t print this help\n");
			printf("  --list \t\t print the list of testable routines\n");
			printf("  --save \t\t save the generated matrices to files\n");
			printf("  --run \t\t run the test(s)\n");
			printf("        \t\t tests are specified by a space-separated list of testable routines\n");
			printf("\n");
			printf("Options are:\n");
			printf("..to be written..\n\n");
			// TODO: help summary
		}
		MAIN_CLEANUP(mpi_rank);
		MPI_Finalize();
		return 0;
	}
	else if ( strcmp(command, "--list" ) == 0 )		// print list of testable functions
	{
		if (mpi_rank==0)
		{
			printf("Testable routines:\n");
			for( j = 0; j < versions_all; j++ )
			{
				printf("     %2d %s\n", j+1 ,versionname_all[j]);
			}
			printf("\n");
		}
		MAIN_CLEANUP(mpi_rank);
		MPI_Finalize();
		return 0;
	}
	else if ( strcmp(command, "--save" ) == 0 )		// save generated matrices
	{
		if (matrix_output_base_name == NULL)
		{
			if (mpi_rank==0)
			{
				printf("ERR: Please specify a base file path to save to\n\n");
			}
			MAIN_CLEANUP(mpi_rank);
			MPI_Finalize();
			return 1;
		}
		else
		{
			if (mpi_rank == 0)
			{
				if (verbose > 0) printf("     Saving reference matrices to binary files..\n");

				matrix_output_file_name = sdsdup(matrix_output_base_name);
				matrix_output_file_name = sdscat(matrix_output_file_name, ".X");
				fp=fopen(matrix_output_file_name,"wb");
				fwrite(x_ref,sizeof(double),n,fp);
				fclose(fp);
				sdsfree(matrix_output_file_name);
				if (verbose > 0) printf("     ..X\n");

				matrix_output_file_name = sdsdup(matrix_output_base_name);
				matrix_output_file_name = sdscat(matrix_output_file_name, ".B");
				fp=fopen(matrix_output_file_name,"wb");
				fwrite(b_ref,sizeof(double),n,fp);
				fclose(fp);
				sdsfree(matrix_output_file_name);
				if (verbose > 0) printf("     ..B\n");

				matrix_output_file_name = sdsdup(matrix_output_base_name);
				matrix_output_file_name = sdscat(matrix_output_file_name, ".A");
				fp=fopen(matrix_output_file_name,"wb");
				fwrite(A_ref,sizeof(double),n*n,fp);
				fclose(fp);
				fp=NULL;
				sdsfree(matrix_output_file_name);
				if (verbose > 0) printf("     ..A\n");
			}
			sdsfree(matrix_output_base_name);

			MAIN_CLEANUP(mpi_rank);
			MPI_Finalize();
			return 0;
		}
	}
	else if ( strcmp(command, "null" ) == 0 )		// no command given
	{
		if (mpi_rank==0)
		{
			printf("ERR: Please specify command: --help|--list|--run|--save\n\n");
		}
		MAIN_CLEANUP(mpi_rank);
		MPI_Finalize();
		return 1;
	}
	else if ( strcmp(command, "--run" ) == 0 )		// run tests
	{
		if (versions_selected == 0)
		{
			if (mpi_rank==0)
			{
				printf("ERR: Please specify at least one test routine\n\n");
			}
			MAIN_CLEANUP(mpi_rank);
			MPI_Finalize();
			return 1;
		}
		else
		{
		/*
		 * *********
		 * run tests
		 * *********
		 */

			/*
			 * print initial summary to file
			 */
			#define fpinfo(string_label,integer_info) if (fp!=NULL && mpi_rank==0) {fprintf(fp,"info,%s,%d\n",string_label,integer_info); \
			}
			#define fpdata(track_num) if (fp!=NULL && mpi_rank==0) {																		\
				fprintf(fp,"data,%s,%d,%d,%.0f,%f\n",versionname_all[track_num], rep+1,	versionrun[track_num][rep].exit_code,				\
																					versionrun[track_num][rep].total_time,					\
																					versionrun[track_num][rep].norm_rel_err);				\
				fprintf(fp,"data,%s%s,%d,%d,%.0f,%f\n",versionname_all[track_num],"(core)", rep+1,	versionrun[track_num][rep].exit_code,	\
																								versionrun[track_num][rep].core_time,		\
																								versionrun[track_num][rep].norm_rel_err);	\
			}

			// get time now
			time_t rawtime;
			struct tm *readtime;
			time ( &rawtime );
			readtime = localtime ( &rawtime );

			// dump summary
			if (output_to_file && mpi_rank == 0)
			{
				fp=fopen(test_output_file_name,"w");

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
				fpinfo("number of OMP threads",omp_threads);
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

			/*
			 * init tests
			 */

			// init totals to 0 before accumulation
			for (i=0; i<versions_selected; i++)
			{
				versiontot[i].total_time   = 0;
				versiontot[i].core_time    = 0;
				versiontot[i].norm_rel_err = 0;
				versiontot[i].exit_code    = -1;
			}

			// init communication channels
			test_dummy(versionname_all[0], verbose, routine_input, mpi_rank, MPI_COMM_WORLD);

			/*
			 * main loop
			 */

			// runs are repeated
			for (rep=0; rep<repetitions; rep++)
			{
				if (mpi_rank==0 && verbose>0) {printf("\n Run #%d:\n",rep+1);}

				// every run calls some selected routines
				for (i=0; i<versions_selected; i++)
				{
					versionrun[i][rep]=tester_routine(0, versionname_selected[i], verbose, routine_env, routine_input, mpi_rank, failing_rank, failing_level, checkpoint_skip_interval);
				}

				if (mpi_rank==0)
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
			if (mpi_rank==0)
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
					if (output_to_file)
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
					if (output_to_file)
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
		}
	}


	/*
	 * *******
	 * cleanup
	 * *******
	 */
	// memory and pointers
	/*
	if (mpi_rank==0)
	{
		if (output_to_file || input_from_file)
		{
			if (fp != NULL) fclose(fp);
		}
		DeallocateMatrix1D(A_ref);
		DeallocateVector(b_ref);
		DeallocateVector(x_ref);
	}
	sdsfree(test_output_file_name);
	*/
	MAIN_CLEANUP(mpi_rank);

	// BLACS
	// TODO: cleanup BLACS
	//Close BLACS environment
	//Cblacs_barrier( context, "All");	// not working! why?
	//MPI_Barrier(MPI_COMM_WORLD);		// working but not needed
	//Cblacs_gridexit( context );		// not needed if calling blacs_exit
	//Cblacs_gridexit( context_global );// not needed if calling blacs_exit
	//Cblacs_exit( one );				// argument not 0: it is assumed the user will continue using the machine after the BLACS are done
										// error, if main function called more tha once, why?
	//Cblacs_barrier( context_all, &order_all);

	// MPI
	MPI_Finalize();
	return 0;
}

