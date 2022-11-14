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
	while (i>=1) // skips the first ("dummy" not callable)
	{
		if( strcmp( all[i], selected ) == 0 )
		{
			break;
		}
		i--;
	}
	return i-1; // returns -1 if no match
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

    int nmat;				// matrix size
    int nrhs;				// num. of r.h.s
    int cnd;				// condition number
    int seed;				// seed for random matrix generation
    int rows;
    int cols;
    int scalapack_iter;		// scalapack total iterations
    int scalapack_nb;		// scalapack blocking factor
    int ime_nb;				// ime blocking factor

    int repetitions;
    char verbose;
    int spare_procs;		// number of processes to allocate for summing (0 = no fault tolerance)
    int calc_procs;			// number of processes for real IMe calc
	int faulty_procs;
	int fault_tolerance;
    int failing_rank;
    int failing_level;
    int failing_level_override;
    int scalapack_checkpoint_interval;
    int scalapack_failing_level;

    sds matrix_gen_type;
    sds matrix_output_base_name;
    sds matrix_output_file_name;
    sds matrix_input_base_name;
    sds matrix_input_file_name;
    sds test_output_file_name;

    char get_nre;
    char set_cnd;
    char get_cnd;
    int output_to_file;
    int input_from_file;

	FILE* fp;

    char* command;

	int cnd_readback;
	double* A_ref;
	double* x_ref;
	double* b_ref;

    /*
     * ************
     * input values
     * ************
     */

		/*
		 * FIXED defaults
		 */
    	matrix_gen_type         = sdsnew("par"); // matrix generation type ( sequential (seq) /  parallel (par) )
		matrix_output_base_name = NULL;
		matrix_output_file_name = NULL;
		matrix_input_base_name  = NULL;
		matrix_input_file_name  = NULL;
		test_output_file_name   = NULL;
		fp = NULL;

		output_to_file			 = 0;		// no output to file
		input_from_file			 = 0;		// no input from file
		command					 = "null";
		nmat					 = -1;		// size (rank) of the input matrix;
		nrhs					 = 1;
		verbose					 = 1;		// minimum output verbosity (0 = none)
		repetitions				 = 1;		// how many calls for each routine
		spare_procs				 = 0;		// no recovery possible in case of fault
		faulty_procs			 = 0;		// no faults will occur
		fault_tolerance			 = -1;		// level of fault tolerance requested (0 = none, -1 = unset)
		failing_rank			 = 2;		// process 2 will fail
		failing_level			 = -1;		// failing level
		failing_level_override	 = -1;		// failing level automatically set
		scalapack_checkpoint_interval = -1;		// -1=never, otherwise do at every (scalapack_checkpoint_interval+1) iteration
		cnd						 = 1;		// condition number for randomly generated matrices
		cnd_readback			 = -1;		// condition number read back from generated or file matrices (-1=value not read back)
		set_cnd					 = 1;		// pre-conditioning (1) of the matrix  or not (0)
		get_cnd					 = 1;		// read back (1) cnd from generated matrix or not (0)
		seed					 = 1;		// seed for random generation
		get_nre					 = 1;		// calc (1) normwise relative error or not (0)
		scalapack_nb			 = SCALAPACKNB;	// scalapack blocking factor, default defined in header
		scalapack_failing_level  = -1;
		ime_nb					 = 1;			// ime blocking factor

		/*
		 * list of testable routines (see tester_labels.h)
		 */

		// TODO: add a description for every routine

		versions_all = 0;
		versionname_all[versions_all++] = "dummy"; // not actually callable

		versionname_all[versions_all++] = IME_DEV;
		// default version removed

		/*
		 * obsolete versions, removed
		 */
		/*
		versionname_all[versions_all++] = IME_PV_SV_WO;
		versionname_all[versions_all++] = IME_PV_SV_CO;

		versionname_all[versions_all++] = IME_PV_SV_WO_OAE;
		versionname_all[versions_all++] = IME_PV_SV_WO_OA;
		versionname_all[versions_all++] = IME_PV_SV_WO_OGE;
		versionname_all[versions_all++] = IME_PV_SV_WO_OG;

		versionname_all[versions_all++] = IME_PV_SV_WO_U1AE;
		versionname_all[versions_all++] = IME_PV_SV_WO_U1A;
		versionname_all[versions_all++] = IME_PV_SV_WO_U1GE;
		versionname_all[versions_all++] = IME_PV_SV_WO_U1G;

		versionname_all[versions_all++] = IME_PV_SV_WO_U2AE;
		versionname_all[versions_all++] = IME_PV_SV_WO_U2A;
		versionname_all[versions_all++] = IME_PV_SV_WO_U2GE;
		versionname_all[versions_all++] = IME_PV_SV_WO_U2G;

		versionname_all[versions_all++] = IME_PV_SV_WO_U3AE;
		versionname_all[versions_all++] = IME_PV_SV_WO_U3A;
		versionname_all[versions_all++] = IME_PV_SV_WO_U3GE;
		versionname_all[versions_all++] = IME_PV_SV_WO_U3G;

		versionname_all[versions_all++] = IME_PV_SV_ICO_G;
		versionname_all[versions_all++] = IME_PV_SV_CO_G_IND;
		versionname_all[versions_all++] = IME_PV_SV_CO_G_2PASS;
		versionname_all[versions_all++] = IME_PV_SV_CO_G_SMALLER;
		versionname_all[versions_all++] = IME_PV_SV_CO_G_SMALLEST;
		versionname_all[versions_all++] = IME_PV_SV_CO_A_SMALL;
		versionname_all[versions_all++] = IME_PV_SV_CO_A_SMALLER;
		versionname_all[versions_all++] = IME_PV_SV_CO_A_SMALLEST;

		versionname_all[versions_all++] = IME_BLACS_SV_CO_1;
		versionname_all[versions_all++] = IME_BLACS_SV_CO_2;
		*/

		versionname_all[versions_all++] = IME_PB_SV_WO_BF1;
		versionname_all[versions_all++] = IME_PB_SV_CO_BF1;
		versionname_all[versions_all++] = IME_PB_SV_CO_BF1_FAULT_0_TOLERANT_X;
		//versionname_all[versions_all++] = IME_PB_SV_CO_BF1_FAULT_X_TOLERANT_0;
		versionname_all[versions_all++] = IME_PB_SV_CO_BF1_FAULT_X_TOLERANT_X;

		versionname_all[versions_all++] = SPK_SV;
		versionname_all[versions_all++] = SPK_SV_NOPIV;
		versionname_all[versions_all++] = SPK_SV_FAULT_0_TOLERANT_1_CP;
		versionname_all[versions_all++] = SPK_SV_FAULT_1_TOLERANT_1_CP;
		//versionname_all[versions_all++] = SPK_SV_FAULT_0_TOLERANT_1_CS;
		versionname_all[versions_all++] = SPK_SV_FAULT_0_TOLERANT_X_CP;

		versionname_all[versions_all++] = SPK_LU;
		versionname_all[versions_all++] = SPK_LU_FAULT_0_TOLERANT_1;
		versionname_all[versions_all++] = SPK_LU_FAULT_1_TOLERANT_1;

		versionname_all[versions_all++] = SPK_QR;
		versionname_all[versions_all++] = SPK_QR_FAULT_0_TOLERANT_1;
		versionname_all[versions_all++] = SPK_QR_FAULT_1_TOLERANT_1;

		versionname_all[versions_all++] = FTLA_LU_FAULT_0_TOLERANT_1;
		versionname_all[versions_all++] = FTLA_LU_FAULT_1_TOLERANT_1;
		versionname_all[versions_all++] = FTLA_LU_FAULT_X_TOLERANT_X;

		versionname_all[versions_all++] = FTLA_LU_SV_FAULT_0_TOLERANT_1;
		versionname_all[versions_all++] = FTLA_LU_SV_FAULT_1_TOLERANT_1;
		versionname_all[versions_all++] = FTLA_LU_SV_FAULT_X_TOLERANT_X;

		versionname_all[versions_all++] = FTLA_QR_FAULT_0_TOLERANT_1;
		versionname_all[versions_all++] = FTLA_QR_FAULT_1_TOLERANT_1;
		versionname_all[versions_all++] = FTLA_QR_FAULT_X_TOLERANT_X;

		versionname_all[versions_all++] = FTLA_QR_SV_FAULT_0_TOLERANT_1;
		versionname_all[versions_all++] = FTLA_QR_SV_FAULT_1_TOLERANT_1;
		versionname_all[versions_all++] = FTLA_QR_SV_FAULT_X_TOLERANT_X;

		versions_selected = 0;

		/*
		 * read command line parameters
		 */
		// TODO: checking for valid arguments
		for( i = 1; i < argc; i++ )
		{
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
			if( strcmp( argv[i], "-v" ) == 0 ) {
				verbose = atoi(argv[i+1]);
				i++;
			}
			if( strcmp( argv[i], "-r" ) == 0 ) {
				repetitions = atoi(argv[i+1]);
				i++;
			}
			if( strcmp( argv[i], "-seed" ) == 0 ) {
				seed = atoi(argv[i+1]);
				i++;
			}
			if( strcmp( argv[i], "-cnd" ) == 0 ) {
				cnd = atoi(argv[i+1]);
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
			if( strcmp( argv[i], "-nm" ) == 0 ) {
				nmat = atoi(argv[i+1]);
				i++;
			}
			if( strcmp( argv[i], "-nrhs" ) == 0 ) {
				nrhs = atoi(argv[i+1]);
				i++;
			}

			// fault tolerance
			if( strcmp( argv[i], "-ft" ) == 0 ) {
				fault_tolerance = atoi(argv[i+1]);
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
			// number of spare procs
			if( strcmp( argv[i], "-nps" ) == 0 ) {
				spare_procs = atoi(argv[i+1]);
				i++;
			}
			// number of procs faulted
			if( strcmp( argv[i], "-npf" ) == 0 ) {
				faulty_procs = atoi(argv[i+1]);
				i++;
			}

			if( strcmp( argv[i], "-cp" ) == 0 ) {
				scalapack_checkpoint_interval = atoi(argv[i+1]);
				i++;
			}

			if( strcmp( argv[i], "-no-cnd-readback" ) == 0 ) {
				get_cnd = 0;
				//i++; // no parameter value, no inc
			}
			if( strcmp( argv[i], "-no-cnd-set" ) == 0 ) {
				set_cnd = 0;
				//i++; // no parameter value, no inc
			}
			if( strcmp( argv[i], "-no-nre-readback" ) == 0 ) {
				get_nre = 0;
				//i++; // no parameter value, no inc
			}
			if( strcmp( argv[i], "-mat-gen" ) == 0 ) {
				matrix_gen_type = sdscpy(matrix_gen_type,argv[i+1]);
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
		rows=nmat;
		cols=nmat;
		scalapack_iter=(int)ceil(rows/scalapack_nb);

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

		calc_procs=mpi_procs-spare_procs;	// number of MPI processes for real calc

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

		if (mpi_rank==0 && verbose>0)
		{
			printf("\n IMe test suite\n================\n");
		}

		/*
		 *******************
		 * preliminary check
		 *******************
		 */
		if ( strcmp(command, "null" ) != 0 && strcmp(command, "--help" ) != 0 && strcmp(command, "--list" ) != 0 ) // if informative commands, skip check
		{
			// checking on some input parameters
			if (verbose < 0)
			{
				//if (mpi_rank==0) DISPLAY_WRN("\b","Verbosity level cannot be less than zero: setting to zero");
				verbose=0;
			}
			if (repetitions <= 0)
			{
				if (mpi_rank==0 && verbose>0) DISPLAY_WRN("\b","Repetitions cannot be less than one: setting to one");
				repetitions=1;
			}
			if (cnd < 1)
			{
				if (mpi_rank==0) DISPLAY_ERR("\b","Condition number invalid: must be at least one");
				MPI_Finalize();
				return ERR_INPUT_ARG;
			}
			if (scalapack_nb <= 0)
			{
				if (mpi_rank==0 && verbose>0) DISPLAY_WRN("\b","ScaLAPACK blocking factor cannot be less than one: setting to one");
				scalapack_nb=1;
			}
			if (ime_nb <= 0)
			{
				if (mpi_rank==0 && verbose>0) DISPLAY_WRN("\b","IMe blocking factor cannot be less than one: setting to one");
				ime_nb=1;
			}
			if (nmat <= 0)
			{
				if (mpi_rank==0) DISPLAY_ERR("\b","Input matrix size not set or invalid: must be grater than zero");
				MPI_Finalize();
				return ERR_INPUT_ARG;
			}
			if (nrhs <= 0)
			{
				if (mpi_rank==0) DISPLAY_ERR("\b","Number of r.h.s invalid: must be greater than zero");
				MPI_Finalize();
				return ERR_INPUT_ARG;
			}
			if (faulty_procs < 0)
			{
				if (mpi_rank==0) DISPLAY_ERR("\b","Number of faulty processes invalid: must be at least zero");
				MPI_Finalize();
				return ERR_INPUT_ARG;
			}
			if (spare_procs < 0)
			{
				if (mpi_rank==0) DISPLAY_ERR("\b","Number of spare processes invalid: must be at least zero");
				MPI_Finalize();
				return ERR_INPUT_ARG;
			}
			if (fault_tolerance < 0)
			{
				if (mpi_rank==0) DISPLAY_ERR("\b","Fault tolerance level not set or invalid: must be at least zero");
				MPI_Finalize();
				return ERR_INPUT_ARG;
			}
			else
			{
				if (fault_tolerance == 0)
				{
					failing_level=-1;
					scalapack_failing_level=-1;
					failing_rank=-1;
				}
				else // (fault_tolerance > 0)
				{
					if (failing_level_override < 0 ) // faulty level NOT set on command line
					{
						if (mpi_rank==0) DISPLAY_ERR("\b","Failing level not set: if fault tolerance is enabled a failing level must be given");
						MPI_Finalize();
						return ERR_INPUT_ARG;
					}
					else // faulty level IS set on command line
					{
						if (failing_level < 0 || failing_level >= nmat ) // out of bounds faulty level
						{
							if (mpi_rank==0 && verbose>0) DISPLAY_WRN("\b","Failing level less than zero or greater than greatest level: never failing!");
							failing_level=-1;
							scalapack_failing_level=-1;
							failing_rank=-1;
						}
						else // valid faulty level
						{
							// calc corresponding level for ScaLAPACK
							scalapack_failing_level=(int)ceil((nmat-failing_level)/scalapack_nb);
						}
					}
					if (faulty_procs == 0)
					{
						if (mpi_rank==0 && verbose>0) DISPLAY_WRN("\b","Fault tolerance is set, but no process will be faulty: never failing!");
						failing_level=-1;
						scalapack_failing_level=-1;
						failing_rank=-1;
					}
					else // (faulty_procs > 0)
					{
						if (failing_level == -1)
						{
							if (mpi_rank==0 && verbose>0) DISPLAY_WRN("\b","At least one faulty process is set but not a faulty level: never failing!");
						}
						else
						{
							if (spare_procs == 0 )
							{
								if (mpi_rank==0 && verbose>0) DISPLAY_WRN("\b","No spare processes: faults cannot be recovered!");
							}
						}
					}
				}
			}

			if (!IS_SQUARE(calc_procs))
			{
				if (mpi_rank==0) DISPLAY_ERR("\b","The number of calc. processes has to be square");
				MPI_Finalize();
				return ERR_INPUT_ARG;
			}
		}

		/*
		 * BLACS
		 */
		int ndims = 2, dims[2];
		dims[0] = (int)sqrt(calc_procs); //square grid as close as possible
		dims[1] = 0;
		MPI_Dims_create(calc_procs, ndims, dims);

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

		int* blacs_ctxt_cp;
		int blacs_row;
		int blacs_col;
		int map_cp[1];

		if (spare_procs>0)
		{
			blacs_ctxt_cp=malloc(spare_procs*sizeof(int));
			for (i=0; i<spare_procs; i++) // fault tolerance enabled
			{
				// context for the checkpointing node
				Cblacs_get( ic, i0, &blacs_ctxt_cp[i] );
				map_cp[0]=calc_procs+i;
				Cblacs_gridmap( &blacs_ctxt_cp[i], map_cp, i1, i1, i1);
			}
		}
		else
		{
			//blacs_ctxt_cp = -1;
			blacs_ctxt_cp = NULL;
		}

		// get coords in general grid context
		Cblacs_gridinfo( blacs_ctxt, &blacs_nprow, &blacs_npcol, &blacs_row, &blacs_col );


	/*
	 * ****************************
	 * common input data structures
	 * ****************************
	 */
		parallel_env routine_env = {
			calc_procs,
			spare_procs,
			mpi_rank,
			MPI_COMM_WORLD,
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
			nmat,
			NULL,
			NULL,
			NULL,
			nrhs,
			ime_nb,
			scalapack_nb,
			get_nre
		};

		fault_env routine_fault	= {
			fault_tolerance,
			faulty_procs,
			failing_rank,
			failing_level,
			scalapack_checkpoint_interval
		};

	/*
	 * ******************************
	 * starting summary to video
	 * ******************************
	 */
	if ( strcmp(command, "null" ) != 0 && strcmp(command, "--help" ) != 0 && strcmp(command, "--list" ) != 0) // if informative commands, skip check
	{
		if (mpi_rank==0 && verbose>0)
		{
			printf("     Total processes:               %d\n",np);
			printf("     OMP threads:                   %d\n",omp_threads);
			printf("     MPI ranks:                     %d\n",mpi_procs);
			printf("     BLACS grid:                    %dx%d\n",blacs_nprow,blacs_npcol);
		}
		/*
		 * check grid ranking
		 */
		if (verbose>2)
		{
			if (mpi_rank==0)
			{
				printf("     MPI ranks of calc. processes:\n");
				printf("       Planned:\n");
				for (i=0; i<blacs_nprow; i++)
				{
					printf("         ");
					for (j=0; j<blacs_npcol; j++)
					{
						printf("%d\t",i*blacs_npcol+j);
					}
					printf("\n");
				}
				printf("       Allocated:\n");
				fflush(stdout);
			}

			for (i=0; i<mpi_procs; i++)
			{
				MPI_Barrier(MPI_COMM_WORLD);
				if (mpi_rank==i)
				{
					printf("         %d @ (%d,%d)\n",i,blacs_row,blacs_col);
				}
				MPI_Barrier(MPI_COMM_WORLD);
				fflush(stdout);
			}
		}
		MPI_Barrier(MPI_COMM_WORLD);

		if (mpi_rank==0 && verbose>0)
		{
			printf("     Calculate n.r.e.:              ");
				if (get_nre)	{printf("yes\n");}
				else			{printf("no\n");}
			printf("     IMe iterations:                %d\n",rows);
			printf("     IMe blocking factor:           %d\n",ime_nb);
			printf("     SPK-like iterations:           %d\n",scalapack_iter);
			printf("     SPK-like blocking factor:      %d\n",scalapack_nb);

			printf("     Fault tolerance:               ");
			if (fault_tolerance > 0)
			{
				printf("enabled = %d\n",fault_tolerance);
				printf("       Calc. processes:             %d\n",calc_procs);
				printf("       Faulty processes:            %d\n",faulty_procs);
				printf("       Spare processes:             %d\n",spare_procs);
				printf("     IMe failing rank:              %d\n",failing_rank);
				printf("     IMe failing level:             ");
					if (failing_level<0) {printf("never = ");}
					printf("%d\n",failing_level);
				printf("     SPK-like failing level:        ");
					if (failing_level<0) {printf("never = -1\n");}
					else {printf("%d\n",nmat-failing_level);}
				printf("     SPK-like failing iteration:    ");
					if (failing_level<0) {printf("never = ");}
					printf("%d\n",scalapack_failing_level);
				printf("     Checkpoint skip interval:      %d\n",scalapack_checkpoint_interval);

				printf("     Checkpoint freq.:              ");
				if (scalapack_checkpoint_interval<0)
				{
					printf("never\n");
				}
				else
				{
					if (scalapack_checkpoint_interval==0)
					{
						printf("always\n");
					}
					else
					{
						printf("every %d iterations\n",scalapack_checkpoint_interval+1);
					}
				}
			}
			else
			{
				printf("disabled\n");
				printf("       Calc. processes:             %d\n",calc_procs);
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
//	}

//	if ( strcmp(command, "--help" ) != 0 && strcmp(command, "--list" ) != 0) // if informative commands, skip checks and preparation
//	{
		/*
		 * ******
		 * checks
		 * ******
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
				return ERR_ROUTINE_UNKNOWN;
			}
		}

		// check specific condition for every selected routine
		j=0; // error accumulation
		for (i=0; i<versions_selected; i++)
		{
			j = j + (tester_routine(1, versionname_selected[i], verbose, routine_env, routine_input, routine_fault)).exit_code;
		}
		MPI_Bcast(&j,1,MPI_INT,0,MPI_COMM_WORLD);
		if (j != 0)
		{
			MPI_Finalize();
			return ERR_ROUTINE_MISCONFIG;
		}
		// continue if no errors

		/*
		 * **********************
		 * prepare input matrices
		 * **********************
		 */
		if (mpi_rank==0)
		{
			A_ref = AllocateMatrix1D(nmat, nmat);
			x_ref = AllocateVector(nmat);
			b_ref = AllocateVector(nmat);
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
					return ERR_INPUT_ARG;
				}
				if (mpi_rank==0)
				{
					if (verbose > 0)
					{
						printf("     Loading reference matrices from binary files..\n");

						if (!get_cnd)
						{
							printf("WRN: Condition number will not read back from loaded matrix\n");
						}
					}
					matrix_input_file_name = sdsdup(matrix_input_base_name);
					matrix_input_file_name = sdscat(matrix_input_file_name, ".X");
					fp=fopen(matrix_input_file_name,"rb");
					fread(x_ref,sizeof(double),nmat,fp);
					fclose(fp);
					sdsfree(matrix_input_file_name);
					if (verbose > 0) printf("     ..X\n");

					matrix_input_file_name = sdsdup(matrix_input_base_name);
					matrix_input_file_name = sdscat(matrix_input_file_name, ".B");
					fp=fopen(matrix_input_file_name,"rb");
					fread(b_ref,sizeof(double),nmat,fp);
					fclose(fp);
					sdsfree(matrix_input_file_name);
					if (verbose > 0) printf("     ..B\n");

					matrix_input_file_name = sdsdup(matrix_input_base_name);
					matrix_input_file_name = sdscat(matrix_input_file_name, ".A");
					fp=fopen(matrix_input_file_name,"rb");
					fread(A_ref,sizeof(double),nmat*nmat,fp);
					fclose(fp);
					fp=NULL;
					sdsfree(matrix_input_file_name);
					if (verbose > 0) printf("     ..A\n");
				}
				if ( get_cnd )
				{
					if (!IS_SQUARE(calc_procs))
					{
						if (mpi_rank==0)
						{
							printf("ERR: To read back the condition number the process grid must be square\n\n");
						}
						MAIN_CLEANUP(mpi_rank);
						MPI_Finalize();
						return ERR_INPUT_ARG;
					}
					else
					{
						cnd_readback = round( pCheckSystemMatrices1D(nmat, A_ref, x_ref, b_ref, scalapack_nb, mpi_rank, calc_procs, blacs_nprow, blacs_npcol, blacs_row, blacs_col, blacs_ctxt, blacs_ctxt_root) );
					}
				}
			}
			/*
			 * generated matrices
			 */
			else
			{
				if (strcmp(matrix_gen_type, "par" ) == 0)
				{
					// init communication channels (generation uses blacs => mpi interference..)
					test_dummy(versionname_all[0], verbose, routine_env, routine_input);

					if (mpi_rank==0 && verbose>0)
					{
						printf("     Generating random input matrices in parallel with ScaLAPACK\n");
						if (!get_cnd)
						{
							printf("WRN: Condition number will not read back from generated matrix\n");
						}
						if (!set_cnd)
						{
							printf("WRN: Matrix will not be pre-conditioned\n");
						}
					}
					if ( get_cnd && !IS_SQUARE(calc_procs) )
					{
						if (mpi_rank==0)
						{
							printf("ERR: To read back the condition number the process grid must be square\n\n");
						}
						MAIN_CLEANUP(mpi_rank);
						MPI_Finalize();
						return ERR_INPUT_ARG;
					}
					cnd_readback = round( pGenSystemMatrices1D(nmat, A_ref, x_ref, b_ref, seed, cnd, set_cnd, get_cnd, scalapack_nb, mpi_rank, calc_procs, blacs_nprow, blacs_npcol, blacs_row, blacs_col, blacs_ctxt, blacs_ctxt_root) );
				}
				else if (strcmp(matrix_gen_type, "seq" ) == 0)
				{
					if (mpi_rank==0 && verbose>0)
					{
						printf("     Generating random input matrices sequentially with LAPACK\n");
						if (!get_cnd)
						{
							printf("WRN: Condition number will not read back from generated matrix\n");
						}
						cnd_readback = round( GenSystemMatrices1D(nmat, A_ref, x_ref, b_ref, seed, cnd, get_cnd) );
					}
				}
				else
				{
					if (mpi_rank==0)
					{
						printf("ERR: Unknown type of matrix generation %s\n\n",matrix_gen_type);
					}
					MAIN_CLEANUP(mpi_rank);
					MPI_Finalize();
					return ERR_INPUT_ARG;
				}
			}
			sdsfree(matrix_input_base_name);

			// matrix for debugging purposes
			/*
			if (mpi_rank==0)
			{
				FillMatrix1D(A_ref, n, n);
				OneMatrix1D(b_ref, n, 1);
			}
			*/

			routine_input.A_ref = A_ref;
			routine_input.x_ref = x_ref;
			routine_input.b_ref = b_ref;

			/*
			 * ******************************
			 * continuing summary to video
			 * ******************************
			 */
			if (mpi_rank==0 && verbose>0)
			{
				printf("     Matrix random generation seed: %d\n",seed);
				printf("     Matrix size:                   %dx%d\n",rows,cols);
				printf("     Matrix condition number set:   %d\n",cnd);
				printf("     Matrix condition number got:   %d\n",cnd_readback);
				if (cnd_readback!=cnd)
				{
					printf("WRN: Condition number (%d) differs from read back (%d)\n",cnd,cnd_readback);
				}
			}
	}
	else
	{
		A_ref = NULL;
		x_ref = NULL;
		b_ref = NULL;
	}


	/*
	 * ********
	 * commands
	 * ********
	 */
	if ( strcmp(command, "--help" ) == 0 )			// print help summary
	{
		if (mpi_rank==0)
		{
			printf("\n");
			printf("Usage: tester [OPTION] command [TEST ROUTINE(S)]\n\n");
			printf("Commands are:\n");
			printf("  --help \t\t print this help\n");
			printf("  --list \t\t print the list of testable routines\n");
			printf("  --save \t\t save the generated matrices to files\n");
			printf("  --run \t\t run the test(s)\n");
			printf("        \t\t tests are specified by a space-separated list of testable routines\n");
			printf("\n");
			printf("Options are:\n");
			printf("  -v    <integer number>     : verbosity level [0-3] (0=quiet)\n" );
			printf("  -nm   <integer number>     : input matrix rank\n" );
			printf("  -nrhs <integer number>     : r.h.s. columns\n" );
			printf("  -seed <integer number>     : seed of the random generation\n" );
			printf("  -cnd  <integer number>     : condition number of the input matrix\n" );
			printf("  -no-cnd-set                : disable matrix pre-conditioning\n" );
			printf("  -no-cnd-readback           : disable condition number checking after generation\n" );
			printf("  -no-nre-readback           : disable normwise relative error checking\n" );
			printf("  -mat-gen <string>          : type of the random generation [par|ser] (parallel or sequential)\n" );
			printf("  -r    <integer number>     : run repetitions\n" );
			printf("  -o    <file path>          : output to CSV file\n" );
			printf("  -i    <file path>          : input matrices base name file path (.A, .X, .B auto appended)\n" );
			printf("  -ft   <integer number>     : fault-tolerance level [0-..] (0=none)\n" );
			printf("  -fr   <integer number>     : simulated faulty mpi rank\n" );
			printf("  -fl   <integer number>     : simulated faulty IMe inhibition level\n" );
			printf("  -npf  <integer number>     : number of simulated faults [0-..] (0=none)\n" );
			printf("  -nps  <integer number>     : number of spare processes [0-..] (0=none)\n" );
			printf("  -cp   <integer number>     : checkpointing interval\n" );
			printf("  -spk-nb <integer number>   : ScaLAPACK blocking factor\n" );
			printf("  -ime-nb <integer number>   : IMe blocking factor\n\n" );
		}
		MAIN_CLEANUP(mpi_rank);
		MPI_Finalize();
		return 0;
	}
	else if ( strcmp(command, "--list" ) == 0 )		// print list of testable functions
	{
		if (mpi_rank==0)
		{
			printf("\nTestable routines:\n");
			for( j = 1; j < versions_all; j++ )
			{
				printf("     %2d %s\n", j ,versionname_all[j]);
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
			return ERR_INPUT_ARG;
		}
		else
		{
			if (mpi_rank == 0)
			{
				if (verbose > 0) printf("     Saving reference matrices to binary files..\n");

				matrix_output_file_name = sdsdup(matrix_output_base_name);
				matrix_output_file_name = sdscat(matrix_output_file_name, ".X");
				fp=fopen(matrix_output_file_name,"wb");
				fwrite(x_ref,sizeof(double),nmat,fp);
				fclose(fp);
				sdsfree(matrix_output_file_name);
				if (verbose > 0) printf("     ..X\n");

				matrix_output_file_name = sdsdup(matrix_output_base_name);
				matrix_output_file_name = sdscat(matrix_output_file_name, ".B");
				fp=fopen(matrix_output_file_name,"wb");
				fwrite(b_ref,sizeof(double),nmat,fp);
				fclose(fp);
				sdsfree(matrix_output_file_name);
				if (verbose > 0) printf("     ..B\n");

				matrix_output_file_name = sdsdup(matrix_output_base_name);
				matrix_output_file_name = sdscat(matrix_output_file_name, ".A");
				fp=fopen(matrix_output_file_name,"wb");
				fwrite(A_ref,sizeof(double),nmat*nmat,fp);
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
		return ERR_INPUT_ARG;
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
			return ERR_INPUT_ARG;
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
				fprintf(fp,"data,%s,%d,%d,%.0f,%.17f\n",versionname_all[track_num], rep+1,	versionrun[track_num][rep].exit_code,				\
																					versionrun[track_num][rep].total_time,					\
																					versionrun[track_num][rep].norm_rel_err);				\
				fprintf(fp,"data,%s%s,%d,%d,%.0f,%.17f\n",versionname_all[track_num],"(core)", rep+1,	versionrun[track_num][rep].exit_code,	\
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
				fpinfo("number of MPI ranks",mpi_procs);
				fpinfo("number of OMP threads",omp_threads);
				fpinfo("number of BLACS rows",blacs_nprow);
				fpinfo("number of BLACS columns",blacs_npcol);
				fpinfo("number of processes",np);
				fpinfo("number of calc processes",calc_procs);
				fpinfo("number of spare processes",spare_procs);
				fpinfo("number of faulty processes",faulty_procs);
				fpinfo("fault tolerance",fault_tolerance);
				fpinfo("failing rank",failing_rank);
				fpinfo("failing level",failing_level);
				fpinfo("checkpoint skip interval",scalapack_checkpoint_interval);
				fpinfo("seed",seed);
				fpinfo("condition number set",cnd);
				fpinfo("condition number got",cnd_readback);
				fpinfo("matrix size",nmat);
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

			// init communication channels if not done during matrix generation
			if (strcmp(matrix_gen_type, "par" ) != 0) test_dummy(versionname_all[0], verbose, routine_env, routine_input);

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
					versionrun[i][rep]=tester_routine(0, versionname_selected[i], verbose, routine_env, routine_input, routine_fault);

					if (output_to_file)
					{
						if (mpi_rank==0)
						{
							fprintf(fp,"data,%s,%d,%d,%.0f,%.17f\n",	versionname_selected[i], rep+1,				\
																		versionrun[i][rep].exit_code,				\
																		versionrun[i][rep].total_time,				\
																		versionrun[i][rep].norm_rel_err				\
							);
							fprintf(fp,"data,%s%s,%d,%d,%.0f,%.17f\n",	versionname_selected[i], "(core)", rep+1,	\
																		versionrun[i][rep].exit_code,				\
																		versionrun[i][rep].core_time,				\
																		versionrun[i][rep].norm_rel_err				\
							);
						}
					}
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
							printf("%-30s    Call    run time: %10.0f (%.0f)\ts\t nre: %.17f\n",	versionname_selected[i],		\
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
			if (mpi_rank==0 && verbose>0)
			{
				printf("\n Summary:\n");

				// total
				for (i=0; i<versions_selected; i++)
				{
					printf("%-30s    Total   run time: %10.0f (%.0f)\ts\n",	versionname_selected[i],	\
																			versiontot[i].total_time,	\
																			versiontot[i].core_time		\
					);
				}
				printf("\n");

				// average
				for (i=0; i<versions_selected; i++)
				{
					printf("%-30s    Average run time: %10.0f (%.0f)\ts\t nre: %.17f\n",	versionname_selected[i],				\
																							versiontot[i].total_time/repetitions,	\
																							versiontot[i].core_time/repetitions,	\
																							versiontot[i].norm_rel_err/repetitions	\
					);
					if (output_to_file)
					{
						fprintf(fp,"data,%s,%d,%d,%.0f,%.17f\n",	versionname_selected[i], 0,				\
																	versionrun[i][0].exit_code,				\
																	versiontot[i].total_time/repetitions,	\
																	versiontot[i].norm_rel_err/repetitions	\
						);
						fprintf(fp,"data,%s%s,%d,%d,%.0f,%.17f\n",	versionname_selected[i], "(core)", 0,	\
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
					printf("%-30s    Median  run time: %10.0f (%.0f)\ts\t nre: %.17f\n",	versionname_selected[i],					\
																							versionrun[i][repetitions/2].total_time,	\
																							versionrun[i][repetitions/2].core_time,		\
																							versionrun[i][repetitions/2].norm_rel_err	\
					);
					if (output_to_file)
					{
						fprintf(fp,"data,%s,%d,%d,%.0f,%.17f\n",	versionname_selected[i], -1,				\
																	versionrun[i][0].exit_code,					\
																	versionrun[i][repetitions/2].total_time,	\
																	versionrun[i][repetitions/2].norm_rel_err	\
						);
						fprintf(fp,"data,%s%s,%d,%d,%.0f,%.17f\n",	versionname_selected[i], "(core)",-1,		\
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
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
	return 0;
}

