/*
 * tester.c
 *
 *  Created on: Apr 25, 2020
 *      Author: marcello
 */

#include <mpi.h>
#include <omp.h>
#include <sched.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

#include "helpers/macros.h"
#include "helpers/lapack.h"
#include "helpers/scalapack.h"
#include "helpers/Cblacs.h"
#include "helpers/simple_dynamic_strings/sds.h"
#include "helpers/matrix_advanced.h"
#include "tester_labels.h"
#include "tester_routine.h"
#include "tester_structures.h"
#include "testers/dummy/test_init.h"
#include "testers/IMe/lib/src/helpers/matrix_basic.h"
#include "testers/IMe/lib/src/helpers/vector_basic.h"


#define MAX_VERSIONS 50
#define MAX_RUNS 10

// TODO: rename auxiliary functions

#ifdef CRESCO_POWERCAP
/*
#define ZONE0 "/home/marcello/rapl/z0/energy_uj"
#define ZONE0SUB0 "/home/marcello/rapl/z0/s0/energy_uj"
#define ZONE1 "/home/marcello/rapl/z1/energy_uj"
#define ZONE1SUB0 "/home/marcello/rapl/z1/s0/energy_uj"
*/
#define ZONE0 "/sys/class/powercap/intel-rapl:0/energy_uj"
#define ZONE0SUB0 "/sys/class/powercap/intel-rapl:0:0/energy_uj"
#define ZONE1 "/sys/class/powercap/intel-rapl:1/energy_uj"
#define ZONE1SUB0 "/sys/class/powercap/intel-rapl:1:0/energy_uj"
#define POWERCAP_ENTRIES 4

unsigned long long int read_energy(char* pathname)
{
	//printf("** %s\n",pathname);
	unsigned long long int value;
	//open and get the file handle
	FILE* fh;
	fh = fopen(pathname, "r");
	//read
	if (!fscanf(fh, "%llu\n", &value)) printf("ERR> %s not read",pathname);
	fclose(fh);
	return(value);
}
#endif

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

	const int real_double = 2;
	const int real_single = 1;

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
    int type;				// precision type
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

    sds matrix_gen_method;
    sds matrix_precision_type;
    sds matrix_output_base_name;
    sds matrix_output_file_name;
    sds matrix_input_base_name;
    sds matrix_input_file_name;
    sds test_output_file_name;

#ifdef CRESCO_POWERCAP
    char get_energy;
    unsigned long long int energy[POWERCAP_ENTRIES];
    unsigned long long int * energy_readings_begin;
    unsigned long long int * energy_readings_end;
#endif
    char get_nre;
    char set_cnd;
    char get_cnd;
    int output_to_file;
    int input_from_file;

	FILE* fp;

    char* command;

	int cnd_readback;

	// double precision
	double* A_ref_d = NULL;
	double* x_ref_d = NULL;
	double* b_ref_d = NULL;

	// single precision
	float* A_ref_s = NULL;
	float* x_ref_s = NULL;
	float* b_ref_s = NULL;


    /*
     * ************
     * input values
     * ************
     */

		/*
		 * FIXED defaults
		 */
    	matrix_gen_method       = sdsnew("par"); // matrix generation type ( sequential (seq) /  parallel (par) )
    	matrix_precision_type   = sdsnew("");
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
		failing_rank			 = 2;		// in case of fault, process 2 will fail
		failing_level			 = -1;		// failing level (-1 = none, never failing)
		failing_level_override	 = -1;		// failing level automatically set
		scalapack_checkpoint_interval = -1;		// -1=never, otherwise do at every (scalapack_checkpoint_interval+1) iteration
		cnd						 = 1;		// condition number for randomly generated matrices
		cnd_readback			 = -1;		// condition number read back from generated or file matrices (-1=value not read back)
		set_cnd					 = 1;		// pre-conditioning (1) of the matrix  or not (0)
		get_cnd					 = 1;		// read back (1) cnd from generated matrix or not (0)
#ifdef CRESCO_POWERCAP
		get_energy				 = 0;		// no energy reading
#endif
		seed					 = 1;		// seed for random generation
		type					 = 0;		// precision type (0=unspecified)
		get_nre					 = 1;		// calc (1) normwise relative error or not (0)
		scalapack_nb			 = SCALAPACKNB;	// scalapack blocking factor, default defined in header
		scalapack_failing_level  = -1;
		ime_nb					 = 1;			// ime blocking factor

		/*
		 * list of testable routines (see tester_labels.h)
		 */

		// TODO: add a description for every routine

		versions_all = 0;
		versionname_all[versions_all++] = TEST_INIT; // not actually callable, used as a workaround to init mpi/blacs
		versionname_all[versions_all++] = DUMMY;     // no actual test (empty routine), but actually callable

		versionname_all[versions_all++] = IME_PSGESV_WO;
		versionname_all[versions_all++] = IME_PDGESV_WO;

		versionname_all[versions_all++] = IME_PSGESV_CO;
		versionname_all[versions_all++] = IME_PDGESV_CO;

		versionname_all[versions_all++] = IME_PSGESV_CO_FT;
		versionname_all[versions_all++] = IME_PDGESV_CO_FT;

		versionname_all[versions_all++] = SPK_PSGESV;
		versionname_all[versions_all++] = SPK_PDGESV;

		versionname_all[versions_all++] = SPK_PSGESV_FT_CP;
		versionname_all[versions_all++] = SPK_PDGESV_FT_CP;

		versionname_all[versions_all++] = SPK_PDGESV_NOPIV;

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
			if( strcmp( argv[i], "-type" ) == 0 ) {
				matrix_precision_type = sdsnew(argv[i+1]);
				i++;
			}
			if( strcmp( argv[i], "-cnd" ) == 0 ) {
				cnd = atoi(argv[i+1]);
				i++;
			}
			if( strcmp( argv[i], "-spk-cp" ) == 0 ) {
				scalapack_checkpoint_interval = atoi(argv[i+1]);
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
#ifdef CRESCO_POWERCAP
			if( strcmp( argv[i], "-energy-reading" ) == 0 ) {
				get_energy = 1;
				//i++; // no parameter value, no inc
			}
#endif
			if( strcmp( argv[i], "-mat-gen" ) == 0 ) {
				matrix_gen_method = sdscpy(matrix_gen_method,argv[i+1]);
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

#ifdef CRESCO_POWERCAP
		int namelen;
	    char mpi_procname[MPI_MAX_PROCESSOR_NAME];
	    char** mpi_hostnames;
		MPI_Comm MPI_COMM_HOST;
		int rank_in_host;
		int size_in_host;
		MPI_Comm MPI_COMM_ROOTS;
		int rank_in_roots;
		int size_in_roots;

		int root_in_host;

		MPI_Comm_split_type(MPI_COMM_WORLD, OMPI_COMM_TYPE_HOST, mpi_rank, MPI_INFO_NULL, &MPI_COMM_HOST);
		MPI_Comm_rank(MPI_COMM_HOST, &rank_in_host);
		MPI_Comm_size(MPI_COMM_HOST, &size_in_host);

		if (rank_in_host==0)
		{
			root_in_host=1;
		}
		else
		{
			root_in_host=0;
		}

		MPI_Comm_split(MPI_COMM_WORLD, root_in_host, mpi_rank, &MPI_COMM_ROOTS);
		MPI_Comm_rank(MPI_COMM_ROOTS, &rank_in_roots);
		MPI_Comm_size(MPI_COMM_ROOTS, &size_in_roots);

		if (mpi_rank==0)
		{
			energy_readings_begin=malloc(size_in_roots * POWERCAP_ENTRIES * sizeof(unsigned long long int));
			energy_readings_end=malloc(size_in_roots * POWERCAP_ENTRIES * sizeof(unsigned long long int));
			mpi_hostnames=malloc(size_in_roots * sizeof(char*));
			mpi_hostnames[0]=malloc(size_in_roots * (MPI_MAX_PROCESSOR_NAME+1) * sizeof(char));
			for (i=1;i<size_in_roots;i++)
			{
				mpi_hostnames[i]=mpi_hostnames[i-1]+(MPI_MAX_PROCESSOR_NAME+1);
			}
		}
		else
		{
			//TODO: delete?
			mpi_hostnames=malloc(sizeof(char*));
			mpi_hostnames[0]=malloc(sizeof(char));
			energy_readings_begin=malloc(sizeof(char));
			energy_readings_end=malloc(sizeof(char));
		}

		if (get_energy)
		{
			if (rank_in_host==0)
			{
				MPI_Get_processor_name(mpi_procname, &namelen);
				MPI_Gather( mpi_procname, MPI_MAX_PROCESSOR_NAME+1, MPI_CHAR, &mpi_hostnames[0][0], MPI_MAX_PROCESSOR_NAME+1, MPI_CHAR, 0, MPI_COMM_ROOTS );
			}
		}
#endif

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
			// checking on some input parameters for "--run" and "--save"
			if (verbose < 0)
			{
				//if (mpi_rank==0) DISPLAY_WRN("\b","Verbosity level cannot be less than zero: setting to zero");
				verbose=0;
			}
			if (cnd < 1)
			{
				if (mpi_rank==0) DISPLAY_ERR_GLB("Condition number invalid: must be at least one");
				MPI_Finalize();
				return ERR_INPUT_ARG;
			}
			if (nmat <= 0)
			{
				if (mpi_rank==0) DISPLAY_ERR_GLB("Input matrix size not set or invalid: must be grater than zero");
				MPI_Finalize();
				return ERR_INPUT_ARG;
			}
			if (nrhs <= 0)
			{
				if (mpi_rank==0) DISPLAY_ERR_GLB("Number of r.h.s invalid: must be greater than zero");
				MPI_Finalize();
				return ERR_INPUT_ARG;
			}
			if (scalapack_nb <= 0)
			{
				if (mpi_rank==0 && verbose>0) DISPLAY_WRN_GLB("ScaLAPACK blocking factor cannot be less than one: setting to one");
				scalapack_nb=1;
			}

			// checking on some input parameters for "--run" only
			if ( strcmp(command, "--run" ) == 0 )
			{
				if (ime_nb <= 0)
				{
					if (mpi_rank==0 && verbose>0) DISPLAY_WRN_GLB("IMe blocking factor cannot be less than one: setting to one");
					ime_nb=1;
				}
				if (repetitions <= 0)
				{
					if (mpi_rank==0 && verbose>0) DISPLAY_WRN_GLB("Repetitions cannot be less than one: setting to one");
					repetitions=1;
				}
				if (spare_procs < 0)
				{
					if (mpi_rank==0) DISPLAY_ERR_GLB("Number of spare processes invalid: must be at least zero");
					MPI_Finalize();
					return ERR_INPUT_ARG;
				}

				if (faulty_procs < 0) // faulty processes invalid
				{
					if (mpi_rank==0) DISPLAY_ERR_GLB("Number of faulty processes invalid: must be at least zero");
					MPI_Finalize();
					return ERR_INPUT_ARG;
				}
				else  // faulty processes valid
				{
					if (faulty_procs > 0) // faults will occur
					{
						if (failing_level_override < 0 ) // faulty level IS NOT set on command line
						{
							if (mpi_rank==0) DISPLAY_ERR_GLB("Failing level not set: if a process will be faulty a failing level must be given");
							MPI_Finalize();
							return ERR_INPUT_ARG;
						}
						else // faulty level IS set on command line
						{
							if (failing_level < 0 || failing_level >= nmat ) // out of bounds faulty level
							{
								if (mpi_rank==0) DISPLAY_ERR_GLB("Failing level invalid: must be at least zero and not greater than greatest level");
								MPI_Finalize();
								return ERR_INPUT_ARG;
							}
							else // valid faulty level
							{
								// calc corresponding level for ScaLAPACK
								scalapack_failing_level=(int)ceil((nmat-failing_level)/scalapack_nb);
							}
						}
					}
				}

				if (fault_tolerance < 0) // fault tolerance level is invalid
				{
					if (mpi_rank==0) DISPLAY_ERR_GLB("Fault tolerance level not set or invalid: must be at least zero");
					MPI_Finalize();
					return ERR_INPUT_ARG;
				}
				else // fault tolerance level is valid
				{
					if (fault_tolerance == 0) // fault tolerance disabled
					{
						if (faulty_procs > 0)
						{
							if (mpi_rank==0 && verbose>0) DISPLAY_WRN_GLB("Fault tolerance level is not set, but at least one process will be faulty: never recovering");
						}
					}
					else // fault tolerance enabled
					{
						if (faulty_procs < 1) // faults will not occur
						{
							if (mpi_rank==0 && verbose>0) DISPLAY_WRN_GLB("Fault tolerance level is set, but no process will be faulty: never failing");
						}
						else // faults will occur
						{
							if (failing_level == -1)
							{
								if (mpi_rank==0 && verbose>0) DISPLAY_WRN_GLB("At least one faulty process is set but not a faulty level: never failing");
							}
							else
							{
								if (spare_procs == 0 )
								{
									if (mpi_rank==0 && verbose>0) DISPLAY_WRN_GLB("No spare processes: never recovering");
								}
							}
						}
					}
				}

				if (!IS_SQUARE(calc_procs))
				{
					if (mpi_rank==0) DISPLAY_ERR_GLB("The number of calc. processes has to be square");
					MPI_Finalize();
					return ERR_INPUT_ARG;
				}
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
		// action commands, check

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
			printf("     SPK iterations:                %d\n",scalapack_iter);
			printf("     SPK blocking factor:           %d\n",scalapack_nb);

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
				printf("     SPK (IMe-like) failing level:  ");
					if (failing_level<0) {printf("never = -1\n");}
					else {printf("%d\n",nmat-failing_level);}
				printf("     SPK (IMe-like) failing iter.:  ");
					if (failing_level<0) {printf("never = ");}
					printf("%d\n",scalapack_failing_level);
				printf("     Checkpoint skip interval:      %d SPK iterations\n",scalapack_checkpoint_interval);

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
						printf("every %d SPK iterations | %d IMe iterations\n",scalapack_checkpoint_interval+1,(scalapack_checkpoint_interval+1)*scalapack_nb);
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
				printf("WRN> No output to file\n");
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
					printf("ERR> Routine '%s' is unknown\n",versionname_selected[i]);
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
		if (strcmp(matrix_precision_type, "double" ) == 0 || strcmp(matrix_precision_type, "d" ) == 0 )
		{
			type = real_double;
		}
		else if (strcmp(matrix_precision_type, "single" ) == 0 || strcmp(matrix_precision_type, "s" ) == 0 )
		{
			type = real_single;
		}
		else
		{
			if (mpi_rank==0)
			{
				printf("ERR> Unknown or unspecified precision type\n\n");

				if (output_to_file || input_from_file)
				{
					if (fp != NULL) fclose(fp);
				}
			}
			sdsfree(test_output_file_name);
			MPI_Finalize();
			return ERR_INPUT_ARG;
		}

		if (type==real_double)
		{
			if (mpi_rank==0)
			{
				A_ref_d = AllocateMatrix1D_double(nmat, nmat);
				x_ref_d = AllocateVector_double(nmat);
				b_ref_d = AllocateVector_double(nmat);
			}
		}
		else if (type==real_single)
		{
			if (mpi_rank==0)
			{
				A_ref_s = AllocateMatrix1D_float(nmat, nmat);
				x_ref_s = AllocateVector_float(nmat);
				b_ref_s = AllocateVector_float(nmat);
			}
		}


		/*
		 * matrices from file
		 */
			if (input_from_file)
			{
				if ( strcmp(command, "--save" ) == 0 )
				{
					if (mpi_rank == 0) printf("ERR> Cannot take input matrices with '--save' command.\n");
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
							printf("WRN> Condition number will not read back from loaded matrix\n");
						}
					}
					matrix_input_file_name = sdsdup(matrix_input_base_name);
					matrix_input_file_name = sdscat(matrix_input_file_name, ".X");
					fp=fopen(matrix_input_file_name,"rb");
							if (type==real_double) {FREAD(x_ref_d, sizeof(x_ref_d[0]), nmat, fp)}
					else	if (type==real_single) {FREAD(x_ref_s, sizeof(x_ref_s[0]), nmat, fp)}
					fclose(fp);
					sdsfree(matrix_input_file_name);
					if (verbose > 0) printf("     ..X\n");

					matrix_input_file_name = sdsdup(matrix_input_base_name);
					matrix_input_file_name = sdscat(matrix_input_file_name, ".B");
					fp=fopen(matrix_input_file_name,"rb");
							if (type==real_double) {FREAD(b_ref_d, sizeof(b_ref_d[0]), nmat, fp)}
					else	if (type==real_single) {FREAD(b_ref_s, sizeof(b_ref_s[0]), nmat, fp)}
					fclose(fp);
					sdsfree(matrix_input_file_name);
					if (verbose > 0) printf("     ..B\n");

					matrix_input_file_name = sdsdup(matrix_input_base_name);
					matrix_input_file_name = sdscat(matrix_input_file_name, ".A");
					fp=fopen(matrix_input_file_name,"rb");
							if (type==real_double) {FREAD(A_ref_d, sizeof(A_ref_d[0]), nmat*nmat, fp)}
					else	if (type==real_single) {FREAD(A_ref_s, sizeof(A_ref_s[0]), nmat*nmat, fp)}
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
							printf("ERR> To read back the condition number the process grid must be square\n\n");
						}
						MAIN_CLEANUP(mpi_rank);
						MPI_Finalize();
						return ERR_INPUT_ARG;
					}
					else
					{
								// init communication channels (generation uses blacs => mpi interference..)
								test_init(versionname_all[0], verbose, routine_env, routine_input);

								if (type==real_double) cnd_readback = round( pCheckSystemMatrices1D_double(nmat, A_ref_d, x_ref_d, b_ref_d, scalapack_nb, mpi_rank, calc_procs, blacs_nprow, blacs_npcol, blacs_row, blacs_col, blacs_ctxt, blacs_ctxt_root) );
						else	if (type==real_single) cnd_readback = round( pCheckSystemMatrices1D_float (nmat, A_ref_s, x_ref_s, b_ref_s, scalapack_nb, mpi_rank, calc_procs, blacs_nprow, blacs_npcol, blacs_row, blacs_col, blacs_ctxt, blacs_ctxt_root) );
					}
				}
			}
			/*
			 * generated matrices
			 */
			else
			{
				if (strcmp(matrix_gen_method, "par" ) == 0)
				{
					// init communication channels (generation uses blacs => mpi interference..)
					test_init(versionname_all[0], verbose, routine_env, routine_input);

					if (mpi_rank==0 && verbose>0)
					{
						printf("     Matrix random generation alg.: in parallel with ScaLAPACK\n");
						if (!get_cnd)
						{
							printf("WRN> Condition number will not read back from generated matrix\n");
						}
						if (!set_cnd)
						{
							printf("WRN> Matrix will not be pre-conditioned\n");
						}
					}
					if ( get_cnd && !IS_SQUARE(calc_procs) )
					{
						if (mpi_rank==0)
						{
							printf("ERR> To read back the condition number the process grid must be square\n\n");
						}
						MAIN_CLEANUP(mpi_rank);
						MPI_Finalize();
						return ERR_INPUT_ARG;
					}
							if (type==real_double) cnd_readback = round( pGenSystemMatrices1D_double(nmat, A_ref_d, x_ref_d, b_ref_d, seed, cnd, set_cnd, get_cnd, scalapack_nb, mpi_rank, calc_procs, blacs_nprow, blacs_npcol, blacs_row, blacs_col, blacs_ctxt, blacs_ctxt_root) );
					else	if (type==real_single) cnd_readback = round( pGenSystemMatrices1D_float (nmat, A_ref_s, x_ref_s, b_ref_s, seed, cnd, set_cnd, get_cnd, scalapack_nb, mpi_rank, calc_procs, blacs_nprow, blacs_npcol, blacs_row, blacs_col, blacs_ctxt, blacs_ctxt_root) );
				}
				else if (strcmp(matrix_gen_method, "seq" ) == 0)
				{
					if (mpi_rank==0 && verbose>0)
					{
						printf("     Matrix random generation alg.: sequentially with LAPACK\n");
						if (!get_cnd)
						{
							printf("WRN> Condition number will not read back from generated matrix\n");
						}
								if (type==real_double) cnd_readback = round( GenSystemMatrices1D_double(nmat, A_ref_d, x_ref_d, b_ref_d, seed, cnd, set_cnd, get_cnd) );
						else	if (type==real_single) cnd_readback = round( GenSystemMatrices1D_float (nmat, A_ref_s, x_ref_s, b_ref_s, seed, cnd, set_cnd, get_cnd) );
					}
				}
				else
				{
					if (mpi_rank==0)
					{
						printf("ERR> Unknown type of matrix generation %s\n\n",matrix_gen_method);
					}
					MAIN_CLEANUP(mpi_rank);
					MPI_Finalize();
					return ERR_INPUT_ARG;
				}
			}
			sdsfree(matrix_input_base_name);


			// debugging
			// if (mpi_rank==0)
			{
				// create matrices for debugging purposes
				/*
				FillMatrix1D_double(A_ref_d, nmat, nmat);
				OneMatrix1D_double(b_ref_d, nmat, 1);
				*/

				// show matrices for debugging purposes
				/*
				PrintMatrix1D_double(A_ref_d,nmat,nmat);
				PrintVector_double(x_ref_d,nmat);
				PrintVector_double(b_ref_d,nmat);
				*/
			}


			if (type==real_double)
			{
				routine_input.A_ref_d = A_ref_d;
				routine_input.x_ref_d = x_ref_d;
				routine_input.b_ref_d = b_ref_d;
			}
			if (type==real_single)
			{
				routine_input.A_ref_s = A_ref_s;
				routine_input.x_ref_s = x_ref_s;
				routine_input.b_ref_s = b_ref_s;
			}

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
					printf("WRN> Condition number (%d) differs from read back (%d)\n",cnd,cnd_readback);
				}
			}
	}
/*	else // informative commands, checks skipped
	{
		if (type==real_double)
		{
			A_ref_d = NULL;
			x_ref_d = NULL;
			b_ref_d = NULL;
		}
		if (type==real_single)
		{
			A_ref_s = NULL;
			x_ref_s = NULL;
			b_ref_s = NULL;
		}
	}
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
			printf("  -type <string>             : precision type [s|d|single|double] of the input matrix\n");
			printf("                             : (routines for different precisions cannot be mixed)\n" );
			printf("  -no-cnd-set                : disable matrix pre-conditioning\n" );
			printf("  -no-cnd-readback           : disable condition number checking after generation\n" );
			printf("  -no-nre-readback           : disable normwise relative error checking\n" );
#ifdef CRESCO_POWERCAP
			printf("  -energy-reading            : read energy counters\n" );
#endif
			printf("  -mat-gen <string>          : type of the random generation [par|ser] (parallel or sequential)\n" );
			printf("  -r    <integer number>     : run repetitions\n" );
			printf("  -o    <file path>          : output to CSV file\n" );
			printf("  -i    <file path>          : input matrices base name file path (.A, .X, .B auto appended)\n" );
			printf("  -ft   <integer number>     : fault-tolerance level [0-..] (0=none)\n" );
			printf("  -fr   <integer number>     : simulated faulty mpi rank\n" );
			printf("  -fl   <integer number>     : simulated faulty IMe inhibition level\n" );
			printf("  -npf  <integer number>     : number of simulated faults [0-..] (0=none)\n" );
			printf("  -nps  <integer number>     : number of spare processes [0-..] (0=none)\n" );
			printf("  -spk-cp <integer number>   : checkpointing interval\n" );
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
				printf("ERR> Please specify a base file path to save to\n\n");
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
						if (type==real_double) fwrite(x_ref_d,sizeof(x_ref_d[0]),nmat,fp);
				else	if (type==real_single) fwrite(x_ref_s,sizeof(x_ref_s[0]),nmat,fp);
				fclose(fp);
				sdsfree(matrix_output_file_name);
				if (verbose > 0) printf("     ..X\n");

				matrix_output_file_name = sdsdup(matrix_output_base_name);
				matrix_output_file_name = sdscat(matrix_output_file_name, ".B");
				fp=fopen(matrix_output_file_name,"wb");
						if (type==real_double) fwrite(b_ref_d,sizeof(b_ref_d[0]),nmat,fp);
				else	if (type==real_single) fwrite(b_ref_s,sizeof(b_ref_s[0]),nmat,fp);
				fclose(fp);
				sdsfree(matrix_output_file_name);
				if (verbose > 0) printf("     ..B\n");

				matrix_output_file_name = sdsdup(matrix_output_base_name);
				matrix_output_file_name = sdscat(matrix_output_file_name, ".A");
				fp=fopen(matrix_output_file_name,"wb");
						if (type==real_double) fwrite(A_ref_d,sizeof(A_ref_d[0]),nmat*nmat,fp);
				else	if (type==real_single) fwrite(A_ref_s,sizeof(A_ref_s[0]),nmat*nmat,fp);
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
			printf("ERR> Please specify command: --help|--list|--run|--save\n\n");
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
				printf("ERR> Please specify at least one test routine\n\n");
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
			if (strcmp(matrix_gen_method, "par" ) != 0 || input_from_file) test_init(versionname_all[0], verbose, routine_env, routine_input);

			/*
			 * main loop
			 */

#ifdef CRESCO_POWERCAP
			if (rank_in_host==0)
			{
				energy[0]=read_energy(ZONE0);
				energy[1]=read_energy(ZONE0SUB0);
				energy[2]=read_energy(ZONE1);
				energy[3]=read_energy(ZONE1SUB0);
				//printf("Hello from root process on host %s (I'm rank %d in roots)\n", mpi_procname, rank_in_roots);
				if (get_energy)
				{
					MPI_Gather( energy, POWERCAP_ENTRIES, MPI_UNSIGNED_LONG_LONG, &energy_readings_begin[rank_in_roots*POWERCAP_ENTRIES], POWERCAP_ENTRIES, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_ROOTS );
				}
			}
			/*
			if (root_in_host==1)
			{
				printf("Hello from root process on host %s (I'm rank %d in roots)\n", mpi_procname, rank_in_roots);
				MPI_Gather( energy, 4, MPI_UNSIGNED_LONG_LONG, &energy_readings_begin[rank_in_roots*4], 4, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_ROOTS );

			}
			*/
#endif

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

#ifdef CRESCO_POWERCAP
				if (rank_in_host==0)
				{
					energy[0]=read_energy(ZONE0);
					energy[1]=read_energy(ZONE0SUB0);
					energy[2]=read_energy(ZONE1);
					energy[3]=read_energy(ZONE1SUB0);
					if (get_energy)
					{
						MPI_Gather( energy, POWERCAP_ENTRIES, MPI_UNSIGNED_LONG_LONG, &energy_readings_end[rank_in_roots*POWERCAP_ENTRIES], POWERCAP_ENTRIES, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_ROOTS );
					}
				}
	/*
				if (root_in_host==1)
				{
					MPI_Gather( energy, 4, MPI_UNSIGNED_LONG_LONG, &energy_readings_end[rank_in_roots*4], 4, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_ROOTS );

				}
	*/
#endif

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

#ifdef CRESCO_POWERCAP
				if (get_energy)
				{
					printf(" Energy counters:\n");
					printf("%-32s %4s %4s %32s %32s\n","host","zone","sub","begin [uJ]","end [uJ]");
					for (i=0;i<size_in_roots;i++)
					{
						printf("%-32s %4d    - %32llu %32llu\n",mpi_hostnames[i],0,energy_readings_begin[i*POWERCAP_ENTRIES+0],energy_readings_end[i*POWERCAP_ENTRIES+0]);
						printf("%-32s %4d %4d %32llu %32llu\n",mpi_hostnames[i],0,0,energy_readings_begin[i*POWERCAP_ENTRIES+1],energy_readings_end[i*POWERCAP_ENTRIES+1]);
						printf("%-32s %4d    - %32llu %32llu\n",mpi_hostnames[i],1,energy_readings_begin[i*POWERCAP_ENTRIES+2],energy_readings_end[i*POWERCAP_ENTRIES+2]);
						printf("%-32s %4d %4d %32llu %32llu\n",mpi_hostnames[i],1,0,energy_readings_begin[i*POWERCAP_ENTRIES+3],energy_readings_end[i*POWERCAP_ENTRIES+3]);
					}
				}
#endif

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
		DeallocateMatrix1D_double(A_ref);
		DeallocateVector_double(b_ref);
		DeallocateVector_double(x_ref);
	}
	sdsfree(test_output_file_name);
	*/
#ifdef CRESCO_POWERCAP
	free(energy_readings_begin);
	free(energy_readings_end);
	if (mpi_rank==0)
	{
		free(mpi_hostnames[0]);
		free(mpi_hostnames);
	}
	MPI_Comm_free(&MPI_COMM_ROOTS);
	MPI_Comm_free(&MPI_COMM_HOST);
#endif
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

