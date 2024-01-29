/*
 * tester_structures.h
 *
 *  Created on: Feb 13, 2020
 *      Author: marcello
 *
 *      https://stackoverflow.com/questions/13156031/measuring-time-in-c
 */

#include <time.h>
#include <mpi.h>

#ifndef __TESTER_STRUCTURES_H__
#define __TESTER_STRUCTURES_H__

#include "IMe/lib/src/helpers/types.h"

typedef struct parallel_env
{
	int calc_procs;
	int spare_procs;
	int mpi_rank;
	MPI_Comm mpi_comm;
	int blacs_nprow;
	int blacs_npcol;
	int blacs_row;
	int blacs_col;
	int blacs_ctxt_onerow;
	int blacs_ctxt_grid;
	int blacs_ctxt_root;
	int* blacs_ctxt_spare;
} parallel_env;


typedef struct fault_env
{
	int fault_tolerance;
	int faulty_procs;
	int failing_rank;
	int failing_level;
	int scalapack_checkpoint_interval;
} fault_env;

typedef struct test_input
{
	int n;
	double* A_ref_d;
	double* x_ref_d;
	double* b_ref_d;
	float* A_ref_s;
	float* x_ref_s;
	float* b_ref_s;
	int nrhs;
	int ime_bf;
	int scalapack_bf;
	char calc_nre;
} test_input;

/*
typedef struct exit_status
{
	time_t total_start_time;
	time_t total_end_time;
	time_t core_start_time;
	time_t core_end_time;
	char   exit_code;
	//double norm_rel_err;
} exit_status;
*/

typedef struct test_result
{
	double total_time;
	double core_time;
	char   exit_code;
	double norm_rel_err;
} test_result;

// as array/struct in #define are not allowed in assignment, a constant has to be defined
const test_result TEST_NOT_IMPLEMENTED	= { -99, -99, -99, -99 };
const test_result TEST_NOT_RUN			= { -1,  -1,  -1,  -1  };
const exit_status EMPTY_OUTPUT			= { -1,  -1,  -1,  -1,  -1};

#endif
