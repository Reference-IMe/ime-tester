/*
 * tester_structures.h
 *
 *  Created on: Feb 13, 2020
 *      Author: marcello
 *
 *      https://stackoverflow.com/questions/13156031/measuring-time-in-c
 */

#include <time.h>

#ifndef __TESTER_STRUCTURES_H__
#define __TESTER_STRUCTURES_H__

typedef struct parallel_env
{
	int mpi_rank;
	int blacs_nprow;
	int blacs_npcol;
	int blacs_row;
	int blacs_col;
	int blacs_ctxt_onerow;
	int blacs_ctxt_grid;
	int blacs_ctxt_root;
	int blacs_ctxt_spare;
} parallel_env;

typedef struct test_input
{
	int n;
	double* A_ref;
	double* x_ref;
	double* b_ref;
	int nrhs;
	int calc_procs;
	int spare_procs;
	int ime_bf;
	int scalapack_bf;
	//int failing_rank;
	//int failing_level;
	//int checkpoint_skip_interval;
} test_input;

typedef struct test_output
{
	time_t total_start_time;
	time_t total_end_time;
	time_t core_start_time;
	time_t core_end_time;
	int	   exit_code;
	double norm_rel_err;
} test_output;

typedef struct test_result
{
	double total_time;
	double core_time;
	int	   exit_code;
	double norm_rel_err;
} test_result;

// as array/struct in #define are not allowed in assignment, a constant has to be defined
const test_result TEST_NOT_IMPLEMENTED	= { -99, -99, -99, -99 };
const test_result TEST_NOT_RUN			= { -1,  -1,  -1,  -1  };
const test_output EMPTY_OUTPUT			= { -1,  -1,  -1,  -1,  -1,  -1};

#endif
