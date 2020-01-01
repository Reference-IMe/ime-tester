/*
 * errors.h
 *
 *  Created on: Dec 29, 2019
 *      Author: marcello
 */

#ifndef __FTLA_COMMONS_H__
#define __FTLA_COMMONS_H__

// code extracted from ftdqr_main.c of FTLA

/*
 * error vector for injection (global var)
 */
int *errors;

/*
 * mpi communicator override
 */
#include <mpi.h>
MPI_Comm ftla_current_comm;

#endif
