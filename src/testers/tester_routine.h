/*
 * tester_routine.c
 *
 *  Created on ) == 0 ) Apr 28, 2020
 *      Author ) == 0 ) marcello
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
#include "tester_structures.h"

#include "test_IMe_pviDGESV.h"
#include "test_ScaLAPACK_pDGESV.h"
#include "test_ScaLAPACK_pDGETRF.h"

/*
#include "test_IMe_pviDGESV_cs.h"
#include "test_IMe_pviDGESV_ft1.h"

#include "test_IMe_pviDGEF.h"
#include "test_IMe_pviDGEF_ft1.h"


#include "test_ScaLAPACK_pDGEQRF.h"

#include "test_FTLA_pDGEQRF.h"
#include "test_FTLA_pDGETRF.h"

#include "test_ScaLAPACK_pDGESV_cp_ft1.h"
#include "test_ScaLAPACK_pDGETRF_cp_ft1.h"
#include "test_ScaLAPACK_pDGEQRF_cp_ft1.h"
*/

test_result tester_routine(const char* routine_name, int verbosity, test_input routine_input, int rank)
{
	test_result info;

	if		( strcmp( routine_name, IME_SV ) 					== 0 )	info = test_IMe_pviDGESV(routine_name, verbosity, routine_input, rank);
//	else if	( strcmp( routine_name, IME_SV_CHECKSUMMED ) 		== 0 )	info = test_IMe_pviDGESV_cs(routine_name, verbosity, rows, cols, 1, rank, sprocs);
/*	else if	( strcmp( routine_name, IME_SV_FAULT_0_TOLERANT_1 ) == 0 )	info = test_IMe_pviDGESV_ft1_sim(routine_name, verbosity, rows, cols, 1, rank, sprocs, -1, -1);
	else if	( strcmp( routine_name, IME_SV_FAULT_1_TOLERANT_1 ) == 0 )	info = test_IMe_pviDGESV_ft1_sim(routine_name, verbosity, rows, cols, 1, rank, sprocs, -1, -1);
	else if	( strcmp( routine_name, IME_XK ) 					== 0 )	info = test_IMe_pviDGEF(routine_name, verbosity, rows, cols, rank, cprocs, sprocs);
	else if	( strcmp( routine_name, IME_XK_FAULT_0_TOLERANT_1 ) == 0 )	info = test_IMe_pviDGEF_ft1_sim(routine_name, verbosity, rows, cols, rank, sprocs, -1, -1);
	else if	( strcmp( routine_name, IME_XK_FAULT_1_TOLERANT_1 ) == 0 )	info = test_IMe_pviDGEF_ft1_sim(routine_name, verbosity, rows, cols, rank, sprocs, failing_rank, failing_level);
*/	else if	( strcmp( routine_name, SPK_SV ) 					== 0 )	info = test_ScaLAPACK_pDGESV(routine_name, verbosity, routine_input, rank);
	/*
	else if	( strcmp( routine_name, SPK_SV_FAULT_0_TOLERANT_1 ) == 0 )	break;
	else if	( strcmp( routine_name, SPK_SV_FAULT_1_TOLERANT_1 ) == 0 )	break;
	*/
	else if	( strcmp( routine_name, SPK_LU ) 					== 0 )	info = test_ScaLAPACK_pDGETRF(routine_name, verbosity, routine_input, rank);
/*	else if	( strcmp( routine_name, SPK_LU_FAULT_0_TOLERANT_1 ) == 0 )	info = test_ScaLAPACK_pDGETRF_cp_ft1_sim(routine_name, verbosity, rows, cols, nb, rank, cprocs, sprocs, -1, checkpoint_skip_interval);
	else if	( strcmp( routine_name, SPK_LU_FAULT_1_TOLERANT_1 ) == 0 )	info = test_ScaLAPACK_pDGETRF_cp_ft1_sim(routine_name, verbosity, rows, cols, nb, rank, cprocs, sprocs, failing_level, checkpoint_skip_interval);
	else if	( strcmp( routine_name, SPK_QR ) 					== 0 )	info = test_ScaLAPACK_pDGEQRF(routine_name, verbosity, rows, cols, nb, rank, cprocs);
	else if	( strcmp( routine_name, SPK_QR_FAULT_0_TOLERANT_1 ) == 0 )	info = test_ScaLAPACK_pDGEQRF_cp_ft1_sim(routine_name, verbosity, rows, cols, nb, rank, cprocs, sprocs, -1, checkpoint_skip_interval);
	else if	( strcmp( routine_name, SPK_QR_FAULT_1_TOLERANT_1 ) == 0 )	info = test_ScaLAPACK_pDGEQRF_cp_ft1_sim(routine_name, verbosity, rows, cols, nb, rank, cprocs, sprocs, failing_level, checkpoint_skip_interval);
*/	/*
	else if	( strcmp( routine_name, FTLA_LU_FAULT_0_TOLERANT_1 ) == 0 )	break;
	else if	( strcmp( routine_name, FTLA_LU_FAULT_1_TOLERANT_1 ) == 0 )	break;
	*/
/*	else if	( strcmp( routine_name, FTLA_QR_FAULT_0_TOLERANT_1 ) == 0 )	info = test_FTLA_pDGEQRF(routine_name, verbosity, rows, cols, nb, rank, cprocs, 0);
	else if	( strcmp( routine_name, FTLA_QR_FAULT_1_TOLERANT_1 ) == 0 )	info = test_FTLA_pDGEQRF(routine_name, verbosity, rows, cols, nb, rank, cprocs, 1);
*/
	else info = TEST_NOT_IMPLEMENTED;

	return info;
}
