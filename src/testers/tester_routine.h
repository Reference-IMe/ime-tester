/*
 * tester_routine.c
 *
 *  Created on: Apr 28, 2020
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
#include "tester_structures.h"

#include "test_ScaLAPACK_pDGESV.h"
#include "test_ScaLAPACK_pDGESV_ft1_cp.h"
#include "test_ScaLAPACK_pDGESV_ft1_cs.h"
#include "test_ScaLAPACK_pDGESV_ftx_cp.h"

#include "test_ScaLAPACK_pDGETRF.h"
#include "test_ScaLAPACK_pDGETRF_ft1_cp.h"

#include "test_ScaLAPACK_pDGEQRF.h"
#include "test_ScaLAPACK_pDGEQRF_ft1_cp.h"

#include "test_FTLA_pDGETRF.h"
#include "test_FTLA_pDGEQRF.h"
#include "test_IMe_pvDGESV.h"
#include "test_IMe_pbDGESV.h"
#include "test_IMe_pbDGESV_ftx.h"

#include "test_IMe_blacsDGESV.h"

/*
#include "test_IMe_pviDGESV_cs.h"
#include "test_IMe_pviDGESV_ft1.h"
#include "test_IMe_pviDGEF.h"
#include "test_IMe_pviDGEF_ft1.h"
*/

test_result tester_routine(const char check, const char* routine_name, int verbosity, parallel_env routine_env, test_input routine_input, int fault_protection, int fault_number, int rank, int failing_rank, int failing_level, int checkpoint_skip_interval)
{
	test_result info;

	/*
	 * obsolete versions, removed
	 */
	/*
	     if	( strcmp( routine_name, IME_PV_SV_WO )						== 0 )	info = test_IMe_pvDGESV(check, routine_name, "PV-WO", verbosity, routine_input, rank);
	else if	( strcmp( routine_name, IME_PV_SV_WO_OAE ) 					== 0 )	info = test_IMe_pvDGESV(check, routine_name, "PV-WO-oae", verbosity, routine_input, rank);
	else if	( strcmp( routine_name, IME_PV_SV_WO_OA ) 					== 0 )	info = test_IMe_pvDGESV(check, routine_name, "PV-WO-oa", verbosity, routine_input, rank);
	else if	( strcmp( routine_name, IME_PV_SV_WO_OGE ) 					== 0 )	info = test_IMe_pvDGESV(check, routine_name, "PV-WO-oge", verbosity, routine_input, rank);
	else if	( strcmp( routine_name, IME_PV_SV_WO_OG ) 					== 0 )	info = test_IMe_pvDGESV(check, routine_name, "PV-WO-og", verbosity, routine_input, rank);

	else if	( strcmp( routine_name, IME_PV_SV_WO_U1AE ) 				== 0 )	info = test_IMe_pvDGESV(check, routine_name, "PV-WO-u1ae", verbosity, routine_input, rank);
	else if	( strcmp( routine_name, IME_PV_SV_WO_U1A ) 					== 0 )	info = test_IMe_pvDGESV(check, routine_name, "PV-WO-u1a", verbosity, routine_input, rank);
	else if	( strcmp( routine_name, IME_PV_SV_WO_U1GE ) 				== 0 )	info = test_IMe_pvDGESV(check, routine_name, "PV-WO-u1ge", verbosity, routine_input, rank);
	else if	( strcmp( routine_name, IME_PV_SV_WO_U1G ) 					== 0 )	info = test_IMe_pvDGESV(check, routine_name, "PV-WO-u1g", verbosity, routine_input, rank);

	else if	( strcmp( routine_name, IME_PV_SV_WO_U2AE ) 				== 0 )	info = test_IMe_pvDGESV(check, routine_name, "PV-WO-u2ae", verbosity, routine_input, rank);
	else if	( strcmp( routine_name, IME_PV_SV_WO_U2A ) 					== 0 )	info = test_IMe_pvDGESV(check, routine_name, "PV-WO-u2a", verbosity, routine_input, rank);
	else if	( strcmp( routine_name, IME_PV_SV_WO_U2GE ) 				== 0 )	info = test_IMe_pvDGESV(check, routine_name, "PV-WO-u2ge", verbosity, routine_input, rank);
	else if	( strcmp( routine_name, IME_PV_SV_WO_U2G ) 					== 0 )	info = test_IMe_pvDGESV(check, routine_name, "PV-WO-u2g", verbosity, routine_input, rank);

	else if	( strcmp( routine_name, IME_PV_SV_WO_U3AE ) 				== 0 )	info = test_IMe_pvDGESV(check, routine_name, "PV-WO-u3ae", verbosity, routine_input, rank);
	else if	( strcmp( routine_name, IME_PV_SV_WO_U3A ) 					== 0 )	info = test_IMe_pvDGESV(check, routine_name, "PV-WO-u3a", verbosity, routine_input, rank);
	else if	( strcmp( routine_name, IME_PV_SV_WO_U3GE ) 				== 0 )	info = test_IMe_pvDGESV(check, routine_name, "PV-WO-u3ge", verbosity, routine_input, rank);
	else if	( strcmp( routine_name, IME_PV_SV_WO_U3G ) 					== 0 )	info = test_IMe_pvDGESV(check, routine_name, "PV-WO-u3g", verbosity, routine_input, rank);

	else if	( strcmp( routine_name, IME_PV_SV_CO ) 						== 0 )	info = test_IMe_pvDGESV(check, routine_name, "PV-CO", verbosity, routine_input, rank);
	else if	( strcmp( routine_name, IME_PV_SV_ICO_G ) 					== 0 )	info = test_IMe_pvDGESV(check, routine_name, "PV-iCO-g", verbosity, routine_input, rank);
	else if	( strcmp( routine_name, IME_PV_SV_CO_G_IND ) 				== 0 )	info = test_IMe_pvDGESV(check, routine_name, "PV-CO-g-ind", verbosity, routine_input, rank);
	else if	( strcmp( routine_name, IME_PV_SV_CO_G_2PASS )				== 0 )	info = test_IMe_pvDGESV(check, routine_name, "PV-CO-g-2pass", verbosity, routine_input, rank);
	else if	( strcmp( routine_name, IME_PV_SV_CO_G_SMALLER )			== 0 )	info = test_IMe_pvDGESV(check, routine_name, "PV-CO-g-smaller", verbosity, routine_input, rank);
	else if	( strcmp( routine_name, IME_PV_SV_CO_G_SMALLEST )			== 0 )	info = test_IMe_pvDGESV(check, routine_name, "PV-CO-g-smallest", verbosity, routine_input, rank);
	else if	( strcmp( routine_name, IME_PV_SV_CO_A_SMALL )				== 0 )	info = test_IMe_pvDGESV(check, routine_name, "PV-CO-a-small", verbosity, routine_input, rank);
	else if	( strcmp( routine_name, IME_PV_SV_CO_A_SMALLER )			== 0 )	info = test_IMe_pvDGESV(check, routine_name, "PV-CO-a-smaller", verbosity, routine_input, rank);
	else if	( strcmp( routine_name, IME_PV_SV_CO_A_SMALLEST )			== 0 )	info = test_IMe_pvDGESV(check, routine_name, "PV-CO-a-smallest", verbosity, routine_input, rank);

	else */

	     if	( strcmp( routine_name, IME_PB_SV_CO_BF1 )					== 0 )	info = test_IMe_pbDGESV    (check, routine_name, "PB-CO-bf1",       verbosity, routine_input, rank);
	else if	( strcmp( routine_name, IME_PB_SV_CO_BF1_FAULT_0_TOLERANT_X)== 0 )	info = test_IMe_pbDGESV_ftx(check, routine_name, "PB-CO-bf1-ftx/0", verbosity, routine_env, routine_input, fault_protection, fault_number, failing_rank, -1,            0);
	else if	( strcmp( routine_name, IME_PB_SV_CO_BF1_FAULT_X_TOLERANT_X)== 0 )	info = test_IMe_pbDGESV_ftx(check, routine_name, "PB-CO-bf1-ftx/x", verbosity, routine_env, routine_input, fault_protection, fault_number, failing_rank, failing_level, 1);
	else if	( strcmp( routine_name, IME_PB_SV_CO_BF1_FAULT_X_TOLERANT_0)== 0 )	info = test_IMe_pbDGESV_ftx(check, routine_name, "PB-CO-bf1-ft0/x", verbosity, routine_env, routine_input, fault_protection, fault_number, failing_rank, failing_level, 0);
	else if	( strcmp( routine_name, IME_DEV )							== 0 )	info = test_IMe_pbDGESV_ftx(check, routine_name, "dev",             verbosity, routine_env, routine_input, fault_protection, fault_number, failing_rank, failing_level, 1);

	else if	( strcmp( routine_name, IME_BLACS_SV_CO_1 )					== 0 )	info = test_IMe_blacsDGESV (check, routine_name, "BLACS-CO-1", verbosity, routine_env, routine_input);
	else if	( strcmp( routine_name, IME_BLACS_SV_CO_2 )					== 0 )	info = test_IMe_blacsDGESV (check, routine_name, "BLACS-CO-2", verbosity, routine_env, routine_input);

	else if	( strcmp( routine_name, SPK_SV ) 						== 0 )	info = test_ScaLAPACK_pDGESV       (check, routine_name, verbosity, routine_env, routine_input);
	else if	( strcmp( routine_name, SPK_SV_FAULT_0_TOLERANT_1_CP ) 	== 0 )	info = test_ScaLAPACK_pDGESV_ft1_cp(check, routine_name, verbosity, routine_env, routine_input, fault_protection, fault_number, -1,            checkpoint_skip_interval);
	else if	( strcmp( routine_name, SPK_SV_FAULT_1_TOLERANT_1_CP ) 	== 0 )	info = test_ScaLAPACK_pDGESV_ft1_cp(check, routine_name, verbosity, routine_env, routine_input, fault_protection, fault_number, failing_level, checkpoint_skip_interval);
	else if	( strcmp( routine_name, SPK_SV_FAULT_0_TOLERANT_1_CS ) 	== 0 )	info = test_ScaLAPACK_pDGESV_ft1_cs(check, routine_name, verbosity, routine_env, routine_input, -1, checkpoint_skip_interval);
	else if	( strcmp( routine_name, SPK_SV_FAULT_0_TOLERANT_X_CP ) 	== 0 )	info = test_ScaLAPACK_pDGESV_ftx_cp(check, routine_name, verbosity, routine_env, routine_input, -1, checkpoint_skip_interval);

	else if	( strcmp( routine_name, SPK_LU ) 					== 0 )	info = test_ScaLAPACK_pDGETRF    (check, routine_name, verbosity, routine_env, routine_input);
	else if	( strcmp( routine_name, SPK_LU_FAULT_0_TOLERANT_1 ) == 0 )	info = test_ScaLAPACK_pDGETRF_ft1(check, routine_name, verbosity, routine_env, routine_input, -1, checkpoint_skip_interval);
	else if	( strcmp( routine_name, SPK_LU_FAULT_1_TOLERANT_1 ) == 0 )	info = test_ScaLAPACK_pDGETRF_ft1(check, routine_name, verbosity, routine_env, routine_input, failing_level, checkpoint_skip_interval);

	else if	( strcmp( routine_name, SPK_QR ) 					== 0 )	info = test_ScaLAPACK_pDGEQRF    (check, routine_name, verbosity, routine_env, routine_input);
	else if	( strcmp( routine_name, SPK_QR_FAULT_0_TOLERANT_1 ) == 0 )	info = test_ScaLAPACK_pDGEQRF_ft1(check, routine_name, verbosity, routine_env, routine_input, -1, checkpoint_skip_interval);
	else if	( strcmp( routine_name, SPK_QR_FAULT_1_TOLERANT_1 ) == 0 )	info = test_ScaLAPACK_pDGEQRF_ft1(check, routine_name, verbosity, routine_env, routine_input, failing_level, checkpoint_skip_interval);

	else if	( strcmp( routine_name, FTLA_LU_FAULT_0_TOLERANT_1 ) == 0 )	info = test_FTLA_pDGETRF(check, routine_name, verbosity, routine_env, routine_input, 0);
	else if	( strcmp( routine_name, FTLA_LU_FAULT_1_TOLERANT_1 ) == 0 )	info = test_FTLA_pDGETRF(check, routine_name, verbosity, routine_env, routine_input, 1);
	else if	( strcmp( routine_name, FTLA_LU_FAULT_X_TOLERANT_X ) == 0 )	info = test_FTLA_pDGETRF(check, routine_name, verbosity, routine_env, routine_input, fault_protection);

	else if	( strcmp( routine_name, FTLA_QR_FAULT_0_TOLERANT_1 ) == 0 )	info = test_FTLA_pDGEQRF(check, routine_name, verbosity, routine_env, routine_input, 0);
	else if	( strcmp( routine_name, FTLA_QR_FAULT_1_TOLERANT_1 ) == 0 )	info = test_FTLA_pDGEQRF(check, routine_name, verbosity, routine_env, routine_input, 1);
	else if	( strcmp( routine_name, FTLA_QR_FAULT_X_TOLERANT_X ) == 0 )	info = test_FTLA_pDGEQRF(check, routine_name, verbosity, routine_env, routine_input, fault_protection);

	else info = TEST_NOT_IMPLEMENTED;

	return info;
}
