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
#include "../helpers/matrix_advanced.h"
#include "../helpers/vector.h"
#include "../helpers/lapack.h"
#include "../helpers/matrix_basic.h"
#include "../helpers/scalapack.h"

#include "tester_labels.h"
#include "tester_structures.h"

#include "test_ScaLAPACK_pDGESV.h"
#include "test_ScaLAPACK_pDGESV_nopivot.h"
//#include "test_ScaLAPACK_pDGESV_ft1_cp.h"
#include "test_ScaLAPACK_pDGESV_ftx_cp.h"

#include "test_ScaLAPACK_pSGESV.h"
#include "test_ScaLAPACK_pSGESV_ftx_cp.h"

#include "test_ScaLAPACK_pDGETRF.h"
#include "test_ScaLAPACK_pDGETRF_ft1_cp.h"

#include "test_ScaLAPACK_pDGEQRF.h"
#include "test_ScaLAPACK_pDGEQRF_ft1_cp.h"

#include "test_FTLA_pDGETRF.h"
#include "test_FTLA_pDGESV_TRF.h"
#include "test_FTLA_pDGEQRF.h"
#include "test_FTLA_pDGESV_QRF.h"

#include "test_IMe_pDGESV-CO.h"
#include "test_IMe_pDGESV-CO-FT.h"
#include "test_IMe_pDGESV-WO.h"

#include "test_IMe_pSGESV-CO.h"
#include "test_IMe_pSGESV-CO-FT.h"
#include "test_IMe_pSGESV-WO.h"

test_result tester_routine(const char check, const char* routine_name, int verbosity, parallel_env routine_env, test_input routine_input, fault_env routine_fault)
{
	test_result info;

		 if	( strcmp( routine_name, IME_PDGESV_CO	)			== 0 )	info = test_IMe_pDGESV_CO		(check, routine_name, "PB-CO-BF1",       verbosity, routine_env, routine_input, routine_fault.fault_tolerance);
	else if	( strcmp( routine_name, IME_PDGESV_CO_FT)			== 0 )	info = test_IMe_pDGESV_CO_FT	(check, routine_name, "PB-CO-BF1-FT",    verbosity, routine_env, routine_input, routine_fault.fault_tolerance, routine_fault.faulty_procs, routine_fault.failing_rank, routine_fault.failing_level, 1);
	else if	( strcmp( routine_name, IME_PDGESV_WO   )			== 0 )	info = test_IMe_pDGESV_WO		(check, routine_name, "PB-WO-BF1",       verbosity, routine_env, routine_input, routine_fault.fault_tolerance);
	else if	( strcmp( routine_name, IME_PSGESV_CO   )			== 0 )	info = test_IMe_pSGESV_CO		(check, routine_name, "PB-CO-BF1",       verbosity, routine_env, routine_input, routine_fault.fault_tolerance);
	else if	( strcmp( routine_name, IME_PSGESV_CO_FT)			== 0 )	info = test_IMe_pSGESV_CO_FT	(check, routine_name, "PB-CO-BF1-FT",    verbosity, routine_env, routine_input, routine_fault.fault_tolerance, routine_fault.faulty_procs, routine_fault.failing_rank, routine_fault.failing_level, 1);
	else if	( strcmp( routine_name, IME_PSGESV_WO   )			== 0 )	info = test_IMe_pSGESV_WO		(check, routine_name, "PB-WO-BF1",       verbosity, routine_env, routine_input, routine_fault.fault_tolerance);

	else if	( strcmp( routine_name, SPK_PDGESV ) 				== 0 )	info = test_ScaLAPACK_pDGESV		(check, routine_name, verbosity, routine_env, routine_input);
	else if	( strcmp( routine_name, SPK_PSGESV ) 				== 0 )	info = test_ScaLAPACK_pSGESV		(check, routine_name, verbosity, routine_env, routine_input);
	else if	( strcmp( routine_name, SPK_PDGESV_NOPIV ) 			== 0 )	info = test_ScaLAPACK_pDGESV_nopivot(check, routine_name, verbosity, routine_env, routine_input);
	else if	( strcmp( routine_name, SPK_PDGESV_FT_CP ) 			== 0 )	info = test_ScaLAPACK_pDGESV_ftx_cp	(check, routine_name, verbosity, routine_env, routine_input, routine_fault.fault_tolerance, routine_fault.faulty_procs, routine_fault.failing_rank, routine_fault.failing_level, routine_fault.scalapack_checkpoint_interval);
	else if	( strcmp( routine_name, SPK_PSGESV_FT_CP ) 			== 0 )	info = test_ScaLAPACK_pSGESV_ftx_cp	(check, routine_name, verbosity, routine_env, routine_input, routine_fault.fault_tolerance, routine_fault.faulty_procs, routine_fault.failing_rank, routine_fault.failing_level, routine_fault.scalapack_checkpoint_interval);

	else if	( strcmp( routine_name, SPK_LU ) 					== 0 )	info = test_ScaLAPACK_pDGETRF    (check, routine_name, verbosity, routine_env, routine_input);
	else if	( strcmp( routine_name, SPK_LU_FAULT_0_TOLERANT_1 ) == 0 )	info = test_ScaLAPACK_pDGETRF_ft1(check, routine_name, verbosity, routine_env, routine_input, 1, 0, -1,                          routine_fault.scalapack_checkpoint_interval);
	else if	( strcmp( routine_name, SPK_LU_FAULT_1_TOLERANT_1 ) == 0 )	info = test_ScaLAPACK_pDGETRF_ft1(check, routine_name, verbosity, routine_env, routine_input, 1, 1, routine_fault.failing_level, routine_fault.scalapack_checkpoint_interval);

	else if	( strcmp( routine_name, SPK_QR ) 					== 0 )	info = test_ScaLAPACK_pDGEQRF    (check, routine_name, verbosity, routine_env, routine_input);
	else if	( strcmp( routine_name, SPK_QR_FAULT_0_TOLERANT_1 ) == 0 )	info = test_ScaLAPACK_pDGEQRF_ft1(check, routine_name, verbosity, routine_env, routine_input, 1, 0, -1,                          routine_fault.scalapack_checkpoint_interval);
	else if	( strcmp( routine_name, SPK_QR_FAULT_1_TOLERANT_1 ) == 0 )	info = test_ScaLAPACK_pDGEQRF_ft1(check, routine_name, verbosity, routine_env, routine_input, 1, 1, routine_fault.failing_level, routine_fault.scalapack_checkpoint_interval);

	else if	( strcmp( routine_name, FTLA_LU_FAULT_0_TOLERANT_1 ) == 0 )	info = test_FTLA_pDGETRF(check, routine_name, verbosity, routine_env, routine_input, 0);
	else if	( strcmp( routine_name, FTLA_LU_FAULT_1_TOLERANT_1 ) == 0 )	info = test_FTLA_pDGETRF(check, routine_name, verbosity, routine_env, routine_input, 1);
	else if	( strcmp( routine_name, FTLA_LU_FAULT_X_TOLERANT_X ) == 0 )	info = test_FTLA_pDGETRF(check, routine_name, verbosity, routine_env, routine_input, routine_fault.faulty_procs);

	else if	( strcmp( routine_name, FTLA_LU_SV_FAULT_0_TOLERANT_1 ) == 0 )	info = test_FTLA_pDGESV_TRF(check, routine_name, verbosity, routine_env, routine_input, 0);
	else if	( strcmp( routine_name, FTLA_LU_SV_FAULT_1_TOLERANT_1 ) == 0 )	info = test_FTLA_pDGESV_TRF(check, routine_name, verbosity, routine_env, routine_input, 1);
	else if	( strcmp( routine_name, FTLA_LU_SV_FAULT_X_TOLERANT_X ) == 0 )	info = test_FTLA_pDGESV_TRF(check, routine_name, verbosity, routine_env, routine_input, routine_fault.faulty_procs);

	else if	( strcmp( routine_name, FTLA_QR_FAULT_0_TOLERANT_1 ) == 0 )	info = test_FTLA_pDGEQRF(check, routine_name, verbosity, routine_env, routine_input, 0);
	else if	( strcmp( routine_name, FTLA_QR_FAULT_1_TOLERANT_1 ) == 0 )	info = test_FTLA_pDGEQRF(check, routine_name, verbosity, routine_env, routine_input, 1);
	else if	( strcmp( routine_name, FTLA_QR_FAULT_X_TOLERANT_X ) == 0 )	info = test_FTLA_pDGEQRF(check, routine_name, verbosity, routine_env, routine_input, routine_fault.faulty_procs);

	else if	( strcmp( routine_name, FTLA_QR_SV_FAULT_0_TOLERANT_1 ) == 0 )	info = test_FTLA_pDGESV_QRF(check, routine_name, verbosity, routine_env, routine_input, 0);
	else if	( strcmp( routine_name, FTLA_QR_SV_FAULT_1_TOLERANT_1 ) == 0 )	info = test_FTLA_pDGESV_QRF(check, routine_name, verbosity, routine_env, routine_input, 1);
	else if	( strcmp( routine_name, FTLA_QR_SV_FAULT_X_TOLERANT_X ) == 0 )	info = test_FTLA_pDGESV_QRF(check, routine_name, verbosity, routine_env, routine_input, routine_fault.faulty_procs);

	else info = TEST_NOT_IMPLEMENTED;

	return info;
}
