
#include "test_IMe_pviDGEF.h"
#include "test_IMe_pviDGEF_ft1.h"

#include "test_ScaLAPACK_pDGETRF.h"
#include "test_ScaLAPACK_pDGEQRF.h"

#include "test_ScaLAPACK_pDGETRF_cp_ft1.h"
#include "test_ScaLAPACK_pDGEQRF_cp_ft1.h"

#include "tester_head_p.c"

	versionname[0]= "SPK-LU          ";
	versionname[1]= "SPK-LU-ft1/0    ";
	versionname[2]= "SPK-LU-ft1/1    ";

	versionname[3]= "SPK-QR          ";
	versionname[4]= "SPK-QR-ft1/0    ";
	versionname[5]= "SPK-QR-ft1/1    ";

	versionname[6]= "IMe-XK          ";

	versions = 7;

#include "tester_shoulder_p.c"

    	//////////////////////////////////////////////////////////////////////////////////

    	versionrun[0][rep]=test_ScaLAPACK_pDGETRF(versionname[0], verbose, rows, cols, nb, main_rank, cprocs); // SPK LU factorization
 		fpdata(0);

		versionrun[1][rep]=test_ScaLAPACK_pDGETRF_cp_ft1_sim(versionname[1], verbose, rows, cols, nb, main_rank, cprocs, sprocs, -1, checkpoint_skip_interval);	// SPKmod LU factorization single FT and 0 faults
		fpdata(1);

		versionrun[2][rep]=test_ScaLAPACK_pDGETRF_cp_ft1_sim(versionname[2], verbose, rows, cols, nb, main_rank, cprocs, sprocs, failing_level, checkpoint_skip_interval);	// SPKmod LU factorization single FT and 1 fault
		fpdata(2);

    	versionrun[3][rep]=test_ScaLAPACK_pDGEQRF(versionname[3], verbose, rows, cols, nb, main_rank, cprocs); // SPK LU factorization
    	fpdata(3);

		versionrun[4][rep]=test_ScaLAPACK_pDGEQRF_cp_ft1_sim(versionname[4], verbose, rows, cols, nb, main_rank, cprocs, sprocs, -1, checkpoint_skip_interval);	// SPKmod QR factorization single FT and 0 faults
		fpdata(4);

		versionrun[5][rep]=test_ScaLAPACK_pDGEQRF_cp_ft1_sim(versionname[5], verbose, rows, cols, nb, main_rank, cprocs, sprocs, failing_level, checkpoint_skip_interval);	// SPKmod QR factorization single FT and 1 fault
		fpdata(5);

		versionrun[6][rep]=test_IMe_pviDGEF(versionname[6], verbose, rows, cols, main_rank, cprocs, sprocs); // IMe factorization only
		fpdata(6);

    	//////////////////////////////////////////////////////////////////////////////////

#include "tester_tail_p.c"
