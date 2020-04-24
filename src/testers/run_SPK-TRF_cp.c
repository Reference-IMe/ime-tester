
#include "test_ScaLAPACK_pDGETRF_cp_ft1.h"
#include "tester_head_p.c"

	versionname[0]= "SPK-TRF-cp  ";

	versions = 1;

#include "tester_shoulder_p.c"

    	//////////////////////////////////////////////////////////////////////////////////

		versionrun[0][rep]=test_ScaLAPACK_pDGETRF_cp_ft1_sim(versionname[0], verbose, n, A_ref, scalapack_nb, main_rank, cprocs, sprocs, failing_level, checkpoint_skip_interval);
		fpdata(0);

    	//////////////////////////////////////////////////////////////////////////////////

#include "tester_tail_p.c"
