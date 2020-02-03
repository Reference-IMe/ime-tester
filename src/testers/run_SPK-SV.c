
#include "test_ScaLAPACK_pDGESV.h"
#include "tester_head_p.c"

	versionname[0]= "SPK-SV      ";

	versions = 1;

#include "tester_shoulder_p.c"

    	//////////////////////////////////////////////////////////////////////////////////

		versionrun[0][rep]=test_ScaLAPACK_pDGESV(versionname[0], verbose, rows, cols, 1, nb, main_rank, cprocs);
		fpdata(0);

    	//////////////////////////////////////////////////////////////////////////////////

#include "tester_tail_p.c"
