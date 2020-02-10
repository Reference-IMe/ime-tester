
#include "test_IMe_pviDGESV.early.h"
#include "test_IMe_pviDGESV.h"
#include "test_ScaLAPACK_pDGESV.h"
#include "tester_head_p.c"

	versionname[0]= "SPK-SV      ";
	versionname[1]= "IMe-SV      ";
	versionname[2]= "IMe-SV-early";


	versions = 3;

#include "tester_shoulder_p.c"

    	//////////////////////////////////////////////////////////////////////////////////

		versionrun[0][rep]=test_ScaLAPACK_pDGESV(versionname[0], verbose, rows, cols, 1, nb, main_rank, cprocs);
		fpdata(0);

    	versionrun[1][rep]=test_IMe_pviDGESV(versionname[1], verbose, rows, cols, 1, main_rank, cprocs, sprocs);
 		fpdata(1);

		versionrun[2][rep]=test_IMe_pviDGESV_early(versionname[2], verbose, rows, cols, 1, main_rank, cprocs, sprocs); // IMe factorization only
		fpdata(2);

    	//////////////////////////////////////////////////////////////////////////////////

#include "tester_tail_p.c"
