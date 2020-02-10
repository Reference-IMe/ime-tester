
#include "test_IMe_pviDGESV.early.h"
#include "tester_head_p.c"

	versionname[0]= "IMe-SV-early";

	versions = 1;

#include "tester_shoulder_p.c"

    	//////////////////////////////////////////////////////////////////////////////////

    	versionrun[0][rep]=test_IMe_pviDGESV_early(versionname[0], verbose, rows, cols, 1, main_rank, cprocs, sprocs);
 		fpdata(0);

    	//////////////////////////////////////////////////////////////////////////////////

#include "tester_tail_p.c"
