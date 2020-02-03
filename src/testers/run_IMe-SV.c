
#include "test_IMe_pviDGESV.h"
#include "tester_head_p.c"

	versionname[0]= "IMe-SV      ";

	versions = 1;

#include "tester_shoulder_p.c"

    	//////////////////////////////////////////////////////////////////////////////////

    	versionrun[0][rep]=test_IMe_pviDGESV(versionname[0], verbose, rows, cols, 1, main_rank, cprocs, sprocs);
 		fpdata(0);

    	//////////////////////////////////////////////////////////////////////////////////

#include "tester_tail_p.c"
