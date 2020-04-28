
#include "test_IMe_pviDGESV.h"
#include "test_IMe_pviDGESV.h"
#include "tester_head_p.c"

	versionname[0]= "IMe-SV-1D  ";

	versions = 1;

#include "tester_shoulder_p.c"

    	//////////////////////////////////////////////////////////////////////////////////

 		versionrun[0][rep]=test_IMe_pviDGESV_1D(versionname[0], verbose, n, A_ref, x_ref, b_ref, 1, ime_nb, main_rank, cprocs, sprocs);
 		fpdata(0);
    	//////////////////////////////////////////////////////////////////////////////////

#include "tester_tail_p.c"
