
#include "test_IMe_pviDGESV.h"
#include "test_IMe_pviDGESV_1D.h"
#include "tester_head_p.c"

	versionname[0]= "IMe-SV     ";
	versionname[1]= "IMe-SV-1D  ";

	versions = 2;

#include "tester_shoulder_p.c"

    	//////////////////////////////////////////////////////////////////////////////////

    	versionrun[0][rep]=test_IMe_pviDGESV(versionname[0], verbose, n, A_ref, x_ref, b_ref, 1, main_rank, cprocs, sprocs);
 		fpdata(0);
 		versionrun[1][rep]=test_IMe_pviDGESV_1D(versionname[1], verbose, rows, cols, 1, ime_nb, main_rank, cprocs, sprocs);
 		fpdata(1);
    	//////////////////////////////////////////////////////////////////////////////////

#include "tester_tail_p.c"
