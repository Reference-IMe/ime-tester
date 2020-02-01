
#include "test_IMe_pviDGESV.h"
#include "tester_head_p.c"

	versionname[0]= "IMe-SV      ";
	versionname[1]= "IMe-SV early";

	versions = 2;

#include "tester_shoulder_p.c"

    	//////////////////////////////////////////////////////////////////////////////////

    	versionrun[ 0][rep]=test_IMe_pviDGESV(versionname[0], verbose, rows, cols, 1, main_rank, cprocs, sprocs); // IMe factorization only
 		fpdata(0);

		versionrun[1][rep]=test_IMe_pviDGESV(versionname[1], verbose, rows, cols, 1, main_rank, cprocs, sprocs); // IMe factorization only
		fpdata(1);

		//versionrun[2][rep]=test_IMe_pviDGESV_gather(versionname[2], verbose, rows, cols, 1, main_rank, cprocs, sprocs); // IMe factorization only
		//fpdata(2);

    	//////////////////////////////////////////////////////////////////////////////////

#include "tester_tail_p.c"
