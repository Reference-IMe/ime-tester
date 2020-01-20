
#include "test_IMe_pviDGESV.h"
#include "test_IMe_pviDGESV.nocacheopt.h"
#include "test_IMe_pviDGESV.allgather.h"
#include "test_IMe_pviDGESV.gather.h"
//#include "test_ScaLAPACK_pDGESV.h"

#include "tester_head_p.c"

	versionname[0]= "IMe-SV no cache opt";
	versionname[1]= "IMe-SV             ";
	versionname[2]= "IMe-SV gather      ";
	versionname[3]= "IMe-SV allgather   ";

	versions = 4;

#include "tester_shoulder_p.c"

    	//////////////////////////////////////////////////////////////////////////////////

    	versionrun[ 0][rep]=test_IMe_pviDGESV_nocacheopt(versionname[0], verbose, rows, cols, 1, main_rank, cprocs, sprocs); // IMe factorization only
 		fpdata(0);

		versionrun[1][rep]=test_IMe_pviDGESV(versionname[1], verbose, rows, cols, 1, main_rank, cprocs, sprocs); // IMe factorization only
		fpdata(1);

		versionrun[2][rep]=test_IMe_pviDGESV_gather(versionname[2], verbose, rows, cols, 1, main_rank, cprocs, sprocs); // IMe factorization only
		fpdata(2);

		versionrun[3][rep]=test_IMe_pviDGESV_allgather(versionname[3], verbose, rows, cols, 1, main_rank, cprocs, sprocs); // IMe factorization only
		fpdata(3);

    	//////////////////////////////////////////////////////////////////////////////////

#include "tester_tail_p.c"
