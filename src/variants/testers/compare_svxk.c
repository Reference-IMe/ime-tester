
#include "test_IMe_pviDGESV.h"
#include "test_IMe_pviDGEF.h"
#include "test_ScaLAPACK_pDGESV.h"

#include "tester_head_p.c"

	versionname[0]= "SPK-SV          ";
	versionname[1]= "IMe-SV          ";
	versionname[2]= "IMe-XK          ";

	versions = 3;

#include "tester_shoulder_p.c"

    	//////////////////////////////////////////////////////////////////////////////////

    	versionrun[ 0][rep]=test_ScaLAPACK_pDGESV(versionname[0], verbose, rows, cols, nb, 1, main_rank, cprocs);		// SPK solve with 1 rhs
 		fpdata(0);

     	versionrun[ 1][rep]=test_IMe_pviDGESV(versionname[1], verbose, rows, cols, 1, main_rank, cprocs, sprocs);		// vanilla IMe solve with 1 rhs
     	fpdata(1);

		versionrun[2][rep]=test_IMe_pviDGEF(versionname[2], verbose, rows, cols, main_rank, cprocs, sprocs); // IMe factorization only
		fpdata(2);

    	//////////////////////////////////////////////////////////////////////////////////

#include "tester_tail_p.c"
