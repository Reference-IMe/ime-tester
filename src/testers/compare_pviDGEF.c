
#include "test_IMe_pviDGEF.h"
#include "test_IMe_pviDGEF.nocacheopt.h"
#include "test_IMe_pviDGEF.allgather.h"
#include "test_IMe_pviDGEF.gather.h"
//#include "test_ScaLAPACK_pDGESV.h"

#include "tester_head_p.c"

	versionname[0]= "IMe-XK no cache opt";
	versionname[1]= "IMe-XK             ";
	versionname[2]= "IMe-XK gather      ";
	versionname[3]= "IMe-XK allgather   ";

	versions = 4;

#include "tester_shoulder_p.c"

    	//////////////////////////////////////////////////////////////////////////////////

    	versionrun[ 0][rep]=test_IMe_pviDGEF_nocacheopt(versionname[0], verbose, rows, cols, main_rank, cprocs, sprocs); // IMe factorization only
 		fpdata(0);

		versionrun[1][rep]=test_IMe_pviDGEF(versionname[1], verbose, rows, cols, main_rank, cprocs, sprocs); // IMe factorization only
		fpdata(1);

		versionrun[2][rep]=test_IMe_pviDGEF_gather(versionname[2], verbose, rows, cols, main_rank, cprocs, sprocs); // IMe factorization only
		fpdata(2);

		versionrun[3][rep]=test_IMe_pviDGEF_allgather(versionname[3], verbose, rows, cols, main_rank, cprocs, sprocs); // IMe factorization only
		fpdata(3);

    	//////////////////////////////////////////////////////////////////////////////////

#include "tester_tail_p.c"
