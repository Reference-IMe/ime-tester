
#include "test_IMe_pviDGESV.h"
#include "test_IMe_pviDGESV_cs.h"
#include "test_IMe_pviDGESV_ft1.h"

#include "test_ScaLAPACK_pDGESV.h"

#include "test_ScaLAPACK_pDGESV_cp_ft1.h"

#include "tester_head_p.c"

	int nRHS=10;

	versionname[ 0]= "IMe-SV        1 ";
	versionname[ 1]= "IMe-SV        10";
	versionname[ 2]= "IMe-SV-cs     1 ";
	versionname[ 3]= "IMe-SV-cs     10";
	versionname[ 4]= "IMe-SV-ft1/0  1 ";
	versionname[ 5]= "IMe-SV-ft1/0  10";
	versionname[ 6]= "IMe-SV-ft1/1  1 ";
	versionname[ 7]= "IMe-SV-ft1/1  10";
	versionname[ 8]= "SPK-SV        1 ";
	versionname[ 9]= "SPK-SV        10";
	versionname[10]= "SPK-SV-ft1    1 ";
	versionname[11]= "SPK-SV-ft1    10";

	versions = 12;

#include "tester_shoulder_p.c"

    	//////////////////////////////////////////////////////////////////////////////////

     	versionrun[ 0][rep]=test_IMe_pviDGESV(versionname[0], verbose, rows, cols, 1, main_rank, cprocs, sprocs);		// vanilla IMe solve with 1 rhs
     	fpdata(0);

     	versionrun[ 1][rep]=test_IMe_pviDGESV(versionname[1], verbose, rows, cols, nRHS, main_rank, cprocs, sprocs);	// vanilla IMe solve with 10 rhs
 		fpdata(1);

 		versionrun[ 2][rep]=test_IMe_pviDGESV_cs(versionname[2], verbose, rows, cols, 1, main_rank, sprocs);	// checksumming IMe solve with 1 rhs
 		fpdata(2);

 		versionrun[ 3][rep]=test_IMe_pviDGESV_cs(versionname[3], verbose, rows, cols, nRHS, main_rank, sprocs);	// checksumming IMe solve with 10 rhs
 		fpdata(3);

 		versionrun[ 4][rep]=test_IMe_pviDGESV_ft1_sim(versionname[4], verbose, rows, cols, 1, main_rank, sprocs, -1, -1);	// IMe single FT solve with 1 rhs and 0 faults
 		fpdata(4);

 		versionrun[ 5][rep]=test_IMe_pviDGESV_ft1_sim(versionname[5], verbose, rows, cols, nRHS, main_rank, sprocs, -1, -1);// IMe single FT solve with 10 rhs and 0 faults
 		fpdata(5);

    	if (sprocs>0)
    	{
    		versionrun[ 6][rep]=test_IMe_pviDGESV_ft1_sim(versionname[6], verbose, rows, cols, 1, main_rank, sprocs, failing_rank, failing_level);   // IMe single FT solve with 1 rhs
     		fpdata(6);

     		versionrun[ 7][rep]=test_IMe_pviDGESV_ft1_sim(versionname[7], verbose, rows, cols, nRHS, main_rank, sprocs, failing_rank, failing_level);// IMe single FT solve with 10 rhs
     		fpdata(7);
    	}
    	else
    	{
    		versionrun[ 6][rep]=not_run; // don't run IMe single FT solve with 1 rhs
     		fpdata(6);

    		versionrun[ 7][rep]=not_run; // don't run IMe single FT solve with 10 rhs
     		fpdata(7);
    	}

    	versionrun[ 8][rep]=test_ScaLAPACK_pDGESV(versionname[8], verbose, A_init, rows, cols, 1, nb, main_rank, cprocs);		// SPK solve with 1 rhs
 		fpdata(8);

    	versionrun[ 9][rep]=test_ScaLAPACK_pDGESV(versionname[9], verbose, A_init, rows, cols, nRHS, nb, main_rank, cprocs);	// SPK solve with 10 rhs
 		fpdata(9);

    	if (sprocs>0)
    	{
    		versionrun[10][rep]=not_implemented; // not yet SPKmod single FT solve with 1 rhs
     		fpdata(10);

    		versionrun[11][rep]=not_implemented; // not yet SPKmod single FT solve with 10 rhs
     		fpdata(11);
    	}
    	else
    	{
    		versionrun[10][rep]=not_run; // don't run SPKmod single FT solve with 1 rhs
     		fpdata(10);

    		versionrun[11][rep]=not_run; // don't run SPKmod single FT solve with 10 rhs
     		fpdata(11);
    	}

    	//////////////////////////////////////////////////////////////////////////////////

#include "tester_tail_p.c"
