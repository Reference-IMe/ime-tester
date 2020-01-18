
#include "test_IMe_pviDGESV.h"
#include "test_IMe_pviDGESV_cs.h"
#include "test_IMe_pviDGESV_ft1.h"

#include "test_IMe_pviDGEF.h"
#include "test_IMe_pviDGEF_ft1.h"

#include "test_ScaLAPACK_pDGESV.h"
#include "test_ScaLAPACK_pDGETRF.h"
#include "test_ScaLAPACK_pDGEQRF.h"

#include "test_FTLA_pDGEQRF.h"
#include "test_FTLA_pDGETRF.h"

#include "test_ScaLAPACK_pDGESV_cp_ft1.h"
#include "test_ScaLAPACK_pDGETRF_cp_ft1.h"
#include "test_ScaLAPACK_pDGEQRF_cp_ft1.h"

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

	versionname[12]= "IMe-XK          ";
	versionname[13]= "IMe-XK-ft1/0    ";
	versionname[14]= "IMe-XK-ft1/1    ";

	versionname[15]= "SPK-LU          ";
	versionname[16]= "SPK-LU-ft1/0    ";
	versionname[17]= "SPK-LU-ft1/1    ";

	versionname[18]= "FTLA-LU-ft1/0   ";
	versionname[19]= "FTLA-LU-ft1/1   ";

	versionname[20]= "SPK-QR          ";
	versionname[21]= "SPK-QR-ft1/0    ";
	versionname[22]= "SPK-QR-ft1/1    ";

	versionname[23]= "FTLA-QR-ft1/0   ";
	versionname[24]= "FTLA-QR-ft1/1   ";

	versions = 25;

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
    		versionrun[ 6][rep]=-1; // don't run IMe single FT solve with 1 rhs
     		fpdata(6);

    		versionrun[ 7][rep]=-1; // don't run IMe single FT solve with 10 rhs
     		fpdata(7);
    	}

    	versionrun[ 8][rep]=test_ScaLAPACK_pDGESV(versionname[8], verbose, rows, cols, 1, nb, main_rank, cprocs);		// SPK solve with 1 rhs
 		fpdata(8);

    	versionrun[ 9][rep]=test_ScaLAPACK_pDGESV(versionname[9], verbose, rows, cols, nRHS, nb, main_rank, cprocs);	// SPK solve with 10 rhs
 		fpdata(9);

    	if (sprocs>0)
    	{
    		versionrun[10][rep]=-99; // not yet SPKmod single FT solve with 1 rhs
     		fpdata(10);

    		versionrun[11][rep]=-99; // not yet SPKmod single FT solve with 10 rhs
     		fpdata(11);
    	}
    	else
    	{
    		versionrun[10][rep]=-1; // don't run SPKmod single FT solve with 1 rhs
     		fpdata(10);

    		versionrun[11][rep]=-1; // don't run SPKmod single FT solve with 10 rhs
     		fpdata(11);
    	}
     	versionrun[12][rep]=test_IMe_pviDGEF(versionname[12], verbose, rows, cols, main_rank, cprocs, sprocs); // IMe factorization only
 		fpdata(12);

 		versionrun[13][rep]=test_IMe_pviDGEF_ft1_sim(versionname[13], verbose, rows, cols, main_rank, sprocs, -1, -1); // IMe factorization only and 0 faults
 		fpdata(13);

    	if (sprocs>0)
    	{
    		versionrun[14][rep]=test_IMe_pviDGEF_ft1_sim(versionname[14], verbose, rows, cols, main_rank, sprocs, failing_rank, failing_level); // IMe factorization only and 1 fault
    		fpdata(14);
    	}
    	else
    	{
    		versionrun[14][rep]=-1; // don't run IMe factorization only
    		fpdata(14);
    	}

    	versionrun[15][rep]=test_ScaLAPACK_pDGETRF(versionname[15], verbose, rows, cols, nb, main_rank, cprocs); // SPK LU factorization
 		fpdata(15);

		if (sprocs>0)
    	{
    		versionrun[16][rep]=test_ScaLAPACK_pDGETRF_cp_ft1_sim(versionname[16], verbose, rows, cols, nb, main_rank, cprocs, sprocs, -1, checkpoint_skip_interval);	// SPKmod LU factorization single FT and 0 faults
     		fpdata(16);

     		versionrun[17][rep]=test_ScaLAPACK_pDGETRF_cp_ft1_sim(versionname[17], verbose, rows, cols, nb, main_rank, cprocs, sprocs, failing_level, checkpoint_skip_interval);	// SPKmod LU factorization single FT and 1 fault
     		fpdata(17);

     		versionrun[18][rep]=test_FTLA_pDGETRF(versionname[18], verbose, rows, cols, nb, main_rank, cprocs, 0);	// FTLA LU with 0 faults
     		fpdata(18);

     		versionrun[19][rep]=test_FTLA_pDGETRF(versionname[19], verbose, rows, cols, nb, main_rank, cprocs, 1);	// FTLA LU with 1 fault
     		fpdata(19);
    	}
    	else
    	{
    		versionrun[16][rep]=-1;	// don't run SPKmod LU factorization single FT and 0 faults
    		fpdata(16);

    		versionrun[17][rep]=-1;	// don't run SPKmod LU factorization single FT and 1 fault
    		fpdata(17);

    		versionrun[18][rep]=-1;	// don't run FTLA LU with 0 faults
    		fpdata(18);

    		versionrun[19][rep]=-1;	// don't run FTLA LU with 1 fault
        	fpdata(19);
    	}

    	versionrun[20][rep]=test_ScaLAPACK_pDGEQRF(versionname[20], verbose, rows, cols, nb, main_rank, cprocs); // SPK LU factorization
    	fpdata(20);

    	if (sprocs>0)
    	{
    		versionrun[21][rep]=test_ScaLAPACK_pDGEQRF_cp_ft1_sim(versionname[21], verbose, rows, cols, nb, main_rank, cprocs, sprocs, -1, checkpoint_skip_interval);	// SPKmod QR factorization single FT and 0 faults
        	fpdata(21);

        	versionrun[22][rep]=test_ScaLAPACK_pDGEQRF_cp_ft1_sim(versionname[22], verbose, rows, cols, nb, main_rank, cprocs, sprocs, failing_level, checkpoint_skip_interval);	// SPKmod QR factorization single FT and 1 fault
        	fpdata(22);

        	versionrun[23][rep]=test_FTLA_pDGEQRF(versionname[23], verbose, rows, cols, nb, main_rank, cprocs, 0);	// FTLA QR with 0 faults
        	fpdata(23);

        	versionrun[24][rep]=test_FTLA_pDGEQRF(versionname[24], verbose, rows, cols, nb, main_rank, cprocs, 1);	// FTLA QR with 1 fault
        	fpdata(24);
    	}
    	else
    	{
    		versionrun[21][rep]=-1;	// don't run SPKmod QR factorization single FT and 0 faults
        	fpdata(21);

    		versionrun[22][rep]=-1;	// don't run SPKmod QR factorization single FT and 1 fault
        	fpdata(22);

        	versionrun[23][rep]=-1;	// don't run FTLA QR with 0 faults
        	fpdata(23);

        	versionrun[24][rep]=-1;	// don't run FTLA QR with 1 fault
        	fpdata(24);
    	}

    	//////////////////////////////////////////////////////////////////////////////////

#include "tester_tail_p.c"
