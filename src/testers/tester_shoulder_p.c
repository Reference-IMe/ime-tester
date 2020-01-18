/*
 * tester_shoulder_p.c
 *
 *  Created on: Jan 17, 2020
 *      Author: marcello
 */

/*
 * mid part of code for verbatim inclusion to create a code tester for some parallel versions
 */

for (i=0; i<versions; i++)
	{
		versiontot[i] = 0;
	}


    for (rep=0; rep<repetitions; rep++) // main loop (see tester_tail_p.c for closing the loop)
    {
    	if (main_rank==0 && verbose>0) {printf("\n Run #%d",rep+1);}

    	//////////////////////////////////////////////////////////////////////////////////
