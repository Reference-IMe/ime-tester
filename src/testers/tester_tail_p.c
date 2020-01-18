/*
 * tester_tail_p.c
 *
 *  Created on: Jan 17, 2020
 *      Author: marcello
 */

/*
 * ending part of code for verbatim inclusion to create a code tester for some parallel versions
 */

    	//////////////////////////////////////////////////////////////////////////////////

    	if (main_rank==0)
		{
			for (i=0; i<versions; i++)
			{
				versiontot[i] += versionrun[i][rep];
				if (verbose>0)
				{
					printf("\n%s    call    run time: %10.0f clk", versionname[i], versionrun[i][rep]);
				}
			}
		}
    } // close main loop (see tester_shoulder_p.c)

	if (main_rank==0)
	{
		printf("\n\n Summary:");
		for (i=0; i<versions; i++)
		{
			printf("\n%s    Total   run time: %10.0f clk", versionname[i], versiontot[i]); // in sec. versiontot[i] / CLOCKS_PER_SEC
		}
		printf("\n");
		for (i=0; i<versions; i++)
		{
			printf("\n%s    Average run time: %10.0f clk", versionname[i], (versiontot[i]/repetitions));

			if (file_name_len>0)
			{
				fprintf(fp,"data,%s,%d,%.0f\n",versionname[i],0,(versiontot[i]/repetitions));
			}
		}
		printf("\n");
	}


	// slow down exit
	//sleep(3);
	//MPI_Barrier(MPI_COMM_WORLD);
	//printf("Done %d.\n",main_rank);

	if (file_name_len>0 && main_rank==0)
	{
		fclose(fp);
	}

	MPI_Finalize();
    return(0);
}
