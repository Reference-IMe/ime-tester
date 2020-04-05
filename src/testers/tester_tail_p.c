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
				versiontot[i].total_time += versionrun[i][rep].total_time;
				versiontot[i].core_time  += versionrun[i][rep].core_time;
				versiontot[i].norm_rel_err  += versionrun[i][rep].norm_rel_err;
				if (verbose>0)
				{
					printf("\n%s    call    run time: %10.0f (%.0f)\ts\t nre: %f", versionname[i], versionrun[i][rep].total_time, versionrun[i][rep].core_time, versionrun[i][rep].norm_rel_err);
				}
			}
		}
    } // close main loop (see tester_shoulder_p.c)

	if (main_rank==0)
	{
		printf("\n\n Summary:");
		for (i=0; i<versions; i++)
		{
			printf("\n%s    Total   run time: %10.0f (%.0f)\ts", versionname[i], versiontot[i].total_time, versiontot[i].core_time); // in sec. versiontot[i] / CLOCKS_PER_SEC
		}
		printf("\n");
		for (i=0; i<versions; i++)
		{
			printf("\n%s    Average run time: %10.0f (%.0f)\ts\t nre: %f", versionname[i], versiontot[i].total_time/repetitions, versiontot[i].core_time/repetitions, versiontot[i].norm_rel_err/repetitions);

			if (file_name_len>0)
			{
				fprintf(fp,"data,%s,%d,%d,%.0f,%f\n",versionname[i], 0, versionrun[i][0].exit_code, versiontot[i].total_time/repetitions,versiontot[i].norm_rel_err/repetitions);
				fprintf(fp,"data,%s%s,%d,%d,%.0f,%f\n",versionname[i], "(core)", 0, versionrun[i][0].exit_code, versiontot[i].core_time/repetitions,versiontot[i].norm_rel_err/repetitions);
			}
		}
		printf("\n");
		for (i=0; i<versions; i++)
		{
			dmedian( versionrun[i], repetitions);
			printf("\n%s    Median  run time: %10.0f (%.0f)\ts\t nre: %f", versionname[i], versionrun[i][repetitions/2].total_time, versionrun[i][repetitions/2].core_time, versionrun[i][repetitions/2].norm_rel_err);

			if (file_name_len>0)
			{
				fprintf(fp,"data,%s,%d,%d,%.0f,%f\n",versionname[i], -1, versionrun[i][0].exit_code, versionrun[i][repetitions/2].total_time, versionrun[i][repetitions/2].norm_rel_err );
				fprintf(fp,"data,%s%s,%d,%d,%.0f,%f\n",versionname[i], "(core)",-1, versionrun[i][0].exit_code, versionrun[i][repetitions/2].core_time, versionrun[i][repetitions/2].norm_rel_err );
			}
		}
		printf("\n");
	}


	// slow down exit
	//sleep(3);
	//MPI_Barrier(MPI_COMM_WORLD);
	//printf("Done %d.\n",main_rank);

	if (main_rank==0)
	{
		if (file_name_len>0)
		{
			fclose(fp);
		}
		DeallocateMatrix1D(A_ref);
		DeallocateVector(b_ref);
		DeallocateVector(x_ref);
	}

	MPI_Finalize();
    return(0);
}
