
for (i=0; i<versions; i++)
	{
		versiontot[i] = 0;
	}


    for (rep=0; rep<repetitions; rep++)
    {
    	if (main_rank==0 && verbose>0) {printf("\n Run #%d",rep+1);}

    	//////////////////////////////////////////////////////////////////////////////////
