#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "../helpers/matrix.h"
#include "../helpers/vector.h"

#include "test_GaussianElimination_GE.h"
#include "test_IMe_DGESV.h"
#include "test_LAPACK_DGESV.h"


int main(int argc, char **argv)
{
    int i,rep;

    double versionrun[20][100];
    const char* versionname[20];
    double versiontot[20];

    int n;		// matrix size
    int cnd;	// condition number
    int seed;	// seed for random matrix generation
    int file_name_len;
    char* file_name;
    int rows;
    int cols;

    int repetitions;
    int verbose;

    /*
     * default values
     */
    n=8;
    cnd=1;
    seed=1;
    verbose=1;			// minimum verbosity
    repetitions=1;
    file_name_len=0;	// no output to file

    /*
     * read command line parameters
     */
	for( i = 1; i < argc; i++ )
	{
		if( strcmp( argv[i], "-n" ) == 0 ) {
			n = atoi(argv[i+1]);
			i++;
		}
		if( strcmp( argv[i], "-cnd" ) == 0 ) {
			cnd = atoi(argv[i+1]);
			i++;
		}
		if( strcmp( argv[i], "-seed" ) == 0 ) {
			seed = atoi(argv[i+1]);
			i++;
		}
		if( strcmp( argv[i], "-v" ) == 0 ) {
			verbose = atoi(argv[i+1]);
			i++;
		}
		if( strcmp( argv[i], "-r" ) == 0 ) {
			repetitions = atoi(argv[i+1]);
			i++;
		}
		if( strcmp( argv[i], "-o" ) == 0 ) {
			file_name_len = strlen(argv[i+1]);
			file_name = malloc(file_name_len+1);
			strcpy(file_name, argv[i+1]);
			i++;
		}
	}

	rows=n;
    cols=n;

    if (verbose>0)
    {
		printf("     Matrix size:                   %dx%d\n",rows,cols);
		printf("     Matrix condition number:       %d\n",cnd);
		printf("     Matrix random generation seed: %d\n",seed);

		printf("     Run repetitions:               %d\n",repetitions);

		if (file_name_len>0)
		{
			printf("     Output file:                   %s\n",file_name);
		}
		else
		{
			printf("WRN: No output to file\n");
		}
    }

    int nRHS=10;
	versionname[0]="LPK   1 ";
	versionname[1]="LPK   10";
	versionname[2]="GE    1 ";
	versionname[3]="GE    10";
	versionname[4]="IMe   1 ";
	versionname[5]="IMe   10";
	int versions = 6;

	for (i=0; i<versions; i++)
	{
		versiontot[i] = 0;
	}


    for (rep=0; rep<repetitions; rep++)
    {
    	if (verbose>0) {printf("\n Run #%d",rep+1);}

    	//////////////////////////////////////////////////////////////////////////////////

    	versionrun[0][rep]=test_Lapack_DGESV(versionname[0], verbose, rows, cols, 1);    			// Lapack with 1 rhs
    	versionrun[1][rep]=test_Lapack_DGESV(versionname[1], verbose, rows, cols, nRHS);			// Lapack with 10 rhs
    	versionrun[2][rep]=test_GaussianElimination(versionname[2], verbose, rows, cols, 1);		// Gaussian Elimination with 1 rhs
    	versionrun[3][rep]=test_GaussianElimination(versionname[3], verbose, rows, cols, nRHS);		// Gaussian Elimination with 10 rhs
    	versionrun[4][rep]=test_IMe_DGESV(versionname[4], verbose, rows, cols, 1);					// IMe with 1 rhs
    	versionrun[5][rep]=test_IMe_DGESV(versionname[5], verbose, rows, cols, nRHS);				// IMe with 10 rhs

		//////////////////////////////////////////////////////////////////////////////////

		for (i=0; i<versions; i++)
		{
			versiontot[i] += versionrun[i][rep];
			if (verbose>0)
			{
				printf("\n%s    call    run time: %10.0f clk", versionname[i], versionrun[i][rep]);
			}
		}
    }

	printf("\n\n Summary:");
	for (i=0; i<versions; i++)
	{
		printf("\n%s    Total   run time: %10.0f clk", versionname[i], versiontot[i]); // in sec. versiontot[i] / CLOCKS_PER_SEC
	}
	printf("\n");
	for (i=0; i<versions; i++)
	{
		printf("\n%s    Average run time: %10.0f clk", versionname[i], (versiontot[i]/repetitions));
	}
	printf("\n");

    return(0);
}
