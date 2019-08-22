#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "../helpers/types.h"
#include "../helpers/selfie.h"
#include "../helpers/matrix.h"
#include "../helpers/vector.h"
#include "../DGESV_WO_1D.h"
#include "../DGESV_WO.h"
#include "GaussJordanElimination/GJE-seq.h"
#include "Lapack/LapackDGESV.h"

double testLapackDGESV(const char* label, cui verbosity, cui rows, cui cols, cui nrhs)
{
	clock_t start, stop;
	double* A;
	double* bb;
	A=AllocateMatrix1D(rows, cols);
	bb=AllocateMatrix1D(rows, nrhs);

	FillMatrixT1D(A, rows, cols);	// input has to be transposed
	//TODO: fill normally and transpose inside the timed section

	OneMatrix1D(bb, rows, nrhs);

	if (verbosity>2)
	{
		printf("\n\n Matrix A:\n");
		PrintMatrix1D(A, rows, cols);
		printf("\n Vector b:\n");
		PrintMatrix1D(bb, rows, nrhs);
	}

	start=clock();

	LapackDGESV(rows, A, nrhs, bb);

	stop=clock();

	if (verbosity>1)
	{
		printf("\nThe %s solution is:\n",label);
		PrintMatrix1D(bb, nrhs, rows);	// result is transposed
		//TODO: transpose the result inside the timed section
	}

	DeallocateMatrix1D(A);
	DeallocateMatrix1D(bb);

	return (double)(stop - start);
}

double testGaussianElimination(const char* label, cui verbosity, cui rows, cui cols, cui nrhs)
{
	clock_t start, stop;
	double** A;
	double** bb;
	double** xx;
	ui i,j;
	A=AllocateMatrix2D(rows,cols,CONTIGUOUS);
	bb=AllocateMatrix2D(rows,nrhs,CONTIGUOUS);
	xx=AllocateMatrix2D(rows,nrhs,CONTIGUOUS);

	FillMatrix2D(A, rows, cols);
	OneMatrix2D(bb, rows, nrhs);

	if (verbosity>2)
	{
		printf("\n\n Matrix A:\n");
		PrintMatrix2D(A, rows, cols);
		printf("\n Vector b:\n");
		PrintMatrix2D(bb, rows, nrhs);
	}

	start =clock();

	GaussianElimination(rows, A, nrhs, bb);
	BackSubstitution(rows, A, nrhs, bb, xx);

	stop = clock();

	if (verbosity>1)
	{
		printf("\nThe %s solution is:\n",label);
		PrintMatrix2D(xx, rows, nrhs);
	}

	DeallocateMatrix2D(A,rows,CONTIGUOUS);
	DeallocateMatrix2D(bb,rows,CONTIGUOUS);
	DeallocateMatrix2D(xx,rows,CONTIGUOUS);

	return (double)(stop - start);
}

double testIMe1D(const char* label, cui verbosity, cui rows, cui cols, cui nrhs)
{
	clock_t start, stop;
	double* A;
	double* bb;
	double* xx;
	A=AllocateMatrix1D(rows, cols);
	bb=AllocateMatrix1D(rows, nrhs);
	xx=AllocateMatrix1D(rows, nrhs);

	FillMatrix1D(A, rows, cols);
	OneMatrix1D(bb, rows, nrhs);

	if (verbosity>2)
	{
		printf("\n\n Matrix A:\n");
		PrintMatrix1D(A, rows, cols);
		printf("\n Vector b:\n");
		PrintMatrix1D(bb, rows, nrhs);
	}

	start = clock();

	DGESV_WO_1D(rows, A, nrhs, bb, xx);

	stop = clock();

	if (verbosity>1)
	{
		printf("\nThe %s solution is:\n",label);
		PrintMatrix1D(xx, rows, nrhs);
	}

	DeallocateMatrix1D(A);
	DeallocateMatrix1D(bb);
	DeallocateVector(xx);

	return (double)(stop - start);
}

double testIMe(const char* label, cui verbosity, cui rows, cui cols, cui nrhs)
{
	clock_t start, stop;
	double** A;
	double** bb;
	double** xx;
	ui i,j;
	A=AllocateMatrix2D(rows, cols, CONTIGUOUS);
	bb=AllocateMatrix2D(rows, nrhs, CONTIGUOUS);
	xx=AllocateMatrix2D(rows, nrhs, CONTIGUOUS);

	FillMatrix2D(A, rows, cols);
	OneMatrix2D(bb, rows, nrhs);

	if (verbosity>2)
	{
		printf("\n\n Matrix A:\n");
		PrintMatrix2D(A, rows, cols);
		printf("\n Vector b:\n");
		PrintMatrix2D(bb, rows, nrhs);
	}

	start = clock();

	DGESV_WO(rows, A, nrhs, bb, xx);

	stop = clock();

	if (verbosity>1)
	{
		printf("\nThe %s solution is:\n",label);
		PrintMatrix2D(xx, rows, nrhs);
	}

	DeallocateMatrix2D(A,rows,CONTIGUOUS);
	DeallocateMatrix2D(bb,rows,CONTIGUOUS);
	DeallocateMatrix2D(xx,rows,CONTIGUOUS);

	return (double)(stop - start);
}

int main(int argc, char **argv)
{
    ui i,rep;

    ui n=atoi(argv[1]);
    ui rows=n;
    ui cols=n;

    ui ft=atoi(argv[2]);

    ui repetitions=atoi(argv[3]);

    ui verbose=atoi(argv[4]);

    ui nRHS=10;

    double versionrun[10][100];
    double versiontot[10];
    const char* versionname[10];
	versionname[0]="LPK   1 ";
	versionname[1]="LPK   10";
	versionname[2]="GJE   1 ";
	versionname[3]="GJE   10";
	versionname[4]="IMe1D 1 ";
	versionname[5]="IMe1D 10";
	versionname[6]="IMe   1 ";
	versionname[7]="IMe   10";
	cui versions = 8;

	for (i=0; i<versions; i++)
	{
		versiontot[i] = 0;
	}

	if (verbose>0)
	{
		printf("\nMatrix size: %dx%d",n,n);
		printf("\nCheckpoint : ");
		if(ft==0)
		{
			printf("no\n");
		}
		else
		{
			printf("yes\n");
		}
	}

    for (rep=0; rep<repetitions; rep++)
    {
    	if (verbose>0) {printf("\n\n Run #%d",rep+1);}

    	//////////////////////////////////////////////////////////////////////////////////

    	versionrun[0][rep]=testLapackDGESV(versionname[0], verbose, rows, cols, 1);    			// Lapack with 1 rhs
    	versionrun[1][rep]=testLapackDGESV(versionname[1], verbose, rows, cols, nRHS);			// Lapack with 10 rhs
    	versionrun[2][rep]=testGaussianElimination(versionname[2], verbose, rows, cols, 1);		// Gaussian Elimination with 1 rhs
    	versionrun[3][rep]=testGaussianElimination(versionname[3], verbose, rows, cols, nRHS);	// Gaussian Elimination with 10 rhs
    	versionrun[4][rep]=testIMe1D(versionname[4], verbose, rows, cols, 1);					// IMe-1D with 1 rhs
    	versionrun[5][rep]=testIMe1D(versionname[5], verbose, rows, cols, nRHS);				// IMe-1D with 10 rhs
    	versionrun[6][rep]=testIMe(versionname[6], verbose, rows, cols, 1);						// IMe with 1 rhs
    	versionrun[7][rep]=testIMe(versionname[7], verbose, rows, cols, nRHS);					// IMe with 10 rhs

		//////////////////////////////////////////////////////////////////////////////////

		for (i=0; i<versions; i++)
		{
			versiontot[i] += versionrun[i][rep];
			if (verbose>0)
			{
				printf("\n%s    call    run time: %f clk", versionname[i], versionrun[i][rep]);
			}
		}
    }
	printf("\n\n Summary:");
	for (i=0; i<versions; i++)
	{
		printf("\n%s    Total   run time: %f clk\t\t%f s", versionname[i], versiontot[i], versiontot[i] / CLOCKS_PER_SEC);
	}
	printf("\n");
	for (i=0; i<versions; i++)
		{
			printf("\n%s    Average run time: %f clk\t\t%f s", versionname[i], versiontot[i]/repetitions, versiontot[i]/repetitions / CLOCKS_PER_SEC);
		}
	printf("\n");
    return(0);
}
