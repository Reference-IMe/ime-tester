#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "../helpers/selfie.h"
#include "../helpers/matrix.h"
#include "../helpers/vector.h"
#include "../SLGEWOS.h"
#include "GaussJordanElimination/GJE-seq.h"
#include "Lapack/LapackDGESV.h"

int main(int argc, char **argv)
{
    int i,j,k,l,rep;

    double** A2;
    double*  b;
    double*  c;
    double*  x;

	double h,hh;

    double** T;
    double* T1;

    double** K;
    double*  H;
    double*  F;
    double*  s;

    double*  A1;
    int*     ipiv;

    double versionrun[10][100];
    const char* versionname[10];
    double versiontot[10];

    int n=atoi(argv[1]);

    int rows=n;
    int cols=n;
    int ft=atoi(argv[2]);
    int repetitions=atoi(argv[3]);
    int verbose=atoi(argv[4]);

    clock_t start, stop;

	versionname[0]="LPK     ";
	versionname[1]="GJE     ";
	versionname[2]="IMe-naif";
	versionname[3]="IMe-uwnd";
	versionname[4]="IMe-iopt";
	versionname[5]="IMe-1D  ";
	int versions = 6;

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
    	// Lapack

		A1=AllocateMatrix1D(rows, cols);
		b=AllocateVector(rows);
		ipiv = malloc(n * sizeof(int));

		FillMatrix1D(A1, rows, cols);
		FillVector(b,rows,1);

		if (verbose>2)
		{
			printf("\n\n Matrix A:\n");
			PrintMatrix1D(A1, rows, cols);
			printf("\n Vector b:\n");
			PrintVector(b, rows);
		}

		start=clock();

		LapackDGESV_calc(A1, b, n, ipiv);

		stop=clock();

		versionrun[0][rep]=(double)(stop - start);

		if (verbose>1)
		{
			printf("\nThe %s solution is:\n",versionname[0]);
			PrintVector(b, rows);
		}

		DeallocateMatrix1D(A1);
		DeallocateVector(b);
		free(ipiv);

		//////////////////////////////////////////////////////////////////////////////////
		// Gaussian Elimination

	    A2=AllocateMatrix2D(rows,cols,CONTIGUOUS);
	    b=AllocateVector(rows);
	    x=AllocateVector(rows);

		FillMatrix2D(A2, rows, cols);
		FillVector(b,rows,1);

		if (verbose>2)
		{
			printf("\n\n Matrix A:\n");
			PrintMatrix2D(A2, rows, cols);
			printf("\n Vector b:\n");
			PrintVector(b, rows);
		}

		start =clock();

		GaussianElimination(A2, b, n);
		BackSubstitution(A2, b, x, n);

		stop = clock();

		versionrun[1][rep]=(double)(stop - start);

		if (verbose>1)
		{
			printf("\nThe %s solution is:\n",versionname[1]);
			PrintVector(x, rows);
		}

	    DeallocateMatrix2D(A2,rows,CONTIGUOUS);
	    DeallocateVector(b);
	    DeallocateVector(x);

		//////////////////////////////////////////////////////////////////////////////////
		// Inhibition Method

	    A2=AllocateMatrix2D(rows,cols,CONTIGUOUS);
	    b=AllocateVector(rows);

		FillMatrix2D(A2, rows, cols);
		FillVector(b,rows,1);

	    K=AllocateMatrix2D(n,n,CONTIGUOUS);
	    H=AllocateVector(n);
	    F=AllocateVector(n);
	    s=AllocateVector(n);

		if (verbose>2)
		{
			printf("\n\n Matrix A:\n");
			PrintMatrix2D(A2, rows, cols);
			printf("\n Vector b:\n");
			PrintVector(b, rows);
		}

		start = clock();

		SLGEWOS_calc(A2, b, s, n, K, H, F);

		stop = clock();

		versionrun[2][rep]=(double)(stop - start);

		if (verbose>1)
		{
			printf("\nThe %s solution is:\n",versionname[2]);
			PrintVector(s, rows);
		}

	    DeallocateMatrix2D(A2,n,CONTIGUOUS);
	    DeallocateMatrix2D(K,n,CONTIGUOUS);
	    DeallocateVector(H);
	    DeallocateVector(F);
	    DeallocateVector(s);
	    DeallocateVector(b);

	    //////////////////////////////////////////////////////////////////////////////////
		// Inhibition Method (loops unwound)

		A2=AllocateMatrix2D(rows,cols,CONTIGUOUS);
		b=AllocateVector(rows);

		FillMatrix2D(A2, rows, cols);
		FillVector(b,rows,1);

		K=AllocateMatrix2D(n,n,CONTIGUOUS);
		H=AllocateVector(n);
		F=AllocateVector(n);
		s=AllocateVector(n);

		if (verbose>2)
		{
			printf("\n\n Matrix A:\n");
			PrintMatrix2D(A2, rows, cols);
			printf("\n Vector b:\n");
			PrintVector(b, rows);
		}

		start = clock();

		SLGEWOS_calc_unwind(A2, b, s, n, K, H, F);

		stop = clock();

		versionrun[3][rep]=(double)(stop - start);

		if (verbose>1)
		{
			printf("\nThe %s solution is:\n",versionname[3]);
			PrintVector(s, rows);
		}

		DeallocateMatrix2D(A2,n,CONTIGUOUS);
		DeallocateMatrix2D(K,n,CONTIGUOUS);
		DeallocateVector(H);
		DeallocateVector(F);
		DeallocateVector(s);
		DeallocateVector(b);

	    //////////////////////////////////////////////////////////////////////////////////
		// Inhibition Method (last, optimized init)

		A2=AllocateMatrix2D(rows,cols,CONTIGUOUS);
		b=AllocateVector(rows);
		T=AllocateMatrix2D(rows,cols*2,CONTIGUOUS);

		FillMatrix2D(A2, rows, cols);
		FillVector(b,rows,1);

		x=AllocateVector(n);

		if (verbose>2)
		{
			printf("\n\n Matrix A:\n");
			PrintMatrix2D(A2, rows, cols);
			printf("\n Vector b:\n");
			PrintVector(b, rows);
		}

		start = clock();

		SLGEWOS_calc_initopt(A2, b, T, x, n, &h, &hh);;

		stop = clock();

		versionrun[4][rep]=(double)(stop - start);

		if (verbose>1)
		{
			printf("\nThe %s solution is:\n",versionname[4]);
			PrintVector(x, rows);
		}

		DeallocateMatrix2D(A2,n,CONTIGUOUS);
		DeallocateMatrix2D(T,n,CONTIGUOUS);
		DeallocateVector(x);
		DeallocateVector(b);

		 //////////////////////////////////////////////////////////////////////////////////
		// Inhibition Method (last, 1D)

		A1=AllocateMatrix1D(rows,cols);
		b=AllocateVector(rows);
		T1=AllocateMatrix1D(rows,cols*2);

		FillMatrix1D(A1, rows, cols);
		FillVector(b,rows,1);

		x=AllocateVector(n);

		if (verbose>2)
		{
			printf("\n\n Matrix A:\n");
			PrintMatrix1D(A1, rows, cols);
			printf("\n Vector b:\n");
			PrintVector(b, rows);
		}

		start = clock();

		SLGEWOS_calc_last(A1, b, T1, x, n, &h, &hh);

		stop = clock();

		versionrun[5][rep]=(double)(stop - start);

		if (verbose>1)
		{
			printf("\nThe %s solution is:\n",versionname[5]);
			PrintVector(x, rows);
		}

		DeallocateMatrix1D(A1);
		DeallocateMatrix1D(T1);
		DeallocateVector(x);
		DeallocateVector(b);

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
