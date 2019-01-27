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
	metrics program;

    int i,j,k,l,rep;

    double** A2;
    double*  b;
    double*  c;
    double*  x;

    double** X;
    double** K;
    double*  H;
    double*  F;
    double*  s;

    double*  A1;
    int*     ipiv;

    /*
    double GErunt[100];
    double LUt[100];
    double BSt[100];
    */
    double GEt[100];
    double IMt[100];
    double LSt[100];
    double GEtotrunt=0.0;
    double IMtotrunt=0.0;
    double LStotrunt=0.0;

    int n=atoi(argv[1]);

    int rows=n;
    int cols=n;
    int ft=atoi(argv[2]);
    int repetitions=atoi(argv[3]);
    int verbose=atoi(argv[4]);

    //clock_t begin1, end1, begin2, end2, begin3, end3;
    time_t start1, stop1, start2, stop2, start3, stop3, start4, stop4;

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
    	if (verbose>0) {printf("\n Run #%d",rep+1);}

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

		//begin1 =clock();
		getmetrics(&program);
		start1=program.wall_clock;

		GaussianElimination(A2, b, n);
		BackSubstitution(A2, b, x, n);

		//end2 = clock();
		getmetrics(&program);
		stop2=program.wall_clock;

	    DeallocateMatrix2D(A2,rows,CONTIGUOUS);
	    DeallocateVector(b);
	    DeallocateVector(x);

		if (verbose>1)
		{
			printf("\nThe GE solution is:\n");
			PrintVector(x, rows);
		}

		//////////////////////////////////////////////////////////////////////////////////
		// Inhibition Method

	    A2=AllocateMatrix2D(rows,cols,CONTIGUOUS);
	    b=AllocateVector(rows);

		FillMatrix2D(A2, rows, cols);
		FillVector(b,rows,1);

	    X=AllocateMatrix2D(n,n,CONTIGUOUS);
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

		//begin3 = clock();
		getmetrics(&program);
		start3=program.wall_clock;

		SLGEWOS_calc(A2, b, s, n, X, K, H, F);

		//end3 = clock();
		getmetrics(&program);
		stop3=program.wall_clock;

		if (verbose>1)
		{
			printf("\nThe IMe solution is:\n");
			PrintVector(s, rows);
		}

	    DeallocateMatrix2D(A2,n,CONTIGUOUS);
	    DeallocateMatrix2D(X,n,CONTIGUOUS);
	    DeallocateMatrix2D(K,n,CONTIGUOUS);
	    DeallocateVector(H);
	    DeallocateVector(F);
	    DeallocateVector(s);
	    DeallocateVector(b);

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

		getmetrics(&program);
		start4=program.wall_clock;

		LapackDGESV_calc(A1, b, n, ipiv);

		getmetrics(&program);
		stop4=program.wall_clock;

		if (verbose>1)
		{
			printf("\nThe LPK solution is:\n");
			PrintVector(b, rows);
		}

		DeallocateMatrix1D(A1);
		DeallocateVector(b);
		free(ipiv);

		//////////////////////////////////////////////////////////////////////////////////

		/*
		LUt[rep]=(double)(end1 - begin1) / CLOCKS_PER_SEC;
		BSt[rep]=(double)(end2 - begin2) / CLOCKS_PER_SEC;
		IMt[rep]=(double)(end3 - begin3) / CLOCKS_PER_SEC;
		*/
		/*
		LUt[rep]=(double)(stop1 - start1);
		BSt[rep]=(double)(stop2 - start2);
		*/
		GEt[rep]=(double)(stop2 - start1);
		IMt[rep]=(double)(stop3 - start3);
		LSt[rep]=(double)(stop4 - start4);

		/*
		GErunt[rep]=LUt[rep]+BSt[rep];
		GEtotrunt += GErunt[rep];
		*/
		GEtotrunt += GEt[rep];
		IMtotrunt += IMt[rep];
		LStotrunt += LSt[rep];

		/*
		printf("\n\nLU  decomposition time: %f", LUt[rep]);
		printf("\nBack substitution time: %f", BSt[rep]);
		printf("\nGaussian elimin.  time: %f\n", GErunt[rep]);
		*/
		if (verbose>0)
		{
			printf("\nGE call    run time: %f", GEt[rep]);
			printf("\nIM call    run time: %f", IMt[rep]);
			printf("\nLS call    run time: %f\n", LSt[rep]);
		}
    }
	printf("\n Summary:");
	printf("\nGE Total   run time: %f", GEtotrunt);
	printf("\nIM Total   run time: %f", IMtotrunt);
	printf("\nLS Total   run time: %f\n", LStotrunt);
	printf("\nGE Average run time: %f", GEtotrunt/repetitions);
	printf("\nIM Average run time: %f", IMtotrunt/repetitions);
	printf("\nLS Average run time: %f\n\n", LStotrunt/repetitions);

    return(0);
}
