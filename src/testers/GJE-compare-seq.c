#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "../helpers/selfie.h"
#include "../helpers/matrix.h"
#include "../helpers/vector.h"
#include "../SLGEWOS.h"
#include "GaussJordanElimination/GJE-seq.h"

int main(int argc, char **argv)
{
	metrics program;

    int i,j,k,l,rep;

    double** A;
    double*  b;
    double*  c;
    double*  x;

    /*
    double GErunt[100];
    double LUt[100];
    double BSt[100];
    */
    double GEt[100];
    double IMt[100];
    double GEtotrunt=0.0;
    double IMtotrunt=0.0;

    int n=atoi(argv[1]);

    int rows=n;
    int cols=n;
    int ft=atoi(argv[2]);
    int repetitions=atoi(argv[3]);

    A=AllocateMatrix2D(rows,cols,CONTIGUOUS);
    b=AllocateVector(rows);
    x=AllocateVector(rows);

    //clock_t begin1, end1, begin2, end2, begin3, end3;
    time_t start1, stop1, start2, stop2, start3, stop3;

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

    for (rep=0; rep<repetitions; rep++)
    {
		printf("\n Run #%d",rep+1);

		FillMatrix2D(A, rows, cols);
		FillVector(b,rows,1);

		/*
		printf("\n\n Matrix A:\n");
		PrintMatrix2D(A, rows, cols);
		printf("\n Vector b:\n");
		PrintVector(b, rows);
		*/

		//////////////////////////////////////////////////////////////////////////////////

		//begin1 =clock();
		getmetrics(&program);
		start1=program.wall_clock;

		GaussianElimination(A, b, n);

		//end1 = clock();
		/*
		getmetrics(&program);
		stop1=program.wall_clock;

		//////////////////////////////////////////////////////////////////////////////////

		//begin2 =clock();
		getmetrics(&program);
		start2=program.wall_clock;
		*/
		BackSubstitution(A, b, x, n);

		//end2 = clock();
		getmetrics(&program);
		stop2=program.wall_clock;

		//////////////////////////////////////////////////////////////////////////////////

		/*
		printf("\nThe GE solution is:\n");
		PrintVector(x, rows);
		*/

		//////////////////////////////////////////////////////////////////////////////////

		FillMatrix2D(A, rows, cols);
		FillVector(b,rows,1);

		/*
		printf("\n\n Matrix A:\n");
		PrintMatrix2D(A, rows, cols);
		printf("\n Vector b:\n");
		PrintVector(b, rows);
		*/

		//begin3 = clock();
		getmetrics(&program);
		start3=program.wall_clock;

		SLGEWOS(A, b, x, n);

		//end3 = clock();
		getmetrics(&program);
		stop3=program.wall_clock;

		//////////////////////////////////////////////////////////////////////////////////

		/*
		printf("\nThe IMe solution is:\n");
		PrintVector(x, rows);
		*/

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

		/*
		GErunt[rep]=LUt[rep]+BSt[rep];
		GEtotrunt += GErunt[rep];
		*/
		GEtotrunt += GEt[rep];
		IMtotrunt += IMt[rep];

		/*
		printf("\n\nLU  decomposition time: %f", LUt[rep]);
		printf("\nBack substitution time: %f", BSt[rep]);
		printf("\nGaussian elimin.  time: %f\n", GErunt[rep]);
		*/
		printf("\nGE call    run time: %f", GEt[rep]);
		printf("\nIM call    run time: %f\n", IMt[rep]);
    }
	printf("\n Summary:");
	printf("\nGE Total   run time: %f", GEtotrunt);
	printf("\nIM Total   run time: %f", IMtotrunt);
	printf("\nGE Average run time: %f", GEtotrunt/repetitions);
	printf("\nIM Average run time: %f\n\n", IMtotrunt/repetitions);

    DeallocateMatrix2D(A,rows,CONTIGUOUS);
    DeallocateVector(b);
    DeallocateVector(x);

    return(0);
}
