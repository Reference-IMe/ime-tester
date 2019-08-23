#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "../helpers/matrix.h"
#include "../helpers/vector.h"

#include <mpi.h>

#include "Scalapack/ScalapackPDGESV.h"
#include "GaussJordanElimination/GJE-par.h"

#include "../pDGESV_WO.h"
#include "../SLGEWOPV-FT.h"

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    int rank, totprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);   /* get current process id */
    MPI_Comm_size(MPI_COMM_WORLD, &totprocs); /* get number of processes */

    int i,j,k,l,rep;

    double** A2;
    double** T;
    double*  b;
    double**  bb;
    double*  c;
    double*  x;
    double**  xx;

    double** K;
    double*  H;
    //double*  F;

    double*  A1;
    int*     ipiv;

    double versionrun[10][100];
    const char* versionname[10];
    double versiontot[10];

    int n=atoi(argv[1]);

    int rows=n;
    int cols=n;
    int sprocs=atoi(argv[2]);		// number of processes to allocate for summing (0 = no fault tolerance)
    int cprocs=totprocs-sprocs;		// number of processes for real IMe calc
    int repetitions=atoi(argv[3]);
    int verbose=atoi(argv[4]);

    int nRHS=10;

    clock_t start, stop;

	versionname[0]="SPK     ";
	versionname[1]="GJE     ";
	versionname[2]="GJE-cp-0";
	versionname[3]="IMe     ";
	versionname[4]="IMe-cs-0";
	int versions = 5;

	for (i=0; i<versions; i++)
	{
		versiontot[i] = 0;
	}

	if (rank==0 && verbose>0)
	{
		printf("\nMatrix size: %dx%d",n,n);
		printf("\nCheckpoint : ");
		if(sprocs==0)
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
    	if (rank==0 && verbose>0) {printf("\n Run #%d",rep+1);}

    	//////////////////////////////////////////////////////////////////////////////////
    	// ScaLapack

		if (rank==0)
		{
			A1=AllocateMatrix1D(rows, cols);
			b=AllocateMatrix1D(rows,nRHS);
			FillMatrixT1D(A1, rows, cols);
			OneMatrix1D(b, rows, nRHS);

			if (verbose>2)
			{
				printf("\n\n Matrix A:\n");
				PrintMatrix1D(A1, rows, cols);
				printf("\n Vector b:\n");
				PrintMatrix1D(b, rows, nRHS);
			}
		}

		start=clock();

		ScalapackPDGESV_calc(n, A1, nRHS, b, rank, cprocs);

		stop=clock();
		versionrun[0][rep]=(double)(stop - start);

		if (rank==0)
		{
			if (verbose>1)
			{
				printf("\nThe %s solution is:\n",versionname[0]);
				PrintMatrix1D(b, nRHS, rows);
			}
			DeallocateMatrix1D(A1);
			DeallocateMatrix1D(b);
		}

		//////////////////////////////////////////////////////////////////////////////////
		// Gaussian Elimination

	    if (rank==0)
	    {
			A2=AllocateMatrix2D(rows,cols,CONTIGUOUS);
		    xx=AllocateMatrix2D(rows,nRHS,CONTIGUOUS);
			bb=AllocateMatrix2D(rows,nRHS,CONTIGUOUS);
			
	    	FillMatrix2D(A2, rows, cols);
	    	OneMatrix2D(bb, rows, nRHS);

			if (verbose>2)
			{
				printf("\n\n Matrix A:\n");
				PrintMatrix2D(A2, rows, cols);
				printf("\n Vector b:\n");
				PrintMatrix2D(bb, rows, nRHS);
			}
	    }
	    else // avoid problem (bug?) in mpi when calling collectives with undefined pointers in buffers even if unused!
	    {
			A2=AllocateMatrix2D(0,0,CONTIGUOUS);
			bb=AllocateMatrix2D(0,0,CONTIGUOUS);
			xx=AllocateMatrix2D(0,0,CONTIGUOUS);
		}

		start=clock();

		pGaussianElimination_partialmatrix(n, A2, nRHS, bb, rank, cprocs);

	    if (rank==0)
		{
			BackSubstitution(n, A2, nRHS, bb, xx);
			stop=clock();
			versionrun[1][rep]=(double)(stop - start);

			if (verbose>1)
			{
				printf("\nThe %s solution is:\n",versionname[1]);
				PrintMatrix2D(xx, rows, nRHS);
			}
			DeallocateMatrix2D(bb,rows,CONTIGUOUS);
			DeallocateMatrix2D(xx,rows,CONTIGUOUS);
			DeallocateMatrix2D(A2,rows,CONTIGUOUS);
		}
	    else
	    {
			DeallocateMatrix2D(bb,0,CONTIGUOUS);
			DeallocateMatrix2D(xx,0,CONTIGUOUS);
	    	DeallocateMatrix2D(A2,0,CONTIGUOUS); // due to mpi problem
	    }

 
		//////////////////////////////////////////////////////////////////////////////////
		// Gaussian Elimination cp-0

		if (rank==0)
		{
			A2=AllocateMatrix2D(rows,cols,CONTIGUOUS);
		    xx=AllocateMatrix2D(rows,nRHS,CONTIGUOUS);
			bb=AllocateMatrix2D(rows,nRHS,CONTIGUOUS);

	    	FillMatrix2D(A2, rows, cols);
	    	OneMatrix2D(bb, rows, nRHS);

			if (verbose>2)
			{
				printf("\n\n Matrix A:\n");
				PrintMatrix2D(A2, rows, cols);
				printf("\n Vector b:\n");
				PrintMatrix2D(bb, rows, nRHS);
			}
		}
	    else // avoid problem (bug?) in mpi when calling collectives with undefined pointers in buffers even if unused!
	    {
			A2=AllocateMatrix2D(0,0,CONTIGUOUS);
			bb=AllocateMatrix2D(0,0,CONTIGUOUS);
			xx=AllocateMatrix2D(0,0,CONTIGUOUS);
		}

		start=clock();

		pGaussianElimination_partialmatrix_cp(n, A2, nRHS, bb, rank, cprocs, sprocs);

		if (rank==0)
		{
			BackSubstitution(n, A2, nRHS, bb, xx);
			stop=clock();
			versionrun[2][rep]=(double)(stop - start);

			if (verbose>1)
			{
				printf("\nThe %s solution is:\n",versionname[2]);
				PrintMatrix2D(xx, rows, nRHS);
			}
			DeallocateMatrix2D(bb,rows,CONTIGUOUS);
			DeallocateMatrix2D(xx,rows,CONTIGUOUS);
			DeallocateMatrix2D(A2,rows,CONTIGUOUS);
		}
	    else
	    {
			DeallocateMatrix2D(bb,0,CONTIGUOUS);
			DeallocateMatrix2D(xx,0,CONTIGUOUS);
	    	DeallocateMatrix2D(A2,0,CONTIGUOUS); // due to mpi problem
	    }

		//////////////////////////////////////////////////////////////////////////////////
		// Inhibition Method (last)

		xx=AllocateMatrix2D(rows,nRHS,CONTIGUOUS);
		bb=AllocateMatrix2D(rows,nRHS,CONTIGUOUS);

		if (rank==0)
		{
			A2=AllocateMatrix2D(rows,cols,CONTIGUOUS);
			FillMatrix2D(A2, rows, cols);

			OneMatrix2D(bb,rows,nRHS);

			if (verbose>2)
			{
				printf("\n\n Matrix A:\n");
				PrintMatrix2D(A2, rows, cols);
				printf("\n Vector b:\n");
				PrintMatrix2D(bb,rows,nRHS);
			}
		}

		start = clock();

		pDGESV_WO(n, A2, nRHS, bb, xx, rank, cprocs);

		stop = clock();
		versionrun[3][rep]=(double)(stop - start);

		if (rank==0 && verbose>1)
		{
			printf("\nThe %s solution is:\n",versionname[3]);
			PrintMatrix2D(xx,rows,nRHS);
		}

		DeallocateMatrix2D(xx,rows,CONTIGUOUS);
		DeallocateMatrix2D(bb,rows,CONTIGUOUS);
		if (rank==0)
		{
			DeallocateMatrix2D(A2,n,CONTIGUOUS);
		}

		//////////////////////////////////////////////////////////////////////////////////
		// Inhibition Method (last + cs-0)

		x=AllocateVector(n);
		b=AllocateVector(rows);

		if (rank==0)
		{
			A2=AllocateMatrix2D(rows,cols,CONTIGUOUS);

			FillMatrix2D(A2, rows, cols);

			FillVector(b,rows,1);

			if (verbose>2)
			{
				printf("\n\n Matrix A:\n");
				PrintMatrix2D(A2, rows, cols);
				printf("\n Vector b:\n");
				PrintVector(b, rows);
			}
		}

		start = clock();

		//SLGEWOPV_calc_last_cs(A2, b, x, n, rank, cprocs, sprocs, T, Tlocal, TlastKc, TlastKr, h, hh);
		SLGEWOPV_calc_last_cs(A2, b, x, n, rank, cprocs, sprocs);

		stop = clock();
		versionrun[4][rep]=(double)(stop - start);

		if (rank==0 && verbose>1)
		{
			printf("\nThe %s solution is:\n",versionname[4]);
			PrintVector(x, rows);
		}

		DeallocateVector(x);
		DeallocateVector(b);
		if (rank==0)
		{
			DeallocateMatrix2D(A2,n,CONTIGUOUS);
		}

		//////////////////////////////////////////////////////////////////////////////////
		if (rank==0)
		{
			for (i=0; i<versions; i++)
			{
				versiontot[i] += versionrun[i][rep];
				if (verbose>0)
				{
					printf("\n%s    call    run time: %f clk", versionname[i], versionrun[i][rep]);
				}
			}
		}
    }

	if (rank==0)
	{
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
	}

	MPI_Finalize();
    return(0);
}
