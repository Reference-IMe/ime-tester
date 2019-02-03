
//#include <stdio.h>
//#include <stdlib.h>
//#include <string.h>
//#include <math.h>
#include <mpi.h>
//#include <cblas.h>
//#include <blas.h>
//#include <scalapack.h>
#include "../../helpers/macros.h"
#include "../../helpers/matrix.h"
#include "../../helpers/vector.h"
#include "../../helpers/Cblacs.h"
#include "../../helpers/scalapack.h"
//#include "../../helpers/selfie.h"

//#define N 10
#define NB 64

void ScalapackPDGESV_calc(double* A_global, double* B_global, int n, int myrank, int nprocs)
{

	// from command line
	int rows;
	int cols;
	rows=cols=n;		// matrix rank
	// general
	int i, j;						//iterators
	int zero = 0, one = 1, six = 6;	//numbers
	// MPI
	int ndims = 2, dims[2] = {0,0};
	// BLACS/SCALAPACK
	int nprow, npcol, info, ic = -1, context, context_global, myrow, mycol;
	int descA_global[9], descB_global[9], descA[9], descB[9];
	char order = 'R', scope = 'A';
	// MATRIX
	//int n = rows;
	int nb = NB, nrhs = 1, nr, nc, lld, lld_global;
	double *A, *B;
	double *work, alpha;
	int *ipiv;
	// MATRIX (global)
	//double *A_global, *B_global;		//defined on every node, but allocated only on root

	// only to equal test conditions
    MPI_Bcast (A_global,n*n,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast (B_global,n,MPI_DOUBLE,0,MPI_COMM_WORLD);

	// Initialize a default BLACS context and the processes grid
	MPI_Dims_create(nprocs, ndims, dims);
	nprow = dims[0];
	npcol = dims[1];
	Cblacs_get( ic, zero, &context );
	Cblacs_gridinit( &context, &order, nprow, npcol );
	Cblacs_get( ic, zero, &context_global );
	Cblacs_gridinit( &context_global, &order, one, one );
	//blacs_gridinfo_( &context_global, &nprow, &npcol, &myrow, &mycol );
	//printf("glob: %d (%d,%d)\n",myrank,myrow,mycol);
	Cblacs_gridinfo( context, &nprow, &npcol, &myrow, &mycol );
	//printf("loca: %d (%d,%d)\n",myrank,myrow,mycol);
	

	// Computation of local matrix size
	nr = numroc_( &n, &nb, &myrow, &zero, &nprow );
	nc = numroc_( &n, &nb, &mycol, &zero, &npcol );
	lld = MAX( 1 , nr );
	A = malloc(nr*nc*sizeof(double));
	B = malloc(nr*sizeof(double));
	work = malloc(nb*sizeof(double));
	ipiv = malloc((lld+nb)*sizeof(int));

	// Descriptors (local)
	descinit_( descB, &n, &one, &nb, &nb, &zero, &zero, &context, &lld, &info );
	descinit_( descA, &n, &n, &nb, &nb, &zero, &zero, &context, &lld, &info );

	if (myrank==0)
	{
		//printf("\nallocating memory..\n");
		//A_global=AllocateMatrix1D(rows,cols);
		//B_global=AllocateVector(rows);

		// Descriptors (global)
		lld_global = n;
		descinit_( descA_global, &n, &n, &one, &one, &zero, &zero, &context_global, &lld_global, &info );
		descinit_( descB_global, &n, &one, &one, &one, &zero, &zero, &context_global, &lld_global, &info );

		//printf("filling matrices..\n\n");

		//FillMatrixT1D(A_global,rows,cols);
		/*
		for( i = 1; i <= n; i++ )
		{
			for( j = 1; j <= n; j++ )
			{
				if( i == j )
					A_global[(i-1)*n+(j-1)] = (double)(10000);
				else
					A_global[(i-1)*n+(j-1)] = ((double)(i)+(double)(j)/(double)(2));
			 }
		}*/
/*
		printf("[A]\n");
#ifdef SCREENDUMP
		Print1DMatrix(A_global,rows,cols);
#endif

		printf("[b]\n"); //the vector of causes will be overwritten by the vector of auxiliary cause
		for(i=0; i<rows; i++)
			//B_global[i]=207-i-1;
			B_global[i]=(double)(1);
		PrintVector(B_global,i);
		*/
	}
	else
	{
		// Descriptors (global, for non-root nodes)
		for (i=0; i<9; i++)
		{
			descA_global[i]=0;
			descB_global[i]=0;
		}
		descA_global[1]=-1;
		descB_global[1]=-1;
	}

	// spread matrices
	pdgemr2d_(&n,&n,A_global,&one,&one,descA_global,A,&one,&one,descA, &context);
	pdgemr2d_(&n,&one,B_global,&one,&one,descB_global,B,&one,&one,descB, &context);

	/*
	printf("I'm %d (%d,%d)\n",myrank,myrow,mycol);
	Print1DMatrix(A,nr,nc);
	PrintVector(B,nr);
	printf("I was %d\n",myrank);
	*/


	// Linear system equations solver
	pdgesv_( &n, &one, A, &one, &one, descA, ipiv, B, &one, &one, descB, &info );
	pdgemr2d_(&n,&one,B,&one,&one,descB,B_global,&one,&one,descB_global, &context);

	/*
	if (myrank==0)
	{


		printf("\n[x]\n");

		PrintVector(B_global,rows);


		DeallocateMatrix1D(A_global);
		DeallocateVector(B_global);
	}
*/
	//pdlaprnt_( &n, &one, B, &one, &one, descB, &zero, &zero, &scope, &six, work );	// Print the results vector

	free(A);
	free(B);
	free(ipiv);
	free(work);



// Close BLACS environment
  Cblacs_gridexit( context );
  //Cblacs_gridexit( context_global );
  //Cblacs_exit( zero );
}
