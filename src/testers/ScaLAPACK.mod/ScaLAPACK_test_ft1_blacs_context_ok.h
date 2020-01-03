/*
 * ScaLAPACK_pDGETRF_cp_ft1_sim.h
 *
 *  Created on: Dec 28, 2019
 *      Author: marcello
 */

#include <mpi.h>
#include "../../helpers/macros.h"
#include "../../helpers/matrix.h"
#include "../../helpers/vector.h"
#include "../../helpers/Cblacs.h"
#include "../../helpers/scalapack.h"

// TODO: checkpointing

void ScaLAPACK_pDGETRF_cp_ft1_sim(int n, double* A_global, int mpi_rank, int cprocs, int sprocs, int failing_rank, int failing_level)
{
	/*
	 * n = system rank (A_global n x n)
	 */

	// general
	int i, j;				//iterators
	int zero = 0, one = 1;	//numbers
	int nprocs = cprocs + sprocs;
	// MPI
	int ndims = 2, dims[2] = {0,0};
	// BLACS/SCALAPACK
	int nprow, npcol, info, ic = -1, context, context_global, context_cp_only, context_all, myrow, mycol;
	int tmprow,tmpcol,tmpmyrow, tmpmycol;
	int descA_global[9], descA[9], descAcp_cp_only[9];
	//int descB_global[9], descB[9];
	char order = 'R', scope = 'A';
	// MATRIX

	int nb, nr, nc, ncrhs, lld, lld_global;
	double *A;
	double *Acp_cp_only;
	//double *B;
	//double *work, alpha;
	int *ipiv;

	// Initialize a default BLACS context and the processes grid
	MPI_Dims_create(cprocs, ndims, dims);
	nprow = dims[0];
	npcol = dims[1];
	Cblacs_get( ic, zero, &context );
	Cblacs_gridinit( &context, &order, nprow, npcol );
	Cblacs_gridinfo( context, &nprow, &npcol, &myrow, &mycol );
	MPI_Barrier(MPI_COMM_WORLD);
	printf("context: %d in %dx%d id %d,%d\n",mpi_rank,nprow,npcol,myrow,mycol);

	/*
	Cblacs_get( ic, zero, &context_global );
	Cblacs_gridinit( &context_global, &order, one, one );
	Cblacs_gridinfo( context_global, &tmprow, &tmpcol, &tmpmyrow, &tmpmycol );
	MPI_Barrier(MPI_COMM_WORLD);
	printf("context_global: %d in %dx%d id %d,%d\n",mpi_rank,tmprow,tmpcol,tmpmyrow,tmpmycol);
	*/

	Cblacs_get( ic, zero, &context_all );
	Cblacs_gridinit( &context_all, &order, one, nprocs );


	Cblacs_get( ic, zero, &context_cp_only );
	//Cblacs_gridinit( &context_cp_only, &order, one, one );
	int map_cp[1][1];
	map_cp[0][0]=cprocs;
	Cblacs_gridmap( &context_cp_only, map_cp, one, one, one);
	Cblacs_gridinfo( context_cp_only, &tmprow, &tmpcol, &tmpmyrow, &tmpmycol );
	MPI_Barrier(MPI_COMM_WORLD);
	printf("context_cp: %d in %dx%d id %d,%d\n",mpi_rank,tmprow,tmpcol,tmpmyrow,tmpmycol);


	if (mpi_rank==cprocs) // spare (checkpointing) node
	{
		// Descriptors (one_row)
		lld_global = n;
		Acp_cp_only = malloc(n*n*sizeof(double));
		descinit_( descAcp_cp_only, &n, &n, &one, &one, &zero, &zero, &context_cp_only, &lld_global, &info );
	}
	else
	{
		// Descriptors (global, for non-spare nodes)
		Acp_cp_only=NULL;
		for (i=0; i<9; i++)
		{
			descAcp_cp_only[i]=0;
		}
		descAcp_cp_only[1]=-1;
	}

	if (mpi_rank < cprocs)
	{
		//printf("\nI %d have to",myrank);

		// Computation of local matrix size
		nb = SCALAPACKNB;
		nr = numroc_( &n, &nb, &myrow, &zero, &nprow );
		nc = numroc_( &n, &nb, &mycol, &zero, &npcol );
		//ncrhs = numroc_( &m, &nb, &mycol, &zero, &npcol );
		lld = MAX( 1 , nr );
		A = malloc(nr*nc*sizeof(double));

		// Descriptors (local)
		descinit_( descA, &n, &n, &nb, &nb, &zero, &zero, &context, &lld, &info );
	}
	else
	{
		A=NULL;
		for (i=0; i<9; i++)
		{
			descA[i]=0;
		}
		descA[1]=-1;
	}


	printf("\ntest\n");
	pdgemr2d_ (&n, &n, A, &one, &one, descA, Acp_cp_only, &one, &one, descAcp_cp_only, &context_all);
	printf("\ntest done\n");

	/*
		if (mpi_rank==0)
		{
			// Descriptors (global)
			lld_global = n;
			descinit_( descA_global, &n, &n, &one, &one, &zero, &zero, &context_global, &lld_global, &info );
			//descinit_( descB_global, &n, &m, &one, &one, &zero, &zero, &context_global, &lld_global, &info );
		}
		else
		{
			// Descriptors (global, for non-root nodes)
			for (i=0; i<9; i++)
			{
				descA_global[i]=0;
				//descB_global[i]=0;
			}
			descA_global[1]=-1;
			//descB_global[1]=-1;
		}
*/


	if (mpi_rank < cprocs)
	{
		// cleanup
		free(A);
	}
	if (mpi_rank==cprocs)
	{
		free(Acp_cp_only);
	}

	MPI_Barrier(MPI_COMM_WORLD);
}
