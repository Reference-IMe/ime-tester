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

void ScaLAPACK_pDGETRF_cp_ft1_sim(int n, double* A_source, int mpi_rank, int cprocs, int sprocs, int failing_rank, int failing_level, int checkpoint_freq)
{
	/*
	 * n = system rank (A_global n x n)
	 */

	// general
	int i, j;				//iterators
	int zero = 0, one = 1;	//numbers
	int nprocs = cprocs + sprocs;
	//int cpfreq = 2; 		// checkpointing frequency
	// MPI
	int ndims = 2, dims[2] = {0,0};
	// BLACS/SCALAPACK
	int ic = -1, info;
	int context_distributed, context_source, context_cp, context_all;
	int nprow, npcol, myrow, mycol;
	int tmprow,tmpcol,tmpmyrow, tmpmycol;
	int descA_source[9], descA[9], descA_cp[9], descIPIV[9], descIPIV_cp[9];
	//int descB_global[9], descB[9];
	char order = 'R', scope = 'A';
	// MATRIX
	int nb, nr, nc, ncrhs, lldA, lldA_source;
	int nIPIV, nrIPIV, ncIPIV, lldIPIV;
	double *A;
	double *A_cp;
	//double *B;
	//double *work, alpha;
	int *IPIV;
	int *IPIV_cp;

	// Initialize a default BLACS context and the processes grid
	MPI_Dims_create(cprocs, ndims, dims);
	nprow = dims[0];
	npcol = dims[1];

	// context for all processes
	Cblacs_get( ic, zero, &context_all );
	Cblacs_gridinit( &context_all, &order, one, nprocs );

	// context for distributed matrix A
	Cblacs_get( ic, zero, &context_distributed );
	Cblacs_gridinit( &context_distributed, &order, nprow, npcol );
	Cblacs_gridinfo( context_distributed, &nprow, &npcol, &myrow, &mycol );
	//MPI_Barrier(MPI_COMM_WORLD);
	//printf("context_distributed: %d in %dx%d id %d,%d\n",mpi_rank,nprow,npcol,myrow,mycol);

	// context for source global matrix A (A_source)
	Cblacs_get( ic, zero, &context_source );
	Cblacs_gridinit( &context_source, &order, one, one );
	Cblacs_gridinfo( context_source, &tmprow, &tmpcol, &tmpmyrow, &tmpmycol );
	//MPI_Barrier(MPI_COMM_WORLD);
	//printf("context_source: %d in %dx%d id %d,%d\n",mpi_rank,tmprow,tmpcol,tmpmyrow,tmpmycol);

	// context for checkpointing global matrix (A_cp)
	Cblacs_get( ic, zero, &context_cp );
	int map_cp[1][1];
	map_cp[0][0]=cprocs;
	Cblacs_gridmap( &context_cp, map_cp, one, one, one);
	Cblacs_gridinfo( context_cp, &tmprow, &tmpcol, &tmpmyrow, &tmpmycol );
	//MPI_Barrier(MPI_COMM_WORLD);
	//printf("context_cp: %d in %dx%d id %d,%d\n",mpi_rank,tmprow,tmpcol,tmpmyrow,tmpmycol);

	// Computation of local matrix size
	nb = SCALAPACKNB;
	nr = numroc_( &n, &nb, &myrow, &zero, &nprow );
	nc = numroc_( &n, &nb, &mycol, &zero, &npcol );
	//ncrhs = numroc_( &m, &nb, &mycol, &zero, &npcol );
	lldA = MAX( 1 , nr );
	MPI_Allreduce( MPI_IN_PLACE, &lldA, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD );

	nIPIV = (lldA+nb);

	lldA_source = n;

	if (mpi_rank==0) // root node
	{
		// Descriptors (global)
		descinit_( descA_source, &n, &n, &one, &one, &zero, &zero, &context_source, &n, &info );
		//descinit_( descB_global, &n, &m, &one, &one, &zero, &zero, &context_global, &lldA_source, &info );
	}
	else
	{
		// Descriptors (global, for non-root nodes)
		A_source=NULL;
		for (i=0; i<9; i++)
		{
			descA_source[i]=0;
			//descB_global[i]=0;
		}
		descA_source[1]=-1;
		//descB_global[1]=-1;
	}

	if (mpi_rank==cprocs) // spare (checkpointing) node
	{
		// Descriptors
		A_cp = malloc(n*n*sizeof(double));
		descinit_( descA_cp, &n, &n, &one, &one, &zero, &zero, &context_cp, &n, &info );
		OneMatrix1D(A_cp, n, n);

		//IPIV_cp=malloc(nIPIV*cprocs*sizeof(int)); // with cprocs is not good because MPI_GATHER wants a buffer for everyone
		IPIV_cp=malloc(nIPIV*nprocs*sizeof(int));
	}
	else
	{
		// Descriptors (global, for non-spare nodes)
		A_cp=NULL;
		IPIV_cp=NULL;
		for (i=0; i<9; i++)
		{
			descA_cp[i]=0;
			descIPIV_cp[i]=0;
		}
		descA_cp[1]=-1;
		descA_cp[4]=nb; // *** important!
		descA_cp[5]=nb; // ***
	}

	if (mpi_rank < cprocs) // non-spare nodes
	{
		A = malloc(nr*nc*sizeof(double));
		IPIV = malloc(nIPIV*sizeof(int));
		//work = malloc(nb*sizeof(double));
		//B = malloc(nr*ncrhs*sizeof(double));

		// Descriptors (local)
		// all processes have to know descA (not only non-spare nodes)
		descinit_( descA, &n, &n, &nb, &nb, &zero, &zero, &context_distributed, &lldA, &info );
	}
	else
	{
		A=NULL;
		//IPIV=NULL;
		IPIV = malloc(nIPIV*sizeof(int)); // also allocated on the spare proc because MPI_GATHER wants a buffer for everyone
		for (i=0; i<9; i++)
		{
			descA[i]=0;
			descIPIV[i]=0;
		}
		// all processes have to know something about descA (not only non-spare nodes)
		// can't use descinint, due to illegal values of spare process not belonging to the right context
		descA[1]=-1;
		descA[4]=nb;
		descA[5]=nb;

		descIPIV[1]=-1;
	}


		// spread matrices
		pdgemr2d_(&n, &n, A_source, &one, &one, descA_source, A, &one, &one, descA, &context_all);
		//if (mpi_rank < cprocs) pdgemr2d_(&n, &n, A_source, &one, &one, descA_source, A, &one, &one, descA, &context_distributed);

		if (failing_level>=0)
		{
			failing_level=n-failing_level;
		}

		// LU factorization
		pdgetrf_cp_  (&n, &n, A, &one, &one, descA, A_cp, &one, &one, descA_cp, IPIV, IPIV_cp, &nIPIV, &checkpoint_freq, &failing_level, &context_all, &info );

		// check factorization
		/*
		pdgetrs_("N",  &n, &m, A, &one, &one, descA, ipiv, B, &one, &one, descB, &info  );
		pdgemr2d_(&n, &m, B, &one, &one, descB, B_global, &one, &one, descB_global, &context);
		*/

		pdgemr2d_ (&n, &n, A, &one, &one, descA, A_source, &one, &one, descA_source, &context_all);
		//if (mpi_rank < cprocs) pdgemr2d_ (&n, &n, A, &one, &one, descA, A_source, &one, &one, descA_source, &context_distributed);

	if (mpi_rank < cprocs)
	{
		// cleanup
		free(A);
		//free(B);
		free(IPIV);
		//free(work);
	}
	else
	{
		free(A_cp);
		free(IPIV_cp);
		free(IPIV);
	}

	MPI_Barrier(MPI_COMM_WORLD);
}
