/*
 * ScaLAPACK_pDGETRF_ft1.h
 *
 *  Created on: Dec 28, 2019
 *      Author: marcello
 */

#include <mpi.h>
#include <time.h>
#include <stdio.h>

#include "../../helpers/macros.h"
#include "../../helpers/matrix.h"
#include "../../helpers/vector.h"
#include "../../helpers/Cblacs.h"
#include "../../helpers/scalapack.h"
#include "../tester_structures.h"

test_output ScaLAPACK_pDGESV_ft1(int n, double* A_source, int m, double* B_source, int nb, int mpi_rank, int cprocs, int sprocs, int failing_level, int checkpoint_freq)
{
	test_output result = EMPTY_OUTPUT;

	result.total_start_time = time(NULL);

	/*
	 * n = system rank (A_global n x n)
	 * m = num. of r.h.s (B_global n x m)
	 */

	// general
	int i;					//iterators
	int zero = 0, one = 1;	//numbers
	double dzero = 0.0, done = 1.0;
	int nprocs = cprocs + sprocs;

	// MPI
	int ndims = 2, dims[2] = {0,0};
	// BLACS/SCALAPACK
	int ic = -1, info;
	int context_distributed, context_source, context_cp, context_all;
	int nprow, npcol, myrow, mycol;
	int tmprow,tmpcol,tmpmyrow, tmpmycol;
	int descA_source[9], descA[9], descAt[9], descA_cp[9];
	int descB_source[9], descB[9], descBt[9];
	char order = 'R';
	// MATRIX
	int nr, nc, lld;
	int ncrhs, nrrhs, lld_global, lldt;
	int ncrhst, nrrhst;
	int nIPIV;
	double *A, *At, *B, *Bt;
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
	//int map_cp[1][1];
	int map_cp[1];
	//map_cp[0][0]=cprocs;
	map_cp[0]=cprocs;
	Cblacs_gridmap( &context_cp, map_cp, one, one, one);
	Cblacs_gridinfo( context_cp, &tmprow, &tmpcol, &tmpmyrow, &tmpmycol );
	//MPI_Barrier(MPI_COMM_WORLD);
	//printf("context_cp: %d in %dx%d id %d,%d\n",mpi_rank,tmprow,tmpcol,tmpmyrow,tmpmycol);

	// Computation of local matrix size
	//nb = SCALAPACKNB;
	nr = numroc_( &n, &nb, &myrow, &zero, &nprow );
	nc = numroc_( &n, &nb, &mycol, &zero, &npcol );

	//ncrhs = numroc_( &m, &nb, &mycol, &zero, &npcol );
	lld = MAX( 1 , nr );
	MPI_Allreduce( MPI_IN_PLACE, &lld, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD );

	nIPIV = (lld+nb);

	ncrhs = numroc_( &m, &nb, &mycol, &zero, &npcol );
	nrrhs = numroc_( &n, &nb, &myrow, &zero, &nprow );

	ncrhst = numroc_( &n, &nb, &mycol, &zero, &npcol );
	nrrhst = numroc_( &m, &nb, &myrow, &zero, &nprow );
	lldt = MAX( 1 , nrrhst );
	MPI_Allreduce( MPI_IN_PLACE, &lldt, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD );



	if (mpi_rank==0) // root node
	{
		// Descriptors (global)
		descinit_( descA_source, &n, &n, &one, &one, &zero, &zero, &context_source, &n, &info );
		descinit_( descB_source, &n, &m, &one, &one, &zero, &zero, &context_source, &n, &info );
	}
	else
	{
		// Descriptors (global, for non-root nodes)
		A_source=NULL;
		B_source=NULL;
		for (i=0; i<9; i++)
		{
			descA_source[i]=0;
			descB_source[i]=0;
		}
		descA_source[1]=-1;
		descB_source[1]=-1;
	}

	if (mpi_rank==cprocs) // spare (checkpointing) node
	{
		// Descriptors
		A_cp = malloc(n*n*sizeof(double));
		descinit_( descA_cp, &n, &n, &one, &one, &zero, &zero, &context_cp, &n, &info );
		//OneMatrix1D(A_cp, n, n);

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
		}
		descA_cp[1]=-1;
		descA_cp[4]=nb; // *** important!
		descA_cp[5]=nb; // ***
	}

	if (mpi_rank < cprocs) // non-spare nodes
	{
		A = malloc(nr*nc*sizeof(double));
		At = malloc(nr*nc*sizeof(double));
		B = malloc(nrrhs*ncrhs*sizeof(double));
		Bt = malloc(nrrhst*ncrhst*sizeof(double));
		IPIV = malloc(nIPIV*sizeof(int));

		// Descriptors (local)
		// all processes have to know descA (not only non-spare nodes)
		descinit_( descA,  &n, &n, &nb, &nb, &zero, &zero, &context_distributed, &lld, &info );
		descinit_( descAt, &n, &n, &nb, &nb, &zero, &zero, &context_distributed, &lld, &info );
		descinit_( descB,  &n, &m, &nb, &nb, &zero, &zero, &context_distributed, &lld, &info );
		descinit_( descBt, &m, &n, &nb, &nb, &zero, &zero, &context_distributed, &lldt, &info );
	}
	else
	{
		A=NULL;
		At=NULL;
		B=NULL;
		Bt=NULL;
		//IPIV=NULL;
		IPIV = malloc(nIPIV*sizeof(int)); // also allocated on the spare proc because MPI_GATHER wants a buffer for everyone

		for (i=0; i<9; i++)
		{
			descA[i]=0;
			descAt[i]=0;
			descB[i]=0;
			descBt[i]=0;
		}
		// all processes have to know something about descA (not only non-spare nodes)
		// can't use descinint, due to illegal values of spare process not belonging to the right context
		descA[1]=-1;
		descA[4]=nb;
		descA[5]=nb;

		descAt[1]=-1;
		descAt[4]=nb;
		descAt[5]=nb;

		descB[1]=-1;
		//descB[4]=nb;
		//descB[5]=nb;
		/*
		descBt[1]=-1;
		descBt[4]=nb;
		descBt[5]=nb;
		*/
	}

	// spread matrices
	pdgemr2d_(&n, &n, A_source, &one, &one, descA_source, A, &one, &one, descA, &context_all);
	pdgemr2d_(&n, &m, B_source, &one, &one, descB_source, B, &one, &one, descB, &context_all);
	MPI_Barrier(MPI_COMM_WORLD);
	printf("%d hic\n",mpi_rank);
	fflush(stdout);

	// transpose system matrix
	if (mpi_rank < cprocs) pdtran_(&n, &n, &done, A, &one, &one, descA, &dzero, At, &one, &one, descAt);
	MPI_Barrier(MPI_COMM_WORLD);
	printf("%d qua\n",mpi_rank);
	fflush(stdout);


	if (failing_level>=0)
	{
		failing_level=n-failing_level;
	}

	// LU factorization (with checkpointing)
	result.core_start_time = time(NULL);
	pdgetrf_cp_  (&n, &n, At, &one, &one, descAt, A_cp, &one, &one, descA_cp, IPIV, IPIV_cp, &nIPIV, &checkpoint_freq, &failing_level, &context_all, &info );

	// solve (without checkpointing)
	sleep(3);
	MPI_Barrier(MPI_COMM_WORLD);
	printf("%d qui\n",mpi_rank);
	fflush(stdout);

	if (mpi_rank < cprocs) // otherwise ERR: {   -1,   -1}:  On entry to PDGETRS parameter number  702 had an illegal value
	{
		descinit_( descA,  &n, &n, &nb, &nb, &zero, &zero, &context_distributed, &lld, &info );
		descinit_( descAt, &n, &n, &nb, &nb, &zero, &zero, &context_distributed, &lld, &info );
		descinit_( descB,  &n, &m, &nb, &nb, &zero, &zero, &context_distributed, &lld, &info );
		descinit_( descBt, &m, &n, &nb, &nb, &zero, &zero, &context_distributed, &lldt, &info );
		pdgetrs_("N", &n, &m, At, &one, &one, descAt, IPIV, B, &one, &one, descB, &info );
	}
	/*
	else
	{
		descinit_( descAt, &n, &n, &nb, &nb, &zero, &zero, &context_distributed, &lld, &info );
		descinit_( descB, &n, &m, &nb, &nb, &zero, &zero, &context_distributed, &lld, &info );
		pdgetrs_("N", &n, &m, At, &one, &one, descAt, IPIV, B, &one, &one, descB, &info );
	}
	 */
	sleep(3);
	MPI_Barrier(MPI_COMM_WORLD);
	printf("%d quo\n",mpi_rank);
	fflush(stdout);

	result.core_end_time = time(NULL);
	result.exit_code = info;

	// transpose system matrix
	MPI_Barrier(MPI_COMM_WORLD);
	//if (mpi_rank < cprocs) pdtran_(&m, &n, &done, B, &one, &one, descB, &dzero, Bt, &one, &one, descBt);

	//pdgemr2d_(&n, &m, B, &one, &one, descB, B_source, &one, &one, descB_source, &context_all);
	//pdgemr2d_(&n, &n, A, &one, &one, descA, A_source, &one, &one, descA_source, &context_all);
	//if (mpi_rank < cprocs) pdgemr2d_ (&n, &n, A, &one, &one, descA, A_source, &one, &one, descA_source, &context_distributed);

	//pdgemr2d_(&m, &n, Bt, &one, &one, descBt, B_source, &one, &one, descB_source, &context_all);



	MPI_Barrier(MPI_COMM_WORLD);
	if (mpi_rank < cprocs)
	{
		// cleanup
		free(A);
		free(At);
		free(B);
		free(Bt);
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

	result.total_end_time = time(NULL);

	return result;
}
