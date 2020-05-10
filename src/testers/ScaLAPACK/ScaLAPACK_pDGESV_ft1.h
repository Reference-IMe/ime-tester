/*
 * ScaLAPACK_pDGESV_ft1.h
 *
 *  Created on: Dec 28, 2019
 *      Author: marcello
 */

#include <mpi.h>
#include <time.h>
#include "../../helpers/macros.h"
#include "../../helpers/matrix.h"
#include "../../helpers/vector.h"
#include "../../helpers/Cblacs.h"
#include "../../helpers/scalapack.h"
#include "../tester_structures.h"


test_output ScaLAPACK_pDGESV_ft1(int n, double* A_global, int m, double* B_global, int nb, int mpi_rank, int cprocs, int sprocs, int failing_level, int checkpoint_freq)
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
	int nprow, npcol, info, ic = -1, context, context_global, context_all, myrow, mycol;
	int descA_global[9], descB_global[9], descA[9], descAt[9], descB[9], descBt[9];
	int context_cp;
	int descA_cp[9];
	char order = 'R'; //, order_all ='A';
	// MATRIX
	int nr, nc, ncrhs, nrrhs, lld, lld_global, lldt;
	int ncrhst, nrrhst;
	double *A, *At, *B, *Bt;
	int *ipiv;
	double *A_cp;
	int *ipiv_cp;
	int nipiv;
	
	// Initialize a default BLACS context and the processes grid
	MPI_Dims_create(cprocs, ndims, dims);
	nprow = dims[0];
	npcol = dims[1];

	// context for checkpointing global matrix (A_cp)
	Cblacs_get( ic, zero, &context_cp );
	int map_cp[1];
	map_cp[0]=cprocs;
	Cblacs_gridmap( &context_cp, map_cp, one, one, one);

	Cblacs_get( ic, zero, &context );
	Cblacs_gridinit( &context, &order, nprow, npcol );
	Cblacs_get( ic, zero, &context_global );
	Cblacs_gridinit( &context_global, &order, one, one );
	Cblacs_get( ic, zero, &context_all );
	Cblacs_gridinit( &context_all, &order, one, nprocs );
	Cblacs_gridinfo( context, &nprow, &npcol, &myrow, &mycol );

	// Computation of local matrix size
	nc = numroc_( &n, &nb, &mycol, &zero, &npcol );
	nr = numroc_( &n, &nb, &myrow, &zero, &nprow );
	lld = MAX( 1 , nr );
	MPI_Allreduce( MPI_IN_PLACE, &lld, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD );

	nipiv = (lld+nb);
	ipiv = malloc(nipiv*sizeof(int));// also allocated on the spare proc because MPI_GATHER wants a buffer for everyone

	ncrhs = numroc_( &m, &nb, &mycol, &zero, &npcol );
	nrrhs = numroc_( &n, &nb, &myrow, &zero, &nprow );

	ncrhst = numroc_( &n, &nb, &mycol, &zero, &npcol );
	nrrhst = numroc_( &m, &nb, &myrow, &zero, &nprow );
	lldt = MAX( 1 , nrrhst );
	MPI_Allreduce( MPI_IN_PLACE, &lldt, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD );

	// global matrices
	if (mpi_rank==0)
	{
		// Descriptors (global)
		lld_global = n;
		descinit_( descA_global, &n, &n, &one, &one, &zero, &zero, &context_global, &lld_global, &info );
		descinit_( descB_global, &n, &m, &one, &one, &zero, &zero, &context_global, &lld_global, &info );
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

	// checkpointing matrices
	if (mpi_rank==cprocs) // spare (checkpointing) node
	{
		// Allocation
		A_cp = malloc(n*n*sizeof(double));
		ipiv_cp=malloc(nipiv*nprocs*sizeof(int)); //ipiv_cp=malloc(nIPIV*cprocs*sizeof(int)); // with cprocs is not good because MPI_GATHER wants a buffer for everyone

		// Descriptors
		descinit_( descA_cp, &n, &n, &one, &one, &zero, &zero, &context_cp, &n, &info );
	}
	else				// non-spare nodes
	{
		// Allocation not needed
		A_cp=NULL;
		ipiv_cp=NULL;

		// Descriptors
		for (i=0; i<9; i++)
		{
			descA_cp[i]=0;
		}
		descA_cp[1]=-1;
		descA_cp[4]=nb; // *** important!
		descA_cp[5]=nb; // ***
	}

	// locally distributed matrices
	if (mpi_rank < cprocs)				// non-spare nodes
	{
		// Allocation
		A = malloc(nr*nc*sizeof(double));
		At = malloc(nr*nc*sizeof(double));
		B = malloc(nrrhs*ncrhs*sizeof(double));
		Bt = malloc(nrrhst*ncrhst*sizeof(double));

		// Descriptors (local)
		descinit_( descA, &n, &n, &nb, &nb, &zero, &zero, &context, &lld, &info );
		descinit_( descAt, &n, &n, &nb, &nb, &zero, &zero, &context, &lld, &info );
		descinit_( descB, &n, &m, &nb, &nb, &zero, &zero, &context, &lld, &info );
		descinit_( descBt, &m, &n, &nb, &nb, &zero, &zero, &context, &lldt, &info );
	}
	else								// spare node
	{
		// Allocation not needed
		A=NULL;
		B=NULL;
		At=NULL;
		Bt=NULL;

		// Descriptors
		for (i=0; i<9; i++)
		{
			descA[i]=0;
			descB[i]=0;
			descAt[i]=0;
			descBt[i]=0;
		}
							// all processes have to know something about locally distributed matrices (not only non-spare nodes)
							// can't use descinint, due to illegal values of spare process not belonging to the right context
		descA[1]=-1;
		descA[4]=nb;
		descA[5]=nb;

		descB[1]=-1;
		descB[4]=nb;
		descB[5]=nb;

		descAt[1]=-1;
		descAt[4]=nb;
		descAt[5]=nb;

		descBt[1]=-1;
		descBt[4]=nb;
		descBt[5]=nb;
	}

	// spread matrices
	pdgemr2d_(&n, &n, A_global, &one, &one, descA_global, A, &one, &one, descA, &context_all);
	pdgemr2d_(&n, &m, B_global, &one, &one, descB_global, B, &one, &one, descB, &context_all);

	// transpose system matrix
	if (mpi_rank<cprocs) pdtran_(&n, &n, &done, A, &one, &one, descA, &dzero, At, &one, &one, descAt);

	result.core_start_time = time(NULL);

	// Linear system equations solver

	if (failing_level>=0)
	{
		failing_level=n-failing_level;
	}

	// split in LU factorization + solve (pdgetrf + pdgetrs) to introduce checkpointing
	// checkpointed factorization called by everyone
	pdgetrf_cp_  (&n, &n, At, &one, &one, descAt, A_cp, &one, &one, descA_cp, ipiv, ipiv_cp, &nipiv, &checkpoint_freq, &failing_level, &context_all, &info );
	// solve called by non-spare nodes only
	if (mpi_rank < cprocs)
	{
		pdgetrs_("N",  &n, &m, At, &one, &one, descAt, ipiv, B, &one, &one, descB, &info  );
	    result.core_end_time = time(NULL);
		result.exit_code = info;

		// re-transpose result
		pdtran_(&m, &n, &done, B, &one, &one, descB, &dzero, Bt, &one, &one, descBt);

		// collect result
		if (mpi_rank==0)
		{
			// Adapt descriptor for B to accept transposed matrix
			descinit_( descB_global, &m, &n, &one, &one, &zero, &zero, &context_global, &m, &info );
		}
		// get matrix back
		pdgemr2d_(&m, &n, Bt, &one, &one, descBt, B_global, &one, &one, descB_global, &context);

		// cleanup
		free(A);
		free(At);
		free(B);
		free(Bt);
		free(ipiv);
	}
	else
	{
		free(A_cp);
		free(ipiv_cp);
		free(ipiv);
	}

		//Close BLACS environment
		//Cblacs_barrier( context, "All");	// not working! why?
		//MPI_Barrier(MPI_COMM_WORLD);		// working but not needed
		//Cblacs_gridexit( context );		// not needed if calling blacs_exit
		//Cblacs_gridexit( context_global );// not needed if calling blacs_exit
		//Cblacs_exit( one );				// argument not 0: it is assumed the user will continue using the machine after the BLACS are done
											// error, if main function called more tha once, why?
	//Cblacs_barrier( context_all, &order_all);
	MPI_Barrier(MPI_COMM_WORLD);

	result.total_end_time = time(NULL);

	return result;
}
