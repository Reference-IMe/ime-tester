#include <mpi.h>
#include "../../helpers/macros.h"
#include "../../helpers/matrix.h"
#include "../../helpers/vector.h"
#include "../../helpers/Cblacs.h"
#include "../../helpers/scalapack.h"

//#define N 10
#define NB 64

void Scalapack_pDGESV_ckp_ft1_sim(int n, double* A_global, int m, double* B_global, int mpi_rank, int cprocs, int sprocs, int failing_rank, int failing_level)
{
	/*
	 * n = system rank (A_global n x n)
	 * m = num. of r.h.s (B_global n x m)
	 */

	// general
	int i, j;				//iterators
	int zero = 0, one = 1;	//numbers
	int nprocs = cprocs + sprocs;
	// MPI
	int ndims = 2, dims[2] = {0,0};
	// BLACS/SCALAPACK
	int nprow, npcol, info, ic = -1, context, context_global, context_all, myrow, mycol;
	int descA_global[9], descB_global[9], descA[9], descB[9];
	char order = 'R', scope = 'A';
	// MATRIX
	int nb = NB, nr, nc, ncrhs, lld, lld_global;
	double *A, *B;
	double *work, alpha;
	int *ipiv;

	// Initialize a default BLACS context and the processes grid
	MPI_Dims_create(cprocs, ndims, dims);
	nprow = dims[0];
	npcol = dims[1];
	Cblacs_get( ic, zero, &context );
	Cblacs_gridinit( &context, &order, nprow, npcol );
	Cblacs_get( ic, zero, &context_global );
	Cblacs_gridinit( &context_global, &order, one, one );
	Cblacs_get( ic, zero, &context_all );
	Cblacs_gridinit( &context_all, &order, one, nprocs );
	Cblacs_gridinfo( context, &nprow, &npcol, &myrow, &mycol );

	/*
	MPI_Barrier(MPI_COMM_WORLD);
	printf("\n mpi-blacs: %d-(%d,%d)",myrank,myrow,mycol);
	MPI_Barrier(MPI_COMM_WORLD);
	//Cblacs_barrier( context, "All");
	*/

if (mpi_rank < cprocs)
{
	//printf("\nI %d have to",myrank);

	// Computation of local matrix size
	nr = numroc_( &n, &nb, &myrow, &zero, &nprow );
	nc = numroc_( &n, &nb, &mycol, &zero, &npcol );
	ncrhs = numroc_( &m, &nb, &mycol, &zero, &npcol );
	lld = MAX( 1 , nr );
	A = malloc(nr*nc*sizeof(double));
	B = malloc(nr*ncrhs*sizeof(double));
	work = malloc(nb*sizeof(double));
	ipiv = malloc((lld+nb)*sizeof(int));

	// Descriptors (local)
	descinit_( descA, &n, &n, &nb, &nb, &zero, &zero, &context, &lld, &info );
	descinit_( descB, &n, &m, &nb, &nb, &zero, &zero, &context, &lld, &info );

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

	// spread matrices
	pdgemr2d_(&n, &n, A_global, &one, &one, descA_global, A, &one, &one, descA, &context);
	pdgemr2d_(&n, &m, B_global, &one, &one, descB_global, B, &one, &one, descB, &context);

	// Linear system equations solver
	// pdgesv_(  &n, &m, A, &one, &one, descA, ipiv, B, &one, &one, descB, &info );
	// split in LU factorization + solve (pdgetrf + pdgetrs) to introduce checkpointing
	pdgetrf_(      &n, &n, A, &one, &one, descA, ipiv, &info );
	pdgetrs_("N",  &n, &m, A, &one, &one, descA, ipiv, B, &one, &one, descB, &info  );

	pdgemr2d_(&n, &m, B, &one, &one, descB, B_global, &one, &one, descB_global, &context);


	free(A);
	free(B);
	free(ipiv);
	free(work);
}
else
{
	//printf("\nI %d NOT have to",myrank);
}

	//Close BLACS environment
		//Cblacs_barrier( context, "All");	// not working! why?
		//MPI_Barrier(MPI_COMM_WORLD);		// working but not needed
		//Cblacs_gridexit( context );		// not needed if calling blacs_exit
		//Cblacs_gridexit( context_global );// not needed if calling blacs_exit
	//Cblacs_exit( one );					// argument not 0: it is assumed the user will continue using the machine after the BLACS are done
											// error, if main function called more tha once, why?
}
