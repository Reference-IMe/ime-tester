#include <mpi.h>
#include "../../helpers/macros.h"
#include "../../helpers/matrix.h"
#include "../../helpers/vector.h"
#include "../../helpers/Cblacs.h"
#include "../../helpers/scalapack.h"
#include "../FTLA/util_matrix.h"
#include "../FTLA/util_inject.h"
#include "../FTLA/ftla_ftwork.h"
#include "../FTLA/ftla_cof.h"
#include "../FTLA/ftla_driver.h"

int *errors;

void FTLA_pDGEQRF_calc(int n, double* A_global, int mpi_rank, int cprocs, int sprocs)
{
	/*
	 * n = system rank (A_global n x n)
	 */
    ftla_work_t ftwork;
    int err = 0;
    int Fstrat='e', F, Fmin=0, Fmax=1, Finc=1;

	// general
	int i, j;						//iterators
	int zero = 0, one = 1;			//numbers
	// MPI
	int ndims = 2, dims[2] = {0,0};
	// BLACS/SCALAPACK
	int nprow, npcol, info, ic = -1, context, context_global, myrow, mycol;
	int descA_global[9], descA[9];
	char order = 'R', scope = 'A';
	// MATRIX
	int nb, nr, nc, lld, lld_global;
	double* A;
	double* work;
	double* tau;

	// Initialize a default BLACS context and the processes grid
	MPI_Dims_create(cprocs, ndims, dims);
	nprow = dims[0];
	npcol = dims[1];
	printf("req procs grid %dx%d\n",nprow,npcol);
	Cblacs_get( ic, zero, &context );
	Cblacs_gridinit( &context, &order, nprow, npcol );
	Cblacs_get( ic, zero, &context_global );
	Cblacs_gridinit( &context_global, &order, one, one );
	Cblacs_gridinfo( context, &nprow, &npcol, &myrow, &mycol );
	printf("conf procs grid %dx%d\n",nprow,npcol);


	int N, M, Ne, NB;
	N = M = n;

if (mpi_rank < cprocs)
{
	// Computation of local matrix size
	nb = NB = SCALAPACKNB;
	nr = numroc_( &n, &nb, &myrow, &zero, &nprow );
	nc = numroc_( &n, &nb, &mycol, &zero, &npcol );
	MPI_Allreduce( MPI_IN_PLACE, &nc, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD );
	lld = MAX( 1 , nr );
	A = malloc(nr*nc*sizeof(double));
	tau = malloc( nc*sizeof(double) );

/*
#ifndef NO_EXTRAFLOPS
            Ne = N + nc*2;// + ((Nc/NB)%Q==0)*NB;
#else
            Ne = N;
#endif
*/
	Ne = N + nc*2;

	// Descriptors (local)
	descinit_( descA, &n, &n, &nb, &nb, &zero, &zero, &context, &lld, &info );

	if (mpi_rank==0)
	{
		// Descriptors (global)
		lld_global = n;
		descinit_( descA_global, &n, &n, &one, &one, &zero, &zero, &context_global, &lld_global, &info );
	}
	else
	{
		// Descriptors (global, for non-root nodes)
		for (i=0; i<9; i++)
		{
			descA_global[i]=0;
		}
		descA_global[1]=-1;
	}

	// spread matrices
	pdgemr2d_(&n, &n, A_global, &one, &one, descA_global, A, &one, &one, descA, &context);

	// init work space
	int lwork=-1;
	double lazywork;
	pdgeqrf_(  &n, &n, A, &one, &one, descA, NULL, &lazywork, &lwork, &info );
	lwork = (int)lazywork;
	work = malloc( lwork*sizeof(double) );

#ifdef INJECT
          for( F = Fmin; F<=Fmax; F+=Finc ) {
            errors = create_error_list( M, NB, F, Fstrat );
#endif

	// QR factorization
    Cftla_work_construct( 0, descA, 0, Ne-N, &ftwork );

    MPI_Barrier( MPI_COMM_WORLD );
    //t1 = MPI_Wtime();
    do { // call ftpdgeqrf until we complete w/o a failure
#ifdef USE_CoF
      if(err) Cftla_cof_dgeqrr( A, descA, tau, work, lwork, &ftwork );
#endif
      err = ftla_pdgeqrf( &M, &Ne, A, &one, &one, descA, tau, work, &lwork, &info, (int*)&ftwork );
      //printf("****** err=%d i1=%d\n",err,i1);

#ifdef USE_CoF
      if(err) Cftla_cof_dgeqrr( A, descA, tau, work, lwork, &ftwork );
#endif
    } while(err);
    //checkerror( info, 0 );

	pdgemr2d_ (&n, &n, A, &one, &one, descA, A_global, &one, &one, descA_global, &context);

    Cftla_cof_cleanup( &ftwork );
    Cftla_work_destruct( &ftwork );
    free( errors );

#ifdef INJECT
          }
#endif
	free(A);
	free(work);
	free(tau);
}
	//Close BLACS environment
		//Cblacs_gridexit( context );
		//Cblacs_gridexit( context_global );
		//Cblacs_exit( zero );
}
