/*
 * matrix_advanced.h
 *
 *  Created on: 12 Mar 2020
 *      Author: Marcello
 */
/*
#include "matrix.h"
#include "vector.h"
#include "Cblacs.h"
#include "lapack.h"
#include "scalapack.h"
#include "types.h"
*/
void OrthogonalizeMatrix1D_double(double* mat, int rows, int cols)
{
	int info;
    int lwork;
    	lwork = MIN(rows,cols);

    double* work;
			work = malloc(lwork*sizeof(double));
    double* tau;
    		tau = malloc(lwork*sizeof(double));

    dgeqrfp_(&rows, &cols, mat, &rows, tau, work, &lwork, &info);
    dorgqr_(&rows, &cols, &cols, mat, &rows, tau, work, &lwork, &info);

    NULLFREE(work);
    NULLFREE(tau);
}

void RandomSquareMatrixCND1D_double(double* mat, int n, int seed, double cnd)
{
	double* mat1;
			mat1=AllocateMatrix1D_double(n, n);
	double* mat2;
			mat2=AllocateMatrix1D_double(n, n);

	double* s;
			s=AllocateVector_double(n);

	double smax = log10(cnd)/2;
	double smin = -smax;

	double gap;
			gap=(smax-smin)/(n-1);

	int r,c;
	for (r=0;r<n;r++)
	{
		s[r]=pow(10,smin + r*gap);
	}

	RandomMatrix1D_double(mat1, n, n, seed);
	OrthogonalizeMatrix1D_double(mat1, n, n);
	RandomMatrix1D_double(mat2, n, n, seed+1);
	OrthogonalizeMatrix1D_double(mat2, n, n);

	// mat <- mat1.S.mat2 = mat1.(S.mat2)
	//
	// mat2 <- S.mat2, where S=diag(s)
	for(r=0;r<n;r++)
	{
		for(c=0;c<n;c++)
		{
			mat2[r*n+c]=s[r]*mat2[r*n+c];
		}
	}
	// mat <- mat1.mat2
	char transA = 'N', transB = 'N';
	double one = 1.0, zero = 0.0;
	dgemm_(&transA, &transB, &n, &n, &n, &one, mat1, &n, mat2, &n, &zero, mat, &n);

	DeallocateVector_double(s);
	DeallocateMatrix1D_double(mat1);
	DeallocateMatrix1D_double(mat2);
}

double ConditionNumber1D_double(double* mat, int rows, int cols)
{
	// condition number calculated by definition:
	// the ratio 'cnd' of the largest to smallest singular value in the singular value decomposition of a matrix 'mat'.
	//
	// example for SVD in:
	// https://software.intel.com/sites/products/documentation/doclib/mkl_sa/11/mkl_lapack_examples/dgesvd_ex.c.htm

    int info, lwork, i;
    double cnd;
    double wkopt;
    double* work;
    double* u=NULL;
//			u=AllocateMatrix1D_double(rows, rows);
	double* vt=NULL;
//			vt=AllocateMatrix1D_double(cols, cols);
	double* s;
			s=AllocateVector_double(cols);
	double* mat_copy;									// SVD algorithm destroys input matrix mat
			mat_copy=AllocateMatrix1D_double(rows, cols);		// http://www.netlib.org/lapack/explore-html/d1/d7e/group__double_g_esing_ga84fdf22a62b12ff364621e4713ce02f2.html

	char nojob='N'; // no need of calculating u and vt

	for (i=0;i<rows*cols;i++)
	{
		mat_copy[i]=mat[i];
	}

	lwork = -1;
	dgesvd_( &nojob, &nojob, &rows, &cols, mat_copy, &rows, s, u, &rows, vt, &cols, &wkopt, &lwork, &info );
	lwork = (int)wkopt;
	work = AllocateVector_double(lwork);
	dgesvd_( &nojob, &nojob, &rows, &cols, mat_copy, &rows, s, u, &rows, vt, &cols, work, &lwork, &info );
    if( info > 0 )
    {
		printf( "The algorithm computing SVD failed to converge.\n" );
		exit( 1 );
    }
    cnd = s[0]/s[cols-1];

    DeallocateMatrix1D_double(mat_copy);
	DeallocateVector_double(s);
	DeallocateVector_double(work);
	DeallocateMatrix1D_double(u);
	DeallocateMatrix1D_double(vt);

	return cnd;
}

double GenSystemMatrices1D_double(int n, double* A, double* x, double* b, int seed, double cnd, char calc_cnd, char cnd_readback)
{
	int i1 = 1;
	double d0 = 0.0;
	double d1 = 1.0;
	char transA = 'T', transx = 'N';
	double read_cnd = -1;

	if (calc_cnd)
	{
		RandomSquareMatrixCND1D_double(A, n, seed, cnd);
	}
	else
	{
		RandomMatrix1D_double(A, n, n, seed);
	}
	FillVector_double(x, n, 1);
	dgemm_(&transA, &transx, &n, &i1, &n, &d1, A, &n, x, &n, &d0, b, &n);
	if (cnd_readback)
	{
		read_cnd = ConditionNumber1D_double(A, n, n);
	}
	return read_cnd;
}

double pGenSystemMatrices1D_double_pdgemr2d(int n, double* A, double* x, double* b, int seed, double cnd, char cnd_readback, int nb, int mpi_rank, int cprocs, int nprow, int npcol, int myrow, int mycol, int context, int context_global)
{
	/*
	 * This version suffers from the pdgemr2d bug: out of memory
	 *
	 * first workaround not completely successful
	 * https://icl.cs.utk.edu/lapack-forum/viewtopic.php?t=465
	 *
	 * final workaround suggested in
	 * https://stackoverflow.com/questions/30167724/how-to-use-pdgemr2d-to-copy-distributed-matrix-in-total-to-all-processes
	 * https://andyspiros.wordpress.com/2011/07/08/an-example-of-blacs-with-c/
	 * see "pGenSystemMatrices1D_double" below
	 *
	 */

	// general
	int i;
	int i0 = 0;
	int i1 = 1;
	double d0 = 0.0;
	double d1 = 1.0;

	double* A1;
	double* A2;
	double* S;
	double* X;
	double* B;

	int descA1[9];
	int descA2[9];
	int descS[9];
	int descA_global[9];
	int descX[9];
	int descX_global[9];
	int descB[9];
	int descB_global[9];

	char trans = 'T', notrans = 'N';
	char nojob = 'N';

	int nr, nc;
	int lld;
	int info;

	double Smax = log10(cnd)/2;
	double Smin = -Smax;

	double gap = (Smax-Smin)/(n-1);

	double* s;
			s=AllocateVector_double(n);

	double read_cnd = -1;

	if (mpi_rank < cprocs)
	{
		// Computation of local matrix size
		nc = numroc_( &n, &nb, &mycol, &i0, &npcol );
		nr = numroc_( &n, &nb, &myrow, &i0, &nprow );
		A1 = malloc(nr*nc*sizeof(double));
		A2 = malloc(nr*nc*sizeof(double));
		S  = malloc(nr*nc*sizeof(double));
		lld = MAX( 1 , nr );

		X = malloc(nr*1*sizeof(double));
		B = malloc(nr*1*sizeof(double));

		// Descriptors (local)
		descinit_( descA1, &n, &n, &nb, &nb, &i0, &i0, &context, &lld, &info );
		descinit_( descA2, &n, &n, &nb, &nb, &i0, &i0, &context, &lld, &info );
		descinit_( descS,  &n, &n, &nb, &nb, &i0, &i0, &context, &lld, &info );
		descinit_( descX,  &n, &i1, &nb, &nb, &i0, &i0, &context, &lld, &info );
		descinit_( descB,  &n, &i1, &nb, &nb, &i0, &i0, &context, &lld, &info );

		// Descriptors (global)
		if (mpi_rank==0)
		{
			descinit_( descA_global, &n, &n, &i1, &i1, &i0, &i0, &context_global, &n, &info );
			descinit_( descX_global, &n, &i1, &i1, &i1, &i0, &i0, &context_global, &n, &info );
			descinit_( descB_global, &n, &i1, &i1, &i1, &i0, &i0, &context_global, &n, &info );
		}
		else
		{
			// Descriptors (global, for non-root nodes)
			for (i=0; i<9; i++)
			{
				descA_global[i]=0;
				descX_global[i]=0;
				descB_global[i]=0;
			}
			descA_global[1]=-1;
			descX_global[1]=-1;
			descB_global[1]=-1;
		}

		int icol,irow,r,c;

		// distributed init to 0 for mat S
		pdlaset_("A", &n, &n, &d0, &d0, S, &i1, &i1, descS);
		// initialize in parallel the local parts of S
		// https://info.gwdg.de/wiki/doku.php?id=wiki:hpc:scalapack
		for (i=1; i<=n; i++)
		{
			r    = indxg2l_(&i,&nb,&i0,&i0,&nprow);
			irow = indxg2p_(&i,&nb,&i0,&i0,&nprow);
			c    = indxg2l_(&i,&nb,&i0,&i0,&npcol);
			icol = indxg2p_(&i,&nb,&i0,&i0,&npcol);
			if (myrow==irow && mycol==icol)
			{
				S[c-1+(r-1)*lld]=pow(10,Smin + (i-1)*gap);
			}
		}

		// A <- mat1.D.mat2 = (mat1.D).mat2

		int local_seed = seed+myrow*npcol+mycol;
		int lwork;
		double lazywork;
		double* work;
		double* tau;

		lwork=-1;
		pdgeqrf_( &n, &n, A1, &i1, &i1, descA1, NULL, &lazywork, &lwork, &info );
		lwork = (int)lazywork;
		work = malloc( lwork*sizeof(double) );
		tau = malloc( nc*sizeof(double) );

		// create mat1 -> A1
		RandomMatrix1D_double(A1, nr, nc, local_seed);

		//OrthogonalizeMatrix1D_double(A1, nr, nc); // to be parallelized
		pdgeqrf_( &n, &n,     A1, &i1, &i1, descA1, tau, work, &lwork, &info );
		pdorgqr_( &n, &n, &n, A1, &i1, &i1, descA1, tau, work, &lwork, &info );

		// mat1.D = A1.S -> A2
		pdgemm_("N", "N", &n, &n, &n, &d1, A1, &i1, &i1, descA1, S, &i1, &i1, descS, &d0, A2, &i1, &i1, descA2);

		// create mat2 -> A1
		RandomMatrix1D_double(A1, nr, nc, local_seed+npcol*nprow);

		//OrthogonalizeMatrix1D_double(A2, nr, nc);// to be parallelized
		pdgeqrf_( &n, &n,     A1, &i1, &i1, descA1, tau, work, &lwork, &info );
		pdorgqr_( &n, &n, &n, A1, &i1, &i1, descA1, tau, work, &lwork, &info );

		// (mat1.D).mat2 = A2.A1 -> S
		pdgemm_("N", "N", &n, &n, &n, &d1, A2, &i1, &i1, descA2, A1, &i1, &i1, descA1, &d0, S, &i1, &i1, descS);

		// transpose result S -> A1
		pdtran_(&n, &n, &d1, S, &i1, &i1, descS, &d0, A1, &i1, &i1, descA1);
		// get back A (A1 -> A)
		pdgemr2d_ (&n, &n, A1, &i1, &i1, descA1, A, &i1, &i1, descA_global, &context);

		// distributed init to 1 for vec X
		pdlaset_("A", &n, &i1, &d1, &d1, X, &i1, &i1, descX);

		// get back X
		pdgemr2d_ (&n, &i1, X, &i1, &i1, descX, x, &i1, &i1, descX_global, &context);

		// A.X -> B  (S.X -> B)
		pdgemm_(&trans, &notrans, &n, &i1, &n, &d1, A1, &i1, &i1, descA1, X, &i1, &i1, descX, &d0, B, &i1, &i1, descB);

		// get back B
		pdgemr2d_ (&n, &i1, B, &i1, &i1, descB, b, &i1, &i1, descB_global, &context);

		// use chunks in A1 to calc the condition number
		NULLFREE(work);

		if (cnd_readback)
		{
			lwork = -1;
			pdgesvd_ ( &nojob, &nojob, &n, &n, A1, &i1, &i1, descA1, s, NULL, &i1, &i1, NULL, NULL, &i1, &i1, NULL, &lazywork, &lwork, &info);
			lwork = (int)lazywork;
			work = AllocateVector_double(lwork);
			pdgesvd_ ( &nojob, &nojob, &n, &n, A1, &i1, &i1, descA1, s, NULL, &i1, &i1, NULL, NULL, &i1, &i1, NULL, work, &lwork, &info);

			read_cnd = s[0]/s[n-1];
		}

		DeallocateMatrix1D_double(A1);
		DeallocateMatrix1D_double(A2);
		DeallocateMatrix1D_double(S);
		DeallocateMatrix1D_double(X);
		DeallocateMatrix1D_double(B);
		NULLFREE(s);
		NULLFREE(work);
		NULLFREE(tau);
	}
	return read_cnd;
}

void MYblacs_scatter_double(int N, int M, double* A_glob, int nrows, int ncols, double* A_loc, int Nb, int Mb, int mpi_rank, int procrows, int proccols, int myrow, int mycol,  int ctxt)
{
	/*
	 * https://andyspiros.wordpress.com/2011/07/08/an-example-of-blacs-with-c/
	 */
	int sendr = 0, sendc = 0, recvr = 0, recvc = 0;
	int r,c;
	for (r = 0; r < N; r += Nb, sendr=(sendr+1)%procrows) {
		sendc = 0;
		// Number of rows to be sent
		// Is this the last row block?
		int nr = Nb;
		if (N-r < Nb)
			nr = N-r;

		for (c = 0; c < M; c += Mb, sendc=(sendc+1)%proccols) {
			// Number of cols to be sent
			// Is this the last col block?
			int nc = Mb;
			if (M-c < Mb)
				nc = M-c;

			if (mpi_rank==0) {
				// Send a nr-by-nc submatrix to process (sendr, sendc)
				Cdgesd2d(ctxt, nr, nc, A_glob+N*c+r, N, sendr, sendc);
			}

			if (myrow == sendr && mycol == sendc) {
				// Receive the same data
				// The leading dimension of the local matrix is nrows!
				Cdgerv2d(ctxt, nr, nc, A_loc+nrows*recvc+recvr, nrows, 0, 0);
				recvc = (recvc+nc)%ncols;
			}
		}
		if (myrow == sendr)
			recvr = (recvr+nr)%nrows;
	}
}

void MYblacs_gather_double(int N, int M, double* A_glob, int nrows, int ncols, double* A_loc, int Nb, int Mb, int mpi_rank, int procrows, int proccols, int myrow, int mycol,  int ctxt)
{
	/*
	 * https://andyspiros.wordpress.com/2011/07/08/an-example-of-blacs-with-c/
	 */
	int sendr = 0, sendc = 0, recvr = 0, recvc = 0;
	int r,c;
	for (r = 0; r < N; r += Nb, sendr=(sendr+1)%procrows) {
		sendc = 0;
		// Number of rows to be sent
		// Is this the last row block?
		int nr = Nb;
		if (N-r < Nb)
			nr = N-r;

		for (c = 0; c < M; c += Mb, sendc=(sendc+1)%proccols) {
			// Number of cols to be sent
			// Is this the last col block?
			int nc = Mb;
			if (M-c < Mb)
				nc = M-c;

			if (myrow == sendr && mycol == sendc) {
				// Send a nr-by-nc submatrix to process (sendr, sendc)
				Cdgesd2d(ctxt, nr, nc, A_loc+nrows*recvc+recvr, nrows, 0, 0);
				recvc = (recvc+nc)%ncols;
			}

			if (mpi_rank==0) {
				// Receive the same data
				// The leading dimension of the local matrix is nrows!
				Cdgerv2d(ctxt, nr, nc, A_glob+N*c+r, N, sendr, sendc);
			}
		}
		if (myrow == sendr)
			recvr = (recvr+nr)%nrows;
	}
}

double pGenSystemMatrices1D_double(int n, double* A, double* x, double* b, int seed, double cnd, char calc_cnd, char cnd_readback, int nb, int mpi_rank, int cprocs, int nprow, int npcol, int myrow, int mycol, int context, int context_global)
{
	// general
	int i;
	int i0 = 0;
	int i1 = 1;
	double d0 = 0.0;
	double d1 = 1.0;

	double* A1;
	double* A2;
	double* S;
	double* X;
	double* B;

	int descA1[9];
	int descA2[9];
	int descS[9];
	int descA_global[9];
	int descX[9];
	int descX_global[9];
	int descB[9];
	int descB_global[9];

	char trans = 'T', notrans = 'N';
	char nojob = 'N';

	int nr, nc;
	int lld;
	int info;

	double Smax = log10(cnd)/2;
	double Smin = -Smax;

	double gap = (Smax-Smin)/(n-1);

	double* s;
			s=AllocateVector_double(n);

	double read_cnd = -1;

	if (mpi_rank < cprocs)
	{
		// Computation of local matrix size
		nc = numroc_( &n, &nb, &mycol, &i0, &npcol );
		nr = numroc_( &n, &nb, &myrow, &i0, &nprow );
		A1 = malloc(nr*nc*sizeof(double));
		A2 = malloc(nr*nc*sizeof(double));
		S  = malloc(nr*nc*sizeof(double));
		lld = MAX( 1 , nr );

		X = malloc(nr*1*sizeof(double));
		B = malloc(nr*1*sizeof(double));

		// Descriptors (local)
		descinit_( descA1, &n, &n, &nb, &nb, &i0, &i0, &context, &lld, &info );
		descinit_( descA2, &n, &n, &nb, &nb, &i0, &i0, &context, &lld, &info );
		descinit_( descS,  &n, &n, &nb, &nb, &i0, &i0, &context, &lld, &info );
		descinit_( descX,  &n, &i1, &nb, &nb, &i0, &i0, &context, &lld, &info );
		descinit_( descB,  &n, &i1, &nb, &nb, &i0, &i0, &context, &lld, &info );

		// Descriptors (global)
		if (mpi_rank==0)
		{
			descinit_( descA_global, &n, &n, &i1, &i1, &i0, &i0, &context_global, &n, &info );
			descinit_( descX_global, &n, &i1, &i1, &i1, &i0, &i0, &context_global, &n, &info );
			descinit_( descB_global, &n, &i1, &i1, &i1, &i0, &i0, &context_global, &n, &info );
		}
		else
		{
			// Descriptors (global, for non-root nodes)
			for (i=0; i<9; i++)
			{
				descA_global[i]=0;
				descX_global[i]=0;
				descB_global[i]=0;
			}
			descA_global[1]=-1;
			descX_global[1]=-1;
			descB_global[1]=-1;
		}

		int icol,irow,r,c;
		int local_seed = seed+myrow*npcol+mycol;
		int lwork;
		double lazywork;
		double* work;
		double* tau;

		if (calc_cnd)
		{
			// distributed init to 0 for mat S
			pdlaset_("A", &n, &n, &d0, &d0, S, &i1, &i1, descS);
			// initialize in parallel the local parts of S
			// https://info.gwdg.de/wiki/doku.php?id=wiki:hpc:scalapack
			for (i=1; i<=n; i++)
			{
				r    = indxg2l_(&i,&nb,&i0,&i0,&nprow);
				irow = indxg2p_(&i,&nb,&i0,&i0,&nprow);
				c    = indxg2l_(&i,&nb,&i0,&i0,&npcol);
				icol = indxg2p_(&i,&nb,&i0,&i0,&npcol);
				if (myrow==irow && mycol==icol)
				{
					S[c-1+(r-1)*lld]=pow(10,Smin + (i-1)*gap);
				}
			}

			// A <- mat1.D.mat2 = (mat1.D).mat2
			lwork=-1;
			pdgeqrf_( &n, &n, A1, &i1, &i1, descA1, NULL, &lazywork, &lwork, &info );
			lwork = (int)lazywork;
			work = malloc( lwork*sizeof(double) );
			tau = malloc( nc*sizeof(double) );

			// create mat1 -> A1
			RandomMatrix1D_double(A1, nr, nc, local_seed);

			//OrthogonalizeMatrix1D_double(A1, nr, nc); // to be parallelized
			pdgeqrf_( &n, &n,     A1, &i1, &i1, descA1, tau, work, &lwork, &info );
			pdorgqr_( &n, &n, &n, A1, &i1, &i1, descA1, tau, work, &lwork, &info );

			// mat1.D = A1.S -> A2
			pdgemm_("N", "N", &n, &n, &n, &d1, A1, &i1, &i1, descA1, S, &i1, &i1, descS, &d0, A2, &i1, &i1, descA2);

			// create mat2 -> A1
			RandomMatrix1D_double(A1, nr, nc, local_seed+npcol*nprow);

			//OrthogonalizeMatrix1D_double(A2, nr, nc);// to be parallelized
			pdgeqrf_( &n, &n,     A1, &i1, &i1, descA1, tau, work, &lwork, &info );
			pdorgqr_( &n, &n, &n, A1, &i1, &i1, descA1, tau, work, &lwork, &info );

			// (mat1.D).mat2 = A2.A1 -> S
			pdgemm_("N", "N", &n, &n, &n, &d1, A2, &i1, &i1, descA2, A1, &i1, &i1, descA1, &d0, S, &i1, &i1, descS);

			// transpose result S -> A1
			pdtran_(&n, &n, &d1, S, &i1, &i1, descS, &d0, A1, &i1, &i1, descA1);
		}
		else
		{
			// create mat1 -> A1
			RandomMatrix1D_double(A1, nr, nc, local_seed);
			tau = NULL;
			work = NULL;
		}
		// get back A (A1 -> A)
		//pdgemr2d_ (&n, &n, A1, &i1, &i1, descA1, A, &i1, &i1, descA_global, &context);
		MYblacs_gather_double(n, n, A, nr, nc, A1, nb, nb, mpi_rank, nprow, npcol, myrow, mycol, context);

		// distributed init to 1 for vec X
		pdlaset_("A", &n, &i1, &d1, &d1, X, &i1, &i1, descX);

		// get back X
		//pdgemr2d_ (&n, &i1, X, &i1, &i1, descX, x, &i1, &i1, descX_global, &context);
		MYblacs_gather_double(n, i1, x, nr, i1, X, nb, nb, mpi_rank, nprow, npcol, myrow, mycol, context);

		// A.X -> B  (S.X -> B)
		pdgemm_(&trans, &notrans, &n, &i1, &n, &d1, A1, &i1, &i1, descA1, X, &i1, &i1, descX, &d0, B, &i1, &i1, descB);

		// get back B
		//pdgemr2d_ (&n, &i1, B, &i1, &i1, descB, b, &i1, &i1, descB_global, &context);
		MYblacs_gather_double(n, i1, b, nr, i1, B, nb, nb, mpi_rank, nprow, npcol, myrow, mycol, context);

		// use chunks in A1 to calc the condition number
		NULLFREE(work);

		if (cnd_readback)
		{
			lwork = -1;
			pdgesvd_ ( &nojob, &nojob, &n, &n, A1, &i1, &i1, descA1, s, NULL, &i1, &i1, NULL, NULL, &i1, &i1, NULL, &lazywork, &lwork, &info);
			lwork = (int)lazywork;
			work = AllocateVector_double(lwork);
			pdgesvd_ ( &nojob, &nojob, &n, &n, A1, &i1, &i1, descA1, s, NULL, &i1, &i1, NULL, NULL, &i1, &i1, NULL, work, &lwork, &info);

			read_cnd = s[0]/s[n-1];
		}

		DeallocateMatrix1D_double(A1);
		DeallocateMatrix1D_double(A2);
		DeallocateMatrix1D_double(S);
		DeallocateMatrix1D_double(X);
		DeallocateMatrix1D_double(B);
		NULLFREE(s);
		NULLFREE(work);
		NULLFREE(tau);
	}
	return read_cnd;
}

double pCheckSystemMatrices1D_double(int n, double* A, double* x, double* b, int nb, int mpi_rank, int cprocs, int nprow, int npcol, int myrow, int mycol, int context, int context_global)
{
	// general
	int i;
	int i0 = 0;
	int i1 = 1;
	double d0 = 0.0;
	double d1 = 1.0;

	double* A1;
	double* A2;

	int descA1[9];
	int descA2[9];
	int descA_global[9];

	char nojob = 'N';

	int nr, nc;
	int lld;
	int info;

	double* s;
			s=AllocateVector_double(n);

	double read_cnd = -1;

	if (mpi_rank < cprocs)
	{
		// Computation of local matrix size
		nc = numroc_( &n, &nb, &mycol, &i0, &npcol );
		nr = numroc_( &n, &nb, &myrow, &i0, &nprow );
		A1 = malloc(nr*nc*sizeof(double));
		A2 = malloc(nr*nc*sizeof(double));
		lld = MAX( 1 , nr );

		// Descriptors (local)
		descinit_( descA1, &n, &n, &nb, &nb, &i0, &i0, &context, &lld, &info );
		descinit_( descA2, &n, &n, &nb, &nb, &i0, &i0, &context, &lld, &info );

		// Descriptors (global)
		if (mpi_rank==0)
		{
			descinit_( descA_global, &n, &n, &i1, &i1, &i0, &i0, &context_global, &n, &info );
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

		int lwork;
		double lazywork;
		double* work;

		pdgemr2d_ (&n, &n, A, &i1, &i1, descA_global, A2, &i1, &i1, descA2, &context);
		pdtran_(&n, &n, &d1, A2, &i1, &i1, descA2, &d0, A1, &i1, &i1, descA1);

		lwork = -1;
		pdgesvd_ ( &nojob, &nojob, &n, &n, A1, &i1, &i1, descA1, s, NULL, &i1, &i1, NULL, NULL, &i1, &i1, NULL, &lazywork, &lwork, &info);
		lwork = (int)lazywork;
		work = AllocateVector_double(lwork);
		pdgesvd_ ( &nojob, &nojob, &n, &n, A1, &i1, &i1, descA1, s, NULL, &i1, &i1, NULL, NULL, &i1, &i1, NULL, work, &lwork, &info);

		read_cnd = s[0]/s[n-1];

		DeallocateMatrix1D_double(A1);
		DeallocateMatrix1D_double(A2);

		NULLFREE(s);
		NULLFREE(work);
	}
	return read_cnd;
}

double NormwiseRelativeError1D_double(double* mat, double* refmat, int rows, int cols)
{
	// errors clearly explained: https://www.cs.cornell.edu/~bindel/class/cs6210-f16/lec/2016-09-02.pdf
	// in LAPACK: https://www.netlib.org/lapack/lug/node78.html
	int i,j;
	double nre = 0;	// used also to signal if input matrix contains NaN (0=no, -1=yes)
	double* diffmat;
			diffmat = AllocateMatrix1D_double(rows, cols);
	double* work;
			work = AllocateVector_double(rows);
	char norm = 'F';
	for (i=0;i<rows;i++)
	{
		for (j=0;j<cols;j++)
		{
			if (isnan(mat[i*cols+j]))
			{
				nre=-1; // err. cannot be calculated
				i=rows; // prepares exit from loop i
				break;  // exits from loop j
			}
			else
			{
				diffmat[i*cols+j] = mat[i*cols+j] - refmat[i*cols+j];
			}
		}
	}
	if (nre==0) // input matrix does not contain NaN, then calc. err.
	{
		nre = dlange_(&norm, &rows, &cols, diffmat, &rows, work) / dlange_(&norm, &rows, &cols, refmat, &rows, work);
	}

	DeallocateMatrix1D_double(diffmat);
	DeallocateVector_double(work);

	return nre;
}
