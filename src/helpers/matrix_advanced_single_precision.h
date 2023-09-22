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

void OrthogonalizeMatrix1D_float(float* mat, int rows, int cols)
{
	int info;
    int lwork;
    	lwork = MIN(rows,cols);

    float* work;
			work = malloc(lwork*sizeof(float));
    float* tau;
    		tau = malloc(lwork*sizeof(float));

    sgeqrfp_(&rows, &cols, mat, &rows, tau, work, &lwork, &info);
    sorgqr_(&rows, &cols, &cols, mat, &rows, tau, work, &lwork, &info);

    NULLFREE(work);
    NULLFREE(tau);
}

void RandomSquareMatrixCND1D_float(float* mat, int n, int seed, float cnd)
{
	float* mat1;
			mat1=AllocateMatrix1D_float(n, n);
	float* mat2;
			mat2=AllocateMatrix1D_float(n, n);

	float* s;
			s=AllocateVector_float(n);

	float smax = log10(cnd)/2;
	float smin = -smax;

	float gap;
			gap=(smax-smin)/(n-1);

	int r,c;
	for (r=0;r<n;r++)
	{
		s[r]=pow(10,smin + r*gap);
	}

	RandomMatrix1D_float(mat1, n, n, seed);
	OrthogonalizeMatrix1D_float(mat1, n, n);
	RandomMatrix1D_float(mat2, n, n, seed+1);
	OrthogonalizeMatrix1D_float(mat2, n, n);

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
	float one = 1.0, zero = 0.0;
	sgemm_(&transA, &transB, &n, &n, &n, &one, mat1, &n, mat2, &n, &zero, mat, &n);

	DeallocateVector_float(s);
	DeallocateMatrix1D_float(mat1);
	DeallocateMatrix1D_float(mat2);
}

float ConditionNumber1D_float(float* mat, int rows, int cols)
{
	// condition number calculated by definition:
	// the ratio 'cnd' of the largest to smallest singular value in the singular value decomposition of a matrix 'mat'.
	//
	// example for SVD in:
	// https://software.intel.com/sites/products/documentation/doclib/mkl_sa/11/mkl_lapack_examples/dgesvd_ex.c.htm

    int info, lwork, i;
    float cnd;
    float wkopt;
    float* work;
    float* u=NULL;
//			u=AllocateMatrix1D_float(rows, rows);
	float* vt=NULL;
//			vt=AllocateMatrix1D_float(cols, cols);
	float* s;
			s=AllocateVector_float(cols);
	float* mat_copy;									// SVD algorithm destroys input matrix mat
			mat_copy=AllocateMatrix1D_float(rows, cols);		// http://www.netlib.org/lapack/explore-html/d1/d7e/group__float_g_esing_ga84fdf22a62b12ff364621e4713ce02f2.html

	char nojob='N'; // no need of calculating u and vt

	for (i=0;i<rows*cols;i++)
	{
		mat_copy[i]=mat[i];
	}

	lwork = -1;
	sgesvd_( &nojob, &nojob, &rows, &cols, mat_copy, &rows, s, u, &rows, vt, &cols, &wkopt, &lwork, &info );
	lwork = (int)wkopt;
	work = AllocateVector_float(lwork);
	sgesvd_( &nojob, &nojob, &rows, &cols, mat_copy, &rows, s, u, &rows, vt, &cols, work, &lwork, &info );
    if( info > 0 )
    {
		printf( "The algorithm computing SVD failed to converge.\n" );
		exit( 1 );
    }
    cnd = s[0]/s[cols-1];

    DeallocateMatrix1D_float(mat_copy);
	DeallocateVector_float(s);
	DeallocateVector_float(work);
	DeallocateMatrix1D_float(u);
	DeallocateMatrix1D_float(vt);

	return cnd;
}

float GenSystemMatrices1D_float(int n, float* A, float* x, float* b, int seed, float cnd, char calc_cnd, char cnd_readback)
{
	int i1 = 1;
	float d0 = 0.0;
	float d1 = 1.0;
	char transA = 'T', transx = 'N';
	float read_cnd = -1;

	if (calc_cnd)
	{
		RandomSquareMatrixCND1D_float(A, n, seed, cnd);
	}
	else
	{
		RandomMatrix1D_float(A, n, n, seed);
	}
	FillVector_float(x, n, 1);
	sgemm_(&transA, &transx, &n, &i1, &n, &d1, A, &n, x, &n, &d0, b, &n);
	if (cnd_readback)
	{
		read_cnd = ConditionNumber1D_float(A, n, n);
	}
	return read_cnd;
}

float pGenSystemMatrices1D_float_psgemr2d(int n, float* A, float* x, float* b, int seed, float cnd, char cnd_readback, int nb, int mpi_rank, int cprocs, int nprow, int npcol, int myrow, int mycol, int context, int context_global)
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
	 * see "pGenSystemMatrices1D_float" below
	 *
	 */

	// general
	int i;
	int i0 = 0;
	int i1 = 1;
	float d0 = 0.0;
	float d1 = 1.0;

	float* A1;
	float* A2;
	float* S;
	float* X;
	float* B;

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

	float Smax = log10(cnd)/2;
	float Smin = -Smax;

	float gap = (Smax-Smin)/(n-1);

	float* s;
			s=AllocateVector_float(n);

	float read_cnd = -1;

	if (mpi_rank < cprocs)
	{
		// Computation of local matrix size
		nc = numroc_( &n, &nb, &mycol, &i0, &npcol );
		nr = numroc_( &n, &nb, &myrow, &i0, &nprow );
		A1 = malloc(nr*nc*sizeof(float));
		A2 = malloc(nr*nc*sizeof(float));
		S  = malloc(nr*nc*sizeof(float));
		lld = MAX( 1 , nr );

		X = malloc(nr*1*sizeof(float));
		B = malloc(nr*1*sizeof(float));

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
		pslaset_("A", &n, &n, &d0, &d0, S, &i1, &i1, descS);
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
		float lazywork;
		float* work;
		float* tau;

		lwork=-1;
		psgeqrf_( &n, &n, A1, &i1, &i1, descA1, NULL, &lazywork, &lwork, &info );
		lwork = (int)lazywork;
		work = malloc( lwork*sizeof(float) );
		tau = malloc( nc*sizeof(float) );

		// create mat1 -> A1
		RandomMatrix1D_float(A1, nr, nc, local_seed);

		//OrthogonalizeMatrix1D_float(A1, nr, nc); // to be parallelized
		psgeqrf_( &n, &n,     A1, &i1, &i1, descA1, tau, work, &lwork, &info );
		psorgqr_( &n, &n, &n, A1, &i1, &i1, descA1, tau, work, &lwork, &info );

		// mat1.D = A1.S -> A2
		psgemm_("N", "N", &n, &n, &n, &d1, A1, &i1, &i1, descA1, S, &i1, &i1, descS, &d0, A2, &i1, &i1, descA2);

		// create mat2 -> A1
		RandomMatrix1D_float(A1, nr, nc, local_seed+npcol*nprow);

		//OrthogonalizeMatrix1D_float(A2, nr, nc);// to be parallelized
		psgeqrf_( &n, &n,     A1, &i1, &i1, descA1, tau, work, &lwork, &info );
		psorgqr_( &n, &n, &n, A1, &i1, &i1, descA1, tau, work, &lwork, &info );

		// (mat1.D).mat2 = A2.A1 -> S
		psgemm_("N", "N", &n, &n, &n, &d1, A2, &i1, &i1, descA2, A1, &i1, &i1, descA1, &d0, S, &i1, &i1, descS);

		// transpose result S -> A1
		pstran_(&n, &n, &d1, S, &i1, &i1, descS, &d0, A1, &i1, &i1, descA1);
		// get back A (A1 -> A)
		psgemr2d_ (&n, &n, A1, &i1, &i1, descA1, A, &i1, &i1, descA_global, &context);

		// distributed init to 1 for vec X
		pslaset_("A", &n, &i1, &d1, &d1, X, &i1, &i1, descX);

		// get back X
		psgemr2d_ (&n, &i1, X, &i1, &i1, descX, x, &i1, &i1, descX_global, &context);

		// A.X -> B  (S.X -> B)
		psgemm_(&trans, &notrans, &n, &i1, &n, &d1, A1, &i1, &i1, descA1, X, &i1, &i1, descX, &d0, B, &i1, &i1, descB);

		// get back B
		psgemr2d_ (&n, &i1, B, &i1, &i1, descB, b, &i1, &i1, descB_global, &context);

		// use chunks in A1 to calc the condition number
		NULLFREE(work);

		if (cnd_readback)
		{
			lwork = -1;
			psgesvd_ ( &nojob, &nojob, &n, &n, A1, &i1, &i1, descA1, s, NULL, &i1, &i1, NULL, NULL, &i1, &i1, NULL, &lazywork, &lwork, &info);
			lwork = (int)lazywork;
			work = AllocateVector_float(lwork);
			psgesvd_ ( &nojob, &nojob, &n, &n, A1, &i1, &i1, descA1, s, NULL, &i1, &i1, NULL, NULL, &i1, &i1, NULL, work, &lwork, &info);

			read_cnd = s[0]/s[n-1];
		}

		DeallocateMatrix1D_float(A1);
		DeallocateMatrix1D_float(A2);
		DeallocateMatrix1D_float(S);
		DeallocateMatrix1D_float(X);
		DeallocateMatrix1D_float(B);
		NULLFREE(s);
		NULLFREE(work);
		NULLFREE(tau);
	}
	return read_cnd;
}

void MYblacs_scatter_float(int N, int M, float* A_glob, int nrows, int ncols, float* A_loc, int Nb, int Mb, int mpi_rank, int procrows, int proccols, int myrow, int mycol,  int ctxt)
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
				Csgesd2d(ctxt, nr, nc, A_glob+N*c+r, N, sendr, sendc);
			}

			if (myrow == sendr && mycol == sendc) {
				// Receive the same data
				// The leading dimension of the local matrix is nrows!
				Csgerv2d(ctxt, nr, nc, A_loc+nrows*recvc+recvr, nrows, 0, 0);
				recvc = (recvc+nc)%ncols;
			}
		}
		if (myrow == sendr)
			recvr = (recvr+nr)%nrows;
	}
}

void MYblacs_gather_float(int N, int M, float* A_glob, int nrows, int ncols, float* A_loc, int Nb, int Mb, int mpi_rank, int procrows, int proccols, int myrow, int mycol,  int ctxt)
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
				Csgesd2d(ctxt, nr, nc, A_loc+nrows*recvc+recvr, nrows, 0, 0);
				recvc = (recvc+nc)%ncols;
			}

			if (mpi_rank==0) {
				// Receive the same data
				// The leading dimension of the local matrix is nrows!
				Csgerv2d(ctxt, nr, nc, A_glob+N*c+r, N, sendr, sendc);
			}
		}
		if (myrow == sendr)
			recvr = (recvr+nr)%nrows;
	}
}

float pGenSystemMatrices1D_float(int n, float* A, float* x, float* b, int seed, float cnd, char calc_cnd, char cnd_readback, int nb, int mpi_rank, int cprocs, int nprow, int npcol, int myrow, int mycol, int context, int context_global)
{
	// general
	int i;
	int i0 = 0;
	int i1 = 1;
	float d0 = 0.0;
	float d1 = 1.0;

	float* A1;
	float* A2;
	float* S;
	float* X;
	float* B;

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

	float Smax = log10(cnd)/2;
	float Smin = -Smax;

	float gap = (Smax-Smin)/(n-1);

	float* s;
			s=AllocateVector_float(n);

	float read_cnd = -1;

	if (mpi_rank < cprocs)
	{
		// Computation of local matrix size
		nc = numroc_( &n, &nb, &mycol, &i0, &npcol );
		nr = numroc_( &n, &nb, &myrow, &i0, &nprow );
		A1 = malloc(nr*nc*sizeof(float));
		A2 = malloc(nr*nc*sizeof(float));
		S  = malloc(nr*nc*sizeof(float));
		lld = MAX( 1 , nr );

		X = malloc(nr*1*sizeof(float));
		B = malloc(nr*1*sizeof(float));

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
		float lazywork;
		float* work;
		float* tau;

		if (calc_cnd)
		{
			// distributed init to 0 for mat S
			pslaset_("A", &n, &n, &d0, &d0, S, &i1, &i1, descS);
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
			psgeqrf_( &n, &n, A1, &i1, &i1, descA1, NULL, &lazywork, &lwork, &info );
			lwork = (int)lazywork;
			work = malloc( lwork*sizeof(float) );
			tau = malloc( nc*sizeof(float) );

			// create mat1 -> A1
			RandomMatrix1D_float(A1, nr, nc, local_seed);

			//OrthogonalizeMatrix1D_float(A1, nr, nc); // to be parallelized
			psgeqrf_( &n, &n,     A1, &i1, &i1, descA1, tau, work, &lwork, &info );
			psorgqr_( &n, &n, &n, A1, &i1, &i1, descA1, tau, work, &lwork, &info );

			// mat1.D = A1.S -> A2
			psgemm_("N", "N", &n, &n, &n, &d1, A1, &i1, &i1, descA1, S, &i1, &i1, descS, &d0, A2, &i1, &i1, descA2);

			// create mat2 -> A1
			RandomMatrix1D_float(A1, nr, nc, local_seed+npcol*nprow);

			//OrthogonalizeMatrix1D_float(A2, nr, nc);// to be parallelized
			psgeqrf_( &n, &n,     A1, &i1, &i1, descA1, tau, work, &lwork, &info );
			psorgqr_( &n, &n, &n, A1, &i1, &i1, descA1, tau, work, &lwork, &info );

			// (mat1.D).mat2 = A2.A1 -> S
			psgemm_("N", "N", &n, &n, &n, &d1, A2, &i1, &i1, descA2, A1, &i1, &i1, descA1, &d0, S, &i1, &i1, descS);

			// transpose result S -> A1
			pstran_(&n, &n, &d1, S, &i1, &i1, descS, &d0, A1, &i1, &i1, descA1);
		}
		else
		{
			// create mat1 -> A1
			RandomMatrix1D_float(A1, nr, nc, local_seed);
			tau = NULL;
			work = NULL;
		}
		// get back A (A1 -> A)
		//pdgemr2d_ (&n, &n, A1, &i1, &i1, descA1, A, &i1, &i1, descA_global, &context);
		MYblacs_gather_float(n, n, A, nr, nc, A1, nb, nb, mpi_rank, nprow, npcol, myrow, mycol, context);

		// distributed init to 1 for vec X
		pslaset_("A", &n, &i1, &d1, &d1, X, &i1, &i1, descX);

		// get back X
		//pdgemr2d_ (&n, &i1, X, &i1, &i1, descX, x, &i1, &i1, descX_global, &context);
		MYblacs_gather_float(n, i1, x, nr, i1, X, nb, nb, mpi_rank, nprow, npcol, myrow, mycol, context);

		// A.X -> B  (S.X -> B)
		psgemm_(&trans, &notrans, &n, &i1, &n, &d1, A1, &i1, &i1, descA1, X, &i1, &i1, descX, &d0, B, &i1, &i1, descB);

		// get back B
		//pdgemr2d_ (&n, &i1, B, &i1, &i1, descB, b, &i1, &i1, descB_global, &context);
		MYblacs_gather_float(n, i1, b, nr, i1, B, nb, nb, mpi_rank, nprow, npcol, myrow, mycol, context);

		// use chunks in A1 to calc the condition number
		NULLFREE(work);

		if (cnd_readback)
		{
			lwork = -1;
			psgesvd_ ( &nojob, &nojob, &n, &n, A1, &i1, &i1, descA1, s, NULL, &i1, &i1, NULL, NULL, &i1, &i1, NULL, &lazywork, &lwork, &info);
			lwork = (int)lazywork;
			work = AllocateVector_float(lwork);
			psgesvd_ ( &nojob, &nojob, &n, &n, A1, &i1, &i1, descA1, s, NULL, &i1, &i1, NULL, NULL, &i1, &i1, NULL, work, &lwork, &info);

			read_cnd = s[0]/s[n-1];
		}

		DeallocateMatrix1D_float(A1);
		DeallocateMatrix1D_float(A2);
		DeallocateMatrix1D_float(S);
		DeallocateMatrix1D_float(X);
		DeallocateMatrix1D_float(B);
		NULLFREE(s);
		NULLFREE(work);
		NULLFREE(tau);
	}
	return read_cnd;
}

float pCheckSystemMatrices1D_float(int n, float* A, float* x, float* b, int nb, int mpi_rank, int cprocs, int nprow, int npcol, int myrow, int mycol, int context, int context_global)
{
	// general
	int i;
	int i0 = 0;
	int i1 = 1;
	float d0 = 0.0;
	float d1 = 1.0;

	float* A1;
	float* A2;

	int descA1[9];
	int descA2[9];
	int descA_global[9];

	char nojob = 'N';

	int nr, nc;
	int lld;
	int info;

	float* s;
			s=AllocateVector_float(n);

	float read_cnd = -1;

	if (mpi_rank < cprocs)
	{
		// Computation of local matrix size
		nc = numroc_( &n, &nb, &mycol, &i0, &npcol );
		nr = numroc_( &n, &nb, &myrow, &i0, &nprow );
		A1 = malloc(nr*nc*sizeof(float));
		A2 = malloc(nr*nc*sizeof(float));
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
		float lazywork;
		float* work;

		psgemr2d_ (&n, &n, A, &i1, &i1, descA_global, A2, &i1, &i1, descA2, &context);
		pstran_(&n, &n, &d1, A2, &i1, &i1, descA2, &d0, A1, &i1, &i1, descA1);

		lwork = -1;
		psgesvd_ ( &nojob, &nojob, &n, &n, A1, &i1, &i1, descA1, s, NULL, &i1, &i1, NULL, NULL, &i1, &i1, NULL, &lazywork, &lwork, &info);
		lwork = (int)lazywork;
		work = AllocateVector_float(lwork);
		psgesvd_ ( &nojob, &nojob, &n, &n, A1, &i1, &i1, descA1, s, NULL, &i1, &i1, NULL, NULL, &i1, &i1, NULL, work, &lwork, &info);

		read_cnd = s[0]/s[n-1];

		DeallocateMatrix1D_float(A1);
		DeallocateMatrix1D_float(A2);

		NULLFREE(s);
		NULLFREE(work);
	}
	return read_cnd;
}

float NormwiseRelativeError1D_float(float* mat, float* refmat, int rows, int cols)
{
	// errors clearly explained: https://www.cs.cornell.edu/~bindel/class/cs6210-f16/lec/2016-09-02.pdf
	// in LAPACK: https://www.netlib.org/lapack/lug/node78.html
	int i,j;
	float nre = 0;	// used also to signal if input matrix contains NaN (0=no, -1=yes)
	float* diffmat;
			diffmat = AllocateMatrix1D_float(rows, cols);
	float* work;
			work = AllocateVector_float(rows);
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
		nre = slange_(&norm, &rows, &cols, diffmat, &rows, work) / slange_(&norm, &rows, &cols, refmat, &rows, work);
	}

	DeallocateMatrix1D_float(diffmat);
	DeallocateVector_float(work);

	return nre;
}
