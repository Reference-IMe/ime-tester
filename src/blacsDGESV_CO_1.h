/*
 * blacsDGESV_CO_1.h
 *
 *  Created on: Jan 4, 2021
 *      Author: marcello
 *
 *      https://www.mathkeisan.com/UsersGuide/ScaLAPACK.cfm
 *
 *      Is there a way to extract the diagonal from a matrix with simple matrix operations
 *      https://math.stackexchange.com/questions/1527670/is-there-a-way-to-extract-the-diagonal-from-a-matrix-with-simple-matrix-operatio
 *      https://stackoverflow.com/questions/59007147/blas-routine-to-compute-diagonal-elements-only-of-a-matrix-product
 *
 *      rank-1 update
 *      http://www.netlib.org/scalapack/explore-html/d9/d15/pdger___8c_source.html
 *
 */

#include <mpi.h>
#include <time.h>
#include "helpers/macros.h"
#include "helpers/matrix.h"
#include "helpers/vector.h"
#include "helpers/Cblacs.h"
#include "helpers/scalapack.h"
#include "testers/tester_structures.h"


test_output blacsDGESV_CO_1(int n, double* A_global, int m, double* B_global, int nb,
								int mpi_rank, int cprocs,
								int nprow, int npcol, int myrow, int mycol,
								int context, int context_global)
{
	test_output result = EMPTY_OUTPUT;

	result.total_start_time = time(NULL);

	/*
	 * n = system rank (A_global n x n)
	 * m = num. of r.h.s (B_global n x m)
	 */

	// general
	int i,j;
	int i0 = 0;
	int i1 = 1;
	double d0 = 0.0;
	double d1 = 1.0;
	int info;
	//int *ipiv;

	// matrix
	int nr, nc;
	double *A;
	//double *At;
	double *C;
	double *T;
	//double *J;
	double *diag;
	double *lcol;

	int ncrhs, nrrhs;
	double *B;
	//int ncrhst, nrrhst;
	//double *Bt;

	int descA_global[9];
	int descB_global[9];
	int descA[9];
	//int descAt[9];
	int descC[9];
	int descT[9];
	//int descJ[9];
	int descB[9];
	//int descBt[9];

	int lld;
	//int lldj;
	//int lldt;

	//int ncj, nrj;

	result.core_start_time = time(NULL);

	if (mpi_rank < cprocs)
	{
		// Computation of local matrix size
		nc = numroc_( &n, &nb, &mycol, &i0, &npcol );
		nr = numroc_( &n, &nb, &myrow, &i0, &nprow );
		lld = MAX( 1 , nr );
		A  = malloc(nr*nc*sizeof(double));
		//At = malloc(nr*nc*sizeof(double));
		C = malloc(nr*nc*sizeof(double));
		T = malloc(nr*nc*sizeof(double));
		diag = malloc(nc*sizeof(double));
		lcol = malloc(nc*sizeof(double));
		/*
		ncj = numroc_( &i1, &nb, &mycol, &i0, &npcol );
		nrj = numroc_( &n,  &nb, &myrow, &i0, &nprow );
		lldj = MAX( 1 , nrj );
		J = malloc(ncj*nrj*sizeof(double));
		*/
		ncrhs = numroc_( &m, &nb, &mycol, &i0, &npcol );
		nrrhs = numroc_( &n, &nb, &myrow, &i0, &nprow );
		B = malloc(nrrhs*ncrhs*sizeof(double));

		//ncrhst = numroc_( &n, &nb, &mycol, &i0, &npcol );
		//nrrhst = numroc_( &m, &nb, &myrow, &i0, &nprow );
		//lldt = MAX( 1 , nrrhst );
		//Bt = malloc(nrrhst*ncrhst*sizeof(double));

		//ipiv = malloc((lld+nb)*sizeof(int));

		// Descriptors (local)
		descinit_( descA, &n, &n, &nb, &nb, &i0, &i0, &context, &lld, &info );
		//descinit_( descAt, &n, &n, &nb, &nb, &i0, &i0, &context, &lld, &info );

		descinit_( descC, &n, &n, &nb, &nb, &i0, &i0, &context, &lld, &info );
		descinit_( descT, &n, &n, &nb, &nb, &i0, &i0, &context, &lld, &info );
		//descinit_( descJ, &n, &i1, &nb, &nb, &i0, &i0, &context, &lldj, &info );

		descinit_( descB, &n, &m, &nb, &nb, &i0, &i0, &context, &lld, &info );
		//descinit_( descBt, &m, &n, &nb, &nb, &i0, &i0, &context, &lldt, &info );

		/*
		for (i=0; i<cprocs; i++)
		{
			MPI_Barrier(MPI_COMM_WORLD);
			if (mpi_rank==i)
			{
				printf("[%d](%d,%d): %d-%d\n",mpi_rank,myrow,mycol,nr,nc);
			}
		}
		*/

		if (mpi_rank==0)
		{
			// Descriptors (global)
			descinit_( descA_global, &n, &n, &i1, &i1, &i0, &i0, &context_global, &n, &info );
			descinit_( descB_global, &n, &m, &i1, &i1, &i0, &i0, &context_global, &n, &info );
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
		pdgemr2d_(&n, &n, A_global, &i1, &i1, descA_global, T, &i1, &i1, descT, &context);
		pdgemr2d_(&n, &m, B_global, &i1, &i1, descB_global, B, &i1, &i1, descB, &context);

		// transpose A, because scalapack expects column major order
		pdtran_(&n, &n, &d1, T, &i1, &i1, descT, &d0, A, &i1, &i1, descA);

		/*
		pdgemr2d_(&n, &n, A, &i1, &i1, descA, A_global, &i1, &i1, descA_global, &context);
		if (mpi_rank==0)
		{
			printf("(A')\n");
			PrintMatrix1D(A_global,n,n);
		}
		MPI_Barrier(MPI_COMM_WORLD);


		if (mpi_rank==1)
		{
			printf("\n[%d](%d,%d)\n",mpi_rank,myrow,mycol);
			printf("(Apart)\n");
			PrintMatrix1D(A,1,nc*nr);
		}
		*/

		// init X and K
		/*
		 * X=Diag(1/A)
		 * K=(A.X)' = X'.A'  ::: factors in reverse order !
		 *
		 * T is the superposition of X and K after having set K[i,i]=0 (i=1..n)
		 *
		 * or, equivalently
		 *
		 * T=(X'.A') after having set A[i,i]=1 (i=1..n)
		 */

		int icol,irow,r,c;

		// distributed init to 0 for mat X (= C)
		pdlaset_("A", &n, &n, &d0, &d0, C, &i1, &i1, descC);

		// initialize in parallel the local parts
		//   of X (diagonal elements to the reciprocal of same elements in A)
		//   of K (diagonal elements to "1")
		for (i=1; i<=n; i++)
		{
			// https://info.gwdg.de/wiki/doku.php?id=wiki:hpc:scalapack
			r    = indxg2l_(&i,&nb,&i0,&i0,&nprow) -1;
			irow = indxg2p_(&i,&nb,&i0,&i0,&nprow);
			c    = indxg2l_(&i,&nb,&i0,&i0,&npcol) -1;
			icol = indxg2p_(&i,&nb,&i0,&i0,&npcol);
			if (myrow==irow && mycol==icol)
			{
				// swap r c, because of column major order
				C[r+c*lld]=1/A[r+c*lld];	// X=Diag(1/A)
				A[r+c*lld]=1;				// A[i,i]=1 (i=1..n)
			}
		}

		// A.C = T'
		// A.X = K'
		pdgemm_("N", "N", &n, &n, &n, &d1, A, &i1, &i1, descA, C, &i1, &i1, descC, &d0, T, &i1, &i1, descT);
		// T actually contains T'

		/*
		pdgemr2d_(&n, &n, T, &i1, &i1, descT, A_global, &i1, &i1, descA_global, &context);
		if (mpi_rank==0)
		{
			printf("\n[%d]\n",n);
			printf("(T')\n");
			PrintMatrix1D(A_global,n,n);
		}
		*/

		//double d1_ = -1;
		int l;
		int l_1;

		for (l=n; l>1; l--)
		{
			l_1=l-1;

			// b*c in A
			pdgemm_("N", "N", &n, &l_1, &i1, &d1, T, &i1, &l, descT, T, &l, &i1, descT, &d0, A, &i1, &i1, descA);

			/*
			pdgemr2d_(&n, &n, A, &i1, &i1, descA, A_global, &i1, &i1, descA_global, &context);
			if (mpi_rank==0)
			{
				printf("\n[%d]\n",l-1);
				printf("(b*c)\n");
				PrintMatrix1D(A_global,n,n);
			}
			*/

			// save elements that will be overwritten
			//
			// !! swapped indices !! transposition
			//
			/*
			 * naive version
			 */
			/*
			for (i=1; i<=l-1; i++) // i should be rows, but treated as cols
			{
				// https://info.gwdg.de/wiki/doku.php?id=wiki:hpc:scalapack
				c    = indxg2l_(&i,&nb,&i0,&i0,&npcol) -1;
				icol = indxg2p_(&i,&nb,&i0,&i0,&npcol);
				if (mycol==icol)
				{
					// diagonal
					r    = indxg2l_(&i,&nb,&i0,&i0,&nprow) -1;
					irow = indxg2p_(&i,&nb,&i0,&i0,&nprow);
					if (myrow==irow)
					{
						diag[c]=T[r+c*lld];
					}
					// l-th column, treated as row
					r    = indxg2l_(&l,&nb,&i0,&i0,&nprow) -1;
					irow = indxg2p_(&l,&nb,&i0,&i0,&nprow);
					if (myrow==irow)
					{
						lcol[c]=T[r+c*lld];
					}
				}
			}
			*/
			/*
			 * faster version
			 */
			int rl;
			int irowl;

			// l-th column, treated as row
			rl    = indxg2l_(&l,&nb,&i0,&i0,&nprow) -1;
			irowl = indxg2p_(&l,&nb,&i0,&i0,&nprow);

			for (i=1; i<=l-1; i++) // i should be rows, but treated as cols
			{
				// https://info.gwdg.de/wiki/doku.php?id=wiki:hpc:scalapack
				c    = indxg2l_(&i,&nb,&i0,&i0,&npcol) -1;
				icol = indxg2p_(&i,&nb,&i0,&i0,&npcol);
				if (mycol==icol)
				{
					// diagonal
					r    = indxg2l_(&i,&nb,&i0,&i0,&nprow) -1;
					irow = indxg2p_(&i,&nb,&i0,&i0,&nprow);
					if (myrow==irow)
					{
						diag[c]=T[r+c*lld];
					}
					// l-th column, treated as row
					//r    = indxg2l_(&l,&nb,&i0,&i0,&nprow) -1;
					//irow = indxg2p_(&l,&nb,&i0,&i0,&nprow);
					if (myrow==irowl)
					{
						lcol[c]=T[rl+c*lld];
					}
				}
			}

			int r1;

			// a-b*c
			/*
			 * faster version ("if" removed)
			 */
			for (i=1; i<=l-1; i++) // "i" should be rows, but treated as cols
			{
				c    = indxg2l_(&i,&nb,&i0,&i0,&npcol) -1;
				icol = indxg2p_(&i,&nb,&i0,&i0,&npcol);
				if (mycol==icol)
				{
					for (r=0; r<nr; r++)
					{
						T[c*lld+r]=T[c*lld+r]-A[c*lld+r]; // difference
						r1=r+1;
						j = indxl2g_( &r1, &nb, &myrow, &i0, &nprow );
						if (i==j)
						{
							A[c*lld+r]=1/(1-A[c*lld+r]);
						}
						else
						{
							A[c*lld+r]=0;
						}
					}
				}
			}

			/*
			pdgemr2d_(&n, &n, T, &i1, &i1, descT, A_global, &i1, &i1, descA_global, &context);
			if (mpi_rank==0)
			{
				printf("(-)\n");
				PrintMatrix1D(A_global,l-1,n);
			}

			pdgemr2d_(&n, &n, A, &i1, &i1, descA, A_global, &i1, &i1, descA_global, &context);
			if (mpi_rank==0)
			{
				printf("(H)\n");
				PrintMatrix1D(A_global,n,n);
			}
			*/

			// restore elements
			//
			// !! swapped indices !! transposition
			//
			/*
			 * naive version
			 */
			/*
			for (i=1; i<=l-1; i++) // "i" should be rows, but treated as cols
			{
				// https://info.gwdg.de/wiki/doku.php?id=wiki:hpc:scalapack
				c    = indxg2l_(&i,&nb,&i0,&i0,&npcol) -1;
				icol = indxg2p_(&i,&nb,&i0,&i0,&npcol);
				if (mycol==icol)
				{
					// diagonal
					r    = indxg2l_(&i,&nb,&i0,&i0,&nprow) -1;
					irow = indxg2p_(&i,&nb,&i0,&i0,&nprow);
					if (myrow==irow)
					{
						T[r+c*lld]=diag[c];
					}
					// l-th column, treated as row
					r    = indxg2l_(&l,&nb,&i0,&i0,&nprow) -1;
					irow = indxg2p_(&l,&nb,&i0,&i0,&nprow);
					if (myrow==irow)
					{
						T[r+c*lld]=T[r+c*lld]-lcol[c];
					}
				}
			}
			*/
			/*
			 * faster version
			 */
			for (i=1; i<=l-1; i++) // i should be rows, but treated as cols
			{
				// https://info.gwdg.de/wiki/doku.php?id=wiki:hpc:scalapack
				c    = indxg2l_(&i,&nb,&i0,&i0,&npcol) -1;
				icol = indxg2p_(&i,&nb,&i0,&i0,&npcol);
				if (mycol==icol)
				{
					// diagonal
					r    = indxg2l_(&i,&nb,&i0,&i0,&nprow) -1;
					irow = indxg2p_(&i,&nb,&i0,&i0,&nprow);
					if (myrow==irow)
					{
						T[r+c*lld]=diag[c];
					}
					// l-th column, treated as row
					//r    = indxg2l_(&l,&nb,&i0,&i0,&nprow) -1;
					//irow = indxg2p_(&l,&nb,&i0,&i0,&nprow);
					if (myrow==irowl)
					{
						T[rl+c*lld]=T[rl+c*lld]-lcol[c];
					}
				}
			}

			/*
			pdgemr2d_(&n, &n, T, &i1, &i1, descT, A_global, &i1, &i1, descA_global, &context);
			if (mpi_rank==0)
			{
				printf("(K pre H)\n");
				PrintMatrix1D(A_global,n,n);
			}
			*/

			// topological formula
			pdgemm_("N", "N", &n, &l_1, &l_1, &d1, T, &i1, &i1, descT, A, &i1, &i1, descA, &d0, C, &i1, &i1, descC);
			pdlacpy_("A", &n, &l_1, C, &i1, &i1, descC, T, &i1, &i1, descT);

			/*
			pdgemr2d_(&n, &n, T, &i1, &i1, descT, A_global, &i1, &i1, descA_global, &context);
			if (mpi_rank==0)
			{
				printf("T(%d)\n",l-1);
				PrintMatrix1D(A_global,n,n);
			}
			*/

		}

		pdtrsm_ ("L", "U", "N", "U", &n, &i1, &d1, T, &i1, &i1, descT, B, &i1, &i1, descB);

		pdtrmm_ ("L", "L", "N", "N", &n, &i1, &d1, T, &i1, &i1, descT, B, &i1, &i1, descB);

		pdgemr2d_(&n, &m, B, &i1, &i1, descB, B_global, &i1, &i1, descB_global, &context);
		/*
		if (mpi_rank==0)
		{
			printf("C(%d)\n",l-1);
			PrintMatrix1D(B_global,n,m);
		}
		*/
		result.core_end_time = time(NULL);
		result.exit_code = 0;
	}
	else
	{
		A  = NULL;
		//At = NULL;
		B  = NULL;
		//Bt = NULL;

		//J = NULL;
		C = NULL;
		T = NULL;

		//ipiv = NULL;

		result.core_start_time = time(NULL);
		result.core_end_time = time(NULL);
		result.exit_code = 0;
	}

	// cleanup
	NULLFREE(A);
	//NULLFREE(At);
	NULLFREE(B);
	//NULLFREE(Bt);

	//NULLFREE(J);
	NULLFREE(C);
	NULLFREE(T);

	//NULLFREE(ipiv);

	MPI_Barrier(MPI_COMM_WORLD);

	result.total_end_time = time(NULL);

	return result;
}
