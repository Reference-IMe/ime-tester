/*
 * matrix_advanced.h
 *
 *  Created on: 12 Mar 2020
 *      Author: Marcello
 */
#include "matrix.h"
#include "vector.h"
#include "lapack.h"
#include "types.h"

#ifndef SRC_HELPERS_MATRIX_ADVANCED_H_
#define SRC_HELPERS_MATRIX_ADVANCED_H_

void OrthogonalizeMatrix1D(double* mat, int rows, int cols)
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

void RandomSquareMatrix1D_cnd(double* mat, int n, int seed, double cnd)
{
	double* mat1;
			mat1=AllocateMatrix1D(n, n);
	double* mat2;
			mat2=AllocateMatrix1D(n, n);

	double* s;
			s=AllocateVector(n);

	double smax = log10(cnd)/2;
	double smin = -smax;

	double gap;
			gap=(smax-smin)/(n-1);

	int r,c;
	for (r=0;r<n;r++)
	{
		s[r]=pow(10,smin + r*gap);
	}

	RandomMatrix1D(mat1, n, n, seed);
	OrthogonalizeMatrix1D(mat1, n, n);
	RandomMatrix1D(mat2, n, n, seed+1);
	OrthogonalizeMatrix1D(mat2, n, n);

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

	DeallocateVector(s);
	DeallocateMatrix1D(mat1);
	DeallocateMatrix1D(mat2);
}

double ConditionNumber1D(double* mat, int rows, int cols)
{
	// condition number calculated by definition:
	// the ratio 'cnd' of the largest to smallest singular value in the singular value decomposition of a matrix 'mat'.
	//
	// example for SVD in:
	// https://software.intel.com/sites/products/documentation/doclib/mkl_sa/11/mkl_lapack_examples/dgesvd_ex.c.htm

    int info, lwork;
    double cnd;
    double wkopt;
    double* work;
    double* u;
			u=AllocateMatrix1D(rows, rows);
	double* vt;
			vt=AllocateMatrix1D(cols, cols);
	double* s;
			s=AllocateVector(cols);

	lwork = -1;
	dgesvd_( "All", "All", &rows, &cols, mat, &rows, s, u, &rows, vt, &cols, &wkopt, &lwork, &info );
	lwork = (int)wkopt;
	work = AllocateVector(lwork);
	dgesvd_( "All", "All", &rows, &cols, mat, &rows, s, u, &rows, vt, &cols, work, &lwork, &info );
    if( info > 0 )
    {
		printf( "The algorithm computing SVD failed to converge.\n" );
		exit( 1 );
    }
    cnd = s[0]/s[cols-1];

	DeallocateVector(s);
	DeallocateVector(work);
	DeallocateMatrix1D(u);
	DeallocateMatrix1D(vt);

	return cnd;
}

double NormwiseRelativeError1D(double* mat, double* refmat, int rows, int cols)
{
	// errors clearly explained: https://www.cs.cornell.edu/~bindel/class/cs6210-f16/lec/2016-09-02.pdf
	// in LAPACK: https://www.netlib.org/lapack/lug/node78.html
	int i,j;
	double nre;
	double* diffmat;
			diffmat = AllocateMatrix1D(rows, cols);
	double* work;
			work = AllocateVector(rows);
	char norm = 'F';
	for (i=0;i<rows;i++)
	{
		for (j=0;j<cols;j++)
		{
			diffmat[i*cols+j] = mat[i*cols+j] - refmat[i*cols+j];
		}
	}
	nre = dlange_(&norm, &rows, &cols, diffmat, &rows, work) / dlange_(&norm, &rows, &cols, refmat, &rows, work);


	DeallocateMatrix1D(diffmat);
	DeallocateVector(work);

	return nre;
}

#endif /* SRC_HELPERS_MATRIX_ADVANCED_H_ */
