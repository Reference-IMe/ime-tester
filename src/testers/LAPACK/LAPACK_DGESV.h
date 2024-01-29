#include "../IMe/lib/src/helpers/matrix_basic.h"

/* DGESV prototype
 * extern void dgesv_( int* n, int* nrhs, double* a, int* lda, int* ipiv, double* b, int* ldb, int* info );
 */
extern void dgesv_( int* n, int* nrhs, double* a, int* lda, int* ipiv, double* b, int* ldb, int* info );

void Lapack_DGESV_calc(double* A, double* b, int n, int* ipiv)
{
	int rows;
	int one=1;
	int info;

	rows=n;

	dgesv_( &rows, &one, A, &rows, ipiv, b, &rows, &info );
}

void Lapack_DGESV(int n, double* A, int m, double* bb)
{
	int* ipiv = malloc(n * sizeof(int));
	int rows;
	int nrhs=m;
	int info;

	rows=n;
	dgesv_( &rows, &nrhs, A, &rows, ipiv, bb, &rows, &info );

	free(ipiv);
}

