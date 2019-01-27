#include <stdio.h>
#include <stdlib.h>

/* DGESV prototype */
//extern void dgesv_( int* n, int* nrhs, double* a, int* lda, int* ipiv, double* b, int* ldb, int* info );
extern void dgesv_( int* n, int* nrhs, double* a, int* lda, int* ipiv, double* b, int* ldb, int* info );

void LapackDGESV_calc(double* A, double* b, cui n, int* ipiv)
{
	int rows;
	int one=1;
	int info;

	rows=n;

	dgesv_( &rows, &one, A, &rows, ipiv, b, &rows, &info );
}

void LapackDGESV(double* A, double* b, cui n)
{
	int* ipiv = malloc(n * sizeof(int));

	LapackDGESV_calc(A, b, n, ipiv);

	free(ipiv);
}

