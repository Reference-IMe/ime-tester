#include "../../helpers/vector.h"

/*
 * Gaussian elimination
 *
 * A with rank n
 * bb multiple r.h.s in m columns
 *
 */
void GaussianElimination(int n, double** A, int m, double** bb)
{
    int i,j,k;
    double* c;

    c=AllocateVector(n);

	for(k=0; k<n; k++)
	{
		for(i= k+1; i<n; i++)
		{
			c[i]=A[i][k]/A[k][k];
		}
		for(i= k+1; i<n; i++)
		{
			for(j=0; j<n; j++)
			{
				A[i][j]=A[i][j]-( c[i]*A[k][j] );
			}
			for(j=0; j<m; j++)
			{
				bb[i][j]=bb[i][j]-( c[i]*bb[k][j] );
			}
		}
	}
	DeallocateVector(c);
}
