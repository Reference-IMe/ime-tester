#include <stdio.h>
#include <stdlib.h>
#include "../../helpers/types.h"

void GaussianElimination(cui n, double** A, cui m, double** b)
{
    ui i,j,k;
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
					b[i][j]=b[i][j]-( c[i]*b[k][j] );
				}
			}
		}
		DeallocateVector(c);
}

void BackSubstitution(cui n, double** A, cui m, double** b, double** x)
{
		int i,j,k;
		double* sum;

		sum=AllocateVector(m);

		for(j=0;j<m;j++)
		{
			x[n-1][j]=b[n-1][j]/A[n-1][n-1];
		}
		for(i=n-2;i>=0;i--)
		{
			for(j=0;j<m;j++)
			{
				sum[j]=0.0;

				for(k=i+1;k<n;k++)
				{
					sum[j]=sum[j]+A[i][k]*x[k][j];
				}
				x[i][j]=(b[i][j]-sum[j])/A[i][i];
			}
		}
}
