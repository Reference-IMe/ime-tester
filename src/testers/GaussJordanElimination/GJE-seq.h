#include <stdio.h>
#include <stdlib.h>
#include "../../helpers/types.h"

void GaussianElimination(double** A, double* b, cui n)
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
				b[i]=b[i]-( c[i]*b[k] );
			}
		}
		DeallocateVector(c);
}

void BackSubstitution(double** A, double* b, double* x, cui n)
{
		double sum;
		int i,j;

		x[n-1]=b[n-1]/A[n-1][n-1];
		for(i=n-2;i>=0;i--)
		{
			sum=0.0;

			for(j=i+1;j<n;j++)
			{
				sum=sum+A[i][j]*x[j];
			}
			x[i]=(b[i]-sum)/A[i][i];
		}
}
