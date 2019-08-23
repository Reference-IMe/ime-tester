#include "../../helpers/types.h"

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
