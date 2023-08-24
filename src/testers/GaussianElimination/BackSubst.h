#include "../../helpers/vector.h"

/*
 * Back substitution
 *
 * A with rank n
 * bb multiple r.h.s in m columns
 * xx corresponding solutions in m columns
 *
 */
#ifndef __BACKSUBST_H__
#define __BACKSUBST_H__

void BackSubstitution(int n, double** A, int m, double** bb, double** xx)
{
		int i,j,k;
		double* sum;

		sum=AllocateVector_double(m);

		for(j=0;j<m;j++)
		{
			xx[n-1][j]=bb[n-1][j]/A[n-1][n-1];
		}
		for(i=n-2;i>=0;i--)
		{
			for(j=0;j<m;j++)
			{
				sum[j]=0.0;

				for(k=i+1;k<n;k++)
				{
					sum[j]=sum[j]+A[i][k]*xx[k][j];
				}
				xx[i][j]=(bb[i][j]-sum[j])/A[i][i];
			}
		}
		DeallocateVector_double(sum);
}

#endif
