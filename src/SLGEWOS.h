#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "helpers/selfie.h"
#include "helpers/matrix.h"
#include "helpers/vector.h"

void SLGEWOS(double** A, double* b, double* s, int n)
{
    int i,j,k,l,rep;

    double** X;
    double** K;
    double*  H;
    double*  F;

    int rows=n;
    int cols=n;

    X=AllocateMatrix2D(rows,cols,CONTIGUOUS);
    K=AllocateMatrix2D(rows,cols,CONTIGUOUS);

    H=AllocateVector(rows);
    F=AllocateVector(rows);

		for (i=0;i<rows;i++)
		{
			for (j=0;j<cols;j++)
			{
				if (i==j)
				{
					X[i][j]=1/A[i][j];
					//K[i][j]=1;
				}
				else
				{
					X[i][j]=0.0;
					K[i][j]=A[j][i]/A[i][i];
				}
			}
			F[i]=b[i];
			s[i]=0.0;
		}

		for (l=rows-1; l>0; l--)
		{
			for (i=0; i<l; i++)
			{
				H[i]=1/(1-K[i][l]*K[l][i]);
				for (j=0; j<cols; j++)
				{
					X[i][j]=H[i]*(X[i][j]-K[i][l]*X[l][j]);
				}
				for (j=0; j<l; j++)
				{
					if (j!=i)
					{
						K[i][j]=H[i]*(K[i][j]-K[i][l]*K[l][j]);
					}
				}
			}
		}

		for (i=rows-2; i>=0; i--)
		{
			for (l=i+1; l<rows; l++)
			{
				F[i]=F[i]-F[l]*K[l][i];
			}
		}

		for (j=0; j<cols; j++)
		{
			for (l=0; l<rows; l++)
			{
				s[j]=F[l]*X[l][j]+s[j];
			}
		}

    DeallocateMatrix2D(X,rows,CONTIGUOUS);
    DeallocateMatrix2D(K,rows,CONTIGUOUS);

    DeallocateVector(H);
    DeallocateVector(F);
}
