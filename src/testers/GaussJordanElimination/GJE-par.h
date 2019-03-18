#include <stdio.h>
#include <stdlib.h>
#include "../../helpers/types.h"
#include <mpi.h>

void pGaussianElimination(double** A, double* b, cui n, int rank, int nprocs)
{
    ui i,j,k;
    double* c;
    int* map;

    c=AllocateVector(n);
    map=malloc(n*sizeof(int));

    //MPI_Bcast (&A[0][0],n*n,MPI_DOUBLE,0,MPI_COMM_WORLD);
	//MPI_Bcast (b,n,MPI_DOUBLE,0,MPI_COMM_WORLD);

	for(i=0; i<n; i++)
	{
		map[i]= i % nprocs;
		if (rank==0)
		{
			MPI_Send (&A[i][0],n,MPI_DOUBLE, map[i],i, MPI_COMM_WORLD);
			MPI_Send (&b[i],1,MPI_DOUBLE, map[i],i+n, MPI_COMM_WORLD);
		}
		else
		{
			if (rank==map[i])
			{
				MPI_Recv (&A[i][0],n,MPI_DOUBLE, 0,i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				MPI_Recv (&b[i],1,MPI_DOUBLE, 0,i+n, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			}
		}
	}

	for(k=0;k<n;k++)
	{
		MPI_Bcast (&A[k][k],n-k,MPI_DOUBLE,map[k],MPI_COMM_WORLD);
		MPI_Bcast (&b[k],1,MPI_DOUBLE,map[k],MPI_COMM_WORLD);
		for(i= k+1; i<n; i++)
		{
			if(map[i] == rank)
			{
				c[i]=A[i][k]/A[k][k];
			}
		}
		for(i= k+1; i<n; i++)
		{
			if(map[i] == rank)
			{
				//printf("I'm %d\n",rank);
				for(j=0;j<n;j++)
				{
					A[i][j]=A[i][j]-( c[i]*A[k][j] );
				}
				b[i]=b[i]-( c[i]*b[k] );
			}
		}
		/*
		if(ft!=0)
		{
			strcpy(filename, "/vagrant/tmp/CHK ");
			sprintf(buffer,"%d",rank);
			strcat(filename, buffer);
			strcat(filename, " ITER ");
			sprintf(buffer,"%d",k);
			strcat(filename, buffer);
			printf("Node %d checkpoint at iteration %d to file %s \n",rank,k,filename);
			fp = fopen(filename, "w");
			fwrite(&A[k][k], sizeof(range), n, fp);
			fwrite(&b[k], sizeof(range), 1, fp);
			fclose(fp);
		}
		*/

	}

	MPI_Barrier(MPI_COMM_WORLD);

	for(i=0; i<n; i++)
	{
		if (rank==0)
		{
			if(map[i]!=0)
			{
				MPI_Recv (&A[i][0], n, MPI_DOUBLE, map[i], i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				MPI_Recv (&b[i], 1, MPI_DOUBLE, map[i], i+n, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			}
		}
		else
		{
			if (rank==map[i])
			{
				MPI_Send (&A[i][0], n, MPI_DOUBLE, 0,i, MPI_COMM_WORLD);
				MPI_Send (&b[i], 1, MPI_DOUBLE, 0, i+n, MPI_COMM_WORLD);
			}
		}
	}

	free(map);
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
