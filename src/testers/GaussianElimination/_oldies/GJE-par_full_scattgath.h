#include <stdio.h>
#include <stdlib.h>
#include "../../helpers/types.h"
#include <mpi.h>

#include "../../helpers/matrix.h"

void pGaussianElimination_localmatrix(double** Alocal, double** A, double* b, int n, int rank, int nprocs)
{
    ui i,j,k;
    double* c;
    int* map;
    int* mymap;
    int mycols;

    c=AllocateVector_double(n);
    map=malloc(n*sizeof(int));
    mycols=n/nprocs;
    mymap=malloc(mycols*sizeof(int));

    for(i=0; i<n; i++)
	{
		map[i]= i % nprocs;
	}
    for(i=0; i<mycols; i++)
	{
    	mymap[i]= i * nprocs + rank;
	}

    /*
     * TODO
     */

    /*
    if (rank==0)
	{
		for(i=0; i<n; i++)
		{
			if (map[i]!=0)
			{
				MPI_Send (&A[i][0],n,MPI_DOUBLE, map[i],i, MPI_COMM_WORLD);
				MPI_Send (&b[i],1,MPI_DOUBLE, map[i],i+n, MPI_COMM_WORLD);
			}
		}
	}
	else
	{
		for(i=0; i<n; i++)
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

	}


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

     */

	free(map);
	free(mymap);
	DeallocateVector_double(c);
}


void pGaussianElimination_fullmatrix(double** A, double* b, int n, int rank, int nprocs)
{
    ui i,j,k;
    double* c;
    int* map;
    int numrows;

    c=AllocateVector_double(n);
    map=malloc(n*sizeof(int));

    numrows=n/nprocs;

    MPI_Datatype multiple_row, multiple_row_type;
	MPI_Datatype multiple_term, multiple_term_type;

	MPI_Type_vector (numrows, n, n * nprocs, MPI_DOUBLE , & multiple_row );
	MPI_Type_commit (& multiple_row);
	MPI_Type_create_resized(multiple_row, 0, n*sizeof(double), &multiple_row_type);
	MPI_Type_commit (& multiple_row_type);

	MPI_Type_vector (numrows, 1, nprocs, MPI_DOUBLE , & multiple_term );
	MPI_Type_commit (& multiple_term);
	MPI_Type_create_resized(multiple_term, 0, 1*sizeof(double), &multiple_term_type);
	MPI_Type_commit (& multiple_term_type);

    for(i=0; i<n; i++)
	{
		map[i]= i % nprocs;
	}

    /*
	if (rank==0)
	{
		for(i=0; i<n; i++)
		{
			if (map[i]!=0)
			{
				MPI_Send (&A[i][0],n,MPI_DOUBLE, map[i],i, MPI_COMM_WORLD);
				MPI_Send (&b[i],1,MPI_DOUBLE, map[i],i+n, MPI_COMM_WORLD);
			}
		}
	}
	else
	{
		for(i=0; i<n; i++)
		{
			if (rank==map[i])
			{
				MPI_Recv (&A[i][0],n,MPI_DOUBLE, 0,i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				MPI_Recv (&b[i],1,MPI_DOUBLE, 0,i+n, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			}
		}
	}
	*/

    MPI_Scatter(&A[0][0], 1, multiple_row_type, &A[rank][0], 1, multiple_row_type, 0, MPI_COMM_WORLD);

    MPI_Scatter(&b[0], 1, multiple_term_type, &b[rank], 1, multiple_term_type, 0, MPI_COMM_WORLD);


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


	/*
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
	*/

    MPI_Gather(&A[rank][0], 1, multiple_row_type, &A[0][0], 1, multiple_row_type, 0, MPI_COMM_WORLD);

    MPI_Gather(&b[rank], 1, multiple_term_type, &b[0], 1, multiple_term_type, 0, MPI_COMM_WORLD);

	free(map);
	DeallocateVector_double(c);
}

void BackSubstitution(double** A, double* b, double* x, int n)
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
