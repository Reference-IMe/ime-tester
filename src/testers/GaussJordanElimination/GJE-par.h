#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../../helpers/types.h"
#include <mpi.h>

#include "../../helpers/matrix.h"

#include <unistd.h>

void pGaussianElimination_partialmatrix(double** A, double* b, cui n, int rank, int nprocs)
{
    ui i,j,k;

    int myrows;
    myrows=n/nprocs;
    int mystart=0;

    double* c;
    c=AllocateVector(n);

	// row for pivoting
    double* Abase;
    Abase=malloc(n*sizeof(double));

    // local storage for a part of the input matrix (continuous rows, not interleaved)
    double** Alocal;
    Alocal=AllocateMatrix2D(myrows, n, CONTIGUOUS);

    int* map;
    map=malloc(n*sizeof(int));

    int* local;
    local=malloc(n*sizeof(int));

    for(i=0; i<n; i++)
	{
		map[i]= i % nprocs;			// who has the row i
		local[i]=floor(i/nprocs);	// position of the row i(global) in the local matrix
	}

    int* global;
    global=malloc(n*sizeof(int));

    for(i=0; i<myrows; i++)
	{
    	global[i]= i * nprocs + rank; // position of the row i(local) in the global matrix
	}

    // derived data types
	MPI_Datatype multiple_row;
	MPI_Type_vector (myrows, n, n * nprocs , MPI_DOUBLE, & multiple_row );
	MPI_Type_commit (& multiple_row);

	MPI_Datatype multiple_row_contiguous;
	MPI_Type_vector (myrows, n, n , MPI_DOUBLE, & multiple_row_contiguous );
	MPI_Type_commit (& multiple_row_contiguous);

	MPI_Datatype multiple_row_type;
	MPI_Type_create_resized (multiple_row, 0, n*sizeof(double), & multiple_row_type);
	MPI_Type_commit (& multiple_row_type);

	MPI_Datatype interleaved_row;
	MPI_Type_vector (myrows, 1, nprocs, MPI_DOUBLE, & interleaved_row );
	MPI_Type_commit (& interleaved_row);

	MPI_Datatype interleaved_row_type;
	MPI_Type_create_resized (interleaved_row, 0, 1*sizeof(double), & interleaved_row_type);
	MPI_Type_commit (& interleaved_row_type);

	// spread input matrix and vector
	MPI_Scatter (&A[0][0], 1, multiple_row_type, &Alocal[0][0], 1, multiple_row_contiguous, 0, MPI_COMM_WORLD);
    MPI_Scatter (&b[0], 1, interleaved_row_type, &b[rank], 1, interleaved_row_type, 0, MPI_COMM_WORLD);

	for(k=0;k<n;k++)
	{

		// prepare (copy) row for broadcast
		// TODO: avoid copy by aliasing buffers
		// TODO: check aliasing buffers with MPI version (MPICH no aliased buffers!)
		if(map[k] == rank)
		{
			for (i=k; i<n; i++)
			{
				Abase[i]=Alocal[local[k]][i];
			}
		}

		MPI_Bcast (&Abase[k],n-k,MPI_DOUBLE,map[k],MPI_COMM_WORLD);
		MPI_Bcast (&b[k],1,MPI_DOUBLE,map[k],MPI_COMM_WORLD);

		/*
		for(i= k+1; i<n; i++)
		{
			if(map[i] == rank)
			{
				c[i]=Alocal[local[i]][k]/Abase[k];
			}
		}
		*/

		/*
		for(i= k+1; i<n; i++)
		{
			if(map[i] == rank)
			{
				for(j=0;j<n;j++)
				{
					Alocal[local[i]][j]=Alocal[local[i]][j]-( c[i]*Abase[j] );
				}
				b[i]=b[i]-( c[i]*b[k] );
			}
		}
		*/

		// at each global index iteration a process drops a row (the first in the set)
		if (map[k] == rank)
		{
			mystart++;
		}

		// to avoid IFs: each process loops on its own set of rows, with indirect addressing
		for (i=mystart; i<myrows; i++)
		{
			c[global[i]]=Alocal[i][k]/Abase[k];
			for(j=0;j<n;j++)
			{
				Alocal[i][j]=Alocal[i][j]-( c[global[i]]*Abase[j] );
				//printf("k=%d: %d=%d,%d -> %d\n", k, map[global[i]], rank, i, global[i]);
			}
			b[global[i]]=b[global[i]]-( c[global[i]]*b[k] );
		}
	}

    MPI_Gather (&Alocal[0][0], 1, multiple_row_contiguous, &A[rank][0], 1, multiple_row_type, 0, MPI_COMM_WORLD);
    MPI_Gather (&b[rank], 1, interleaved_row_type, &b[0], 1, interleaved_row_type, 0, MPI_COMM_WORLD);

	free(map);
	free(local);
	free(global);
	DeallocateVector(c);
    DeallocateMatrix2D(Alocal,myrows,CONTIGUOUS);
}


void pGaussianElimination_fullmatrix(double** A, double* b, cui n, int rank, int nprocs)
{
    ui i,j,k;
    double* c;
    int* map;
    int numrows;

    c=AllocateVector(n);
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
