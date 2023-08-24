#include <mpi.h>
//#include <stdio.h>
//#include <stdlib.h>
//#include <math.h>
//#include <unistd.h>

#include "../../helpers/matrix.h"
#include "../GaussianElimination/BackSubst.h"

/*
 * Gaussian elimination, parallel version
 *
 * A with rank n
 * bb multiple r.h.s in m columns
 *
 */
void pGaussianElimination(int n, double** A, int m, double** bb, int rank, int nprocs)
{
    int i,j,k,rhs;

    int myrows;
    myrows=n/nprocs;
    int mystart=0;

    double* c;
    c=AllocateVector(n);
		
	double** blocal;
    blocal=AllocateMatrix2D_double(n,m,CONTIGUOUS);

	// row for pivoting
    double* Abase;
    Abase=malloc(n*sizeof(double));

    // local storage for a part of the input matrix (continuous rows, not interleaved)
    double** Alocal;
    Alocal=AllocateMatrix2D_double(myrows, n, CONTIGUOUS);

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
    global=malloc(myrows*sizeof(int));

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
	MPI_Type_vector (myrows, m*1, m*nprocs, MPI_DOUBLE, & interleaved_row );
	MPI_Type_commit (& interleaved_row);

	MPI_Datatype interleaved_row_type;
	MPI_Type_create_resized (interleaved_row, 0, m*sizeof(double), & interleaved_row_type);
	MPI_Type_commit (& interleaved_row_type);

	// spread input matrix and vector
	//
	MPI_Scatter (&A[0][0], 1, multiple_row_type, &Alocal[0][0], 1, multiple_row_contiguous, 0, MPI_COMM_WORLD);

	if (rank == 0)
	{
		for (i=0; i<n; i++)
		{
			for(rhs=0;rhs<m;rhs++)
			{
				blocal[i][rhs]=bb[i][rhs];
			}
		}
	}
	MPI_Bcast (&blocal[0][0],n*m,MPI_DOUBLE,0,MPI_COMM_WORLD);

	//
	// old workaround for openmpi + mpich - superseded by introducing "blocal"
	/*
	if (rank!=0)
	{
		//MPI_Scatter (&b[0], 1, interleaved_row_type, &b[rank], 1, interleaved_row_type, 0, MPI_COMM_WORLD);		// openmpi ok
	}
	else
	{
		//MPI_Scatter (&b[0], 1, interleaved_row_type, MPI_IN_PLACE, 1, interleaved_row_type, 0, MPI_COMM_WORLD);	// mpich ok
	}
	*/

	for(k=0;k<n;k++)
	{
		// prepare (copy) row for broadcast
		if(map[k] == rank)
		{
			for (i=k; i<n; i++)
			{
				Abase[i]=Alocal[local[k]][i];
			}
		}

		MPI_Bcast (&Abase[k],n-k,MPI_DOUBLE,map[k],MPI_COMM_WORLD);
		MPI_Bcast (&blocal[k][0],m,MPI_DOUBLE,map[k],MPI_COMM_WORLD);

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
			}
			for(rhs=0;rhs<m;rhs++)
			{
				blocal[global[i]][rhs]=blocal[global[i]][rhs]-( c[global[i]]*blocal[k][rhs] );
			}
		}
	}

	//collect final solution
	MPI_Gather (&Alocal[0][0], 1, multiple_row_contiguous, &A[0][0], 1, multiple_row_type, 0, MPI_COMM_WORLD);
	MPI_Gather (&blocal[rank][0], 1, interleaved_row_type, &bb[0][0], 1, interleaved_row_type, 0, MPI_COMM_WORLD);
	//REMEMBER: MPI_Gather is picky: buffer pointers MUST HAVE a (not NULL) value

	// old workaround for openmpi
	/*
	for(i=0;i<n;i++)
	{
		if(rank==map[i])
		{
			MPI_Send (&blocal[i][0], m, MPI_DOUBLE, 0, i, MPI_COMM_WORLD);
		}
		if(rank==0)
		{
			MPI_Recv (&b[i][0], m, MPI_DOUBLE, map[i], i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
	}
	*/
	//
	// old workaround for openmpi + mpich - superseded by introducing "blocal"
	/*
	if (rank!=0)
	{
		MPI_Gather (&b[rank], 1, interleaved_row_type, &b[0], 1, interleaved_row_type, 0, MPI_COMM_WORLD);		// openmpi ok
	}
	else
	{
		MPI_Gather (MPI_IN_PLACE, 1, interleaved_row_type, &b[0], 1, interleaved_row_type, 0, MPI_COMM_WORLD);	// mpich ok
	}
	*/

	free(map);
	free(local);
	free(global);
	DeallocateVector(c);
	DeallocateMatrix2D_double(blocal,n,CONTIGUOUS);
	DeallocateVector(Abase);
    DeallocateMatrix2D_double(Alocal,myrows,CONTIGUOUS);
}
