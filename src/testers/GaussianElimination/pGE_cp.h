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

/*
 * checkpointing version (cp)
 *
 * if sprocs=0, then no checkpoints
 */
void pGaussianElimination_cp(int n, double** A, int m, double** b, int rank, int cprocs, int sprocs)
{
	int i,j,k,rhs;

    int nprocs=cprocs+sprocs;

    int myrows;
    myrows=n/cprocs;
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

    int Srows;
    	Srows=myrows*sprocs;

    int Trows=n+Srows;

    int* map;
    map=malloc(Trows*sizeof(int));

    int* local;
    local=malloc(Trows*sizeof(int));

    for(i=0; i<n; i++)
	{
		map[i]= i % cprocs;			// who has the row i
		local[i]=floor(i/cprocs);	// position of the row i(global) in the local matrix
	}
    for(i=n; i<Trows; i++)
    {
    	map[i]=nprocs-1;
    	local[i]=i-n;
    }


    int* global;
    global=malloc(myrows*sizeof(int));

    if (rank>=cprocs)
    {
    	for(i=0; i<myrows; i++)
    		{
    			global[i]=n+i;
    		}
    }
    else
    {
        for(i=0; i<myrows; i++)
    	{
        	global[i]= i * cprocs + rank; // position of the row i(local) in the global matrix
    	}
    }


    // derived data types
	MPI_Datatype multiple_row;
	MPI_Type_vector (myrows, n, n * cprocs , MPI_DOUBLE, & multiple_row );
	MPI_Type_commit (& multiple_row);

	MPI_Datatype multiple_row_contiguous;
	MPI_Type_vector (myrows, n, n , MPI_DOUBLE, & multiple_row_contiguous );
	MPI_Type_commit (& multiple_row_contiguous);

	MPI_Datatype multiple_row_type;
	MPI_Type_create_resized (multiple_row, 0, n*sizeof(double), & multiple_row_type);
	MPI_Type_commit (& multiple_row_type);

	MPI_Datatype interleaved_row;
	MPI_Type_vector (myrows, m*1, m*cprocs, MPI_DOUBLE, & interleaved_row );
	MPI_Type_commit (& interleaved_row);

	MPI_Datatype interleaved_row_type;
	MPI_Type_create_resized (interleaved_row, 0, m*sizeof(double), & interleaved_row_type);
	MPI_Type_commit (& interleaved_row_type);

	int i_am_calc;
		if (rank<cprocs)
		{
			i_am_calc=1;
		}
		else
		{
			i_am_calc=0;
		}

	MPI_Comm comm_calc;
	int rank_calc;
	if (sprocs>0)
	{
		MPI_Comm_split(MPI_COMM_WORLD, i_am_calc, rank, &comm_calc);
		MPI_Comm_rank(comm_calc, &rank_calc);
	}
	else
	{
		comm_calc=MPI_COMM_WORLD;
		rank_calc=rank;
	}

	// spread input matrix and vector
	if (i_am_calc)
	{
		MPI_Scatter (&A[0][0], 1, multiple_row_type, &Alocal[0][0], 1, multiple_row_contiguous, 0, comm_calc);

		if (rank == 0)
		{
			for (i=0; i<n; i++)
			{
				for(rhs=0;rhs<m;rhs++)
				{
					blocal[i][rhs]=b[i][rhs];
				}
			}
		}
		MPI_Bcast (&blocal[0][0],n*m,MPI_DOUBLE,0,MPI_COMM_WORLD);

		//
		// old workaround for openmpi + mpich - superseded by introducing "blocal"
		/*
		if (rank!=0)
		{
			MPI_Scatter (&b[0], 1, interleaved_row_type, &b[rank], 1, interleaved_row_type, 0, comm_calc);		// openmpi ok
		}
		else
		{
			MPI_Scatter (&b[0], 1, interleaved_row_type, MPI_IN_PLACE, 1, interleaved_row_type, 0, comm_calc);	// mpich ok
		}
		*/
	}
	for(k=0;k<n;k++)
	{
		if (i_am_calc)
		{
			// prepare (copy) row for broadcast
			if(map[k] == rank_calc)
			{
				for (i=k; i<n; i++)
				{
					Abase[i]=Alocal[local[k]][i];
				}
			}

			MPI_Bcast (&Abase[k],n-k,MPI_DOUBLE,map[k],comm_calc);
			MPI_Bcast (&blocal[k][0],m,MPI_DOUBLE,map[k],MPI_COMM_WORLD);

			// at each global index iteration a process drops a row (the first in the set)
			if (map[k] == rank_calc)
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
				}			}

			if (sprocs>0)
			{
				// checkpointing
				MPI_Send (&Alocal[0][0], 1, multiple_row_contiguous, nprocs-1, rank, MPI_COMM_WORLD);
			}
		}
		else
		{
			for (j=0; j<cprocs; j++)
			{
				MPI_Recv (&A[j][0], 1, multiple_row, j, j, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			}
		}
	}

	//collect final solution
	if (i_am_calc)
	{
	    MPI_Gather (&Alocal[0][0], 1, multiple_row_contiguous, &A[rank][0], 1, multiple_row_type, 0, comm_calc);
		MPI_Gather (&blocal[rank][0], 1, interleaved_row_type, &b[0][0], 1, interleaved_row_type, 0, MPI_COMM_WORLD);
		//REMEMBER: MPI_Gather is picky: buffer pointers MUST HAVE a (not NULL) value

		//
	    // old workaround for openmpi + mpich - superseded by introducing "blocal"
	    /*
		if (rank!=0)
		{
	    		MPI_Gather (&b[rank], 1, interleaved_row_type, &b[0], 1, interleaved_row_type, 0, comm_calc);		// openmpi ok
		}
		else
		{
	    		MPI_Gather (MPI_IN_PLACE, 1, interleaved_row_type, &b[0], 1, interleaved_row_type, 0, comm_calc);	// mpich ok
		}
		*/
	}


	free(map);
	free(local);
	free(global);
	DeallocateVector(c);
	DeallocateMatrix2D_double(blocal,n,CONTIGUOUS);
	DeallocateVector(Abase);
    DeallocateMatrix2D_double(Alocal,myrows,CONTIGUOUS);
}
