#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include "BackSubst.h"
#include "../../helpers/types.h"
#include "../../helpers/matrix.h"



/*
 * checkpointing version (cp)
 *
 * if sprocs=0, then no checkpoints
 */
void pGaussianElimination_partialmatrix_cp(cui n, double** A, cui m, double** b, int rank, int cprocs, int sprocs)
{
    ui i,j,k,rhs;

    int nprocs=cprocs+sprocs;

    int myrows;
    myrows=n/cprocs;
    int mystart=0;

    double* c;
    c=AllocateVector(n);
    
    double** blocal;
    blocal=AllocateMatrix2D(n,m,CONTIGUOUS);

	// row for pivoting
    double* Abase;
    Abase=malloc(n*sizeof(double));

    // local storage for a part of the input matrix (continuous rows, not interleaved)
    double** Alocal;
    Alocal=AllocateMatrix2D(myrows, n, CONTIGUOUS);

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
	DeallocateMatrix2D(blocal,n,CONTIGUOUS);
	DeallocateVector(Abase);
    DeallocateMatrix2D(Alocal,myrows,CONTIGUOUS);
}

void pGaussianElimination_partialmatrix(cui n, double** A, cui m, double** b, int rank, int nprocs)
{
    ui i,j,k,rhs;

    int myrows;
    myrows=n/nprocs;
    int mystart=0;

    double* c;
    c=AllocateVector(n);
		
	double** blocal;
    blocal=AllocateMatrix2D(n,m,CONTIGUOUS);

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
	MPI_Gather (&blocal[rank][0], 1, interleaved_row_type, &b[0][0], 1, interleaved_row_type, 0, MPI_COMM_WORLD);
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
	DeallocateMatrix2D(blocal,n,CONTIGUOUS);
	DeallocateVector(Abase);
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
