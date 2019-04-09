#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "helpers/matrix.h"
#include "helpers/vector.h"
#include <mpi.h>

#include <unistd.h>

void SLGEWOPV_calc_last_cs(double** A, double* b, double* x, int n, int rank, int cprocs, int sprocs, double** T, double** Tlocal, double* TlastKc, double* TlastKr, double* h, double* hh)
{
	int i,j,l;						// indexes
	int XKcols=2*n;					// num of cols X + K
    int myTcols;					// num of cols per process
    	myTcols=XKcols/cprocs;
    int Scols;
    	Scols=myTcols*sprocs;
    int Tcols=XKcols+Scols;			// num of cols X + K + S
    int myend;						// loop boundaries on local cols
    int mystart;

    int nprocs=cprocs+sprocs;
    /*
     * local storage for a part of the input matrix (continuous columns, not interleaved)
     */

    /*
    double** Tlocal;
    Tlocal=AllocateMatrix2D(n, myTcols, CONTIGUOUS);

    double** T;
			if (rank==0)
			{
				T=AllocateMatrix2D(n,Tcols,CONTIGUOUS);
			}
			else
			{
				//T=AllocateMatrix2D(n,Tcols,CONTIGUOUS);
				T=Tlocal;				// dummy assignment to avoid segfault (i.e. in  next scatter)
			}

    double* TlastKc;					// last col of T (K)
    		TlastKc=AllocateVector(n);
    double* TlastKr;					// last row of T (K part)
    		TlastKr=AllocateVector(n);
    double* h;							// helper vectors
    		h=AllocateVector(n);
    double* hh;
			hh=AllocateVector(n);
*/
	/*
	 * map columns to process
	 */

	int*	local;
    		local=malloc(Tcols*sizeof(int));
    int*	map;
    		map=malloc(Tcols*sizeof(int));
    		if (sprocs>0)						// with checksum cols
    		{
				for (i=0; i<myTcols; i++)
				{
					map[i]= 0;					// root has the first cols of T (S)
					local[i]=i;					// position of the column i(global) in the local matrix
				}
				for (i=myTcols; i<Tcols; i++)
				{
					map[i]= (i-myTcols) % cprocs +sprocs;	// who has the other cols i (from rank 1 onwards)
					local[i]=floor((i-myTcols)/cprocs);	// position of the column i(global) in the local matrix
				}
    		}
    		else								// without checksum cols
    		{
				for (i=0; i<Tcols; i++)
				{
					map[i]= i % cprocs;			// who has the col i
					local[i]=floor(i/cprocs);	// position of the column i(global) in the local matrix
				}
    		}
    int*	global;
    		global=malloc(myTcols*sizeof(int));
    		if (sprocs>0)						// with checksum cols
    		{
    			if (rank==0)
    			{
					for(i=0; i<myTcols; i++)
					{
						global[i]= i ; 			// root has the checksum cols (in the first positions of the column i(local) in the global matrix)
					}
    			}
    			else
    			{
    				for(i=0; i<myTcols; i++)
					{
						global[i]= myTcols + i * cprocs + rank -1; // position of the column i(local) in the global matrix
					}
    			}
    		}
    		else									// without checksum cols
    		{
    			for(i=0; i<myTcols; i++)
    			{
    				global[i]= i * cprocs + rank; 	// position of the column i(local) in the global matrix
    			}
    		}

    /*
     * MPI derived types
     */

	MPI_Datatype single_column;
	MPI_Type_vector (n, 1, Tcols, MPI_DOUBLE, & single_column );
	MPI_Type_commit (& single_column);

	MPI_Datatype local_single_column;
	MPI_Type_vector (n, 1, myTcols, MPI_DOUBLE, & local_single_column );
	MPI_Type_commit (& local_single_column);

	MPI_Datatype interleaved_row;
	MPI_Type_vector (myTcols, 1, cprocs, MPI_DOUBLE, & interleaved_row );
	MPI_Type_commit (& interleaved_row);

	MPI_Datatype interleaved_row_resized;
	MPI_Type_create_resized (interleaved_row, 0, 1*sizeof(double), & interleaved_row_resized);
	MPI_Type_commit (& interleaved_row_resized);

	MPI_Datatype multiple_column;
	MPI_Type_vector (n * myTcols, 1, cprocs , MPI_DOUBLE, & multiple_column );
	MPI_Type_commit (& multiple_column);

	MPI_Datatype local_multiple_column_contiguous;
	MPI_Type_vector (n, myTcols, myTcols , MPI_DOUBLE, & local_multiple_column_contiguous );
	MPI_Type_commit (& local_multiple_column_contiguous);

	MPI_Datatype multiple_column_contiguous;
	MPI_Type_vector (n, myTcols, Tcols , MPI_DOUBLE, & multiple_column_contiguous );
	MPI_Type_commit (& multiple_column_contiguous);

	MPI_Datatype multiple_column_resized;
	MPI_Type_create_resized (multiple_column, 0, 1*sizeof(double), & multiple_column_resized);
	MPI_Type_commit (& multiple_column_resized);


	/*
	 *  init inhibition table
	 */

    MPI_Bcast (&b[0], n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    for (i=0; i<n; i++)
	{
		x[i]=0.0;
	}

	if (rank==0)
	{
		for (j=0; j<n; j++)
		{
			for (i=0; i<n; i++)
			{
				T[i][Scols+j+n]=A[j][i]/A[i][i];
				T[i][Scols+j]=0;
			}
			T[j][Scols+j]=1/A[j][j];

			MPI_Send (&T[0][Scols+j], 1, single_column, map[Scols+j], Scols+j, MPI_COMM_WORLD);
			MPI_Send (&T[0][Scols+j+n], 1, single_column, map[Scols+j+n], Scols+j+n, MPI_COMM_WORLD);
		}
		for (j=0; j<Scols; j++)
		{
			for (i=0; i<n; i++)
			{
				T[i][j]=0;
				for (l=0; l<cprocs; l++)
				{
					T[i][j]=T[i][j]+T[i][Scols+j*cprocs+l];
				}
			}
			MPI_Send (&T[0][j], 1, single_column, map[j], j, MPI_COMM_WORLD);
		}

		if (rank==0)
		{
			printf("\ninited:\n");
			PrintMatrix2D(T,n,Tcols);
		}

		// copy data into local buffer before broadcast
		for (i=0; i<n; i++)
		{
			TlastKc[i]=T[i][Scols+n*2-1];
			TlastKr[i]=T[n-1][Scols+n+i];
		}
	}

	// scatter columns to nodes
	// MPI_Scatter (&T[0][0], 1, multiple_column_resized, &Tlocal[0][0], 1, multiple_column_contiguous, 0, MPI_COMM_WORLD);
	// TODO: define new communicator for scattering, to avoid multiple p2p

	// receive
	for (j=0; j<myTcols; j++)
	{
		MPI_Recv (&Tlocal[0][j], 1, local_single_column, 0, global[j], MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}

	// broadcast of the last col and the last row of T (K part)
	MPI_Bcast (&TlastKc[0], n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast (&TlastKr[0], n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	//TODO: combine previous broadcasts in a single one

	/*
	 *  calc inhibition sequence
	 */


	// all levels but last one (l=0)
	for (l=n-1; l>0; l--)
	{
		MPI_Barrier(MPI_COMM_WORLD);
		// update helpers
		// every process
		for (i=0; i<=l-1; i++)
		{
			h[i]   = 1/(1-TlastKc[i]*TlastKr[i]);
			hh[i]  = TlastKc[i]*h[i];
			b[i]   = b[i]-TlastKr[i]*b[l];
		}

		// only to debug
		for (j=0; j<Tcols; j++)
		{
			if (rank==map[j])
			{
				//printf("%d#%d (%d) sends %d with tag %d\n",l,rank,map[j],j,j);
				MPI_Send (&Tlocal[0][local[j]], 1, local_single_column, 0, j, MPI_COMM_WORLD);
			}
		}
		if (rank==0)
		{
			for (j=0; j<Tcols; j++)
			{
				//printf("%d#%d recvs %d from %d with tag %d\n",l,rank,j,map[j],j);
				MPI_Recv (&T[0][j], 1, single_column, map[j], j, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			}
			printf("*%d::\n",l);
			PrintMatrix2D(T,n,Tcols);
			printf("lastC:\n");
			PrintVector(TlastKc,n);
			printf("lastR:\n");
			PrintVector(TlastKr,n);
			printf("h:\n");
			PrintVector(h,n);
			printf("\nb:\n");
			PrintVector(b,n);
		}
		// end debug


		if (rank>0)
		{
			// update solutions
			// l .. n-1

			/*
			mystart=local[Scols+l];
			if (rank<map[Scols+l])
			{
				mystart++;
			}
			for (i=mystart; i<=local[Scols+n-1]; i++)
			{
				printf("\n%d@%d: %d",rank,l,global[i]-Scols);
				x[global[i]-Scols]=x[global[i]-Scols]+Tlocal[l][i]*b[l];
			}
			*/
			for (j=Scols; j<n+Scols; j++)
			{
				if (rank==map[j])
				{
					x[j-Scols]=x[j-Scols]+Tlocal[l][local[j]]*b[l];
				}
			}

			// update T
			// to avoid IFs: each process loops on its own set of cols, with indirect addressing

			// 0 .. l-1
			// processes with diagonal elements not null
			mystart=local[Scols+0];
			myend=local[Scols+l-1];
			if (rank>map[Scols+l-1])
			{
				myend--;
			}
			for (i=mystart; i<=myend; i++)
			{
				Tlocal[global[i]-Scols][i]=Tlocal[global[i]-Scols][i]*h[global[i]-Scols];
			}


			// l
			// process with full not null column (column l)
			if (rank==map[Scols+l])
			{
				printf("\nsinglel %d (local %d) in %d\n",l,local[Scols+l],rank);
				//PrintMatrix2D(Tlocal,n,myTcols);
				for (i=0; i<=l-1; i++)
				{
					Tlocal[i][local[Scols+l]]= -Tlocal[l][local[Scols+l]]*hh[i];
					printf("\n> %d:%g,%g\n",i,Tlocal[i][local[Scols+l]],Tlocal[l][local[Scols+l]]);
				}
			}


			// l+1 .. n+l-1
			// all other cases
			mystart=local[Scols+l+1];
			if (rank<map[Scols+l+1])
			{
				mystart++;
			}
			myend=local[Scols+n+l-1];
			if (rank>map[Scols+n+l-1])
			{
				myend--;
			}

			//printf("%d@%d: %d(%d)-%d(%d)\n",l,rank,mystart,global[mystart],myend,global[myend]);
			for (j=mystart; j<=myend; j++)
			{
				for (i=0; i<=l-1; i++)
				{
					//Tlocal[i][j]=(Tlocal[i][j]-Tlocal[l][j]*TlastKc[i])*h[i];
					Tlocal[i][j]=Tlocal[i][j]*h[i]-Tlocal[l][j]*hh[i];
				}
			}

		}
		else // root (containg S)
		{
			for (j=0; j<myTcols; j++)
			{
				for (i=0; i<=l-1; i++)
				{
					Tlocal[i][j]=Tlocal[i][j]*h[i]-Tlocal[l][j]*hh[i];
				}
			}

		}

		// collect chunks of last row of K to "future" last node
		//MPI_Gather (&Tlocal[l-1][local[n]], myTcols/2, MPI_DOUBLE, &TlastKr[0], 1, interleaved_row_resized, map[l-1], MPI_COMM_WORLD);

		if (rank>0)
		{
			//printf("%d sending to %d tag %d\n",rank,map[Scols+l-1],Scols+l-1);
			MPI_Send (&Tlocal[l-1][local[Scols+n]], myTcols/2, MPI_DOUBLE, map[Scols+n+l-1], rank, MPI_COMM_WORLD);
		}
		if (rank==map[Scols+n+l-1])
		{
			for (j=1; j<nprocs; j++)
			{
				//if (j!=map[Scols+n+l-1])
				{
					//printf("%d wants to receive from %d tag %d\n",rank,j,Scols+l-1);
					MPI_Recv (&TlastKr[j-1],  1, interleaved_row_resized, j, j, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				}
			}
		}

		//future last node broadcasts last row and col of K
		if (rank==map[Scols+n+l-1])
		{
			printf("**%d %d %d\n",rank,l, local[Scols+n+l-1]);
			// copy data into local buffer before broadcast
			for (i=0; i<n; i++)
			{
				TlastKc[i]=Tlocal[i][local[Scols+n+l-1]];
			}
		}
		MPI_Bcast (&TlastKc[0], n, MPI_DOUBLE, map[Scols+n+l-1], MPI_COMM_WORLD);
		MPI_Bcast (&TlastKr[0], n, MPI_DOUBLE, map[Scols+n+l-1], MPI_COMM_WORLD);
		//TODO: substitute Gather with an All-to-All
		//TODO: or combine broadcasts in a single one

	}

	// only to debug
	for (j=0; j<Tcols; j++)
	{
		if (rank==map[j])
		{
			//printf("%d#%d (%d) sends %d with tag %d\n",l,rank,map[j],j,j);
			MPI_Send (&Tlocal[0][local[j]], 1, local_single_column, 0, j, MPI_COMM_WORLD);
		}
	}
	if (rank==0)
	{
		for (j=0; j<Tcols; j++)
		{
			//printf("%d#%d recvs %d from %d with tag %d\n",l,rank,j,map[j],j);
			MPI_Recv (&T[0][j], 1, single_column, map[j], j, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
		printf("*%d::\n",l);
		PrintMatrix2D(T,n,Tcols);
		printf("lastC:\n");
		PrintVector(TlastKc,n);
		printf("lastR:\n");
		PrintVector(TlastKr,n);
		printf("h:\n");
		PrintVector(h,n);
		printf("\nb:\n");
		PrintVector(b,n);
	}
	// end debug

	// last level (l=0)
	/*
	if (rank>0)
	{
		for (i=0; i<myTcols/2; i++)
		{
			printf("%d:%d\n",rank,global[i]-Scols);
			x[global[i]-Scols]=x[global[i]-Scols]+Tlocal[0][i]*b[0];
		}
		if (rank==2) 		PrintVector(x,n);
	}
	*/

	//MPI_Barrier(MPI_COMM_WORLD);
	//printf("\nfine x\n");
	// collect solution
	//MPI_Gather (&x[rank], 1, interleaved_row_resized, &x[0], 1, interleaved_row_resized, 0, MPI_COMM_WORLD);


	for (j=Scols; j<n+Scols; j++)
	{
		if (rank==map[j])
		{
			x[j-Scols]=x[j-Scols]+Tlocal[0][local[j]]*b[0];
			MPI_Send (&x[j-Scols], 1, MPI_DOUBLE, 0, j, MPI_COMM_WORLD);
		}
		if (rank==0)
		{
			MPI_Recv (&x[j-Scols],  1, MPI_DOUBLE, map[j], j, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
	}

	// cleanup
	free(local);
	free(global);
	free(map);

	/*
	DeallocateVector(TlastKc);
	DeallocateVector(TlastKr);
	DeallocateVector(h);
	DeallocateVector(hh);
	DeallocateMatrix2D(Tlocal,n,CONTIGUOUS);
    if (rank==0)
    {
    	DeallocateMatrix2D(T,n,CONTIGUOUS);
    }
	*/
}


void SLGEWOPV_calc_last(double** A, double* b, double* x, int n, int rank, int cprocs)
{
	int i,j,l;						// indexes
    int Tcols=2*n;					// total num of cols (X + K)
    int myTcols;					// num of cols per process
    	myTcols=Tcols/cprocs;
    int myend;						// loop boundaries on local cols =myTcols/2;
    int mystart;

    /*
     * local storage for a part of the input matrix (continuous columns, not interleaved)
     */
    double** Tlocal;
    		Tlocal=AllocateMatrix2D(n, myTcols, CONTIGUOUS);

    double** T;
			if (rank==0)
			{
				T=AllocateMatrix2D(n,Tcols,CONTIGUOUS);
			}
			else
			{
				T=Tlocal;				// dummy assignment to avoid segfault (i.e. in  next scatter)
			}

    double* TlastKc;					// last col of T (K)
    		TlastKc=AllocateVector(n);
    double* TlastKr;					// last row of T (K part)
    		TlastKr=AllocateVector(n);
    double* h;							// helper vectors
    		h=AllocateVector(n);
    double* hh;
			hh=AllocateVector(n);

	/*
	 * map columns to process
	 */

	int*	local;
    		local=malloc(Tcols*sizeof(int));
    int*	map;
    		map=malloc(Tcols*sizeof(int));
			for (i=0; i<Tcols; i++)
			{
				map[i]= i % cprocs;			// who has the col i
				local[i]=floor(i/cprocs);	// position of the column i(global) in the local matrix
			}
    int*	global;
    		global=malloc(myTcols*sizeof(int));
			for(i=0; i<myTcols; i++)
			{
				global[i]= i * cprocs + rank; // position of the column i(local) in the global matrix
			}

    /*
     * MPI derived types
     */

	MPI_Datatype single_column;
	MPI_Type_vector (n, 1, Tcols, MPI_DOUBLE, & single_column );
	MPI_Type_commit (& single_column);

	MPI_Datatype interleaved_row;
	MPI_Type_vector (n/cprocs, 1, cprocs, MPI_DOUBLE, & interleaved_row );
	MPI_Type_commit (& interleaved_row);

	MPI_Datatype interleaved_row_resized;
	MPI_Type_create_resized (interleaved_row, 0, 1*sizeof(double), & interleaved_row_resized);
	MPI_Type_commit (& interleaved_row_resized);

	MPI_Datatype multiple_column;
	MPI_Type_vector (n * myTcols, 1, cprocs , MPI_DOUBLE, & multiple_column );
	MPI_Type_commit (& multiple_column);

	MPI_Datatype multiple_column_contiguous;
	MPI_Type_vector (n, myTcols, myTcols , MPI_DOUBLE, & multiple_column_contiguous );
	MPI_Type_commit (& multiple_column_contiguous);

	MPI_Datatype multiple_column_resized;
	MPI_Type_create_resized (multiple_column, 0, 1*sizeof(double), & multiple_column_resized);
	MPI_Type_commit (& multiple_column_resized);


	/*
	 *  init inhibition table
	 */

    MPI_Bcast (&b[0], n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    for (i=0; i<n; i++)
	{
		x[i]=0.0;
	}

	if (rank==0)
	{
		for (j=0; j<n; j++)
		{
			for (i=0; i<n; i++)
			{
				T[i][j+n]=A[j][i]/A[i][i];
				T[i][j]=0;
			}
			T[j][j]=1/A[j][j];
		}
		// copy data into local buffer before broadcast
		for (i=0; i<n; i++)
		{
			TlastKc[i]=T[i][n*2-1];
			TlastKr[i]=T[n-1][n+i];
		}
	}

	// scatter columns to nodes
	MPI_Scatter (&T[0][0], 1, multiple_column_resized, &Tlocal[0][0], 1, multiple_column_contiguous, 0, MPI_COMM_WORLD);

	// broadcast of the last col and the last row of T (K part)
	MPI_Bcast (&TlastKc[0], n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast (&TlastKr[0], n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	//TODO: combine previous broadcasts in a single one

	/*
	 *  calc inhibition sequence
	 */

	// all levels but last one (l=0)
	for (l=n-1; l>0; l--)
	{
		// update solutions
		// l .. n-1
		mystart=local[l];
		if (rank<map[l])
		{
			mystart++;
		}
		for (i=mystart; i<=local[n-1]; i++)
		{
			x[global[i]]=x[global[i]]+Tlocal[l][i]*b[l];
		}

		// update helpers
		// every process
		for (i=0; i<=l-1; i++)
		{
			h[i]   = 1/(1-TlastKc[i]*TlastKr[i]);
			hh[i]  = TlastKc[i]*h[i];
			b[i]   = b[i]-TlastKr[i]*b[l];
		}

		// update T
		// to avoid IFs: each process loops on its own set of cols, with indirect addressing

		// 0 .. l-1
		// processes with diagonal elements not null
		mystart=0;
		myend=local[l-1];
		if (rank>map[l-1])
		{
			myend--;
		}
		for (i=mystart; i<=myend; i++)
		{
			Tlocal[global[i]][i]=Tlocal[global[i]][i]*h[global[i]];
		}

		// l
		// process with full not null column (column l)
		if (rank==map[l])
		{
			for (i=0; i<=l-1; i++)
			{
				Tlocal[i][local[l]]= -Tlocal[l][local[l]]*hh[i];
			}
		}

		// l+1 .. n+l-1
		// all other cases
		mystart=local[l+1];
		if (rank<map[l+1])
		{
			mystart++;
		}
		myend=local[n+l-1];
		if (rank>map[n+l-1])
		{
			myend--;
		}
		for (j=mystart; j<=myend; j++)
		{
			for (i=0; i<=l-1; i++)
			{
				Tlocal[i][j]=Tlocal[i][j]*h[i]-Tlocal[l][j]*hh[i];
			}
		}

		// collect chunks of last row of K to "future" last node
		MPI_Gather (&Tlocal[l-1][local[n]], myTcols/2, MPI_DOUBLE, &TlastKr[0], 1, interleaved_row_resized, map[l-1], MPI_COMM_WORLD);

		//future last node broadcasts last row and col of K
		if (rank==map[l-1])
		{
			// copy data into local buffer before broadcast
			for (i=0; i<n; i++)
			{
				TlastKc[i]=Tlocal[i][local[n+l-1]];
			}
		}
		MPI_Bcast (&TlastKc[0], n, MPI_DOUBLE, map[l-1], MPI_COMM_WORLD);
		MPI_Bcast (&TlastKr[0], n, MPI_DOUBLE, map[l-1], MPI_COMM_WORLD);
		//TODO: substitute Gather with an All-to-All
		//TODO: or combine broadcasts in a single one
	}

	// last level (l=0)
	for (i=0; i<myTcols/2; i++)
	{
		x[global[i]]=x[global[i]]+Tlocal[0][i]*b[0];
	}

	// collect solution
	MPI_Gather (&x[rank], 1, interleaved_row_resized, &x[0], 1, interleaved_row_resized, 0, MPI_COMM_WORLD);

	// cleanup
	free(local);
	free(global);
	free(map);
	DeallocateVector(TlastKc);
	DeallocateVector(TlastKr);
	DeallocateVector(h);
	DeallocateVector(hh);
	DeallocateMatrix2D(Tlocal,n,CONTIGUOUS);
    if (rank==0)
    {
    	DeallocateMatrix2D(T,n,CONTIGUOUS);
    }
}


/*
 * IF removed from init
 * IF removed from calc
 * init with p2p
 * calc with collectives
 */
void SLGEWOPV_calc_sendopt(double** A, double* b, double* s, int n, double** K, double* H, int rank, int cprocs)
{
	int i,j,l;				// indexes
    int mycols;				// num of cols per process
    	mycols=n/cprocs;
    int myend;				// loop boundaries on local cols =mycols;
    int mystart;

	double** X=A;			// aliases
    double*  F=b;

    /*
     * local storage for a part of the input matrix (continuous columns, not interleaved)
     */

    double**	localK;
    			localK=AllocateMatrix2D(n,mycols,CONTIGUOUS);
    double**	localX;
    			localX=AllocateMatrix2D(n,mycols,CONTIGUOUS);

    double* 	lastKc;						// last col of K
    			lastKc=AllocateVector(n);
    double*		lastKr;						// last row of K
				lastKr=AllocateVector(n);
	double		tmpDiag;					// to hold a diagonal element during init

	/*
	 *  map columns to process
	 */

    int*	local;
    		local=malloc(n*sizeof(int));
	int*	map;
    		map=malloc(n*sizeof(int));
			for (i=0; i<n; i++)
			{
				map[i]= i % cprocs;			// who has the col i
				local[i]=floor(i/cprocs);	// position of the column i(global) in the local matrix

				s[i]=0.0;					// and init solution vector
			}
	int* 	global;
    		global=malloc(n*sizeof(int));
			for(i=0; i<mycols; i++)
			{
				global[i]= i * cprocs + rank; // position of the column i(local) in the global matrix
			}

	/*
	 * MPI derived types
	 */

	MPI_Datatype single_column;
	MPI_Type_vector (n, 1, n, MPI_DOUBLE, & single_column );
	MPI_Type_commit (& single_column);

	MPI_Datatype local_single_column;
	MPI_Type_vector (n, 1, mycols, MPI_DOUBLE, & local_single_column );
	MPI_Type_commit (& local_single_column);

	MPI_Datatype multiple_column;
	MPI_Type_vector (n * mycols, 1, cprocs , MPI_DOUBLE, & multiple_column );
	MPI_Type_commit (& multiple_column);

	MPI_Datatype multiple_column_resized;
	MPI_Type_create_resized (multiple_column, 0, 1*sizeof(double), & multiple_column_resized);
	MPI_Type_commit (& multiple_column_resized);

	MPI_Datatype interleaved_row;
	MPI_Type_vector (n/cprocs, 1, cprocs, MPI_DOUBLE, & interleaved_row );
	MPI_Type_commit (& interleaved_row);

	MPI_Datatype interleaved_row_resized;
	MPI_Type_create_resized (interleaved_row, 0, 1*sizeof(double), & interleaved_row_resized);
	MPI_Type_commit (& interleaved_row_resized);


	/*
	 * init inhibition table
	 */

    MPI_Bcast (&b[0], n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	//send and recv inside loops, trying to overlap with calculations
	if (rank==0)
	{
		for (j=0; j<n; j++)
		{
			// init col of K
			for (i=0; i<n; i++)
			{
				K[i][j]=A[j][i]/A[i][i];
				// when i==j -> K[i][j]=1
			}

			MPI_Send (&K[0][j], 1, single_column, map[j], j+n, MPI_COMM_WORLD);
		}
		for (j=0; j<n; j++)
		{
			// init col of X
			tmpDiag=1/A[j][j];
			for (i=0; i<n; i++)
			{
				X[i][j]=0.0;
			}
			X[j][j]=tmpDiag;

			MPI_Send (&X[0][j], 1, single_column, map[j], j, MPI_COMM_WORLD);
		}
		// copy data into local buffer before broadcast
		for (i=0; i<n; i++)
		{
			lastKc[i]=K[i][n-1];
			lastKr[i]=K[n-1][i];
		}
	}

	// receive
	for (j=0; j<mycols; j++)
	{
		MPI_Recv (&localX[0][j], 1, local_single_column, 0, global[j], MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Recv (&localK[0][j], 1, local_single_column, 0, global[j]+n, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}

	// broadcast last column and last row of K
	MPI_Bcast (&lastKc[0], n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast (&lastKr[0], n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	//TODO: combine previous broadcasts in a single one

	/*
	 *  calc inhibition sequence
	 */

	for (l=n-1; l>0; l--)
	{

		for (i=0; i<l; i++)
		{
			H[i]=1/(1-lastKr[i]*lastKc[i]);
			F[i]=F[i]-F[l]*lastKr[i];
		}
		for (j=rank; j<l; j+=cprocs)
		{
			localX[j][local[j]]=H[j]*(localX[j][local[j]]-lastKc[j]*localX[l][local[j]]);
		}
		for (i=0; i<l; i++)
		{
			for (j=n-(cprocs-rank-1)-1; j>=l; j-=cprocs)
			{
				localX[i][local[j]]=H[i]*(localX[i][local[j]]-lastKc[i]*localX[l][local[j]]);
			}
			for (j=rank; j<l; j+=cprocs)
			{
				localK[i][local[j]]=H[i]*(localK[i][local[j]]-lastKc[i]*localK[l][local[j]]);
			}
		}

		// collect chunks of last row of K to "future" last node
		MPI_Gather (&localK[l-1][0], mycols, MPI_DOUBLE, &lastKr[0], 1, interleaved_row_resized, map[l-1], MPI_COMM_WORLD);

		//future last node broadcasts last row and col of K
		if (rank==map[l-1])
		{
			// copy data into local buffer before broadcast
			for (i=0; i<n; i++)
			{
				lastKc[i]=localK[i][local[l-1]];
			}
		}
		MPI_Bcast (&lastKc[0], n, MPI_DOUBLE, map[l-1], MPI_COMM_WORLD);
		MPI_Bcast (&lastKr[0], n, MPI_DOUBLE, map[l-1], MPI_COMM_WORLD);
		//TODO: combine previous broadcasts in a single one
	}

	// calc solution on nodes
	for (j=rank; j<n; j+=cprocs)
	{
		for (l=0; l<n; l++)
		{
			s[j]=F[l]*localX[l][local[j]]+s[j];
		}
	}

	// collect solution
	MPI_Gather (&s[rank], 1, interleaved_row, &s[0], 1, interleaved_row_resized, 0, MPI_COMM_WORLD);

	// cleanup
	free(local);
	free(global);
	free(map);
	DeallocateMatrix2D(localK,n,CONTIGUOUS);
	DeallocateMatrix2D(localX,n,CONTIGUOUS);
	DeallocateVector(lastKr);
	DeallocateVector(lastKc);
}


/*
 * IF _not_ removed from init
 * IF removed from calc
 * every node broadcasts chunks of last row
 * last node broadcast last col
 * master calcs and broadcasts auxiliary causes
 */
void SLGEWOPV_calc_unswitch(double** A, double* b, double* s, int n, double** K, double* H, int rank, int cprocs)
{
    double** X=A;
    double*  F=b;
    double tmpAdiag;

	// indexes
    int i,j,l;

    // map columns to process
    int* map;
    map=malloc(n*sizeof(int));
	for (i=0; i<n; i++)
	{
		map[i]= i % cprocs;
		s[i]=0.0;			// and init solution vector
	}

    // num of cols per process
    int numcols;
    numcols=n/cprocs;

    // MPI derived types
	MPI_Datatype single_column;
	MPI_Type_vector (n, 1, n, MPI_DOUBLE, & single_column );
	MPI_Type_commit (& single_column);

	MPI_Datatype interleaved_row;
	MPI_Type_vector (numcols, 1, cprocs, MPI_DOUBLE, & interleaved_row );
	MPI_Type_commit (& interleaved_row);

	MPI_Datatype interleaved_row_type;
	MPI_Type_create_resized (interleaved_row, 0, 1*sizeof(double), & interleaved_row_type);
	MPI_Type_commit (& interleaved_row_type);

	MPI_Datatype multiple_column;
	MPI_Type_vector (n * numcols, 1, cprocs, MPI_DOUBLE, & multiple_column );
	MPI_Type_commit (& multiple_column);

	MPI_Datatype multiple_column_type;
	MPI_Type_create_resized (multiple_column, 0, 1*sizeof(double), &multiple_column_type);
	MPI_Type_commit (& multiple_column_type);


	/*
	 *  init inhibition table
	 */

	MPI_Bcast (&b[0], n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	if (rank==0)
	{
		for (i=0; i<n; i++)
		{
			tmpAdiag=1/A[i][i];
			for (j=0; j<n; j++)
			{
				if (i==j)
				{
					X[i][j]=tmpAdiag;
					K[i][j]=1;
				}
				else
				{
					K[i][j]=A[j][i]*tmpAdiag;

					// ATTENTION : transposed
					X[j][i]=0.0;
				}
			}
		}
	}

	// spread inhibition table
	MPI_Scatter (&X[0][0], 1, multiple_column_type, &X[0][rank], 1, multiple_column_type, 0, MPI_COMM_WORLD);
	MPI_Scatter (&K[0][0], 1, multiple_column_type, &K[0][rank], 1, multiple_column_type, 0, MPI_COMM_WORLD);

	// broadcast last column of K
	MPI_Bcast (&K[0][n-1], 1, single_column, 0, MPI_COMM_WORLD);

	// broadcast of the last row of K
	MPI_Bcast (&K[n-1][0], n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	/*
	 *  calc inhibition sequence
	 */
	for (l=n-1; l>0; l--)
	{
		for (i=0; i<l; i++)
		{
			H[i]=1/(1-K[i][l]*K[l][i]);
		}
		for (j=rank; j<l; j+=cprocs)
		{
			X[j][j]=H[j]*(X[j][j]-K[j][l]*X[l][j]);
		}
		for (i=0; i<l; i++)
		{

			for (j=n-(cprocs-rank-1)-1; j>=l; j-=cprocs)
			{
				X[i][j]=H[i]*(X[i][j]-K[i][l]*X[l][j]);
			}
			for (j=rank; j<l; j+=cprocs)
			{
				K[i][j]=H[i]*(K[i][j]-K[i][l]*K[l][j]);
			}
		}

		// every node broadcasts its chunks of the last row of K
		MPI_Allgather (&K[l-1][rank], 1, interleaved_row_type, &K[l-1][0], 1, interleaved_row_type, MPI_COMM_WORLD);

		// broadcast last column of K
		MPI_Bcast (&K[0][l-1], 1, single_column, map[l-1], MPI_COMM_WORLD);
	}

	// calc auxiliary causes on root and then broadcast
    if (rank==0)
    {
		for (i=n-2; i>=0; i--)
		{
			for (l=i+1; l<n; l++)
			{
				F[i]=F[i]-F[l]*K[l][i];
			}
		}
    }
	MPI_Bcast (F, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	// final solution on nodes
	for (j=rank; j<n; j+=cprocs)
	{
		for (l=0; l<n; l++)
		{
			s[j]=F[l]*X[l][j]+s[j];
		}
	}

	// collect solution
	MPI_Gather (&s[rank], 1, interleaved_row_type, &s[0], 1, interleaved_row_type, 0, MPI_COMM_WORLD);
}


/*
 * IF _not_ removed from init
 * IF _not_ removed from calc
 * send chunks of last row to last node
 * broadcast last row and last col
 * every node calcs auxiliary causes
 * calc and broadcast solution
 */
void SLGEWOPV_calc_base(double** A, double* b, double* s, int n, double** K, double* H, int rank, int cprocs)
{
	double** X=A;
	double*  F=b;
	double tmpAdiag;

	// indexes
    int i,j,l;

    // map columns to process
    int* map;
    map=malloc(n*sizeof(int));
	for (i=0; i<n; i++)
	{
		map[i]= i % cprocs;
		s[i]=0.0;			// and init solution vector
	}

    // MPI derived types
	MPI_Datatype single_column;
	MPI_Type_vector (n, 1, n, MPI_DOUBLE, & single_column );
	MPI_Type_commit (& single_column);

	MPI_Datatype interleaved_row;
	MPI_Type_vector (n/cprocs, 1, cprocs, MPI_DOUBLE, & interleaved_row );
	MPI_Type_commit (& interleaved_row);

	MPI_Datatype interleaved_row_type;
	MPI_Type_create_resized (interleaved_row, 0, 1*sizeof(double), &interleaved_row_type);
	MPI_Type_commit (& interleaved_row_type);


	/*
	 * init inhibition table
	 */
    MPI_Bcast (&b[0], n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	if (rank==0)
	{
		for (i=0; i<n; i++)
		{
			tmpAdiag=1/A[i][i];
			for (j=0; j<n; j++)
			{
				if (i==j)
				{
					X[i][j]=tmpAdiag;
					K[i][j]=1;
				}
				else
				{
					K[i][j]=A[j][i]*tmpAdiag;

					// ATTENTION : transposed
					X[j][i]=0.0;
				}
			}
		}
	}

	// spread inhibition table
	if (rank==0)
	{
		for (j=0; j<n-1; j++) // to last but one
		{
			// send columns to single nodes
			if (map[j]!=0)
			{
				MPI_Send (&X[0][j], 1, single_column, map[j], j, MPI_COMM_WORLD);
				MPI_Send (&K[0][j], 1, single_column, map[j], j+n, MPI_COMM_WORLD);
			}
		}
		// last one of K is broadcasted to all nodes later
		MPI_Send (&X[0][n-1], 1, single_column, map[n-1], n-1, MPI_COMM_WORLD);
	}
	else
	{
		for (j=0; j<n-1; j++)
		{
			// receive
			if (rank==map[j])
			{
				MPI_Recv (&X[0][j], 1, single_column, 0, j, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				MPI_Recv (&K[0][j], 1, single_column, 0, j+n, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			}
		}
		if (rank==map[n-1])
			{
			MPI_Recv (&X[0][n-1], 1, single_column, 0, n-1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			}
	}

	// broadcast last column of K
	MPI_Bcast (&K[0][n-1], 1, single_column, 0, MPI_COMM_WORLD);

	// broadcast of the last row of K
	MPI_Bcast (&K[n-1][0], n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	/*
	 *  calc inhibition sequence
	 */
	for (l=n-1; l>0; l--)
	{
		for (i=0; i<l; i++)
		{
			H[i]=1/(1-K[i][l]*K[l][i]);
			if (map[i]==rank)
			{
				X[i][i]=H[i]*(X[i][i]);
			}
			for (j=l; j<n; j++)
			{
				if (map[j]==rank)
				{
					X[i][j]=H[i]*(X[i][j]-K[i][l]*X[l][j]);
				}
			}
			for (j=0; j<l; j++)
			{
				if (map[j]==rank)
				{
					K[i][j]=H[i]*(K[i][j]-K[i][l]*K[l][j]);
				}
			}
		}

		// every node broadcasts its chunks of the last row of K
		for (j=0; j<cprocs; j++)
		{
			MPI_Bcast (&K[l-1][j], 1, interleaved_row, j, MPI_COMM_WORLD);
		}

		// broadcast last column of K
		MPI_Bcast (&K[0][l-1], 1, single_column, map[l-1], MPI_COMM_WORLD);

		// calc auxiliary causes on every node
		for (i=l-1; i>=0; i--)
		{
				F[i]=F[i]-F[l]*K[l][i];
		}
	}

	// final solution on nodes
	for (j=0; j<n; j++)
	{
		if (map[j]==rank)
		{
			for (l=0; l<n; l++)
			{
				s[j]=F[l]*X[l][j]+s[j];
			}
		}
	}

	// collect solution
	MPI_Gather (&s[rank], 1, interleaved_row_type, &s[0], 1, interleaved_row_type, 0, MPI_COMM_WORLD);
}
