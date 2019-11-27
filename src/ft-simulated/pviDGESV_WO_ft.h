#include <mpi.h>
#include "helpers/matrix.h"
#include "helpers/vector.h"
#include "DGEZR.h"
#include <signal.h>
#include "../ft-real/pDGEIT_Wx.h"

/*
 *	solve (SV) system with general (GE) matrix A of doubles (D)
 *	of order n and with m r.h.s in matrix bb[n,m] and solutions in xx[n,m]
 *	with:
 *	wide overwrite (WO) memory model
 *	parallelized in interleaved columns (pvi) over cprocs calculating processors
 *	parallelized initialization
 *	optimized loops
 *	ifs removed
 *
 */

extern int failing_rank;
extern int failing_level;
extern int faulty;
extern int faulty_location;
extern MPI_Comm current_comm;
extern void common_error_handler(int s);

void pviDGESV_WO_ft(int n, double** A, int m, double** bb, double** xx, MPI_Comm comm, int sprocs)
{
    int rank;
	int nprocs, cprocs;
	MPI_Comm_size(comm, &nprocs);
	MPI_Comm_rank(comm, &rank);
	cprocs=nprocs-sprocs;

	int i,j,l;						// indexes
    int Tcols=2*n;					// total num of cols (X + K)
    int myTcols;					// num of T cols per process
    	myTcols=Tcols/cprocs;
    int myAchunks;					// num of A rows/cols per process
    	myAchunks=n/cprocs;
    int myxxrows=myAchunks;
    int myend;						// loop boundaries on local cols =myTcols/2;
    int mystart;
    int rhs;

    int avoidif;					// for boolean --> int conversion

    int myAcols=n/cprocs;
    /*
     * local storage for a part of the input matrix (continuous columns, not interleaved)
     */
    double** Tlocal;
    		 Tlocal=AllocateMatrix2D(n, myTcols, CONTIGUOUS);

	double** TlastK;
			 TlastK=AllocateMatrix2D(2,n, CONTIGUOUS);	// last col [0] and row [1] of T (K part)
	double*  TlastKc=&TlastK[0][0];						// alias for last col
	double*  TlastKr=&TlastK[1][0];						// alias for last row

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
			for (i=0; i<myTcols; i++)
			{
				global[i]= i * cprocs + rank; // position of the column i(local) in the global matrix
			}

    /*
     * MPI derived types
     */

	MPI_Datatype TlastKr_chunks;
	MPI_Type_vector (n/cprocs, 1, cprocs, MPI_DOUBLE, & TlastKr_chunks );
	MPI_Type_commit (& TlastKr_chunks);

	MPI_Datatype TlastKr_chunks_resized;
	MPI_Type_create_resized (TlastKr_chunks, 0, 1*sizeof(double), & TlastKr_chunks_resized);
	MPI_Type_commit (& TlastKr_chunks_resized);

	// rows of xx to be extracted
	MPI_Datatype xx_rows_interleaved;
	MPI_Type_vector (myxxrows, m, m*cprocs, MPI_DOUBLE, & xx_rows_interleaved );
	MPI_Type_commit (& xx_rows_interleaved);

	// rows of xx to be extracted, properly resized for gathering
	MPI_Datatype xx_rows_interleaved_resized;
	MPI_Type_create_resized (xx_rows_interleaved, 0, m*sizeof(double), & xx_rows_interleaved_resized);
	MPI_Type_commit (& xx_rows_interleaved_resized);


	// FT



    /*
	 *  init inhibition table
	 */
	DGEZR(xx, n, m);													// init (zero) solution vectors
	pDGEIT_W(n, A, Tlocal, TlastK, comm, map, global, local);	// init inhibition table

	faulty_location=1;
	// send all r.h.s to all procs
    MPI_Bcast (&bb[0][0], n*m, MPI_DOUBLE, 0, comm);
    faulty_location=10;
	/*
	 *  calc inhibition sequence
	 */
if (!faulty)
{
	// all levels but last one (l=0)
	for (l=n-1; l>0; l--)
	{

		// ALL procs
		// update solutions
		// l .. n-1

		/*
		mystart=local[l];
		if (rank<map[l])
		{
			mystart++;
		}
		*/
		avoidif=rank<map[l];
		mystart = local[l] + avoidif;

		for (i=mystart; i<=local[n-1]; i++)
		{
			for (rhs=0;rhs<m;rhs++)
			{
				xx[global[i]][rhs]=xx[global[i]][rhs]+Tlocal[l][i]*bb[l][rhs];
			}
		}

		// ALL procs
		// update helpers
		for (i=0; i<=l-1; i++)
		{
			h[i]   = 1/(1-TlastKc[i]*TlastKr[i]);
			hh[i]  = TlastKc[i]*h[i];
			for (rhs=0;rhs<m;rhs++)
			{
				bb[i][rhs] = bb[i][rhs]-TlastKr[i]*bb[l][rhs];
			}
		}

		// update T
		// to avoid IFs: each process loops on its own set of cols, with indirect addressing

		// 0 .. l-1
		// ALL procs
		// processes with diagonal elements not null
		mystart=0;
		/*
		myend=local[l-1];
		if (rank>map[l-1])
		{
			myend--;
		}
		*/
		avoidif = rank>map[l-1];
		myend = local[l-1] - avoidif;

		for (i=mystart; i<=myend; i++)
		{
			Tlocal[global[i]][i]=Tlocal[global[i]][i]*h[global[i]];
		}

		// l .. n+l-1
		// ALL procs
		/*
		mystart=local[l];
		if (rank<map[l])
		{
			mystart++;
		}
		*/
		avoidif=rank<map[l];
		mystart=local[l]+avoidif;
		/*
		myend=local[n+l-1];
		if (rank>map[n+l-1])
		{
			myend--;
		}
		*/
		avoidif=rank>map[n+l-1];
		myend=local[n+l-1]-avoidif;

		for (i=0; i<=l-1; i++)
		{
			for (j=mystart; j<=myend; j++)
			{
				Tlocal[i][j]=Tlocal[i][j]*h[i] - Tlocal[l][j]*hh[i];
			}
		}

		faulty_location=2;
		// collect chunks of last row of K to "future" last node
		MPI_Gather (&Tlocal[l-1][local[n]], myTcols/2, MPI_DOUBLE, &TlastKr[0], 1, TlastKr_chunks_resized, map[l-1], comm);
		faulty_location=20;

		//future last node broadcasts last row and col of K
		if (rank==map[l-1])
		{
			// copy data into local buffer before broadcast
			for (i=0; i<n; i++)
			{
				TlastKc[i]=Tlocal[i][local[n+l-1]];
			}
		}

		// inject ERROR
		if (l==failing_level)
		{
			//MPI_Barrier(current_comm);
			//sleep(5);
			if (rank==failing_rank)
			{
				printf("\n===============================\n");
				printf("rank %d failing at level %d",failing_rank,failing_level);
				printf("\n===============================\n");
				raise(SIGKILL);
			}
			//else
			{
				faulty=1;
				sleep(1);
				break;
			}
		}

		faulty_location=3;
		//TODO: substitute Gather with an All-to-All
		MPI_Bcast (&TlastK[0][0], Tcols, MPI_DOUBLE, map[l-1], comm);
		faulty_location=30;
	}

	// last level (l=0)
	for (i=0; i<myTcols/2; i++)
	{
		for(rhs=0;rhs<m;rhs++)
		{
			xx[global[i]][rhs]=xx[global[i]][rhs]+Tlocal[0][i]*bb[0][rhs];
		}
	}

	faulty_location=4;
	// collect solution
	// MPI_IN_PLACE required for MPICH based versions
	if (rank==0)
	{
		MPI_Gather (MPI_IN_PLACE, 1, xx_rows_interleaved_resized, &xx[0][0], 1, xx_rows_interleaved_resized, 0, comm);
	}
	else
	{
		MPI_Gather (&xx[rank][0], 1, xx_rows_interleaved_resized, &xx[0][0], 1, xx_rows_interleaved_resized, 0, comm);
	}
}
	// cleanup
	free(local);
	free(global);
	free(map);

	DeallocateMatrix2D(TlastK,2,CONTIGUOUS);
	DeallocateVector(h);
	DeallocateVector(hh);
	DeallocateMatrix2D(Tlocal,n,CONTIGUOUS);
}
