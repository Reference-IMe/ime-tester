/*  
 * Copyright (c) 2011-2013 The University of Tennessee and The University                                                                          
 *                         of Tennessee Research Foundation.  All rights
 *                         reserved.
 *
 *
 * $COPYRIGHT$
 * 
 * Additional copyrights may follow
 * 
 * $HEADER$
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "slp.h"
#include "util_matrix.h"
#include "util_inject.h"
#include "ftla_ftwork.h"
#include "ftla_cof.h"
#include "ftla_driver.h"

static int i0=0, i1=1;
static double p0=0e0, p1=1e0, m1=-1e0;
static double eps;
static int rank, np;

#define print0( args ) do { if( 0 == rank ) { printf args ; fflush( stdout ); } } while( 0 )

#define checkerror(a,b) do { if( (a) != (b) ) { \
   printf( "error on line %d in %s\n", __LINE__, __FILE__ );\
   exit( 0 );} \
} while( 0 )

static void create_matrix (int ctxt, int seed, double **A, int *descA, int M, int N, int NB, int *np_A, int *nq_A);
static double condA( int M, int N, double* A0, int* descA );
static double verifyQRm( int M, int N, double *A0, double *Aqr, int *descA, double *tau );
static double verifyQRb( int M, int N, int S, double *A0, double *Aqr, int *descA, double *tau, double *B0, double *X, int *descX );
static double verifyX( int M, int S, double *X0, double *X, int *descX );


int *errors;

int main( int argc, char **argv ) 
{
    int ictxt, info;
    int P=1, Q=1, PxQ=1;
    int nprow, npcol, myrow, mycol;
    ftla_work_t ftwork;

    int NB=100;
    int M, N, Nc, Ne, S=1;
    int Fstrat='e', Fmin=0, Fmax=0, Finc=1, F;
    double *A0=NULL, *A=NULL;    
    int descA[9]; 
    double *B0=NULL, *X0=NULL, *Xo=NULL, *Xf=NULL;
    int descX[9];

    double t1, t2, To, Tf;
    double Rom, Rfm, Rob, Rfb, Rox, Rfx, KA;
    int err = 0;
    int i, start, end, step;
    start = end = step = 640;

    {/* read command line parameter */
        for( i = 1; i < argc; i++ ) {
            if( strcmp( argv[i], "-p" ) == 0 ) {
                P = atoi(argv[i+1]);
                i++;
            }
            if( strcmp( argv[i], "-q" ) == 0 ) {
                Q = atoi(argv[i+1]);
                i++;
            }
            if( strcmp( argv[i], "-b" ) == 0 ) {
                NB = atoi(argv[i+1]);
                i++;
            }
            if( strcmp( argv[i], "-s" ) == 0 ) {
                S = atoi(argv[i+1]);
                i++;
            }
            if( strcmp( argv[i], "-r" ) == 0 ) {
                start  = atoi(argv[i+1]);
                end = atoi(argv[i+2]);
                step = atoi(argv[i+3]);
                i+=3;
            }
            if( strcmp( argv[i], "-fl" ) == 0 ) {
                Fstrat = 'l';
                Fmin = atoi(argv[++i]);
                Fmax = atoi(argv[++i]);
                Finc = atoi(argv[++i]);
            }
            if( strcmp( argv[i], "-fr" ) == 0 ) {
                Fstrat = 'r';
                Fmin = atoi(argv[++i]);
                Fmax = atoi(argv[++i]);
                Finc = atoi(argv[++i]);
            }
            if( strcmp( argv[i], "-fe" ) == 0 ) {
                Fstrat = 'e';
                Fmin = atoi(argv[++i]);
                Fmax = atoi(argv[++i]);
                Finc = atoi(argv[++i]);
            }
            if( strcmp( argv[i], "-cof" ) == 0 ) {
                err = 1; // force cof reload
            }
        }
        PxQ = P * Q;
    }
    {/* init BLACS */
        Cblacs_pinfo( &rank, &np );
        print0( ("Setting up MPI and BLACS...") );
        Cblacs_get( -1, 0, &ictxt );
        Cblacs_gridinit( &ictxt, "Row", P, Q );
        Cblacs_gridinfo( ictxt, &nprow, &npcol, &myrow, &mycol );
        print0( ("done\tPxQ=%dx%d, NB=%d\n", nprow, npcol, NB) );
        eps = pdlamch_( &ictxt, "Epsilon" );
    }

    print0( ("M\tNe\t : time(s)\tgflops/proc\t(realflops)\t: |A-A~|\t|B-AX~|\t\t|X-X~|\n") );

    for( i = start; i <= end; i += step) {
        {/* allocate matrices */
            /* determine checksum size, generate A matrix */
            N = M = i;
            Nc = numroc_( &N, &NB, &mycol, &i0, &npcol ); //LOCc(N_A) 
            MPI_Allreduce( MPI_IN_PLACE, &Nc, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD );
#ifndef NO_EXTRAFLOPS
            Ne = N + Nc*2;// + ((Nc/NB)%Q==0)*NB;
#else
            Ne = N;
#endif
            print0( ("%-7d\t%-7d\t", M, Ne) );

            create_matrix( ictxt, 800, &A0, descA, M, Ne, NB, NULL, NULL );
            create_matrix( ictxt, 0,   &A,  NULL, M, Ne, NB, NULL, NULL );
            pdlacpy_( "All", &M, &N, A0, &i1, &i1, descA, A, &i1, &i1, descA );
/*            KA = condA( M, N, A, descA );
            print0( ("%-7.2g\t", KA) );
            pdlacpy_( "All", &M, &N, A0, &i1, &i1, descA, A, &i1, &i1, descA );*/
            /* allocate local buffer for the Q-wide local panel copy */
            create_matrix( ictxt, 0, (typeof(&A))&(ftwork.pcopy.Pc), ftwork.pcopy.descPc, M, (Q+2)*NB, NB, &(ftwork.pcopy.nrPc), &(ftwork.pcopy.ncPc) );
            /* generate X0, compute B from A.X0, set X=B0 */
            create_matrix( ictxt, 500, &X0, descX, M, S, NB, NULL, NULL );
            create_matrix( ictxt, 0,   &B0,  NULL, M, S, NB, NULL, NULL );
            pdgemm_( "N", "N", &M, &S, &N, &p1, A0, &i1, &i1, descA,
                                                X0, &i1, &i1, descX,
                                           &p0, B0, &i1, &i1, descX );
            create_matrix( ictxt, 0,   &Xo,  NULL, M, S, NB, NULL, NULL );
            create_matrix( ictxt, 0,   &Xf,  NULL, M, S, NB, NULL, NULL );
            pdlacpy_( "All", &M, &S, B0, &i1, &i1, descX, Xo, &i1, &i1, descX );
            pdlacpy_( "All", &M, &S, B0, &i1, &i1, descX, Xf, &i1, &i1, descX );
        }
        {/*    call ScaLAPACK QR without FT */
            int lwork=-1;
            double lazywork;
            pdgeqrf_( &M, &N, NULL, &i1, &i1, descA, NULL, &lazywork, &lwork, &info );
            lwork = (int)lazywork;
            double *work = (double*)malloc( lwork*sizeof(double) );
            double *tau  = (double*)malloc( Nc*sizeof(double) );

            MPI_Barrier( MPI_COMM_WORLD );
            t1 = MPI_Wtime();
            pdgeqrf_( &M, &N, A, &i1, &i1, descA, tau, work, &lwork, &info );
            checkerror( info, 0 );
            t2 = MPI_Wtime();
            t2 -= t1;
            MPI_Reduce( &t2, &To, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD );
            print0( ("o: %f\t%-11g\t\t\t: ", To, 4.0/3.0*M*M*M/1e9/To/PxQ) );
            dprintmatrix(A, descA, "Ao", 2);
            
            Rom = verifyQRm( M, N, A0, A, descA, tau );
            print0( ("%-10e\t", Rom) );
            Rob = verifyQRb( M, N, S, A0, A, descA, tau, B0, Xo, descX );
            print0( ("%-10e\t", Rob) );
            Rox = verifyX( M, S, X0, Xo, descX );
            print0( ("%-10e\n", Rox) );

            free( work );
            free( tau );
        }        

        {/* call resilient QR */
            int err;
            int lwork=-1;
            double lazywork;
            pdgeqrf_( &M, &Ne, NULL, &i1, &i1, descA, NULL, &lazywork, &lwork, &info );
            lwork = (int)lazywork;
            double *work = (double*)malloc( lwork*sizeof(double) );
            double *tau  = (double*)malloc( Nc*sizeof(double) );
            
#ifdef INJECT        
          for( F = Fmin; F<=Fmax; F+=Finc ) {
            errors = create_error_list( M, NB, F, Fstrat );
#endif
            /* Reset A for the second run */
            pdlacpy_( "All", &M, &N, A0, &i1, &i1, descA, A, &i1, &i1, descA );
            pdlacpy_( "All", &M, &S, B0, &i1, &i1, descX, Xf, &i1, &i1, descX );
            
            Cftla_work_construct( 0, descA, 0, Ne-N, &ftwork );
            
            MPI_Barrier( MPI_COMM_WORLD );
            t1 = MPI_Wtime();
            do { // call ftpdgeqrf until we complete w/o a failure
#ifdef USE_CoF
              if(err) Cftla_cof_dgeqrr( A, descA, tau, work, lwork, &ftwork );
#endif
              err = ftla_pdgeqrf( &M, &Ne, A, &i1, &i1, descA, tau, work, &lwork, &info, (int*)&ftwork );
              //printf("****** err=%d i1=%d\n",err,i1);
            
#ifdef USE_CoF
              if(err) Cftla_cof_dgeqrr( A, descA, tau, work, lwork, &ftwork );
#endif
            } while(err);
            checkerror( info, 0 );
            t2 = MPI_Wtime();
            t2 -= t1;
            MPI_Reduce( &t2, &Tf, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD );
            print0( ("\t% 8df: %f\t%-11g\t%-11g\t: ", F, Tf, 4.0/3.0*M*M*M/1e9/Tf/PxQ, 4.0/3.0*M*M*Ne/1e9/Tf/PxQ) );
            dprintmatrix(A, descA, "Af", 2);
                        
            Rfm = verifyQRm( M, N, A0, A, descA, tau );
            print0( ("%-10e\t", Rfm) );
            Rfb = verifyQRb( M, N, S, A0, A, descA, tau, B0, Xf, descX );
            print0( ("%-10e\t", Rfb) );
            Rfx = verifyX( M, S, X0, Xf, descX );
            print0( ("%-10e\n", Rfx) );

            Cftla_cof_cleanup( &ftwork );
            Cftla_work_destruct( &ftwork );
            free( errors );
        
        {/* Final prints */
            char who[3]={' ', ' ', '\0'};
            double Rofx;
            if( isnan(Rom) || Rom>1e-10 )
                who[1]='o';
            if( isnan(Rfm) || Rfm>1e-10 ) {
                if( who[1]=='o' )
                    who[1]='a';
                else
                    who[1]='f';
            }
            if( ' ' != who[1] ) who[0]='!';
            Rofx = verifyX( M, S, Xo, Xf, descX );
            print0( ("%s\t\t :%12coverhead%% = %-7.2f\t\t: \t\t%16s%-10e\n", who, ' ', (Tf-To)/To*100, "|Xo-Xf| = ", Rofx) );
        }
#ifdef INJECT
          }
#endif
            free( work );
            free( tau );
        }
        {/* Cleanup */
            if( NULL != A0 ) free( A0 );
            if( NULL != A  ) free( A );
            A0 = A = NULL;
            if( NULL != B0 ) free( B0 );
            if( NULL != X0 ) free( X0 );
            if( NULL != Xo ) free( Xo );
            if( NULL != Xf ) free( Xf );
            B0 = X0 = Xo = Xf = NULL;

            if( NULL != ftwork.pcopy.Pc ) free( ftwork.pcopy.Pc);
            ftwork.pcopy.Pc = NULL;
        }
    }

    fflush( stdout );
    Cblacs_gridexit( ictxt );
    MPI_Finalize();
    return 0;
}

#define MIN(a,b) ((a>b)?b:a)
#define MAX(a,b) ((a>b)?a:b)

/*
 * produce distributed matrix,  
 */
void create_matrix( int ctxt, int seed, double **A, int *descA, 
        int M, int N, int NB, int *np_A, int *nq_A) 
{
    int info;
    int nprow, npcol, myrow, mycol;
    Cblacs_gridinfo( ctxt, &nprow, &npcol, &myrow, &mycol );

    // allocate the generator matrix and check matrix
    int np_iA = numroc_( &M, &NB, &myrow, &i0, &nprow );
    int nq_iA = numroc_( &N, &NB, &mycol, &i0, &npcol );

    if (np_iA*nq_iA!=0)
    {
        *A = malloc(np_iA*nq_iA*sizeof(**A)) ;
        if (*A == NULL) Cblacs_abort( ctxt, 10 );
        memset (*A, 0, np_iA*nq_iA*sizeof(**A));
    }
    else *A = NULL;

    if (descA != NULL)
    {
        int itemp = MAX( 1, np_iA );
        descinit_( descA, &M, &N, &NB, &NB, &i0, &i0, &ctxt, &itemp, &info );
        if (info != 0) Cblacs_abort( ctxt, 12 );
    }

    if (seed)
    {
        // fill in random numbers
        pdmatgen_ (&ctxt, "N", "N", &M, &N, &NB, &NB, *A, 
                descA+8, descA+6, descA+7, 
                &seed, &i0, &np_iA, &i0, &nq_iA, 
                &myrow, &mycol, &nprow, &npcol);
    }

    /* set np and nq */
    if (np_A != NULL)
        *np_A = np_iA;
    if (nq_A != NULL)
        *nq_A = nq_iA;
}

/* Verify that Q*R=A */
static double verifyQRm( int M, int N, double *A0, double *Aqr, int *descA, double *tau ) {
    double resid;
    double *work=NULL;
    double *A=NULL;
    //----- non-ft -----//
    int Mp0, Nq0, iarow, iacol, lwork;
    // grid parameters
    int ctxt=descA[CTXT_];
    int nb = descA[NB_];
    int nprow, npcol, myrow, mycol;
    Cblacs_gridinfo( ctxt, &nprow, &npcol, &myrow, &mycol );

    //----------- allocate workspace ----------// -_-||
    //LWORK = NB_A * ( 2*Mp0 + Nq0 + NB_A ), where
    //Mp0   = NUMROC( M+IROFF, MB_A, MYROW, IAROW, NPROW ) * NB_A,
    //Nq0   = NUMROC( N+ICOFF, NB_A, MYCOL, IACOL, NPCOL ) * MB_A,
    //IROFF = MOD( IA-1, MB_A ), ICOFF = MOD( JA-1, NB_A ),
    //IAROW = INDXG2P( IA, DESCA( MB_ ), MYROW, DESCA( RSRC_ ),NPROW ),
    //IACOL = INDXG2P( JA, DESCA( NB_ ), MYCOL, DESCA( CSRC_ ),NPCOL )
    iarow = indxg2p_( &i1, &nb, &myrow, &i0, &nprow );
    iacol = indxg2p_( &i1, &nb, &mycol, &i0, &npcol );
    Mp0 = nb * numroc_( &M , &nb, &myrow, &iarow, &nprow );
    Nq0 = nb * numroc_( &N , &nb, &mycol, &iacol, &npcol );
    lwork = nb * ( 2*Mp0 + Nq0 + nb );
    work = (double *)malloc( lwork*sizeof(double) );

    /* Work on a copy of A */
    create_matrix( ctxt,  0, &A, NULL, M, N, nb, NULL, NULL );
    pdlacpy_( "All", &M, &N, Aqr, &i1, &i1, descA, A, &i1, &i1, descA );

    pdgeqrrv_( &M, &N, A, &i1, &i1, descA, tau, work );
    pdmatadd_( &M, &N, &p1, A0, &i1, &i1, descA, &m1, A, &i1, &i1, descA );
    resid = pdlange_( "F", &M, &N, A, &i1, &i1, descA, NULL ) / 
            pdlange_( "F", &M, &N, A0, &i1, &i1, descA, NULL ) / 
            M;

    free( work );
    if( NULL != A ) free( A );
    return resid;
}

/* Verify the backward error of the solution of Ax=B obtained from the 
 * QR decomposition of A */
static double verifyQRb( int M, int N, int S, double *A0, double *Aqr, int *descA, double *tau, double* B0, double *X, int *descX ) {
    double  *A=NULL, *R=NULL;
    int     ctxt=descA[CTXT_];
    int     nb = descA[NB_];
    int     lwork; double lazywork; double *work=NULL;
    int     info;
    double  Anorm, Xnorm, Rnorm, resid;
    
    create_matrix( ctxt,  0,   &A,   NULL, M, N, nb, NULL, NULL );
    pdlacpy_( "All", &M, &N, Aqr, &i1, &i1, descA, A, &i1, &i1, descA );
    create_matrix( ctxt,  0,   &R,   NULL, M, S, nb, NULL, NULL );
    pdlacpy_( "All", &M, &S, B0,  &i1, &i1, descX, R, &i1, &i1, descX );

    /* Compute the solution of QRx=B */
    lwork = -1;
    pdormqr_( "L", "T", &M, &S, &M, A, &i1, &i1, descA, tau, 
              X, &i1, &i1, descX, &lazywork, &lwork, &info );
    lwork = (int) lazywork;
    work = (double*) malloc( lwork * sizeof(double) );
    //multiply X = Q**T * X , where Q is found in A
    //since X is passed Xo or Xf, which are init to B --> X=Q**T * B
    pdormqr_( "L", "T", &M, &S, &M, A, &i1, &i1, descA, tau, 
              X, &i1, &i1, descX, work, &lwork, &info );
    //solves A*X = alpha*X,
    //where A is upper triangular (U), assumed not unit triangular, alpha=p1=1.0
    //since X=Q**T * B --> solves A*X=Q**T * B
    pdtrsm_( "L", "U", "N", "N", &M, &S, &p1, A, &i1, &i1, descA,
             X, &i1, &i1, descX );
    
    /* Compute norms and backward error residual */
    Anorm = pdlange_( "F", &M, &N, A0, &i1, &i1, descA, NULL );
    Xnorm = pdlange_( "F", &M, &S,  X, &i1, &i1, descX, NULL );
    //multiply R = A0 * X
    pdgemm_( "N", "N", &M, &S, &N, &p1, A0, &i1, &i1, descA, 
                                         X, &i1, &i1, descX, 
                                   &m1,  R, &i1, &i1, descX );
    Rnorm = pdlange_( "F", &M, &S,  R, &i1, &i1, descX, NULL );
    resid = Rnorm / ( Anorm * Xnorm * eps );
    free( work );
    if( NULL != A ) free( A ); if( NULL != R ) free( R );
    return resid;
}

/* Compute the norm of X-X~ */
static double verifyX( int M, int S, double *X0, double *X, int *descX ) {
    int ctxt=descX[CTXT_];
    int nb = descX[NB_];
    double *Xs=NULL;
    double resid;
    create_matrix( ctxt,  0, &Xs,  NULL, M, S, nb, NULL, NULL );
    pdlacpy_( "All", &M, &S, X, &i1, &i1, descX, Xs, &i1, &i1, descX );

    pdmatadd_( &M, &S, &p1, X0, &i1, &i1, descX, 
                       &m1, Xs, &i1, &i1, descX );
    resid = pdlange_( "F", &M, &S, Xs, &i1, &i1, descX, NULL ) / 
            pdlange_( "F", &M, &S, X0, &i1, &i1, descX, NULL );
    
    if( NULL != Xs ) free( Xs );
    return resid;
}

/* Compute Condition number of A */
static double condA( int M, int N, double* A0, int* descA ) {
    double lazywork;
    int info, lwork, liwork;
    double* work; int *iwork;
    double rcond;
    double anorm;
    int ictxt=descA[CTXT_];
    int nb = descA[NB_];
    int nprow, npcol, myrow, mycol;
    Cblacs_gridinfo( ictxt, &nprow, &npcol, &myrow, &mycol );
        
    lwork = numroc_( &M, &nb, &myrow, &i1, &myrow );
    work = malloc( lwork*sizeof(double) );
    anorm = pdlange_( "I", &M, &N, A0, &i1, &i1, descA, work );
    free( work );
    iwork = malloc( (M+nb)*sizeof(int) );
    pdgetrf_( &M, &N, A0, &i1, &i1, descA, iwork, &info );
    free( iwork );
    lwork = liwork = -1;
    pdgecon_( "I", &M, A0, &i1, &i1, descA, &anorm, &rcond, &lazywork, &lwork, &liwork, &liwork, &info );
    lwork = (int)lazywork;
    work = malloc( lwork*sizeof(double) );
    iwork = malloc( liwork*sizeof(int) );
    pdgecon_( "I", &N, A0, &i1, &i1, descA, &anorm, &rcond, work, &lwork, iwork, &liwork, &info );
    free( work ); free( iwork );
    return rcond;
}
    
#if 0
/* Verify orthogonality of QR decomposition */
double verifyQRo() {
    
}
#endif
