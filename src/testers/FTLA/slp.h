#ifndef __INLCUDE_SLP_H__
#define __INLCUDE_SLP_H__
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

#define DTYPE_      0
#define CTXT_       1
#define M_          2
#define N_          3
#define MB_         4
#define NB_         5
#define RSRC_       6
#define CSRC_       7
#define LLD_        8
#define TAG_L    1021
#define TAG_LL    904    

#include <math.h>
#include <complex.h>

typedef float complex complexfloat_t;
typedef double complex complexdouble_t;

/********************************************************************/
extern double get_cur_time();

extern int iceil_( int * A, int * B );

extern void Cblacs_pinfo( int * iam, int * nprocs );
extern void Cblacs_get( int, int, int * ictxt );
extern void Cblacs_set( int ictxt, int what, int * val);
extern void Cblacs_gridinfo( int ictxt, int * nprow, int * npcol, int * myrow, int * mycol );
extern void Cblacs_gridinit( int *, char *, int nprow, int npcol );
extern void Cblacs_gridexit( int ictxt );
extern void Cblacs_abort( int ictxt, int ErrNo);
extern int  Cblacs_pnum( int ictxt, int prow, int pcol );
extern void Cblacs_pcoord( int ictxt, int pnum, int *prow, int *pcol );

extern void Cblacs_barrier(int ConTxt, char *scope);

extern void Csgebs2d ( int ConTxt, char *scope, char *top, int m, int n, float *A, int lda);
extern void Csgebr2d ( int ConTxt, char *scope, char *top, int m, int n, float *A, int lda, int rsrc, int csrc);
extern void Cpsgemr2d( int m, int n, float *A, int IA, int JA, int *descA, float *B, int IB, int JB, int *descB, int gcontext);
extern void Cstrbs2d ( int ConTxt, char *scope, char *top, char *uplo, char *diag, int m, int n, float *A, int lda);
extern void Cstrbr2d ( int ConTxt, char *scope, char *top, char *uplo, char *diag, int m, int n, float *A, int lda, int rsrc, int csrc);

extern void Cdgebs2d ( int ConTxt, char *scope, char *top, int m, int n, double *A, int lda);
extern void Cdgebr2d ( int ConTxt, char *scope, char *top, int m, int n, double *A, int lda, int rsrc, int csrc);
extern void Cpdgemr2d( int m, int n, double *A, int IA, int JA, int *descA, double *B, int IB, int JB, int *descB, int gcontext);
extern void Cdtrbs2d ( int ConTxt, char *scope, char *top, char *uplo, char *diag, int m, int n, double *A, int lda);
extern void Cdtrbr2d ( int ConTxt, char *scope, char *top, char *uplo, char *diag, int m, int n, double *A, int lda, int rsrc, int csrc);

extern void Ccgebs2d ( int ConTxt, char *scope, char *top, int m, int n, float complex *A, int lda);
extern void Ccgebr2d ( int ConTxt, char *scope, char *top, int m, int n, float complex *A, int lda, int rsrc, int csrc);
extern void Cpcgemr2d( int m, int n, float complex *A, int IA, int JA, int *descA, float complex *B, int IB, int JB, int *descB, int gcontext);
extern void Cctrbs2d ( int ConTxt, char *scope, char *top, char *uplo, char *diag, int m, int n, float complex *A, int lda);
extern void Cctrbr2d ( int ConTxt, char *scope, char *top, char *uplo, char *diag, int m, int n, float complex *A, int lda, int rsrc, int csrc);

extern void Czgebs2d ( int ConTxt, char *scope, char *top, int m, int n, double complex *A, int lda);
extern void Czgebr2d ( int ConTxt, char *scope, char *top, int m, int n, double complex *A, int lda, int rsrc, int csrc);
extern void Cpzgemr2d( int m, int n, double complex *A, int IA, int JA, int *descA, double complex *B, int IB, int JB, int *descB, int gcontext);
extern void Cztrbs2d ( int ConTxt, char *scope, char *top, char *uplo, char *diag, int m, int n, double complex *A, int lda);
extern void Cztrbr2d ( int ConTxt, char *scope, char *top, char *uplo, char *diag, int m, int n, double complex *A, int lda, int rsrc, int csrc);

extern void descinit_( int * desc, int * M, int * N, int * MB, int * NB, int * IRSRC, int * ICSRC, int * ICTXT, int * LLD, int * INFO );
extern int  numroc_  ( int * N, int * NB, int * IPROC, int * ISRCPROC, int * NPROCS );
extern int infog2l_ ( int * GRINDX, int * GCINDX, int * DESC, int * NPROW, int * NPCOL, int * MYROW, int * MYCOL, 
                       int * LRINDX, int * LCINDX, int * RSRC, int * CSRC );
extern int infog1l_ ( int * GRINDX, int * NB, int * NPROC, int * MYPROC, int * ISRCPROC, int * LINDX, int * ROCSRC ); 
extern int  indxl2g_ ( int * INDXLOC, int * NB, int * IPROC, int * ISRCPROC, int * NPROCS );
extern int  indxg2p_ ( int * INDXGLOB, int * NB, int * IPROC, int * ISRCPROC, int * NPROCS );

//extern void pb_topget_( int * ICTXT, char * OP, char * SCOPE, char * TOP );
//extern void pb_topset_( int * ICTXT, char * OP, char * SCOPE, char * TOP );
//extern void igamn2d_   ( int * ICTXT, char * SCOPE, char * TOP, int * M, int * N, int * A, int * LDA, int * RA, int * CA, 
//                         int * RCFLAG, int * RDEST, int * CDEST );

//extern void  chk1mat_ (int * MA, int * MAPOSO, int * NA, int * NAPOS0, int * IA, int * JA, int * DESCA, int * DESCAPOS0, int * INFO);
//extern void pchk1mat_ (int * MA, int * MAPOSO, int * NA, int * NAPOS0, int * IA, int * JA, int * DESCA, 
//                       int * DESCAPOS0, int * NEXTRA, int * EX, int * EXPOS, int * INFO);

// LAPACK
//extern void   dgemm_ ( char *transa,char *transb,int *m,int *n,int *k,double *alpha,double *a,int *lda,double *b,int *ldb,double *beta,double *c,int *ldc );
extern void   dtrsm_ ( const char *side, const char *uplo, const char *transa, const char *diag, const int *m, const int *n, const double *alpha, const double *a, const int *lda, double *b, const int *ldb );
extern void   dsyrk_ ( const char *uplo, const char *trans, const int *n, const int *k, const double *alpha, const double *a, const int *lda, const double *beta, double *c, const int *ldc );
extern void   dscal_ ( int *, double *, double *, int * );
extern float slamch_( char * type);
extern double dlamch_( char * type);

extern float slange_( char *norm, int *m, int *n, float *A, int *lda, float *work );
extern double dlange_( char *norm, int *m, int *n, double *A, int *lda, double *work );
extern double clange_( char *norm, int *m, int *n, float complex *A, int *lda, float *work );
extern double zlange_( char *norm, int *m, int *n, double complex *A, int *lda, double *work );

// SCALAPACK

extern void pdgeqrf_ ( int *m, int *n, double *A, int *ia, int *ja, int *descA, double *tau, double *work, int *lwork, int *info );
extern void pdgetrf_ ( int *M, int *N, double *A, int *ia, int *ja, int *descA, int *ipiv, int *info );
extern void pdgesv_  ( int *n, int *nrhs, double *A, int *ia, int *ja, int *descA, int *ipiv, double *B, int *ib, int *jb, int *descB, int *info );

extern void pdgetrs_( char *trans, int*n, int *nrhs, double *A, int *ia, int *ja, int *descA, int *ipiv, double *B, int *ib, int *jb, int *descB, int *info );


//extern int pdtrsm_ ( char * SIDE, char * UPLO, char * TRANS, char * DIAG, int * M, int * N, 
//                        double * ALPHA, double * A, int * IA, int * JA, int * DESCA,
//                                        double * B, int * IB, int * JB, int * DESCB );
//extern int pdgemm_ ( char * TRANSA, char * TRANSB,    int * M, int * N, int * K,
//                        double * ALPHA, double * A, int * IA, int * JA, int * DESCA,
//                                          double * B, int * IB, int * JB, int * DESCB, 
//                          double * BETA,    double * C, int * IC, int * JC, int * DESCC );

extern void pdormqr_( char *side, char *tran, int *m, int *n, int *k, 
                                        double *A, int *ia, int *ja, int *descA,
                           double *tau, double *C, int *ic, int *jc, int *descC, 
                           double *work, int *lwork, int *info );
extern void pdtrmm_ ( char *side, char *uplo, char *transa, char *diag, int *m, int *n, double *alpha, double *A, int *ia, int *ja, int *desca, double *B, int *ib, int *jb, int *descb );

extern void pdamax_ ( int * N, double * AMAX, int * INDX, double * X, int * IX, int * JX, int * DESCX, int * INCX );
extern void pdlapiv_( char * direc, char * rowcol, char * pivroc, int * m, int * n, double * A, int * ia, int * ja, int * descA, int * ipiv, int * ip,int * jp, int * descIP, int * iwork );
//extern int pdlaswp_( char * direc, char * rwocol, int * n, double * A, int * ia, int * ja, int * descA, int * k1, int * k2, int * ipiv );

extern float pslamch_( int * ictxt, char * type);
extern double pdlamch_( int * ictxt, char * type);
extern float pslange_( char *norm, int *m, int *n, float *a, int *ia, int *ja, int *desca, float *work);
extern double pdlange_( char *norm, int *m, int *n, double *a, int *ia, int *ja, int *desca, double *work);
extern float pclange_( char *norm, int *m, int *n, float complex *a, int *ia, int *ja, int *desca, double *work);
extern double pzlange_( char *norm, int *m, int *n, double complex *a, int *ia, int *ja, int *desca, double *work);
extern void   pdgecon_( char * NORM, int * N, double * A, int * IA, int * JA, int * DESCA, double * ANORM, double * RCOND, double * WORK, int * LWORK, int * IWORK, int * LIWORK, int * INFO ); 

extern void pdelset_ (double *A, int *IA, int *JA, int *descA, double *alpha);
extern void pdelget_ (char *scope, char *top, double *alpha, double *A, int *IA, int *JA, int *descA);
extern void pdlacpy_ (char *, int *, int *, double *, int *, int *, int *, double *, int *, int *, int *);
extern void pdcopy_( int * N, double * X, int * IX, int * JX, int * DESCX, int * INCX,double * Y, int * IY, int * JY, int * DESCY, int * INCY );

extern void pddot_( int * N, double * DOT, double * X, int * IX, int * JX, int * DESCX, int * INCX, double * Y, int * IY, int * JY, int * DESCY, int * INCY );
extern void pdlafchk_ (char *aform, char *diag, int *M, int *N, double *A, int *ia, int *ja, int *descA, int *iaseed, double *anorm, double *fresid, double *work);

extern void pdgeadd_( char *TRANS, int * M, int * N,  double * ALPHA, double * A, int * IA, int * JA, int * DESCA, 
                                                      double * BETA, double * C, int * IC, int * JC, int * DESCC ); //sub( C ) := beta*sub( C ) + alpha*op( sub( A ) )

extern void psmatadd_ ( int *, int *, float *, float *, int *, int *, int *, float *, float *, int *, int *, int *);
extern void psmatgen_ (int *, char *, char *, int *, int *, int *, int *, float *, int *, int *, int *, int *, int *, int *, int *, int *, int *, int *, int *, int *);
extern void pslaprnt_ (int *m, int * n, float * A, int * ia, int * ja, int * descA, int * prow, int * pcol, 
                        char * NAME, int * NOUT, float * work, int len);
extern void pdmatadd_ ( int *, int *, double *, double *, int *, int *, int *, double *, double *, int *, int *, int *);
extern void pdmatgen_ (int *, char *, char *, int *, int *, int *, int *, double *, int *, int *, int *, int *, int *, int *, int *, int *, int *, int *, int *, int *);
extern void pdlaprnt_ (int *m, int * n, double * A, int * ia, int * ja, int * descA, int * prow, int * pcol, char * NAME, int * NOUT, double * work, int len);
extern void pcmatadd_ ( int *, int *, float complex *, float complex *, int *, int *, int *, float complex *, float complex *, int *, int *, int *);
extern void pcmatgen_ (int *, char *, char *, int *, int *, int *, int *, float complex *, int *, int *, int *, int *, int *, int *, int *, int *, int *, int *, int *, int *);
extern void pclaprnt_ (int *m, int * n, float complex * A, int * ia, int * ja, int * descA, int * prow, int * pcol, 
                        char * NAME, int * NOUT, double * work, int len);
extern void pzmatadd_ ( int *, int *, double complex *, double complex *, int *, int *, int *, double complex *, double complex *, int *, int *, int *);
extern void pzmatgen_ (int *, char *, char *, int *, int *, int *, int *, double complex *, int *, int *, int *, int *, int *, int *, int *, int *, int *, int *, int *, int *);
extern void pzlaprnt_ (int *m, int * n, double complex * A, int * ia, int * ja, int * descA, int * prow, int * pcol, 
                        char * NAME, int * NOUT, double complex * work, int len);

extern void pdpotrrv_ (char * uplo, int * m, double * A, int * i, int * j, int * descA, double *work);
extern void pdgetrrv_ (int *m, int *n, double *A, int *ia, int *ja, int *descA, int *ipiv, double *work);
extern void pdgeqrrv_ (int *m, int *n, double *A, int *ia, int *ja, int *descA, double*tau, double *work);

#endif
