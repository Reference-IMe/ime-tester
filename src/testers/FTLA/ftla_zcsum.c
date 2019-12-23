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

/*
 * @precisions normal z -> s d c
 */


#include "ftla_csum.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "util_matrix.h"
#include "slp.h"


static void recover_checksum( char side, complexdouble_t *A, int m, int n, int *descA, ftla_work_t *ftw );
static void recover_abft( complexdouble_t *A, int m, int n, int *descA, ftla_work_t *ftw );
static void hor_checkpointing( complexdouble_t *A, int *descA, int IA, int JA, int iscopy, ftla_work_t *ftw );
static void hor_checkpointing_local (complexdouble_t *A, int *descA, int IA, int JA, int iscopy, ftla_work_t *ftw);
static void recover_local (char side, complexdouble_t *A, int MA, int NA, int *descA, ftla_work_t *ftw);
static void recover_local_fillin (complexdouble_t *A, int MA, int NA, int *descA, ftla_work_t *ftw);
static void create_matrix( int ctxt, complexdouble_t **A, int *descA, int M, int N, int NB, int *np_A, int *nq_A);

static int i0=0, i1=1;

#define MIN(a,b) ((a>b)?b:a)
#define MAX(a,b) ((a>b)?a:b)

#if 0
#define ROOTSAY(s) do { if ((myrow+mycol)==0) printf ("%s\n", s); MPI_Barrier (MPI_COMM_WORLD); } while(0)
#else 
#define ROOTSAY(s)
#endif

/*********************************************************************
 ***** RECOVERY 
 */
void Cftla_pzcsumr ( char zone, complexdouble_t *A, int m, int n, int *descA, ftla_work_t *ftw)
{
    int ctxt = descA[CTXT_];
    int nprow, npcol, myrow, mycol, myrank, nranks;
    Cblacs_gridinfo( ctxt, &nprow, &npcol, &myrow, &mycol );
    Cblacs_pinfo( &myrank, &nranks );
            
    // broadcast victim's id
    MPI_Allreduce ( MPI_IN_PLACE, &ftw->errproc, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

    zprintmatrix(A, descA, "Af_hole", 7);

    // recover checksum
    if( zone == 'C' || zone == 'A' ) {
        recover_checksum ('b', A, m, n, descA, ftw );
        ROOTSAY ("checksum recovered");
        zprintmatrix(A, descA, "Af_csum", 7);
    }

    // recover non-Q panels 
    if( zone == 'I' || zone == 'A' ) {
        recover_abft (A, m, n, descA, ftw );
        ROOTSAY ("ABFT region recovered");
        zprintmatrix(A, descA, "Af_abft", 7);
    }

    // repair localcopy Pc
    if( zone == 'L' || zone == 'A' ) {
        recover_local ('b', ftw->pcopy.Pc, m, n, ftw->pcopy.descPc, ftw );
        ROOTSAY ("localcopy recovered");
        zprintmatrix(A, descA, "Af_lcpy", 7);
    }
}
extern void ftla_pzcsumr_( char *zone, complexdouble_t *A, int *m, int *n, int *descA, int *ftw ) {
    Cftla_pzcsumr( zone[1], A, *m, *n, descA, (ftla_work_t*)ftw );
}

/*
 * Use the duplicate checksum to restore missing checksum blocks
 */
static void recover_checksum (char side, complexdouble_t *A, int MA, int NA, int *descA, ftla_work_t *ftw)
{
    int i, j, k;
    int tt[2], g[2];

    // grid parameters
    int ctxt = descA[CTXT_];
    int nb = descA[NB_];
    int nprow, npcol, myrow, mycol, myrank, nranks;
    Cblacs_gridinfo( ctxt, &nprow, &npcol, &myrow, &mycol );
    Cblacs_pinfo( &myrank, &nranks );
    int iamvic = (myrank==ftw->errproc);

    // get nchkc
    int np_A = numroc_( &MA, &nb, &myrow, &i0, &nprow );

    int NAA = NA+ftw->Nc;
    int nq_A = numroc_( &NA, &nb, &mycol, &i0, &npcol );
    int nqq_A = numroc_( &NAA, &nb, &mycol, &i0, &npcol );

    // horizontal, and only horizontal
    if (side == 'b' || side == 'r')
    {
        /*
        if (iamvic)
            printf ("(%d,%d): recovering checksum\n",myrow, mycol);
            */

        if (nqq_A<nq_A)
        {
            printf ("(%d,%d): quick return from recovering checksum\n",myrow, mycol);
            return;        // quick return is victim doesn't carry checksum
        }
        //int tc = ((NA-JA+1)/(nb*npcol)+((NA-JA+1)%(nb*npcol)!=0)-1)*2*nb+NA+1;
        //Cpzgemr2d (MA-IA+1, nb, A, IA, tc, descA, A, IA, tc+nb, descA, ctxt);
        int i_c, j_c;
        if (iamvic)
        {
            i_c = np_A;
            j_c = (nqq_A-nq_A);
        }
        MPI_Bcast( &i_c, 1, MPI_INT, ftw->errproc, MPI_COMM_WORLD );
        MPI_Bcast( &j_c, 1, MPI_INT, ftw->errproc, MPI_COMM_WORLD );

        for (i=0; i<i_c; i+=nb)
            for (j=0, k=0; j<j_c; j+=nb, k++)
            {
                if (iamvic)
                {
                    int li = i+1;
                    g[0] = indxl2g_( &li, &nb, &myrow, &i0, &nprow );
                    int j_1 = nq_A + j + 1;
                    g[1] = indxl2g_( &j_1, &nb, &mycol, &i0, &npcol );

                    tt[1]=(((g[1]-NA-1)/nb)%2==0)?(g[1]+nb):(g[1]-nb);
                    tt[0]=g[0];
                    //printf ("(%d,%d): block (%d,%d) <------ (%d,%d)\n", myrow, mycol, g[0], g[1], tt[0],tt[1]);
                }
                MPI_Bcast( tt, 2, MPI_INT, ftw->errproc, MPI_COMM_WORLD );
                MPI_Bcast( g, 2, MPI_INT, ftw->errproc, MPI_COMM_WORLD );

                Cpzgemr2d (nb, nb, A, tt[0], tt[1], descA, A, g[0], g[1], descA, ctxt);
            }
    }
}

/*
 * recover the checksum protected area (update+final)
 */
static void recover_abft (complexdouble_t *A, int MA, int NA, int *descA, ftla_work_t *ftw)
{
    int j;
    int lrindx, lcindx, rsrc, csrc;
    complexdouble_t *Achk_old, *Achk_new;
    int descAchk[9], npp, nqq;

    int ctxt = descA[CTXT_];
    int nb = descA[NB_];
    int nprow, npcol, myrow, mycol, myrank, nranks;
    Cblacs_gridinfo( ctxt, &nprow, &npcol, &myrow, &mycol );
    Cblacs_pinfo( &myrank, &nranks );
    int iamvic = (myrank==ftw->errproc);

    // determine local data size (non-checksum)
    int np_A = numroc_( &MA, &nb, &myrow, &i0, &nprow );

    int NAA = NA+ftw->Nc;
    int nq_A = numroc_( &NA, &nb, &mycol, &i0, &npcol );
    //int nqq_A = numroc_( &NAA, &nb, &mycol, &i0, &npcol );



    create_matrix( ctxt, &Achk_old, descAchk, MA, NAA-NA, nb, &npp, &nqq);
    create_matrix( ctxt, &Achk_new, descAchk, MA, NAA-NA, nb, &npp, &nqq);

    // copy the old checksum into buffer
    Cpzgemr2d (MA, NAA-NA, A, 1, NA+1, descA, Achk_old, 1, 1, descAchk, ctxt);

    // re-calculate the checksum
    Cftla_pzcsum( 'r', A, MA, NA, 1, 1, descA, A, 1, NA+1, descA, ftw );

    // copy the new checksum into buffer
    Cpzgemr2d (MA, NAA-NA, A, 1, NA+1, descA, Achk_new, 1, 1, descAchk, ctxt);

    int N_M = NAA-NA, i1=1;
    complexdouble_t done = 1, mone = -1;
    pzmatadd_ ( &MA, &N_M, &done, Achk_old, &i1, &i1, descAchk, &mone, Achk_new, &i1, &i1, descAchk);
    Cpzgemr2d (MA, NAA-NA, Achk_new, 1, 1, descAchk, A, 1, NA+1, descA, ctxt);

    // copy resut back in place
    int j_c;
    if (iamvic)
        j_c = nq_A;
    MPI_Bcast( &j_c, 1, MPI_INT, ftw->errproc, MPI_COMM_WORLD );

    for (j=0; j<j_c; j+=nb)
    {
        int gi, gj, g;
        if (iamvic)
        {
            // source matrix column (global)
            int j_1 = j+1;
            g = indxl2g_( &j_1, &nb, &mycol, &i0, &npcol );
            gj = NAA-((int)(g/(npcol*nb))+1)*nb*2+1;
            //gj = NA+1+2*nb*((NA-g+1)/(nb*npcol)+(((NA-g+1)%(nb*npcol)!=0)-1));

            // source matrix row (global)
            gi = indxl2g_( &i1, &nb, &myrow, &i0, &nprow );
        }
        MPI_Bcast( &gi, 1, MPI_INT, ftw->errproc, MPI_COMM_WORLD );
        MPI_Bcast( &gj, 1, MPI_INT, ftw->errproc, MPI_COMM_WORLD );
            
        infog2l_ (&gi, &gj, descA, &nprow, &npcol, &myrow, &mycol, &lrindx, &lcindx, &rsrc, &csrc );
        //printf ("(%d,%d) <--- (%d,%d)\n", rsrc, ftw->errproc%npcol, rsrc, csrc);

        // copy the new checksum into buffer
        if (csrc==(ftw->errproc%npcol))
        {
            // local copy
            if (iamvic)
            {
                //printf ("(%d,%d) local copy when j=%d\n", myrow, mycol, j);
                memcpy (A+j*np_A, A+(lcindx-1)*np_A, np_A*nb*sizeof(complexdouble_t));
            }
        }
        else
        {
            // MPI copy
            if (myrow == rsrc)
            {
                MPI_Status status;
                if (iamvic)
                {
                    //printf ("(%d,%d) reciving from (%d,%d) for MPI_COPY when j=%d\n", myrow, mycol, myrow, csrc, j);
                    MPI_Recv(A+j*np_A, np_A*nb, MPI_DOUBLE_COMPLEX, csrc, 1, ftw->comms.rowcomm, &status);
                }
                else if (mycol==csrc)
                {
                    //printf ("(%d,%d) sending to (%d,%d) for MPI_COPY when j=%d\n", myrow, mycol, myrow, ftw->errproc%npcol, j);
                    MPI_Send(A+(lcindx-1)*np_A, np_A*nb, MPI_DOUBLE_COMPLEX, ftw->errproc%npcol, 1, ftw->comms.rowcomm);
                }
            }
        }
    }

    // copy the correct checksum back //
    Cpzgemr2d (MA, NAA-NA, Achk_old, 1, 1, descAchk, A, 1, NA+1, descA, ctxt);

    if (npp*nqq>0)
    {
        free (Achk_new);
        free (Achk_old);
    }
}

/*
 * Restore local panel copy from checksum for the failed rank
 */
static void recover_local (char side, complexdouble_t *A, int MA, int NA, int *descA, ftla_work_t *ftw)
{
    int i, j, k;
    int tt[2], g[2];

    int ctxt = descA[CTXT_];
    int nb = descA[NB_];
    int nprow, npcol, myrow, mycol, myrank, nranks;
    Cblacs_gridinfo( ctxt, &nprow, &npcol, &myrow, &mycol );
    Cblacs_pinfo( &myrank, &nranks );
    int iamvic = (myrank==ftw->errproc);

    // get nchkc
    int np_A = ftw->pcopy.nrPc; 

    NA = npcol*nb;
    int nq_A = numroc_( &NA, &nb, &mycol, &i0, &npcol );

    int nqq_A = ftw->pcopy.ncPc; 

    // horizontal, and only horizontal
    if (side == 'b' || side == 'r')
    {
        if (nqq_A<nq_A)
        {
            printf ("(%d,%d): quick return from recovering checksum\n",myrow, mycol);
            return;        // quick return is victim doesn't carry checksum
        }

        int i_c, j_c;
        if (iamvic)
        {
            i_c = np_A;
            j_c = (nqq_A-nq_A);
        }
        MPI_Bcast( &i_c, 1, MPI_INT, ftw->errproc, MPI_COMM_WORLD );
        MPI_Bcast( &j_c, 1, MPI_INT, ftw->errproc, MPI_COMM_WORLD );

        for (i=0; i<i_c; i+=nb)
            for (j=0, k=0; j<j_c; j+=nb, k++)
            {
                if (iamvic)
                {
                    int i_1 = i+1;
                    g[0] = indxl2g_( &i_1, &nb, &myrow, &i0, &nprow );
                    int j_1 = nq_A + j + 1;
                    g[1] = indxl2g_( &j_1, &nb, &mycol, &i0, &npcol );

                    tt[1]=(((g[1]-NA-1)/nb)%2==0)?(g[1]+nb):(g[1]-nb);
                    tt[0]=g[0];
                    //printf ("(%d,%d): block (%d,%d) <------ (%d,%d)\n", myrow, mycol, g[0], g[1], tt[0],tt[1]);
                }
                MPI_Bcast( tt, 2, MPI_INT, ftw->errproc, MPI_COMM_WORLD );
                MPI_Bcast( g, 2, MPI_INT, ftw->errproc, MPI_COMM_WORLD );

                Cpzgemr2d (nb, nb, A, tt[0], tt[1], descA, A, g[0], g[1], descA, ctxt);
            }

        recover_local_fillin (A, MA, NA, descA, ftw );
    }
}

static void recover_local_fillin (complexdouble_t *A, int MA, int NA, int *descA, ftla_work_t *ftw)
{
    int j;
    int lrindx, lcindx, rsrc, csrc;
    complexdouble_t *Achk_old, *Achk_new;
    int descAchk[9], npp, nqq;

    int ctxt = descA[CTXT_];
    int nb = descA[NB_];
    int nprow, npcol, myrow, mycol, myrank, nranks;
    Cblacs_gridinfo( ctxt, &nprow, &npcol, &myrow, &mycol );
    Cblacs_pinfo( &myrank, &nranks );
    int iamvic = (myrank==ftw->errproc);
    
    // determine local data size (non-checksum)
    int np_A = numroc_( &MA, &nb, &myrow, &i0, &nprow );

    //int NAA = NA+nb*2;
    int nq_A = numroc_( &NA, &nb, &mycol, &i0, &npcol );
    //int nqq_A = numroc_( &NAA, &nb, &mycol, &i0, &npcol );

    create_matrix( ctxt, &Achk_old, descAchk, MA, nb*2, nb, &npp, &nqq);
    create_matrix( ctxt, &Achk_new, descAchk, MA, nb*2, nb, &npp, &nqq);

    // copy the old checksum into buffer
    Cpzgemr2d (MA, nb*2, A, 1, NA+1, descA, Achk_old, 1, 1, descAchk, ctxt);

    // re-calculate the checksum
    hor_checkpointing_local (A, descA, 1, 1, 1, ftw );

    // copy the new checksum into buffer
    Cpzgemr2d (MA, nb*2, A, 1, NA+1, descA, Achk_new, 1, 1, descAchk, ctxt);

    int N_M = nb*2, i1=1;
    complexdouble_t done = 1, mone = -1;
    pzmatadd_ ( &MA, &N_M, &done, Achk_old, &i1, &i1, descAchk, &mone, Achk_new, &i1, &i1, descAchk);
    Cpzgemr2d (MA, nb*2, Achk_new, 1, 1, descAchk, A, 1, NA+1, descA, ctxt);

    // copy resut back in place
    int j_c;
    if (iamvic)
        j_c = nq_A;
    MPI_Bcast( &j_c, 1, MPI_INT, ftw->errproc, MPI_COMM_WORLD );

    for (j=0; j<j_c; j+=nb)
    {
        int gi, gj;
        if (iamvic)
        {
            // source matrix column (global)
            gj = NA+1;

            // source matrix row (global)
            gi = indxl2g_( &i1, &nb, &myrow, &i0, &nprow );
        }
        MPI_Bcast( &gi, 1, MPI_INT, ftw->errproc, MPI_COMM_WORLD );
        MPI_Bcast( &gj, 1, MPI_INT, ftw->errproc, MPI_COMM_WORLD );
            
        infog2l_ (&gi, &gj, descA, &nprow, &npcol, &myrow, &mycol, &lrindx, &lcindx, &rsrc, &csrc );
        //printf ("(%d,%d) <--- (%d,%d)\n", rsrc, ftw->errproc%npcol, rsrc, csrc);

        // copy the new checksum into buffer
        if (csrc==(ftw->errproc%npcol))
        {
            // local copy
            if (iamvic)
            {
                //printf ("(%d,%d) local copy when j=%d\n", myrow, mycol, j);
                memcpy (A+j*np_A, A+(lcindx-1)*np_A, np_A*nb*sizeof(complexdouble_t));
            }
        }
        else
        {
            // MPI copy
            if (myrow == rsrc)
            {
                MPI_Status status;
                if (iamvic)
                {
                    //printf ("(%d,%d) reciving from (%d,%d) for MPI_COPY when j=%d\n", myrow, mycol, myrow, csrc, j);
                    MPI_Recv(A+j*np_A, np_A*nb, MPI_DOUBLE_COMPLEX, csrc, 1, ftw->comms.rowcomm, &status);
                }
                else if (mycol==csrc)
                {
                    //printf ("(%d,%d) sending to (%d,%d) for MPI_COPY when j=%d\n", myrow, mycol, myrow, ftw->errproc%npcol, j);
                    MPI_Send(A+(lcindx-1)*np_A, np_A*nb, MPI_DOUBLE_COMPLEX, ftw->errproc%npcol, 1, ftw->comms.rowcomm);
                }
            }
        }
    }

    // copy the correct checksum back //
    Cpzgemr2d (MA, nb*2, Achk_old, 1, 1, descAchk, A, 1, NA+1, descA, ctxt);

    if (npp*nqq>0)
    {
        free (Achk_new);
        free (Achk_old);
    }
}


/*********************************************************************
 **** LOCAL COPY OF THE PANEL
 */

/**
 * create a copy of the local panel in the current Q-wide section
 */
void Cftla_pzqplcpy( complexdouble_t *A, int *descA, int IA, int JA, ftla_work_t *ftw ) {
#ifndef NO_QLOCALCOPY
    int i;
    complexdouble_t *Aoff;
    int ctxt = descA[CTXT_];
    int nb = descA[NB_];
    int nprow, npcol, myrow, mycol;
    Cblacs_gridinfo( ctxt, &nprow, &npcol, &myrow, &mycol );
    int MA = descA[2]; 
    int NA = descA[3]-ftw->Nc; 
    int np_A = numroc_( &MA, &nb, &myrow, &i0, &nprow );
    int nq_A = numroc_( &NA, &nb, &mycol, &i0, &npcol );
                        
    // determine number of rows to checkpoint locally 
    int lrindx, lcindx, rsrc, csrc;
    infog2l_ (&IA, &JA, descA, &nprow, &npcol, &myrow, &mycol, &lrindx, &lcindx, &rsrc, &csrc);

    int len = (np_A-lrindx+1)*sizeof(complexdouble_t);
    if (len>0 && lcindx<nq_A && lrindx<np_A)
    {
        // make a panel_copy if demanded
//        printf ("(%d,%d) copying %d row %d col at (%d,%d)\n", myrow, mycol, (np_A-lrindx+1), nb, IA, JA);
        Aoff = A+(lcindx-1)*np_A+(lrindx-1);
        for (i=0; i<nb; i++)
            memcpy (((typeof(A))ftw->pcopy.Pc)+i*np_A, Aoff+i*np_A, len);
    }

    // checkpoint local copy
    hor_checkpointing_local (ftw->pcopy.Pc, ftw->pcopy.descPc, 1, 1, 1, ftw );
#endif
}
extern void ftla_pzqplcpy_(complexdouble_t *A, int *descA, int *IA, int *JA, int *ftw ) {
    Cftla_pzqplcpy( A, descA, *IA, *JA, (ftla_work_t*)ftw );
}

/**
 * Reset the local copy of the panel in the current Q-wide section
 */
void Cftla_pzqplrst( complexdouble_t *A, int *descA, int IA, int JA, ftla_work_t *ftw )
{
#ifndef NO_QLOCALCOPY
    int i;
    complexdouble_t *Aoff;
    int ctxt = descA[CTXT_];
    int nb = descA[NB_];
    int nprow, npcol, myrow, mycol;
    Cblacs_gridinfo( ctxt, &nprow, &npcol, &myrow, &mycol );
    int MA = descA[2];
    int NA = descA[3]-ftw->Nc; 
    int np_A = numroc_( &MA, &nb, &myrow, &i0, &nprow );
    int nq_A = numroc_( &NA, &nb, &mycol, &i0, &npcol );
    
    // determine number of rows to checkpoint locally 
    int lrindx, lcindx, rsrc, csrc;
    infog2l_ (&IA, &JA, descA, &nprow, &npcol, &myrow, &mycol, &lrindx, &lcindx, &rsrc, &csrc);

    int len = (np_A-lrindx+1)*sizeof(complexdouble_t);
    if (len>0 && lcindx<nq_A && lrindx<np_A)
    {
        // make a panel_copy if demanded
        //printf ("(%d,%d) copying %d row %d col at (%d,%d)\n", myrow, mycol, (np_A-lrindx+1), nb, IA, JA);
        Aoff = A+(lcindx-1)*np_A+(lrindx-1);
        for (i=0; i<nb; i++)
            memcpy (Aoff+i*np_A, ((typeof(A))ftw->pcopy.Pc)+i*np_A, len);
    }
#endif
}
extern void ftla_pzqplrst_(complexdouble_t *A, int *descA, int *IA, int *JA, int *ftw ) {
    Cftla_pzqplrst( A, descA, *IA, *JA, (ftla_work_t*)ftw );
}


/*
 * perform checkpointing for the local copy
 * assuming JA is on process column 0
 */
static void hor_checkpointing_local( complexdouble_t *A, int *descA, int IA, int JA, int iscopy, ftla_work_t *ftw )
{
    int i, k;
    complexdouble_t *Aoff;
    int ctxt = descA[CTXT_];
    int nb = descA[NB_];
    int nprow, npcol, myrow, mycol;
    Cblacs_gridinfo( ctxt, &nprow, &npcol, &myrow, &mycol );
    int MA = descA[2]; 
    int np_A = ftw->pcopy.nrPc; 
    
    //int NA = descA[3]-2*nb; 
    //int nq_A = numroc_( &NA, &nb, &mycol, &i0, &npcol );
    
    // temp buffers
    complexdouble_t *S = ftw->comms.S;
    complexdouble_t *R = ftw->comms.R;

    // determine number of rows to checkpoint locally 
    int lrindx, lcindx, rsrc, csrc;
    infog2l_ (&IA, &JA, descA, &nprow, &npcol, &myrow, &mycol, &lrindx, &lcindx, &rsrc, &csrc);
    int rows = np_A-lrindx+1;

    // determine target locatoin
    int tc = nb*npcol+1;
    //printf ("(%d,%d) JA=%d, tc=%d\n", myrow, mycol, JA, tc);    

    // find root
    int root, rootcr, rootcc;
    infog2l_ (&IA, &tc, descA, &nprow, &npcol, &myrow, &mycol, &rootcr, &rootcc, &rsrc, &csrc);
    root = (mycol == csrc);
//    if (root)
//        printf ("(%d,%d) is root (csrc=%d) when JA=%d\n", myrow, mycol, csrc, JA);

    // checkpointing
    for (i=0; i<rows; i+=nb)
    {
        // determine participants 
        int in = ((lrindx+i)<np_A);
        if (in)
        {
            Aoff = A+(lcindx-1)*np_A+(lrindx-1)+i;
            for (k=0; k<nb; k++)
                memcpy (S+k*nb, Aoff+k*np_A, nb*sizeof(complexdouble_t));
        }
        else
            memset (S, 0, nb*nb*sizeof(complexdouble_t));

        if ((lrindx+i)<np_A)
        {
            MPI_Reduce (S, R, nb*nb, MPI_DOUBLE_COMPLEX, MPI_SUM, csrc, ftw->comms.rowcomm);

            if (root)
            {
                Aoff = A+(rootcc-1)*np_A+(rootcr-1)+i;  
                for (k=0; k<nb; k++)
                    memcpy (Aoff+k*np_A, R+k*nb, nb*sizeof(complexdouble_t));
            }
        }
    }
        
    if (iscopy)
    {
        //make a copy
        Cpzgemr2d (MA-IA+1, nb, A, IA, tc, descA, A, IA, tc+nb, descA, ctxt);
    }
}


/*********************************************************************
 ***** CHECKPOINT / CHECKSUM
 */
#define COPY 1
#define NOCOPY 0

/*
 * Row-wise or column-wise checksum of A into Acd
 */
void Cftla_pzcsum (char side,
        complexdouble_t *A, int MA, int NA, int IA, int JA, int *descA, 
        complexdouble_t *Acd, int IAcd, int JAcd, int *descAcd,
        ftla_work_t *ftw)
{
#ifndef NO_CHECKPOINT
    // quick return if possible
    if (side != 'r' && side != 'b' && side != 'd') {
        pxerbla_(&descA[CTXT_], "PDCSUM", &i1, 6);
        return;
    }

    int i;
    int ctxt = descA[CTXT_];
    int nb = descA[NB_];
    int nprow, npcol, myrow, mycol;
    Cblacs_gridinfo( ctxt, &nprow, &npcol, &myrow, &mycol );

    if (side=='d' || side=='b') {
        hor_checkpointing (A, descA, IA, JA, COPY, ftw );
    }
    else if (side=='r' || side=='b') {
        // determine local data size (non-checksum)
        int MAA = MA+ftw->Mc;
        int np_A = numroc_( &MAA, &nb, &myrow, &i0, &nprow );
        //int npp_A = numroc_( &MA, &nb, &myrow, &i0, &nprow ); // no checksum np_A 

        int NAA = NA+ftw->Nc;
        int nq_A = numroc_( &NAA, &nb, &mycol, &i0, &npcol );
        int nqq_A = numroc_( &NA, &nb, &mycol, &i0, &npcol );

        // zero out the initial checksum
        if (nq_A > nqq_A)
            memset (A+np_A*nqq_A, 0, (nq_A-nqq_A)*np_A*sizeof(complexdouble_t));

        // checkpointing
        int Qnb=npcol*nb;
        for (i=0; i<NA; i+=Qnb)
        {
            hor_checkpointing (A, descA, IA, JA+i, COPY, ftw );
        }
    }
#endif
}
extern void ftla_pzcsum_( char *side, complexdouble_t *A, int *MA, int *NA, int *IA, int *JA, int *descA, complexdouble_t *Acd, int *IAcd, int *JAcd, int *descAcd, int *ftw ) {
    Cftla_pzcsum( side[1], A, *MA, *NA, *IA, *JA, descA, Acd, *IAcd, *JAcd, descAcd, (ftla_work_t*)ftw );
}

/*
 * perform one step of checkpointing
 * assuming JA is on process column 0
 */
static void hor_checkpointing (complexdouble_t *A, int *descA, int IA, int JA, int iscopy, ftla_work_t *ftw)
{
    int i, k;
    complexdouble_t *Aoff;
    int ctxt = descA[CTXT_];
    int nb = descA[NB_];
    int nprow, npcol, myrow, mycol;
    Cblacs_gridinfo( ctxt, &nprow, &npcol, &myrow, &mycol );
    
    if( 0 == ftw->Nc ) return;
    int MA = descA[2]; 
    int NA = descA[3]-ftw->Nc; 
    int np_A = numroc_( &MA, &nb, &myrow, &i0, &nprow );
    int nq_A = numroc_( &NA, &nb, &mycol, &i0, &npcol );

    // temp buffers
    complexdouble_t *S = ftw->comms.S;
    complexdouble_t *R = ftw->comms.R;

    // determine number of rows to checkpoint locally 
    int lrindx, lcindx, rsrc, csrc;
    infog2l_ (&IA, &JA, descA, &nprow, &npcol, &myrow, &mycol, &lrindx, &lcindx, &rsrc, &csrc);
    int rows = np_A-lrindx+1;

    /*
    // make a panel_copy if demanded
    if (iscopy && (lcindx<nq_A) && (lrindx+i)<np_A)
    {
        printf ("(%d,%d) copying %d row %d col at (%d,%d)\n", myrow, mycol, (np_A-lrindx), nb, IA, JA);
        int len = (np_A-lrindx)*sizeof(complexdouble_t);
        Aoff = A+(lcindx-1)*np_A+(lrindx-1);
        for (i=0; i<nb; i++)
            memcpy (((typeof(A))ftw->pcopy.Pc)+i*np_A, Aoff+i*np_A, len);
    }
    */

    // determine target location
    int tc = ((NA-JA+1)/(nb*npcol)+((NA-JA+1)%(nb*npcol)!=0)-1)*2*nb+NA+1;
    //printf ("(%d,%d) JA=%d, tc=%d\n", myrow, mycol, JA, tc);    

    // find root
    int root, rootcr, rootcc;
    infog2l_ (&IA, &tc, descA, &nprow, &npcol, &myrow, &mycol, &rootcr, &rootcc, &rsrc, &csrc);
    root = (mycol == csrc);
    //if (root)
    //    printf ("(%d,%d) is root (csrc=%d) when JA=%d\n", myrow, mycol, csrc, JA);

    // checkpointing
    for (i=0; i<rows; i+=nb)
    {
        // determine participants 
        int in = ((lcindx<nq_A) && (lrindx+i)<np_A);
        //printf ("(%d,%d) IA=%d, JA=%d, lrindx=%d, lcindx=%d, rows=%d, i=%d, IN=%d\n", myrow, mycol, IA, JA, lrindx, lcindx, rows, i, in);    
        if (in)
        {
            Aoff = A+(lcindx-1)*np_A+(lrindx-1)+i;
            for (k=0; k<nb; k++)
                memcpy (S+k*nb, Aoff+k*np_A, nb*sizeof(complexdouble_t));
        }
        else
            memset (S, 0, nb*nb*sizeof(complexdouble_t));

        if ((lrindx+i)<np_A)
        {
            MPI_Reduce (S, R, nb*nb, MPI_DOUBLE_COMPLEX, MPI_SUM, csrc, ftw->comms.rowcomm);

            if (root)
            {
                Aoff = A+(rootcc-1)*np_A+(rootcr-1)+i;  
                for (k=0; k<nb; k++)
                    memcpy (Aoff+k*np_A, R+k*nb, nb*sizeof(complexdouble_t));
            }
        }
    }
        
    if (iscopy) {
        //make a copy
        Cpzgemr2d (MA-IA+1, nb, A, IA, tc, descA, A, IA, tc+nb, descA, ctxt);
    }
}

/*
 * verify checkpointing after computing
 */
double Cftla_pzcsumv( int M, int N, complexdouble_t *A, int *descA, ftla_work_t *ftw) {
    complexdouble_t *Achk_old, *Achk_new;
    int descAchk[9], npp, nqq;
    int nb = descA[NB_];
    int ctxt = descA[CTXT_];
    int nprow, npcol, myrow, mycol;
    Cblacs_gridinfo( ctxt, &nprow, &npcol, &myrow, &mycol );

    create_matrix( ctxt, &Achk_old, descAchk, M, N-M, nb, &npp, &nqq);
    create_matrix( ctxt, &Achk_new, descAchk, M, N-M, nb, &npp, &nqq);

    // copy the old checksum into buffer
    Cpzgemr2d (M, N-M, A, 1, M+1, descA, Achk_old, 1, 1, descAchk, ctxt);

    // re-calculate the checksum
    Cftla_pzcsum ('r', A, M, M, 1, 1, descA, A, 1, M, descA, ftw );

    // copy the new checksum into buffer
    Cpzgemr2d (M, N-M, A, 1, M+1, descA, Achk_new, 1, 1, descAchk, ctxt);

    // cross check
    int N_M = N-M, i1=1;
    complexdouble_t done = 1, mone = -1;
    pzmatadd_ ( &M, &N_M, &done, Achk_new, &i1, &i1, descAchk, &mone, Achk_old, &i1, &i1, descAchk);
    double resid = pzlange_("F", &M, &N_M, Achk_old, &i1, &i1, descAchk, NULL)/pzlange_("F", &M, &N_M, Achk_new, &i1, &i1, descAchk, NULL)/M;

    if (npp*nqq>0) {
        free (Achk_new);
        free (Achk_old);
    }
    return resid;
}
extern double ftla_pzcsumv_( int *M, int *N, complexdouble_t *A, int *descA, int *ftw ) {
    return Cftla_pzcsumv( *M, *N, A, descA, (ftla_work_t*)ftw );
}


static void create_matrix( int ctxt, complexdouble_t **A, int *descA, 
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

    /* set np and nq */
    if (np_A != NULL)
        *np_A = np_iA;
    if (nq_A != NULL)
        *nq_A = nq_iA;
}
