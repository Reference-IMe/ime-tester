/*
 * create_matrix.h
 *
 *  Created on: Dec 29, 2019
 *      Author: marcello
 */

#include "../../helpers/Cblacs.h"
#include "../../helpers/scalapack.h"

#ifndef __FTLA_CREATE_MATRIX_H__
#define __FTLA_CREATE_MATRIX_H__

// code extracted from ftdqr_main.c of FTLA

/*
 * produce distributed matrix,  
 */
void create_matrix( int ctxt, int seed, double **A, int *descA, int M, int N, int nb, int *np_A, int *nq_A)
{
	int i0=0;
    int info;
    int nprow, npcol, myrow, mycol;
    Cblacs_gridinfo( ctxt, &nprow, &npcol, &myrow, &mycol );

    // allocate the generator matrix and check matrix
    int np_iA = numroc_( &M, &nb, &myrow, &i0, &nprow );
    int nq_iA = numroc_( &N, &nb, &mycol, &i0, &npcol );

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
        descinit_( descA, &M, &N, &nb, &nb, &i0, &i0, &ctxt, &itemp, &info );
        if (info != 0) Cblacs_abort( ctxt, 12 );
    }

    if (seed)
    {
        // fill in random numbers
        pdmatgen_ (&ctxt, "N", "N", &M, &N, &nb, &nb, *A,
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

#endif
