#ifndef __INCLUDE_FTLA_CSUM_H__
#define __INCLUDE_FTLA_CSUM_H__
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

#include "ftla_ftwork.h"
#include <complex.h>

extern void Cftla_pscsum( char side, float *A, int MA, int NA, int IA, int JA, int *descA, float *Acd, int IAcd, int JAcd, int *descAcd, ftla_work_t *ftwork );
extern void ftla_pscsum_( char *side, float *A, int *MA, int *NA, int *IA, int *JA, int *descA, float *Acd, int *IAcd, int *JAcd, int *descAcd, int *ftwork );
extern void Cftla_pscsumr( char zone, float *A, int m, int n, int *descA, ftla_work_t *ftwork );
extern void ftla_pscsumr_( char *zone, float *A, int *m, int *n, int *descA, int *ftwork );
extern float Cftla_pscsumv( int M, int N, float *A, int *descA, ftla_work_t *ftwork );
extern float ftla_pscsumv_( int *M, int *N, float *A, int *descA, int *ftwork );

extern void Cftla_psqplcpy(float *A, int *descA, int IA, int JA, ftla_work_t *ftwork );
extern void ftla_psqplcpy_(float *A, int *descA, int *IA, int *JA, int *ftwork );
extern void Cftla_psqplrst(float *A, int *descA, int IA, int JA, ftla_work_t *ftwork );
extern void ftla_psqplrst_(float *A, int *descA, int *IA, int *JA, int *ftwork );


extern void Cftla_pdcsum( char side, double *A, int MA, int NA, int IA, int JA, int *descA, double *Acd, int IAcd, int JAcd, int *descAcd, ftla_work_t *ftwork );
extern void ftla_pdcsum_( char *side, double *A, int *MA, int *NA, int *IA, int *JA, int *descA, double *Acd, int *IAcd, int *JAcd, int *descAcd, int *ftwork );
extern void Cftla_pdcsumr( char zone, double *A, int m, int n, int *descA, ftla_work_t *ftwork );
extern void ftla_pdcsumr_( char *zone, double *A, int *m, int *n, int *descA, int *ftwork );
extern double Cftla_pdcsumv( int M, int N, double *A, int *descA, ftla_work_t *ftwork );
extern double ftla_pdcsumv_( int *M, int *N, double *A, int *descA, int *ftwork );

extern void Cftla_pdqplcpy(double *A, int *descA, int IA, int JA, ftla_work_t *ftwork );
extern void ftla_pdqplcpy_(double *A, int *descA, int *IA, int *JA, int *ftwork );
extern void Cftla_pdqplrst(double *A, int *descA, int IA, int JA, ftla_work_t *ftwork );
extern void ftla_pdqplrst_(double *A, int *descA, int *IA, int *JA, int *ftwork );


extern void Cftla_pccsum( char side, float complex *A, int MA, int NA, int IA, int JA, int *descA, float complex *Acd, int IAcd, int JAcd, int *descAcd, ftla_work_t *ftwork );
extern void ftla_pccsum_( char *side, float complex *A, int *MA, int *NA, int *IA, int *JA, int *descA, float complex *Acd, int *IAcd, int *JAcd, int *descAcd, int *ftwork );
extern void Cftla_pccsumr( char zone, float complex *A, int m, int n, int *descA, ftla_work_t *ftwork );
extern void ftla_pccsumr_( char *zone, float complex *A, int *m, int *n, int *descA, int *ftwork );
extern float Cftla_pccsumv( int M, int N, float complex *A, int *descA, ftla_work_t *ftwork );
extern float ftla_pccsumv_( int *M, int *N, float complex *A, int *descA, int *ftwork );

extern void Cftla_pcqplcpy(float complex *A, int *descA, int IA, int JA, ftla_work_t *ftwork );
extern void ftla_pcqplcpy_(float complex *A, int *descA, int *IA, int *JA, int *ftwork );
extern void Cftla_pcqplrst(float complex *A, int *descA, int IA, int JA, ftla_work_t *ftwork );
extern void ftla_pcqplrst_(float complex *A, int *descA, int *IA, int *JA, int *ftwork );


extern void Cftla_pzcsum( char side, double complex *A, int MA, int NA, int IA, int JA, int *descA, double complex *Acd, int IAcd, int JAcd, int *descAcd, ftla_work_t *ftwork );
extern void ftla_pzcsum_( char *side, double complex *A, int *MA, int *NA, int *IA, int *JA, int *descA, double complex *Acd, int *IAcd, int *JAcd, int *descAcd, int *ftwork );
extern void Cftla_pzcsumr( char zone, double complex *A, int m, int n, int *descA, ftla_work_t *ftwork );
extern void ftla_pzcsumr_( char *zone, double complex *A, int *m, int *n, int *descA, int *ftwork );
extern double Cftla_pzcsumv( int M, int N, double complex *A, int *descA, ftla_work_t *ftwork );
extern double ftla_pzcsumv_( int *M, int *N, double complex *A, int *descA, int *ftwork );

extern void Cftla_pzqplcpy(double complex *A, int *descA, int IA, int JA, ftla_work_t *ftwork );
extern void ftla_pzqplcpy_(double complex *A, int *descA, int *IA, int *JA, int *ftwork );
extern void Cftla_pzqplrst(double complex *A, int *descA, int IA, int JA, ftla_work_t *ftwork );
extern void ftla_pzqplrst_(double complex *A, int *descA, int *IA, int *JA, int *ftwork );

#endif

