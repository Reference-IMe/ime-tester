/*
 * Cblacs.h
 *
 *  Created on: Feb 3, 2019
 *      Author: marcello
 */

#ifndef __CBLACS_H__
#define __CBLACS_H__

extern void Cblacs_pinfo( int* mypnum, int* nprocs);
extern void Cblacs_get( int context, int request, int* value);
extern void Cblacs_gridinit( int* context, char* order, int np_row, int np_col);
extern void Cblacs_gridmap( int* context, int* usermap, int ldumap, int  np_row, int np_col);
extern void Cblacs_gridinfo( int context, int*  np_row, int* np_col, int* my_row, int* my_col);
extern void Cblacs_gridexit( int context);
extern void Cblacs_exit( int error_code);
extern void Cblacs_abort( int ictxt, int ErrNo);
extern void Cblacs_pcoord(int, int, int*, int*);
extern void Cblacs_barrier(int, const char*);

extern void Cdgerv2d(int, int, int, double*, int, int, int);
extern void Cdgesd2d(int, int, int, double*, int, int, int);

extern void Csgerv2d(int, int, int, float*, int, int, int);
extern void Csgesd2d(int, int, int, float*, int, int, int);

#endif /* __CBLACS_H__ */
