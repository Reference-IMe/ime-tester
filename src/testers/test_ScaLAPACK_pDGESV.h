/*
 * test_ScaLAPACK_pDGESV.h
 *
 *  Created on: Dec 27, 2019
 *      Author: marcello
 */

#include <mpi.h>
#include <time.h>
#include "../helpers/info.h"
#include "../helpers/macros.h"
#include "../helpers/matrix.h"
#include "../helpers/matrix_advanced.h"
#include "ScaLAPACK/ScaLAPACK_pDGESV.h"

duration_t test_ScaLAPACK_pDGESV(const char* label, int verbosity, int n, double* A_ref, double* x_ref, double* b_ref, int nrhs, int nb, int rank, int cprocs)
{
	duration_t timing, timing_max;
	result_info info;

	int i,j;

	double* A;
	double* bb;
	double* xx_ref;

	if (rank==0)
	{
		A=AllocateMatrix1D(n, n);
		bb=AllocateMatrix1D(n, nrhs);
		xx_ref=AllocateMatrix1D(n, nrhs);

		CopyMatrix1D(A_ref, A, n, n);
		for (i=0;i<n;i++)
		{
			for (j=0;j<nrhs;j++)
			{
				bb[i*nrhs+j] = b_ref[i];
				xx_ref[i*nrhs+j] = x_ref[i];
			}
		}


		if (verbosity>2)
		{
			printf("\n\n Matrix A:\n");
			PrintMatrix1D(A, n, n);
			printf("\n Vector b:\n");
			PrintMatrix1D(bb, n, nrhs);
		}
	}

	info = ScaLAPACK_pDGESV_calc(n, A, nrhs, bb, nb, rank, cprocs);

	if (rank==0)
	{
		// check exit condition
		if (info.exit_code!=0)
		{
			printf("\n** Dangerous exit code.. (%d)**\n",info.exit_code);
		}
		// calc error
		info.norm_rel_err = NormwiseRelativeError1D(bb, xx_ref, n, nrhs);

		if (verbosity>1)
		{
			printf("\nThe %s solution is:\n",label);
			PrintMatrix1D(bb, n, nrhs);
			printf("\n with exit code     %d\n",info.exit_code);
			printf("      norm.rel.err. %f\n",info.norm_rel_err);
		}
		DeallocateMatrix1D(A);
		DeallocateMatrix1D(bb);
		DeallocateMatrix1D(xx_ref);
	}

	TEST_END(info, timing, timing_max);
}
