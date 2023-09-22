
#if TYPE == REAL_DOUBLE
	#define PRECISION double
#endif

#if TYPE == REAL_SINGLE
	#define PRECISION float
#endif

/*
 * test_ScaLAPACK_show-solution.inc
 *
 *  Created on: Nov 13, 2022
 *      Author: marcello
 */

			if (verbosity>1)
			{
				printf("\nThe %s solution is:\n",label);
				FUNCNAME(PrintMatrix1D_, PRECISION)(bb, input.n, input.nrhs);
				printf("\n with exit code     %d\n",output.exit_code);
				printf("      norm.rel.err. %.17f\n",rank_result.norm_rel_err);
			}

#undef PRECISION
