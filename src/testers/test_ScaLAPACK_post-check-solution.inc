
#if TYPE == REAL_DOUBLE
	#define PRECISION double
#endif

#if TYPE == REAL_SINGLE
	#define PRECISION float
#endif

/*
 * test_ScaLAPACK_post-check.h
 *
 *  Created on: Nov 13, 2022
 *      Author: marcello
 */

			// check exit condition
			if (output.exit_code!=0)
			{
				if (verbosity>0) printf("\n** Dangerous exit code.. (%d)**\n",output.exit_code);
			}
			// calc error
			if (input.calc_nre) rank_result.norm_rel_err = FUNCNAME(NormwiseRelativeError1D_, PRECISION)(bb, xx_ref, input.n, input.nrhs);

#undef PRECISION
