/*
 * test_ScaLAPACK_show-factors.inc
 *
 *  Created on: Nov 13, 2022
 *      Author: marcello
 */
			if (verbosity>1)
			{
				printf("\nThe %s factorization is:\n",label);
				PrintMatrix1D_double(A, input.n, input.n);
				printf("\n with exit code     %d\n",output.exit_code);
				printf("      norm.rel.err. %.17f\n",rank_result.norm_rel_err);
			}
