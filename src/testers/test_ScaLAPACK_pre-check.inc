/*
 * test_ScaLAPACK_check.inc
 *
 *  Created on: Nov 13, 2022
 *      Author: marcello
 */


		if (env.mpi_rank==0)
		{
			if (verbosity>0)
			{
				if (env.spare_procs > 0)
				{
					DISPLAY_WRN(label,"Can run also with fault tolerance enabled, but calc. processes differ from total processes")
				}
				if (input.scalapack_bf < 64)
				{
					DISPLAY_WRN(label,"Blocking factor < 64")
				}
				DISPLAY_MSG(label,"OK");
			}
			output.exit_code = 0;
		}
