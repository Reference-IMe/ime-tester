/*
 * test_ScaLAPACK_check_ft.inc
 *
 *  Created on: Nov 13, 2022
 *      Author: marcello
 */


		if (env.mpi_rank==0)
		{
			if (checkpoint_freq < 0)
			{
				if (verbosity>0) DISPLAY_WRN(label,"No checkpoint will be done: (never recovering)");
			}
			if (faulty_procs == 0)
			{
				if (verbosity>0) DISPLAY_WRN(label,"No fault will be actually injected: (never failing)");
			}
			if (faulty_procs > fault_tolerance)
			{
				DISPLAY_ERR(label,"Requested fault occurrences greater than fault tolerance (single fault routine)");
			}
			else if (fault_tolerance < 1)
			{
				DISPLAY_ERR(label,"Fault tolerance disabled ('-ft 0')");
			}
			else if (fault_tolerance > 1)
			{
				DISPLAY_ERR(label,"Requested fault tolerance too high (single fault routine)");
			}
			else
			{
				if (input.scalapack_bf < 64)
				{
					if (verbosity>0) DISPLAY_WRN(label,"Blocking factor < 64");
				}
				if (env.spare_procs > 1)
				{
					DISPLAY_ERR(label,"Too many spare processes allocated (single fault routine)");
				}
				else if (env.spare_procs > 0)
				{
					if (IS_SQUARE(env.calc_procs))
					{
						if (IS_MULT(input.n, rank_calc_procs))
						{
							if (IS_MULT(input.n / rank_calc_procs, input.scalapack_bf))
							{
								if (verbosity>0) DISPLAY_MSG(label,"OK");
								output.exit_code = 0;
							}
							else
							{
								DISPLAY_ERR(label,"The number of columns (rows) per calc. process has to be a multiple of the blocking factor");
							}
						}
						else
						{
							DISPLAY_ERR(label,"The matrix size has to be a multiple of the calc. processes per rank");
						}
					}
					else
					{
						DISPLAY_ERR(label,"The number of the calc. processes must be square");
					}
				}
				else
				{
					DISPLAY_ERR(label,"Fault tolerance enabled, but no requested spare processes");
				}
			}
		}
