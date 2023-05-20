/*
 * test_ScaLAPACK_check_ft.inc
 *
 *  Created on: Nov 13, 2022
 *      Author: marcello
 */


		if (env.mpi_rank==0)
		{
			if (input.scalapack_bf < 64)
			{
				if (verbosity>0) DISPLAY_WRN(label,"blocking factor < 64");
			}

			if (checkpoint_freq <= 0)
			{
				if (verbosity>0) DISPLAY_WRN(label,"no checkpoint will be done: never recovering");
			}

			if (fault_tolerance < 1)
			{
				if (verbosity>0) DISPLAY_WRN(label,"fault tolerance disabled: never recovering");
			}
			else
			{
				if (faulty_procs == 0)
				{
					if (verbosity>0) DISPLAY_WRN(label,"no fault will be actually injected: never failing");
				}
			}

			if (fault_tolerance < 1 && env.spare_procs < 1)
			{
				DISPLAY_ERR(label,"spare processes must be allocated even if fault tolerance is disabled");
			}
			else if (fault_tolerance > 0 && env.spare_procs < fault_tolerance)
			{
				DISPLAY_ERR(label,"fault tolerance set but not enough spare processes for the requested fault tolerance level");
			}
			else if (fault_tolerance > 0 && env.spare_procs > fault_tolerance )
			{
				DISPLAY_ERR(label,"fault tolerance set but too many spare processes for the requested fault tolerance level");
			}
			else if (fault_tolerance > 0 && faulty_procs > fault_tolerance)
			{
				DISPLAY_ERR(label,"fault tolerance set but requested fault occurrences greater than fault tolerance level");
			}
			else if (input.scalapack_bf < 1)
			{
				DISPLAY_ERR(label,"the blocking factor has to be greater than 0");
			}
			else
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
							DISPLAY_ERR(label,"the number of columns (rows) per calc. process has to be a multiple of the blocking factor");
						}
					}
					else
					{
						DISPLAY_ERR(label,"the matrix size has to be a multiple of the calc. processes per rank");
					}
				}
				else
				{
					DISPLAY_ERR(label,"the number of the calc. processes must be square");
				}
			}
		}
