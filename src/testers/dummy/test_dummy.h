#include <mpi.h>
#include <time.h>
#include <unistd.h>
#include "../../tester_structures.h"
#include "../../helpers/macros.h"
#include "../../helpers/matrix_advanced.h"
#include "../../helpers/simple_dynamic_strings/sds.h"
#include "../IMe/lib/src/helpers/matrix_basic.h"

test_result test_dummy (	const char check,
							const char* tag,
							int verbosity,
							parallel_env env,
							test_input input,
							int fault_tolerance	)
{
	test_result rank_result = TEST_NOT_RUN;
	test_result team_result = TEST_NOT_RUN;
	exit_status output      = EMPTY_OUTPUT;

	sds label=sdsempty();
	TAG2LABEL(tag,label);

	MPI_Comm comm_calc;

	int i_calc; // participating in ime calc = 1, checksumming = 0

	if (check)
	{
		if (env.mpi_rank==0)
		{
			if (verbosity>0) DISPLAY_MSG(label,"OK");
		}
		output.exit_code = 0;
	}
	else
	{
		if (env.mpi_rank >= env.calc_procs)
		{
			i_calc=0;
			MPI_Comm_split(MPI_COMM_WORLD, MPI_UNDEFINED, MPI_UNDEFINED, &comm_calc); // checksumming procs don't belong to calc communicator
		}
		else
		{
			i_calc=1;
			MPI_Comm_split(MPI_COMM_WORLD, i_calc, env.mpi_rank, &comm_calc); // calc procs belong to calc communicator
		}

		if (i_calc)
		{
			output.exit_code=0;

			if (env.mpi_rank==0)
			{
				// calc error
				if (input.calc_nre) rank_result.norm_rel_err = 0;

				if (verbosity>1)
				{
					printf("\nThe %s solution is void\n",tag);
					printf("\n with exit code     %d\n",output.exit_code);
					printf("      norm.rel.err. %.17f\n",rank_result.norm_rel_err);
				}
			}

		}
		else
		{
			rank_result.total_time=0;
			rank_result.core_time=0;
		}

		if (env.spare_procs>0)
		{
			if (comm_calc != MPI_COMM_NULL)
			{
				MPI_Comm_free(&comm_calc);
			}
		}
	}
	sdsfree(label);
	TEST_END(output, rank_result, team_result);
}
