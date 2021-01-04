#include <mpi.h>
#include <time.h>
#include <unistd.h>
#include "../helpers/macros.h"
#include "../helpers/matrix.h"
#include "../helpers/matrix_advanced.h"
#include "tester_structures.h"

#include "../pviDGESV_WO.h"
#include "../pviDGESV_WO.ordered.allgather.early.h"
#include "../pviDGESV_WO.ordered.allgather.h"
#include "../pviDGESV_WO.ordered.gather.early.h"
#include "../pviDGESV_WO.ordered.gather.h"

#include "../pviDGESV_WO.unordered1.allgather.early.h"
#include "../pviDGESV_WO.unordered1.allgather.h"
#include "../pviDGESV_WO.unordered1.gather.early.h"
#include "../pviDGESV_WO.unordered1.gather.h"

#include "../pviDGESV_WO.unordered2.allgather.early.h"
#include "../pviDGESV_WO.unordered2.allgather.h"
#include "../pviDGESV_WO.unordered2.gather.early.h"
#include "../pviDGESV_WO.unordered2.gather.h"

#include "../pviDGESV_WO.unordered3.allgather.early.h"
#include "../pviDGESV_WO.unordered3.allgather.h"
#include "../pviDGESV_WO.unordered3.gather.early.h"
#include "../pviDGESV_WO.unordered3.gather.h"

#include "../pviDGESV_CO.gather.h"
#include "../pvDGESV_CO.h"
#include "../pvDGESV_CO.gather.2pass.h"
#include "../pvDGESV_CO.gather.ind.h"
#include "../pvDGESV_CO.allgather.small.h"
#include "../pvDGESV_CO.allgather.smaller.h"
#include "../pvDGESV_CO.allgather.smallest.h"
#include "../pvDGESV_CO.gather.smaller.h"
#include "../pvDGESV_CO.gather.smallest.h"


test_result test_IMe_pvDGESV(const char check, const char* label, const char* variant, int verbosity, test_input input, int rank)
{
	test_result rank_result = TEST_NOT_RUN;
	test_result team_result = TEST_NOT_RUN;
	test_output output      = EMPTY_OUTPUT;

	int i,j;

	double** A2;
	double*  A2_1D;
	double** bb;
	double** xx;
	double*  xx_ref;

	MPI_Comm comm_calc;

	int i_calc; // participating in ime calc = 1, checksumming = 0

	if (check)
	{
		if (rank==0)
		{
			if (input.ime_bf < 1)
			{
				DISPLAY_ERR(label,"the blocking factor has to be greater than 0");
			}
			else
			{
				if (IS_MULT(input.n, input.calc_procs))
				{
					if (input.n / input.calc_procs > 0)
					{
						if (input.spare_procs > 0)
						{
							DISPLAY_WRN(label,"can run also with FT enabled, but calc. processes differ from total processes")
						}
						if (IS_MULT(input.n / input.calc_procs, input.ime_bf))
						{
							DISPLAY_MSG(label,"OK");
							output.exit_code = 0;
						}
						else
						{
							DISPLAY_ERR(label,"the number of columns per calc. process has to be a multiple of the blocking factor");
						}
					}
					else DISPLAY_ERR(label,"the number of columns per calc. process has to be greater than 0");
				}
				else DISPLAY_ERR(label,"the number of columns has to be a multiple of the calc. processes");
			}
		}
	}
	else
	{
		if (rank >= input.calc_procs)
		{
			i_calc=0;
			MPI_Comm_split(MPI_COMM_WORLD, MPI_UNDEFINED, MPI_UNDEFINED, &comm_calc); // checksumming procs don't belong to calc communicator
		}
		else
		{
			i_calc=1;
			MPI_Comm_split(MPI_COMM_WORLD, i_calc, rank, &comm_calc); // calc procs belong to calc communicator
		}

		if (i_calc)
		{
			xx=AllocateMatrix2D(input.n, input.nrhs, CONTIGUOUS);
			bb=AllocateMatrix2D(input.n, input.nrhs, CONTIGUOUS);

			if (rank==0)
			{
				xx_ref=AllocateMatrix1D(input.n, input.nrhs);
				A2=AllocateMatrix2D(input.n, input.n, CONTIGUOUS);
				A2_1D=&A2[0][0];
				CopyMatrix1D(input.A_ref, A2_1D, input.n, input.n);

				for (i=0;i<input.n;i++)
				{
					for (j=0;j<input.nrhs;j++)
					{
						bb[i][j] = input.b_ref[i];
						xx_ref[i*input.nrhs+j] = input.x_ref[i];
					}
				}

				if (verbosity>2)
				{
					printf("\n\n Matrix A:\n");
					PrintMatrix2D(A2, input.n, input.n);
					printf("\n Vector b:\n");
					PrintMatrix2D(bb, input.n, input.nrhs);
				}
			}
			else
			{
				A2=AllocateMatrix2D(1, 1, CONTIGUOUS); // to avoid segmentation fault in mpi collectives with 2D arrays
				xx_ref=NULL;
			}
			     if ( strcmp( variant, "PV-WO"      )     == 0) output = pviDGESV_WO_default (input.ime_bf, input.n, A2, input.nrhs, bb, xx, comm_calc);
			else if ( strcmp( variant, "PV-WO-oae"  )     == 0) output = pviDGESV_WO_oae (input.ime_bf, input.n, A2, input.nrhs, bb, xx, comm_calc);
			else if ( strcmp( variant, "PV-WO-oa"   )     == 0) output = pviDGESV_WO_oa  (input.ime_bf, input.n, A2, input.nrhs, bb, xx, comm_calc);
			else if ( strcmp( variant, "PV-WO-oge"  )     == 0) output = pviDGESV_WO_oge (input.ime_bf, input.n, A2, input.nrhs, bb, xx, comm_calc);
			else if ( strcmp( variant, "PV-WO-og"   )     == 0) output = pviDGESV_WO_og  (input.ime_bf, input.n, A2, input.nrhs, bb, xx, comm_calc);

			else if ( strcmp( variant, "PV-WO-u1ae" )     == 0) output = pviDGESV_WO_u1ae (input.ime_bf, input.n, A2, input.nrhs, bb, xx, comm_calc);
			else if ( strcmp( variant, "PV-WO-u1a"  )     == 0) output = pviDGESV_WO_u1a  (input.ime_bf, input.n, A2, input.nrhs, bb, xx, comm_calc);
			else if ( strcmp( variant, "PV-WO-u1ge" )     == 0) output = pviDGESV_WO_u1ge (input.ime_bf, input.n, A2, input.nrhs, bb, xx, comm_calc);
			else if ( strcmp( variant, "PV-WO-u1g"  )     == 0) output = pviDGESV_WO_u1g  (input.ime_bf, input.n, A2, input.nrhs, bb, xx, comm_calc);

			else if ( strcmp( variant, "PV-WO-u2ae" )     == 0) output = pviDGESV_WO_u2ae (input.ime_bf, input.n, A2, input.nrhs, bb, xx, comm_calc);
			else if ( strcmp( variant, "PV-WO-u2a"  )     == 0) output = pviDGESV_WO_u2a  (input.ime_bf, input.n, A2, input.nrhs, bb, xx, comm_calc);
			else if ( strcmp( variant, "PV-WO-u2ge" )     == 0) output = pviDGESV_WO_u2ge (input.ime_bf, input.n, A2, input.nrhs, bb, xx, comm_calc);
			else if ( strcmp( variant, "PV-WO-u2g"  )     == 0) output = pviDGESV_WO_u2g  (input.ime_bf, input.n, A2, input.nrhs, bb, xx, comm_calc);

			else if ( strcmp( variant, "PV-WO-u3ae" )     == 0) output = pviDGESV_WO_u3ae (input.ime_bf, input.n, A2, input.nrhs, bb, xx, comm_calc);
			else if ( strcmp( variant, "PV-WO-u3a"  )     == 0) output = pviDGESV_WO_u3a  (input.ime_bf, input.n, A2, input.nrhs, bb, xx, comm_calc);
			else if ( strcmp( variant, "PV-WO-u3ge" )     == 0) output = pviDGESV_WO_u3ge (input.ime_bf, input.n, A2, input.nrhs, bb, xx, comm_calc);
			else if ( strcmp( variant, "PV-WO-u3g"  )     == 0) output = pviDGESV_WO_u3g  (input.ime_bf, input.n, A2, input.nrhs, bb, xx, comm_calc);

			else if ( strcmp( variant, "PV-iCO-g"        ) == 0) output = pviDGESV_CO_g (input.ime_bf, input.n, A2, input.nrhs, bb, xx, comm_calc);
			//TODO: create a separate test launcher for non-interleaved routines or rename this launcher
			else if ( strcmp( variant, "PV-CO"           ) == 0) output = pvDGESV_CO_default    (input.ime_bf, input.n, A2, input.nrhs, bb, xx, comm_calc);
			else if ( strcmp( variant, "PV-CO-g-ind"     ) == 0) output = pvDGESV_CO_g_ind      (input.ime_bf, input.n, A2, input.nrhs, bb, xx, comm_calc);
			else if ( strcmp( variant, "PV-CO-g-2pass"   ) == 0) output = pvDGESV_CO_g_2pass    (input.ime_bf, input.n, A2, input.nrhs, bb, xx, comm_calc);
			else if ( strcmp( variant, "PV-CO-g-smaller" ) == 0) output = pvDGESV_CO_g_smaller  (input.ime_bf, input.n, A2, input.nrhs, bb, xx, comm_calc);
			else if ( strcmp( variant, "PV-CO-g-smallest") == 0) output = pvDGESV_CO_g_smallest (input.ime_bf, input.n, A2, input.nrhs, bb, xx, comm_calc);
			else if ( strcmp( variant, "PV-CO-a-small"   ) == 0) output = pvDGESV_CO_a_small    (input.ime_bf, input.n, A2, input.nrhs, bb, xx, comm_calc);
			else if ( strcmp( variant, "PV-CO-a-smaller" ) == 0) output = pvDGESV_CO_a_smaller  (input.ime_bf, input.n, A2, input.nrhs, bb, xx, comm_calc);
			else if ( strcmp( variant, "PV-CO-a-smallest") == 0) output = pvDGESV_CO_a_smallest (input.ime_bf, input.n, A2, input.nrhs, bb, xx, comm_calc);
			else
			{
				DISPLAY_ERR(label,"not yet implemented! UNDEFINED BEHAVIOUR!");
			}

			if (rank==0)
			{
				// check exit condition
				if (output.exit_code!=0)
				{
					printf("\n** Dangerous exit code.. (%d)**\n",output.exit_code);
				}
				// calc error
				if (input.calc_nre) rank_result.norm_rel_err = NormwiseRelativeError1D(&xx[0][0], xx_ref, input.n, input.nrhs);

				if (verbosity>1)
				{
					printf("\nThe %s solution is:\n",label);
					PrintMatrix2D(xx, input.n, input.nrhs);
					printf("\n with exit code     %d\n",output.exit_code);
					printf("      norm.rel.err. %f\n",rank_result.norm_rel_err);
				}
			}
			//cleanup
			if (rank==0)
			{
				DeallocateMatrix2D(A2, input.n, CONTIGUOUS);
			}
			else
			{
				DeallocateMatrix2D(A2, 1, CONTIGUOUS);
			}
			DeallocateMatrix1D(xx_ref);
			DeallocateMatrix2D(xx, input.n, CONTIGUOUS);
			DeallocateMatrix2D(bb, input.n, CONTIGUOUS);

		}
		else
		{
			xx=NULL;
			bb=NULL;

			rank_result.total_time=0;
			rank_result.core_time=0;
		}

		if (input.spare_procs>0)
		{
			if (comm_calc != MPI_COMM_NULL)
			{
				MPI_Comm_free(&comm_calc);
			}
		}
	}
	TEST_END(output, rank_result, team_result);
}
