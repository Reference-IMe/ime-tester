
MACHINELIST = marconi|galileo|ubuntu|cresco6
MPIFLAVOURS = intel|ch|open
PROJECT_DIR = $(CURDIR)
BIN_DIR = $(PROJECT_DIR)/bin
SRC_DIR = $(PROJECT_DIR)/src
TST_DIR = $(SRC_DIR)/testers

NO_WARN_UNUSED = -Wno-unused-but-set-variable -Wno-unused-variable

CC = icc

DEBUG = -g -O3
CFLAGS = $(DEBUG) -DINJECT -Wall -w3 -wd1418 -wd2259 # $(NO_WARN_UNUSED)
FFLAGS = $(DEBUG)

ifeq ($(mpi),)
  $(error Unspecified mpi flavour, please specify 'mpi = $(MPIFLAVOURS)')
endif
ifeq ($(machine),)
  $(error Unspecified machine, please specify 'machine = $(MACHINELIST)')
endif

#ifeq ($(machine),marconi) # if marconi
ifeq ($(machine),$(filter $(machine),marconi galileo))
  ifeq ($(mpi),intel)
    # marconi intelmpi
    MPICC = mpiicc
    MPIFC = mpiifort
    MACHINEFLAGS_ser = -I$(MKL_INC) -L$(MKL_LIB) -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lm
    MACHINEFLAGS = -I$(MKL_INC) -L$(MKL_LIB) -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_blacs_intelmpi_lp64 -lmkl_core -lmkl_sequential -lmkl_gf_lp64 -lm
  else
    ifeq ($(mpi),ch)
      # marconi mpich
      MPICC = mpicc
      MACHINEFLAGS = -I$(MKL_INC) -L$(MKL_LIB) -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_blacs_intelmpi_lp64 -lmkl_core -lmkl_sequential -lm
    else
      ifeq ($(mpi),open)
        # marconi openmpi
        MPICC = mpicc
        MACHINEFLAGS= -L${SCALAPACK_LIB} -lscalapack -L${LAPACK_LIB} -llapack -L${BLAS_LIB} -lblas -lifcore -lm
      else
        $(error Unknown mpi flavour '$(mpi)', please specify 'mpi = $(MPIFLAVOURS)')
      endif
    endif
  endif # end marconi
else
  ifeq ($(machine),ubuntu) # if ubuntu
    MACHINEFLAGS=
  else # end ubuntu
    ifeq ($(machine),cresco6) # if cresco6
      ifeq ($(mpi),intel) # cresco6 intelmpi
        # cresco6 intelmpi
        MPICC = mpiicc
        MPIFC = mpiifort
        MACHINEFLAGS_ser= -I${MKLROOT}/include -L${MKLROOT}/lib -mkl -ldl -lm
        MACHINEFLAGS= -I${MKLROOT}/include -L${MKLROOT}/lib -mkl -lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64 -ldl -lm
      else # end cresco6 intel
        ifeq ($(mpi),ch) # cresco6 mpich
          # marconi mpich
          MPICC = mpicc
          MACHINEFLAGS= 
        else # end cresco6 mpich
          ifeq ($(mpi),open) # cresco6 openmpi
            # marconi openmpi
            MPICC = mpicc
            MACHINEFLAGS=
          else # end cresco6 openmpi
            $(error Unknown mpi flavour '$(mpi)', please specify 'mpi = $(MPIFLAVOURS)')
          endif # end cresco6 openmpi (and all)
        endif # end cresco6 mpich (and all)
      endif # end cresco6 intelmpi (and all)
    else # end cresco 6
      $(error Unknown machine '$(machine)', please specify 'machine = $(MACHINELIST)')
    endif # end cresco6 (and all)
  endif # end ubuntu
endif # end marconi (and all)

SEQ_EXE = $(BIN_DIR)/compare_DGESV
PAR_EXE = $(BIN_DIR)/compare_all_versions $(BIN_DIR)/compare_checkpointing $(BIN_DIR)/compare_svxk # $(BIN_DIR)/compare_pDGESV $(BIN_DIR)/compare_pDGESV_versions
EXE = $(SEQ_EXE) $(PAR_EXE)

PAR_STD_DEP = $(SRC_DIR)/*.h $(TST_DIR)/test_*.h $(TST_DIR)/tester_head_p.c $(TST_DIR)/tester_shoulder_p.c $(TST_DIR)/tester_tail_p.c $(SRC_DIR)/helpers/*.h
SEQ_STD_DEP = $(SRC_DIR)/*.h $(TST_DIR)/tester_head.c $(TST_DIR)/tester_shoulder.c $(TST_DIR)/tester_tail.c $(SRC_DIR)/helpers/*.h

all: $(BIN_DIR) $(TST_DIR)/FTLA/libftla.a $(EXE) # $(SRC_DIR)/FTLA/pdmatgen.o 
.PHONY: all

$(BIN_DIR):
	mkdir -p $(BIN_DIR)

#$(SRC_DIR)/FTLA/pdmatgen.o: $(SRC_DIR)/FTLA/pdmatgen.f
#	$(MPIFC) $(FFLAGS) -c $< -o $@ $(MACHINEFLAGS)
	
$(TST_DIR)/FTLA/libftla.a:
	cd $(TST_DIR)/FTLA && $(MAKE)

$(TST_DIR)/ScaLAPACK.mod/%.o: $(TST_DIR)/ScaLAPACK.mod/%.f
	$(MPIFC) $(FFLAGS) -c $< -o $@ $(MACHINEFLAGS)
	
$(BIN_DIR)/compare_DGESV : $(TST_DIR)/compare_DGESV.c $(TST_DIR)/test_IMe_DGESV.h $(TST_DIR)/test_GaussianElimination_GE.h $(TST_DIR)/test_LAPACK_DGESV.h
	$(CC) $(CFLAGS) $(TST_DIR)/compare_DGESV.c -o $(BIN_DIR)/compare_DGESV $(MACHINEFLAGS_ser)

$(BIN_DIR)/compare_all_versions: $(TST_DIR)/compare_all_versions.c $(PAR_STD_DEP) $(SRC_DIR)/ft-simulated/*.h $(TST_DIR)/ScaLAPACK.mod/*.h $(TST_DIR)/ScaLAPACK.mod/pdgetrf_cp.o $(TST_DIR)/ScaLAPACK.mod/pdgeqrf_cp.o $(TST_DIR)/ScaLAPACK/*.h $(TST_DIR)/FTLA.mod/*.h $(TST_DIR)/FTLA/libftla.a
	$(MPICC) $(CFLAGS) -lifcore $< -o $(BIN_DIR)/compare_all_versions $(TST_DIR)/ScaLAPACK.mod/pdgetrf_cp.o $(TST_DIR)/ScaLAPACK.mod/pdgeqrf_cp.o $(TST_DIR)/FTLA/libftla.a $(TST_DIR)/FTLA/helpersftla.a -L$(TST_DIR)/ScaLAPACK $(MACHINEFLAGS)

$(BIN_DIR)/compare_checkpointing: $(TST_DIR)/compare_checkpointing.c $(PAR_STD_DEP) $(TST_DIR)/ScaLAPACK.mod/*.h $(TST_DIR)/ScaLAPACK.mod/pdgetrf_cp.o $(TST_DIR)/ScaLAPACK.mod/pdgeqrf_cp.o $(TST_DIR)/ScaLAPACK/*.h
	$(MPICC) $(CFLAGS) -lifcore $< -o $(BIN_DIR)/compare_checkpointing $(TST_DIR)/ScaLAPACK.mod/pdgetrf_cp.o $(TST_DIR)/ScaLAPACK.mod/pdgeqrf_cp.o -L$(TST_DIR)/ScaLAPACK $(MACHINEFLAGS)

$(BIN_DIR)/compare_svxk: $(TST_DIR)/compare_svxk.c $(PAR_STD_DEP) $(TST_DIR)/../*.h $(TST_DIR)/*.h $(TST_DIR)/ScaLAPACK/*.h
	$(MPICC) $(CFLAGS) -lifcore $< -o $(BIN_DIR)/compare_svxk -L$(TST_DIR)/ScaLAPACK $(MACHINEFLAGS)

#$(BIN_DIR)/%: $(SRC_DIR)/%.c $(SRC_DIR)/*.h $(SRC_DIR)/../ft-simulated/*.h $(SRC_DIR)/ScaLAPACK/*.o $(SRC_DIR)/ScaLAPACK/*.f $(SRC_DIR)/ScaLAPACK/*.h $(SRC_DIR)/FTLA/libftla.a $(SRC_DIR)/FTLA/helpersftla.a
#	$(MPICC) $(CFLAGS)  -lifcore $< -o $@ $(SRC_DIR)/FTLA/helpersftla.a $(SRC_DIR)/FTLA/libftla.a -L$(SRC_DIR)/ScaLAPACK  $(MACHINEFLAGS)

clean:
	rm -f $(EXE)
	rm -f $(TST_DIR)/ScaLAPACK.mod/*.o
	cd $(TST_DIR)/FTLA && $(MAKE) clean