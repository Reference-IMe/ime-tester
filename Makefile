# machine	= marconi|zavai|cresco6
# mpi		= intel|ch|open

PROJECT_DIR = $(CURDIR)
BIN_DIR = $(PROJECT_DIR)/bin
SRC_DIR = $(PROJECT_DIR)/src/testers
NO_WARN_UNUSED = -Wno-unused-but-set-variable -Wno-unused-variable

CC = icc

DEBUG = -g -O3
CFLAGS = -Wall $(NO_WARN_UNUSED) $(DEBUG) -DINJECT
FFLAGS = $(DEBUG)

SEQ_EXE = $(BIN_DIR)/compare_DGESV
PAR_EXE = $(BIN_DIR)/compare_all_versions # $(BIN_DIR)/compare_pDGESV $(BIN_DIR)/compare_pDGESV_versions
EXE = $(SEQ_EXE) $(PAR_EXE)


ifeq ($(mpi),)
  $(error Unspecified mpi flavour, please specify 'mpi = intel|ch|open')
endif
ifeq ($(machine),)
  $(error Unspecified machine, please specify 'machine = marconi|zavai|cresco6')
endif

ifeq ($(machine),marconi) # if marconi
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
        $(error Unknown mpi flavour '$(mpi)', please specify 'mpi = intel|ch|open')
      endif
    endif
  endif # end marconi
else
  ifeq ($(machine),zavai) # if zavai
    MACHINEFLAGS=
  else # end zavai
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
            $(error Unknown mpi flavour '$(mpi)', please specify 'mpi = intel|ch|open')
          endif # end cresco6 openmpi (and all)
        endif # end cresco6 mpich (and all)
      endif # end cresco6 intelmpi (and all)
    else # end cresco 6
      $(error Unknown machine '$(machine)', please specify 'machine = marconi|zavai|cresco6')
    endif # end cresco6 (and all)
  endif # end zavai
endif # end marconi (and all)

all: $(BIN_DIR) $(SRC_DIR)/FTLA/libftla.a $(EXE) # $(SRC_DIR)/FTLA/pdmatgen.o 
.PHONY: all

#$(SRC_DIR)/FTLA/pdmatgen.o: $(SRC_DIR)/FTLA/pdmatgen.f
#	$(MPIFC) $(FFLAGS) -c $< -o $@ $(MACHINEFLAGS)
	
$(SRC_DIR)/FTLA/libftla.a:
	cd $(SRC_DIR)/FTLA && $(MAKE)

$(SRC_DIR)/ScaLAPACK/%.o: $(SRC_DIR)/ScaLAPACK/%.f
	$(MPIFC) $(FFLAGS) -c $< -o $@ $(MACHINEFLAGS)
	
$(BIN_DIR):
	mkdir -p $(BIN_DIR)
	
$(BIN_DIR)/compare_DGESV : $(SRC_DIR)/compare_DGESV.c
	$(CC) $(CFLAGS) $(SRC_DIR)/compare_DGESV.c -o $(BIN_DIR)/compare_DGESV $(MACHINEFLAGS_ser)

$(BIN_DIR)/%: $(SRC_DIR)/%.c $(SRC_DIR)/*.h $(SRC_DIR)/../ft-simulated/*.h $(SRC_DIR)/ScaLAPACK/*.o $(SRC_DIR)/ScaLAPACK/*.f $(SRC_DIR)/ScaLAPACK/*.h $(SRC_DIR)/FTLA/libftla.a $(SRC_DIR)/FTLA/helpersftla.a
	$(MPICC) $(CFLAGS)  -lifcore $< -o $@ $(SRC_DIR)/FTLA/helpersftla.a $(SRC_DIR)/FTLA/libftla.a -L$(SRC_DIR)/ScaLAPACK  $(MACHINEFLAGS)

	
clean:
	rm -f $(EXE)
	cd $(SRC_DIR)/FTLA && $(MAKE) clean