CC = icc
MPICC = mpiicc
DEBUG = -g -O3
CFLAGS = -Wall $(NO_WARN_UNUSED) $(DEBUG)

all: compare_DGESV compare_pDGESV compare_pDGESV_versions
.PHONY: all
  
compare_DGESV : compare_DGESV.c
	$(CC) $(CFLAGS) compare_DGESV.c -o $(BIN_DIR)/compare_DGESV -I$(MKL_INC) -L$(MKL_LIB) -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lm

compare_pDGESV : compare_pDGESV.c
	$(MPICC) $(CFLAGS) compare_pDGESV.c -o $(BIN_DIR)/compare_pDGESV  -I$(MKL_INC) -L$(MKL_LIB) -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_blacs_intelmpi_lp64 -lmkl_core -lmkl_sequential -lmkl_gf_lp64 -lm

compare_pDGESV_versions : compare_pDGESV_versions.c
	$(MPICC) $(CFLAGS) compare_pDGESV_versions.c -o $(BIN_DIR)/compare_pDGESV_versions  -I$(MKL_INC) -L$(MKL_LIB) -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_blacs_intelmpi_lp64 -lmkl_core -lmkl_sequential -lm

clean:
	\rm *.o *.a *~ *.so
