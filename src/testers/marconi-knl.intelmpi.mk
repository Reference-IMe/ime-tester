CC = icc
MPICC = mpiicc
DEBUG = -g -O3
CFLAGS = -Wall $(NO_WARN_UNUSED) $(DEBUG)

all : testers

testers : SLGEWOS_compare SLGEWOPV_compare SLGEWOPV-FT_compare
  
SLGEWOS_compare : SLGEWOS_compare.c
	$(CC) $(CFLAGS) SLGEWOS_compare.c -o $(BIN_DIR)/SLGEWOS_compare -I$(MKL_INC) -L$(MKL_LIB) -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lm

SLGEWOPV_compare : SLGEWOPV_compare.c
	$(MPICC) $(CFLAGS) SLGEWOPV_compare.c -o $(BIN_DIR)/SLGEWOPV_compare  -I$(MKL_INC) -L$(MKL_LIB) -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_blacs_intelmpi_lp64 -lmkl_core -lmkl_sequential -lm

SLGEWOPV-FT_compare : SLGEWOPV-FT_compare.c
	$(MPICC) $(CFLAGS) SLGEWOPV-FT_compare.c -o $(BIN_DIR)/SLGEWOPV-FT_compare  -I$(MKL_INC) -L$(MKL_LIB) -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_blacs_intelmpi_lp64 -lmkl_core -lmkl_sequential -lm

clean:
	\rm *.o *.a *~ *.so
