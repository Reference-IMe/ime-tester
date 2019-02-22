
CC = icc
MPICC = mpicc
DEBUG = -g -O3
CFLAGS = -Wall $(NO_WARN_UNUSED) $(DEBUG)

all : testers

testers : SLGEWOS-compare SLGEWOPV-compare
  
SLGEWOS-compare : SLGEWOS-compare.c
        $(CC) $(CFLAGS) SLGEWOS-compare.c -o $(BIN_DIR)/SLGEWOS-compare -lifcore -I$(LAPACK_INC) -L$(LAPACK_LIB) -llapack -I$(BLAS_INC) -L$(BLAS_LIB) -lblas

SLGEWOPV-compare : SLGEWOPV-compare.c
        $(MPICC) $(CFLAGS) SLGEWOPV-compare.c -o $(BIN_DIR)/SLGEWOPV-compare  -I$(MKL_INC) -L$(MKL_LIB) -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_blacs_intelmpi_lp64 -lmkl_core -lmkl_sequential -lm

clean:
        \rm *.o *.a *~ *.so
