CC = icc
MPICC = mpiicc
DEBUG = -g -O3
CFLAGS = -Wall $(NO_WARN_UNUSED) $(DEBUG)

all : testers

testers : SLGEWOS-compare SLGEWOPV-compare
  
SLGEWOS-compare : SLGEWOS-compare.c
	$(CC) $(CFLAGS) SLGEWOS-compare.c -o $(BIN_DIR)/SLGEWOS-compare -I$(MKL_INC) -L$(MKL_LIB) -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lm


SLGEWOPV-compare : SLGEWOPV-compare.c
	$(MPICC) $(CFLAGS) SLGEWOPV-compare.c -o $(BIN_DIR)/SLGEWOPV-compare  -I$(MKL_INC) -L$(MKL_LIB) -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_blacs_intelmpi_lp64 -lmkl_core -lmkl_sequential -lm

clean:
	\rm *.o *.a *~ *.so
