CC = icc
MPICC = mpicc
DEBUG = -g -O3
CFLAGS = -Wall $(NO_WARN_UNUSED) $(DEBUG)

all : testers

testers : SLGEWOS-compare SLGEWOPV-compare
  
SLGEWOS-compare : SLGEWOS-compare.c
	$(CC) $(CFLAGS) SLGEWOS-compare.c -o $(BIN_DIR)/SLGEWOS-compare -mkl -lmkl_sequential -ldl

SLGEWOPV-compare : SLGEWOPV-compare.c
	$(MPICC) $(CFLAGS) SLGEWOPV-compare.c -o $(BIN_DIR)/SLGEWOPV-compare -mkl -lmkl_scalapack_lp64 -lmkl_blacs_openmpi_lp64 -ldl

all : testers



clean:
	\rm *.o *.a *~ *.so
