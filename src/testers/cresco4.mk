CC = icc
MPICC = mpicc
DEBUG = -g -O3
CFLAGS = -Wall $(NO_WARN_UNUSED) $(DEBUG)

all : testers

testers : SLGEWOS_compare SLGEWOPV_compare
  
SLGEWOS_compare : SLGEWOS_compare.c
	$(CC) $(CFLAGS) SLGEWOS_compare.c -o $(BIN_DIR)/SLGEWOS_compare -mkl -lmkl_sequential -ldl

SLGEWOPV_compare : SLGEWOPV_compare.c
	$(MPICC) $(CFLAGS) SLGEWOPV_compare.c -o $(BIN_DIR)/SLGEWOPV_compare -mkl -lmkl_scalapack_lp64 -lmkl_blacs_openmpi_lp64 -ldl

clean:
	\rm *.o *.a *~ *.so
