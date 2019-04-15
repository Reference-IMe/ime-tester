CC = icc
MPICC = mpicc
DEBUG = -g -O3
CFLAGS = -Wall $(NO_WARN_UNUSED) $(DEBUG)

all : testers

testers : SLGEWOS_compare SLGEWOPV_compare SLGEWOPV-FT_compare
  
SLGEWOS_compare : SLGEWOS_compare.c
	$(CC) $(CFLAGS) SLGEWOS_compare.c -o $(BIN_DIR)/SLGEWOS_compare -L${SCALAPACK_LIB} -lscalapack -L${LAPACK_LIB} -llapack -L${BLAS_LIB} -lblas -lifcore -lm

SLGEWOPV_compare : SLGEWOPV_compare.c
	$(MPICC) $(CFLAGS) SLGEWOPV_compare.c -o $(BIN_DIR)/SLGEWOPV_compare  -L${SCALAPACK_LIB} -lscalapack -L${LAPACK_LIB} -llapack -L${BLAS_LIB} -lblas -lifcore -lm

SLGEWOPV-FT_compare : SLGEWOPV-FT_compare.c
	$(MPICC) $(CFLAGS) SLGEWOPV-FT_compare.c -o $(BIN_DIR)/SLGEWOPV-FT_compare  -L${SCALAPACK_LIB} -lscalapack -L${LAPACK_LIB} -llapack -L${BLAS_LIB} -lblas -lifcore -lm

clean:
	\rm *.o *.a *~ *.so
