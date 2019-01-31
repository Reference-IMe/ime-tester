CC = gcc
MPICC = mpicc
DEBUG = -g -O3
CFLAGS = -Wall $(NO_WARN_UNUSED) $(DEBUG)

all : testers

testers : SLGEWOS-compare SLGEWOPV-compare
  
SLGEWOS-compare : SLGEWOS-compare.c
	$(CC) $(CFLAGS) SLGEWOS-compare.c -o $(BIN_DIR)/SLGEWOS-compare -llapack -lblas

SLGEWOPV-compare : SLGEWOPV-compare.c
	$(MPICC) $(CFLAGS) SLGEWOPV-compare.c -o $(BIN_DIR)/SLGEWOPV-compare

clean:
	\rm *.o *.a *~ *.so
