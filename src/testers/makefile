CC = gcc
MPICC = mpicc
DEBUG = -g -O3
CFLAGS = -Wall $(NO_WARN_UNUSED) $(DEBUG)

all : testers

testers : SLGEWOS-compare SLGEWOPV-compare
  
SLGEWOS-compare : SLGEWOS-compare.c
	$(CC) $(CFLAGS) SLGEWOS-compare.c -o $(BIN_DIR)/SLGEWOS-compare -llapack -lblas

SLGEWOPV-compare : SLGEWOPV-compare.c
	$(MPICC) $(CFLAGS) SLGEWOPV-compare.c -o $(BIN_DIR)/SLGEWOPV-compare /usr/lib/x86_64-linux-gnu/libmpi.so /usr/lib/x86_64-linux-gnu/libscalapack-openmpi.so.2.0 /usr/lib/x86_64-linux-gnu/libblacsCinit-openmpi.a /usr/lib/x86_64-linux-gnu/libblacs-openmpi.a

clean:
	\rm *.o *.a *~ *.so
